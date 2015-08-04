#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "StEvent/StDcaGeometry.h"
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "StRoot/StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoD0EventMaker/StPicoD0Event.h"
#include "StPicoD0EventMaker/StKaonPion.h"
#include "StPicoDstarAnaMaker.h"
#include "StPicoDstarAnaHists.h"
#include "StAnaCuts.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StEventPlane/StEventPlane.h"

#include "StMemStat.h"

ClassImp(StPicoDstarAnaMaker)

StPicoDstarAnaMaker::StPicoDstarAnaMaker(char const * name, TString const inputFilesList,
					 TString const outFileBaseName, StPicoDstMaker* picoDstMaker, StRefMultCorr* grefmultCorrUtil, 
					  StEventPlane* eventPlaneMaker):
   StMaker(name), mPicoDstMaker(picoDstMaker), mPicoD0Event(NULL), mGRefMultCorrUtil(grefmultCorrUtil),
  mEventPlane(eventPlaneMaker), mInputFilesList(inputFilesList), mOutFileBaseName(outFileBaseName),
   mChain(NULL), mEventCounter(0),mHists(NULL)
{}

Int_t StPicoDstarAnaMaker::Init()
{
   mPicoD0Event = new StPicoD0Event();

   mChain = new TChain("T");
   std::ifstream listOfFiles(mInputFilesList.Data());
   if (listOfFiles.is_open())
   {
      std::string file;
      while (getline(listOfFiles, file))
      {
         LOG_INFO << "StPicoDstarAnaMaker - Adding :" << file << endm;
         mChain->Add(file.c_str());
      }
   }
   else
   {
      LOG_ERROR << "StPicoDstarAnaMaker - Could not open list of files. ABORT!" << endm;
      return kStErr;
   }

   mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
   mChain->SetBranchAddress("dEvent", &mPicoD0Event);

   // -------------- USER VARIABLES -------------------------
   mHists = new StPicoDstarAnaHists(mOutFileBaseName);

   return kStOK;
}
//-----------------------------------------------------------------------------
StPicoDstarAnaMaker::~StPicoDstarAnaMaker()
{
   delete mGRefMultCorrUtil;
}
//-----------------------------------------------------------------------------
Int_t StPicoDstarAnaMaker::Finish()
{
   mHists->closeFile();
   return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoDstarAnaMaker::Make()
{
   StMemStat mem;
   readNextEvent();

   if (!mPicoDstMaker)
   {
      LOG_WARN << " StPicoDstarAnaMaker - No PicoDstMaker! Skip! " << endm;
      return kStWarn;
   }

   StPicoDst const* picoDst = mPicoDstMaker->picoDst();

   if (!picoDst)
   {
      LOG_WARN << "StPicoDstarAnaMaker - No PicoDst! Skip! " << endm;
      return kStWarn;
   }

   if (mPicoD0Event->runId() != picoDst->event()->runId() ||
         mPicoD0Event->eventId() != picoDst->event()->eventId())
   {
      LOG_ERROR << " StPicoDstarAnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!" << endm;
      LOG_ERROR << " StPicoDstarAnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync." << endm;
      exit(1);
   }
   // -------------- USER ANALYSIS -------------------------

   mHists->addEventBeforeCut(picoDst->event());
   if (isGoodEvent(picoDst->event()))
   {
     TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();
     if (aKaonPion->GetEntries()) mHists->addEvent(picoDst->event());
     
     StThreeVectorF const pVtx = picoDst->event()->primaryVertex();
     StThreeVectorF const pVtxErr = picoDst->event()->primaryVertexError();
     

      if (!mGRefMultCorrUtil)
      {
         LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
         return kStWarn;
      }
      mGRefMultCorrUtil->init(picoDst->event()->runId());
      mGRefMultCorrUtil->initEvent(picoDst->event()->grefMult(), pVtx.z(), picoDst->event()->ZDCx()) ;
      if (mGRefMultCorrUtil->isBadRun(picoDst->event()->runId()))
      {
         //cout<<"This is a bad run from mGRefMultCorrUtil! Skip! " << endl;
         return kStOK;
      }

      int centrality  = mGRefMultCorrUtil->getCentralityBin9();
      const double reweight = mGRefMultCorrUtil->getWeight();
      const double refmultCor = mGRefMultCorrUtil->getRefMultCorr();
      mHists->addCent(refmultCor, centrality, reweight, pVtx.z());

      if (!loadEventPlaneCorr(mEventPlane)){
	LOG_WARN << "Event plane calculations unavalable! Skipping" << endm;
	//mFailedRunnumber = picoDst->event()->runId();
	return kStOK;
      }
     float const eventPlane = mEventPlane->getEventPlane();
     int const eventPlane_bin = (int)((eventPlane) / 0.3141592) ;
     
     mHists->addEventPlane(centrality, eventPlane, mEventPlane->getResolutionRandom(), mEventPlane->getResolutionEta());
     //if (eventPlane_bin < 0  ||  eventPlane_bin > 9) return kStOk;
      //Get array of soft pions and basic QA for production
      std::vector<int> softPions;

      UInt_t nTracks = picoDst->numberOfTracks();
      for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack)
      {
         StPicoTrack const* trk = picoDst->track(iTrack);
         if (!trk) continue;
	 //Soft pion array
	 if( isSoftPion(trk,pVtx) )
	   softPions.push_back(iTrack);
         StPhysicalHelixD helix = trk->helix();
         float dca = float(helix.geometricSignedDistance(pVtx));
         StThreeVectorF momentum = trk->gMom(pVtx, picoDst->event()->bField());

         bool tofMatch = getTofBeta(trk, &pVtx) > 0;

         int tofMatchFlag =  tofMatch ? 1 : 0 ;
         int hftMatchFlag =  trk->isHFTTrack() ? 1 : 0 ;

         if (!isGoodQaTrack(trk, momentum, dca)) continue;

         StThreeVectorF dcaPoint = helix.at(helix.pathLength(pVtx.x(), pVtx.y()));
         float dcaZ = dcaPoint.z() - pVtx.z();
         double dcaXy = helix.geometricSignedDistance(pVtx.x(), pVtx.y());

         bool isPion = kFALSE;
         bool isKaon = kFALSE;
         if (fabs(trk->nSigmaPion()) < anaCuts::nSigmaPion)  isPion = kTRUE;
         if (fabs(trk->nSigmaKaon()) < anaCuts::nSigmaKaon)  isKaon = kTRUE;
         if (trk && tofMatch && fabs(dca) < 1.0 && trk->isHFTTrack() && (isPion || isKaon))
         {
            mHists->addDcaPtCent(dca, dcaXy, dcaZ, isPion, isKaon, momentum.perp(), centrality, momentum.pseudoRapidity(), momentum.phi(), pVtx.z(), picoDst->event()->ZDCx() / 1000.); //add Dca distribution
         }
         if (trk && tofMatch && fabs(dca) < 1.5 && (isPion || isKaon))
         {
            mHists->addTpcDenom1(isPion, isKaon, momentum.perp(), centrality, momentum.pseudoRapidity(), momentum.phi(), pVtx.z(), picoDst->event()->ZDCx() / 1000.); //Dca cut on 1.5cm, add Tpc Denominator
         }
         if (trk && tofMatch && fabs(dca) < 1.5 && trk->isHFTTrack() && (isPion || isKaon))
         {
            mHists->addHFTNumer1(isPion, isKaon, momentum.perp(), centrality,  momentum.pseudoRapidity(), momentum.phi(), pVtx.z(), picoDst->event()->ZDCx() / 1000.); //Dca cut on 1.5cm, add HFT Numerator
         }
      } // .. end tracks loop

      for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
      {
         StKaonPion const* kp = (StKaonPion*)aKaonPion->UncheckedAt(idx);

         if (!isGoodPair(kp)) continue;

	 int trkIndex_d0[2] = { kp->kaonIdx(), kp->pionIdx()};
	 float psi_d0 = mEventPlane->getEventPlane(2, trkIndex_d0);
	 float dPhi_d0 = kp->phi() - psi_d0;
	 
	 while(dPhi_d0 < 0) dPhi_d0 += TMath::Pi();
	 while(dPhi_d0 >= TMath::Pi()) dPhi_d0 -= TMath::Pi();	 
	 
	 mHists->addDzero(kp,dPhi_d0, centrality, reweight);
         StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
         StPicoTrack const* pion = picoDst->track(kp->pionIdx());

         if (!isGoodTrack(kaon) || !isGoodTrack(pion)) continue;

         bool tpcPion = isTpcPion(pion);
         bool tpcKaon = isTpcKaon(kaon);
         float pBeta = getTofBeta(pion, &pVtx);
         float kBeta = getTofBeta(kaon, &pVtx);
         bool pTofAvailable = pBeta > 0;
         bool kTofAvailable = kBeta > 0;
         bool tofPion = isTofPion(pion, pBeta);
         bool tofKaon = isTofKaon(kaon, kBeta);

         bool goodPion = (pTofAvailable && tofPion) || (!pTofAvailable && tpcPion);
         bool goodKaon = (kTofAvailable && tofKaon) || (!kTofAvailable && tpcKaon);
         bool tof = goodPion && goodKaon;
         bool tpc = tpcPion && tpcKaon;

         if (tpc || tof)
         {
            bool unlike = kaon->charge() * pion->charge() < 0 ? true : false;
            mHists->addKaonPion(kp, unlike, tpc, tof, centrality, reweight);
         }
	 //End of D0 calculaitons
	 for(int iPion = 0; iPion < softPions.size(); ++iPion)
	 {
	   StPicoTrack const* softpion = picoDst->track( softPions.at(iPion) );
	   if( isGoodDstar(softpion,kp,pVtx) )
	   {
	     StPhysicalHelixD helix = softpion->dcaGeometry().helix();
	     helix.moveOrigin(helix.pathLength(pVtx));
	     StThreeVectorF const softPionMom(softpion->pMom().x(), softpion->pMom().y(), softpion->pMom().z());
	     StLorentzVectorF const pionFourMom(softPionMom,softPionMom.massHypothesis(M_PION_PLUS));
	     StLorentzVectorF const Dstar = pionFourMom + (kp->lorentzVector());
	     
	     //Event plane stuff
	     int trkIndex[3] = { kp->kaonIdx(), kp->pionIdx(), softPions.at(iPion)};
	     float psi = mEventPlane->getEventPlane(3, trkIndex);
	     float dPhi = Dstar.phi() - psi;
      
	     while(dPhi < 0) dPhi += TMath::Pi();
	     while(dPhi >= TMath::Pi()) dPhi -= TMath::Pi();

	     const bool isHft =  softpion->isHFTTrack();

	     mHists->addDstar(Dstar, *kp, dPhi, centrality, reweight, isHft);
	     //test
	     mHists->addDstarNtuple(Dstar, *softpion,  *kp);
	   }
	 }  
      } // end of kaonPion loop
   } // end of isGoodEvent

   return kStOK;
}
//-----------------------------------------------------------------------------
int StPicoDstarAnaMaker::getD0PtIndex(StKaonPion const * kp) const
{
   for (int i = 0; i < anaCuts::nPtBins; i++)
   {
      if ((kp->pt() >= anaCuts::PtBinsEdge[i]) && (kp->pt() < anaCuts::PtBinsEdge[i + 1]))
         return i;
   }
   return anaCuts::nPtBins - 1;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isGoodEvent(StPicoEvent const * const picoEvent) const
{
   return (picoEvent->triggerWord() & anaCuts::triggerWord) &&
          fabs(picoEvent->primaryVertex().z()) < anaCuts::vz &&
          fabs(picoEvent->primaryVertex().z() - picoEvent->vzVpd()) < anaCuts::vzVpdVz &&
          !(fabs(picoEvent->primaryVertex().x()) < anaCuts::Verror && fabs(picoEvent->primaryVertex().y()) < anaCuts::Verror && fabs(picoEvent->primaryVertex().z()) < anaCuts::Verror) &&
          sqrt(TMath::Power(picoEvent->primaryVertex().x(), 2) + TMath::Power(picoEvent->primaryVertex().y(), 2)) <=  anaCuts::Vrcut;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isGoodQaTrack(StPicoTrack const * const trk, StThreeVectorF const momentum, const double dca) const
{
   return trk->gPt() > anaCuts::qaGPt && trk->nHitsFit() >= anaCuts::qaNHitsFit;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isGoodTrack(StPicoTrack const * const trk) const
{
   StThreeVectorF mom = trk->gMom(mPicoDstMaker->picoDst()->event()->primaryVertex(), mPicoDstMaker->picoDst()->event()->bField());

   return trk->gPt() > anaCuts::minPt &&
          trk->nHitsFit() >= anaCuts::nHitsFit &&
          fabs(mom.pseudoRapidity()) <= anaCuts::Eta;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isTpcPion(StPicoTrack const * const trk) const
{
   return fabs(trk->nSigmaPion()) < anaCuts::nSigmaPion;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isTpcKaon(StPicoTrack const * const trk) const
{
   return fabs(trk->nSigmaKaon()) < anaCuts::nSigmaKaon;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isGoodPair(StKaonPion const* const kp) const
{
   int tmpIndex = getD0PtIndex(kp);
   return cos(kp->pointingAngle()) > anaCuts::cosTheta[tmpIndex] &&
          kp->pionDca() > anaCuts::pDca[tmpIndex] && kp->kaonDca() > anaCuts::kDca[tmpIndex] &&
          kp->dcaDaughters() < anaCuts::dcaDaughters[tmpIndex] &&
          kp->decayLength() > anaCuts::decayLength[tmpIndex] &&
          fabs(kp->lorentzVector().rapidity()) < anaCuts::RapidityCut &&
          ((kp->decayLength()) * sin(kp->pointingAngle())) < anaCuts::dcaV0ToPv[tmpIndex];
}
bool StPicoDstarAnaMaker::isGoodDstar(StPicoTrack const * const softPion, StKaonPion const * const Dzero, StThreeVectorF const &pVtx) const
{
  StPhysicalHelixD helix = softPion->dcaGeometry().helix();
  helix.moveOrigin(helix.pathLength(pVtx));
  StThreeVectorF const softPionMom(softPion->pMom().x(), softPion->pMom().y(), softPion->pMom().z());
  StLorentzVectorF const pionFourMom(softPionMom,softPionMom.massHypothesis(M_PION_PLUS));
  StLorentzVectorF const Dstar = pionFourMom + (Dzero->lorentzVector());
  //float const pT = sqrt(pow(Dstar.px(), 2.0) + pow(Dstar.py(), 2.0) );
  //if( (Dstar.m()-Dzero->m() ) >0.175 || (Dstar.m()-Dzero->m() ) < 0.115) continue;
  return (Dstar.m()-Dzero->m() )>anaCuts::dstarMinMass && 
         (Dstar.m()-Dzero->m() )<anaCuts::dstarMaxMass;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isTofKaon(StPicoTrack const * const trk, float beta) const
{
   bool tofKaon = false;

   if (beta > 0)
   {
      double ptot = trk->dcaGeometry().momentum().mag();
      float beta_k = ptot / sqrt(ptot * ptot + M_KAON_PLUS * M_KAON_PLUS);
      tofKaon = fabs(1 / beta - 1 / beta_k) < anaCuts::kTofBetaDiff ? true : false;
   }

   return tofKaon;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isTofPion(StPicoTrack const * const trk, float beta) const
{
   bool tofPion = false;

   if (beta > 0)
   {
      double ptot = trk->dcaGeometry().momentum().mag();
      float beta_pi = ptot / sqrt(ptot * ptot + M_PION_PLUS * M_PION_PLUS);
      tofPion = fabs(1 / beta - 1 / beta_pi) < anaCuts::pTofBetaDiff ? true : false;
   }

   return tofPion;
}
//-----------------------------------------------------------------------------
float StPicoDstarAnaMaker::getTofBeta(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{
   int index2tof = trk->bTofPidTraitsIndex();

   float beta = std::numeric_limits<float>::quiet_NaN();

   if (index2tof >= 0)
   {
      StPicoBTofPidTraits *tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);

      if (tofPid)
      {
         beta = tofPid->btofBeta();

         if (beta < 1e-4)
         {
            StThreeVectorF const btofHitPos = tofPid->btofHitPos();

            StPhysicalHelixD helix = trk->helix();
            float L = tofPathLength(pVtx, &btofHitPos, helix.curvature());
            float tof = tofPid->btof();
            if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
            else beta = std::numeric_limits<float>::quiet_NaN();
         }
      }
   }

   return beta;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isSoftPion(StPicoTrack const * const trk, StThreeVectorF const &pVtx) const 
{
  //Primary momentum!
  StThreeVectorF mom = trk->pMom();

  //track
  if( trk->nHitsFit() < anaCuts::nHitsFit && 
      fabs(mom.pseudoRapidity()) > anaCuts::Eta ) return false;

  // -- good pion
  if( trk->pMom().perp()<anaCuts::softPionMinPt || !isTpcPion(trk) ) return false;
  
  //From vertex
  StPhysicalHelixD helix = trk->dcaGeometry().helix();
  helix.moveOrigin(helix.pathLength(pVtx));
  //Should be safe, already tpc primary
  if( (helix.origin() - pVtx).mag() > 3.0 ) return false;
  
  return true;
}
// _________________________________________________________
bool StPicoDstarAnaMaker::loadEventPlaneCorr(StEventPlane const * mEventPlane)
{
  //needs to implement, will currently break maker
  if (!mEventPlane)
    {
      LOG_WARN << "No EventPlane ! Skipping! " << endm;
      return kFALSE;
    }
  if (!mEventPlane->getAcceptEvent())
    {
      // LOG_WARN << "StPicoMixedEvent::THistograms and TProiles NOT found! shoudl check the input Qvector files From HaoQiu ! Skipping this run! " << endm;
      return kFALSE;
    }
  return kTRUE;
}
