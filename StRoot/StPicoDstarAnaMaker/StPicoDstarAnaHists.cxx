#include <cmath>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TString.h"
#include "../StPicoDstMaker/StPicoEvent.h"
#include "../StPicoDstMaker/StPicoTrack.h"
#include "../StPicoPrescales/StPicoPrescales.h"
#include "../StPicoD0EventMaker/StKaonPion.h"
#include "StAnaCuts.h"
#include "TNtuple.h"
#include "THnSparse.h"
#include "THn.h"

#include "StPicoDstarAnaHists.h"

//-----------------------------------------------------------------------
StPicoDstarAnaHists::StPicoDstarAnaHists(TString fileBaseName) : 
  mPrescales(NULL), mOutFile(NULL),
  mh2InvariantMassVsPt(NULL), mh2InvariantMassVsPtLike(NULL), mh2InvariantMassVsPtTof(NULL), mh2InvariantMassVsPtTofLike(NULL),
  mh1Cent(NULL), mh1CentWg(NULL), mh1gRefmultCor(NULL), mh1gRefmultCorWg(NULL), mh2CentVz(NULL), mh2CentVzWg(NULL), mh3InvariantMassVsPtVsCent(NULL), mh3InvariantMassVsPtVsCentLike(NULL), mh3InvariantMassVsPtVsCentTof(NULL), mh3InvariantMassVsPtVsCentTofLike(NULL), mh2Tpc1PtCent(NULL),  mh2Tpc1PhiVz(NULL), mh2HFT1PtCent(NULL),  mh2HFT1PhiVz(NULL),  mh3DcaXyPtCent(NULL), mh3DcaZPtCent(NULL),
  mh2EPresCentRandom(NULL), mh2EPresCentEta(NULL), mh2EventPlaneCent(NULL),
  mh3D0InvariantMassVsPtVsCent(NULL), mh3D0InvMassPtdPhi(NULL),mh3DstarInvariantMassVsPtVsCent(NULL), mh3DstarInvMassPtdPhi(NULL)
{
   mPrescales = new StPicoPrescales(anaCuts::prescalesFilesDirectoryName);

   mOutFile = new TFile(Form("%s.hists.root", fileBaseName.Data()), "RECREATE");

   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < anaCuts::nEtas; iEta++)
      {
         for (int iVz = 0; iVz < anaCuts::nVzs; iVz++)
         {
            mh3DcaXyPtCentPartEtaVz[iParticle][iEta][iVz] = NULL;
            mh3DcaZPtCentPartEtaVz[iParticle][iEta][iVz] = NULL;
            mh2Tpc1PtCentPartEtaVz[iParticle][iEta][iVz] = NULL;
            mh2HFT1PtCentPartEtaVz[iParticle][iEta][iVz] = NULL;
         }
      }
   }
   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < anaCuts::nEtas; iEta++)
      {
         for (int iVz = 0; iVz < anaCuts::nVzs; iVz++)
         {
            for (int iPhi = 0; iPhi < anaCuts::nPhis; iPhi++)
            {
               mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi] = NULL;
               mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi] = NULL;
            }
         }
      }
   }
   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iPhi = 0; iPhi < anaCuts::nPhis; iPhi++)
      {
         mh3DcaXyPtCentPartPhi[iParticle][iPhi] = NULL;
         mh3DcaZPtCentPartPhi[iParticle][iPhi] = NULL;
         mh2Tpc1PtCentPartPhi[iParticle][iPhi] = NULL;
         mh2HFT1PtCentPartPhi[iParticle][iPhi] = NULL;
      }
   }
   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iZdcx = 0; iZdcx < anaCuts::nZdcxs; iZdcx++)
      {
         mh3DcaXyPtCentPartZdcx[iParticle][iZdcx] = NULL;
         mh3DcaZPtCentPartZdcx[iParticle][iZdcx] = NULL;
         mh2Tpc1PtCentPartZdcx[iParticle][iZdcx] = NULL;
         mh2HFT1PtCentPartZdcx[iParticle][iZdcx] = NULL;
      }
   }

   int nRuns = mPrescales->numberOfRuns();
   TH1::SetDefaultSumw2();
   mh1TotalEventsInRun         = new TH1F("mh1TotalEventsInRun", "totalEventsInRun;runIndex;totalEventsInRun", nRuns + 1, 0, nRuns + 1);
   mh1TotalEventsInRunBeforeCut = new TH1F("mh1TotalEventsInRunBeforeCut", "totalEventsInRun;runIndex;totalEventsInRun", nRuns + 1, 0, nRuns + 1);
   mh2InvariantMassVsPt        = new TH2F("mh2InvariantMassVsPt", "invariantMassVsPt;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})", 120, 0, 12, 50, 1.6, 2.1);
   mh2InvariantMassVsPtLike    = new TH2F("mh2InvariantMassVsPtLike", "invariantMassVsPtLike;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})", 120, 0, 12, 50, 1.6, 2.1);
   mh2InvariantMassVsPtTof     = new TH2F("mh2InvariantMassVsPtTof", "invariantMassVsPtTof;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})", 120, 0, 12, 50, 1.6, 2.1);
   mh2InvariantMassVsPtTofLike = new TH2F("mh2InvariantMassVsPtTofLike", "invariantMassVsPtTofLike;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})", 120, 0, 12, 50, 1.6, 2.1);
   //add centrality
   mh1Cent         = new TH1F("mh1Cent", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5);
   mh1CentWg         = new TH1F("mh1CentWg", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5);
   mh1gRefmultCor  = new TH1F("mh1gRefmultCor", "gRefmultCor;gRefmult;Counts", 700, 0, 700);
   mh1gRefmultCorWg  = new TH1F("mh1gRefmultCorWg", "gRefmultCorWg;gRefmultCorWg;Counts", 700, 0, 700);
   mh2CentVz         = new TH2F("mh2CentVz", "CentralityVsVz;cent;Vz", 10, -1.5, 8.5, 200, -10, 10);
   mh2CentVzWg         = new TH2F("mh2CentVzWg", "CentralityVsVzWg;cent;Vz", 10, -1.5, 8.5, 200, -10, 10);
   mh3InvariantMassVsPtVsCent        = new TH3F("mh3InvariantMassVsPtVsCent", "invariantMassVsPtVsCent;p_{T}(K#pi)(GeV/c);Cent;m_{K#pi}(GeV/c^{2})", 120, 0, 12, 10, -1.5, 8.5, 50, 1.6, 2.1);
   mh3InvariantMassVsPtVsCentLike    = new TH3F("mh3InvariantMassVsPtVsCentLike", "invariantMassVsPtVsCentLike;p_{T}(K#pi)(GeV/c);Cent;m_{K#pi}(GeV/c^{2})", 120, 0, 12, 10, -1.5, 8.5, 50, 1.6, 2.1);
   mh3InvariantMassVsPtVsCentTof     = new TH3F("mh3InvariantMassVsPtVsCentTof", "invariantMassVsPtVsCentTof;p_{T}(K#pi)(GeV/c);Cent;m_{K#pi}(GeV/c^{2})", 120, 0, 12, 10, -1.5, 8.5, 50, 1.6, 2.1);
   mh3InvariantMassVsPtVsCentTofLike = new TH3F("mh3InvariantMassVsPtVsCentTofLike", "invariantMassVsPtVsCentTofLike;p_{T}(K#pi)(GeV/c);Cent;m_{K#pi}(GeV/c^{2})", 120, 0, 12, 10, -1.5, 8.5, 50, 1.6, 2.1);
   //Add some HFT ratio plots
   mh2Tpc1PtCent  = new TH2F("mh2Tpc1PtCent", "Tpc tacks;p_{T}(GeV/c);cent", 120, 0, 12, 10, -1.5, 8.5); //Dca 1.5cm
   mh2HFT1PtCent  = new TH2F("mh2HFT1PtCent", "HFT tacks;p_{T}(GeV/c);cent", 120, 0, 12, 10, -1.5, 8.5); //Dca 1.5cm
   mh2Tpc1PhiVz  = new TH2F("mh2Tpc1PhiVz", "Tpc tacks;#Phi;Vz", 100, -3.1415, 3.1415, 20, -10, 10); //Dca 1.5cm
   mh2HFT1PhiVz  = new TH2F("mh2HFT1PhiVz", "HFT tacks;#Phi;Vz", 100, -3.1415, 3.1415, 20, -10, 10); //Dca 1.5cm
   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < anaCuts::nEtas; iEta++)
      {
         for (int iVz = 0; iVz < anaCuts::nVzs; iVz++)
         {
            mh2Tpc1PtCentPartEtaVz[iParticle][iEta][iVz]  = new TH2F(Form("mh2Tpc1PtCentPartEtaVz_%d_%d_%d", iParticle, iEta, iVz), Form("mh2Tpc1PtCent_%s_Eta%2.1f_Vz%2.1f;p_{T}(GeV/c);cent", anaCuts::ParticleName[iParticle], anaCuts::EtaEdge[iEta], anaCuts::VzEdge[iVz]), anaCuts::nPts, anaCuts::PtEdge, anaCuts::nCents, anaCuts::CentEdge); //Dca 1.cm
            mh2HFT1PtCentPartEtaVz[iParticle][iEta][iVz]  = new TH2F(Form("mh2HFT1PtCentPartEtaVz_%d_%d_%d", iParticle, iEta, iVz), Form("mh2HFT1PtCent_%s_Eta%2.1f_Vz%2.1f;p_{T}(GeV/c);cent", anaCuts::ParticleName[iParticle], anaCuts::EtaEdge[iEta], anaCuts::VzEdge[iVz]), anaCuts::nPts, anaCuts::PtEdge, anaCuts::nCents, anaCuts::CentEdge); //Dca 1.cm
         }
      }
   }

   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < anaCuts::nEtas; iEta++)
      {
         for (int iVz = 0; iVz < anaCuts::nVzs; iVz++)
         {
            for (int iPhi = 0; iPhi < anaCuts::nPhis; iPhi++)
            {
               mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]  = new TH2F(Form("mh2Tpc1PtCentPartEtaVz_%d_%d_%d_%d", iParticle, iEta, iVz, iPhi), Form("mh2Tpc1PtCent_%s_Eta%2.1f_Vz%2.1f_Phi%2.1f;p_{T}(GeV/c);cent", anaCuts::ParticleName[iParticle], anaCuts::EtaEdge[iEta], anaCuts::VzEdge[iVz], anaCuts::PhiEdge[iPhi]), anaCuts::nPts, anaCuts::PtEdge, anaCuts::nCents, anaCuts::CentEdge); //Dca 1.cm
               mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]  = new TH2F(Form("mh2HFT1PtCentPartEtaVz_%d_%d_%d_%d", iParticle, iEta, iVz, iPhi), Form("mh2HFT1PtCent_%s_Eta%2.1f_Vz%2.1f_Phi%2.1f;p_{T}(GeV/c);cent", anaCuts::ParticleName[iParticle], anaCuts::EtaEdge[iEta], anaCuts::VzEdge[iVz], anaCuts::PhiEdge[iPhi]), anaCuts::nPts, anaCuts::PtEdge, anaCuts::nCents, anaCuts::CentEdge); //Dca 1.cm
            }
         }
      }
   }

   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iPhi = 0; iPhi < anaCuts::nPhis; iPhi++)
      {
         mh2Tpc1PtCentPartPhi[iParticle][iPhi]  = new TH2F(Form("mh2Tpc1PtCentPartPhi_%d_%d", iParticle, iPhi), Form("mh2Tpc1PtCent_%s_Phi%2.1f;p_{T}(GeV/c);cent", anaCuts::ParticleName[iParticle], anaCuts::PhiEdge[iPhi]), anaCuts::nPts, anaCuts::PtEdge, anaCuts::nCents, anaCuts::CentEdge); //Dca 1.cm
         mh2HFT1PtCentPartPhi[iParticle][iPhi]  = new TH2F(Form("mh2HFT1PtCentPartPhi_%d_%d", iParticle, iPhi), Form("mh2HFT1PtCent_%s_Phi%2.1f;p_{T}(GeV/c);cent", anaCuts::ParticleName[iParticle], anaCuts::PhiEdge[iPhi]), anaCuts::nPts, anaCuts::PtEdge, anaCuts::nCents, anaCuts::CentEdge); //Dca 1.cm
      }
   }

   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iZdcx = 0; iZdcx < anaCuts::nZdcxs; iZdcx++)
      {
         mh2Tpc1PtCentPartZdcx[iParticle][iZdcx]  = new TH2F(Form("mh2Tpc1PtCentPartZdcx_%d_%d", iParticle, iZdcx), Form("mh2Tpc1PtCent_%s_Zdcx%2.1f;p_{T}(GeV/c);cent", anaCuts::ParticleName[iParticle], anaCuts::ZdcxEdge[iZdcx]), anaCuts::nPts, anaCuts::PtEdge, anaCuts::nCents, anaCuts::CentEdge); //Dca 1.cm
         mh2HFT1PtCentPartZdcx[iParticle][iZdcx]  = new TH2F(Form("mh2HFT1PtCentPartZdcx_%d_%d", iParticle, iZdcx), Form("mh2HFT1PtCent_%s_Zdcx%2.1f;p_{T}(GeV/c);cent", anaCuts::ParticleName[iParticle], anaCuts::ZdcxEdge[iZdcx]), anaCuts::nPts, anaCuts::PtEdge, anaCuts::nCents, anaCuts::CentEdge); //Dca 1.cm
      }
   }

   // Add some QA HFT , Dca, resolution

   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < anaCuts::nEtas; iEta++)
      {
         for (int iVz = 0; iVz < anaCuts::nVzs; iVz++)
         {
            mh3DcaXyPtCentPartEtaVz[iParticle][iEta][iVz]  = new TH3F(Form("mh3DcaXyPtCentPartEtaVz_%d_%d_%d", iParticle, iEta, iVz), Form("mh3DcaXyPtCent_%s_Eta%2.1f_Vz%2.1f;p_{T}(GeV/c);cent;DcaXy(cm)", anaCuts::ParticleName[iParticle], anaCuts::EtaEdge[iEta], anaCuts::VzEdge[iVz]), anaCuts::nPts, anaCuts::PtEdge, anaCuts::nCents, anaCuts::CentEdge, anaCuts::nDcas, anaCuts::DcaEdge); //Dca 1.cm
            mh3DcaZPtCentPartEtaVz[iParticle][iEta][iVz]  = new TH3F(Form("mh3DcaZPtCentPartEtaVz_%d_%d_%d", iParticle, iEta, iVz), Form("mh3DcaZPtCent_%s_Eta%2.1f_Vz%2.1f;p_{T}(GeV/c);cent;DcaZ(cm)", anaCuts::ParticleName[iParticle], anaCuts::EtaEdge[iEta], anaCuts::VzEdge[iVz]), anaCuts::nPts, anaCuts::PtEdge, anaCuts::nCents, anaCuts::CentEdge, anaCuts::nDcas, anaCuts::DcaEdge); //Dca 1.cm
         }
      }
   }

   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iPhi = 0; iPhi < anaCuts::nPhis; iPhi++)
      {
         mh3DcaXyPtCentPartPhi[iParticle][iPhi]  = new TH3F(Form("mh3DcaXyPtCentPartPhi_%d_%d", iParticle, iPhi), Form("mh3DcaXyPtCent_%s_Phi%2.1f;p_{T}(GeV/c);cent;DcaXy(cm)", anaCuts::ParticleName[iParticle], anaCuts::PhiEdge[iPhi]), anaCuts::nPts, anaCuts::PtEdge, anaCuts::nCents, anaCuts::CentEdge, anaCuts::nDcas, anaCuts::DcaEdge); //Dca 1.cm
         mh3DcaZPtCentPartPhi[iParticle][iPhi]  = new TH3F(Form("mh3DcaZPtCentPartPhi_%d_%d", iParticle, iPhi), Form("mh3DcaZPtCent_%s_Phi%2.1f;p_{T}(GeV/c);cent;DcaZ(cm)", anaCuts::ParticleName[iParticle], anaCuts::PhiEdge[iPhi]), anaCuts::nPts, anaCuts::PtEdge, anaCuts::nCents, anaCuts::CentEdge, anaCuts::nDcas, anaCuts::DcaEdge); //Dca 1.cm
      }
   }

   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iZdcx = 0; iZdcx < anaCuts::nZdcxs; iZdcx++)
      {
         mh3DcaXyPtCentPartZdcx[iParticle][iZdcx]  = new TH3F(Form("mh3DcaXyPtCentPartZdcx_%d_%d", iParticle, iZdcx), Form("mh3DcaXyPtCent_%s_Zdcx%2.1f;p_{T}(GeV/c);cent;DcaXy(cm)", anaCuts::ParticleName[iParticle], anaCuts::ZdcxEdge[iZdcx]), anaCuts::nPts, anaCuts::PtEdge, anaCuts::nCents, anaCuts::CentEdge, anaCuts::nDcas, anaCuts::DcaEdge); //Dca 1.cm
         mh3DcaZPtCentPartZdcx[iParticle][iZdcx]  = new TH3F(Form("mh3DcaZPtCentPartZdcx_%d_%d", iParticle, iZdcx), Form("mh3DcaZPtCent_%s_Zdcx%2.1f;p_{T}(GeV/c);cent;DcaZ(cm)", anaCuts::ParticleName[iParticle], anaCuts::ZdcxEdge[iZdcx]), anaCuts::nPts, anaCuts::PtEdge, anaCuts::nCents, anaCuts::CentEdge, anaCuts::nDcas, anaCuts::DcaEdge); //Dca 1.cm
      }
   }

   mh3DcaXyPtCent  = new TH3F("mh3DcaXyPtCent", "mh3DcaXyPtCent;p_{T}(GeV/c);cent;DcaXy(cm)", 120, 0, 12, 10, -1.5, 8.5, 1000, -1, 1); //Dca 1.cm
   mh3DcaZPtCent  = new TH3F("mh3DcaZPtCent", "mh3DcaZPtCent;p_{T}(GeV/c);cent;DcaZ(cm)", 120, 0, 12, 10, -1.5, 8.5, 1000, -1, 1); //Dca 1.cm
   //DSTAR
   mh2EPresCentRandom = new TH2F("hEPresCentRandom","EPresVsCentEta",10,-1.5,8.5,100,-1,1);
   mh2EPresCentEta = new TH2F("hEPresCentEta","EPresVsCentEta",10,-1.5,8.5,100,-1,1);
   mh2EventPlaneCent = new TH2F("hEventPlaneCent","EventPlanevsCent",10,-1.5,8.5,100,0,M_PI);
   //all
   mh3D0InvariantMassVsPtVsCent = new TH3F("mh3D0InvariantMassVsPtVsCent", "invariantMassVsPtVsCent;p_{T}(K#pi)(GeV/c);Cent;m_{K#pi}(GeV/c^{2})", 120, 0, 12, 10, -1.5, 8.5, 50, 1.6, 2.1);
   mh3D0InvMassPtdPhi = new TH3F("hD0InvMassPtdPhi","D0VsInvMassVsPtVsdPhi;p_{T}(K#pi)(GeV/c);dPhi;m_{K#pi}(GeV/c^{2})", 120, 0, 12,10,0,M_PI,50, 1.6, 2.1);
   mh3DstarInvariantMassVsPtVsCent= new TH3F("mh3DstarInvariantMassVsPtVsCent", "invariantMassVsPtVsCent;p_{T}(K#pi#pi)(GeV/c);Cent;m_{K#pi#pi - K#pi}(GeV/c^{2})", 120, 0, 12, 10, -1.5, 8.5, 50, 0.138, 0.158);
   mh3DstarInvMassPtdPhi = new TH3F("hDstarInvMassPtdPhi","DstarVsInvMassVsPtVsdPhi;p_{T}(K#pi#pi)(GeV/c);dPhi;m_{K#pi#pi-K#pi}(GeV/c^{2})", 120, 0, 12,10,0,M_PI,50,0.138,0.158 );
   //HFT
   mh3DstarInvariantMassVsPtVsCent_hft= new TH3F("mh3DstarInvariantMassVsPtVsCent_hft", "invariantMassVsPtVsCent_hft;p_{T}(K#pi#pi)(GeV/c);Cent;m_{K#pi#pi - K#pi}(GeV/c^{2})", 120, 0, 12, 10, -1.5, 8.5, 50, 0.138, 0.158);
   mh3DstarInvMassPtdPhi_hft = new TH3F("hDstarInvMassPtdPhi_hft","DstarVsInvMassVsPtVsdPhi_hft;p_{T}(K#pi#pi)(GeV/c);dPhi;m_{K#pi#pi-K#pi}(GeV/c^{2})", 120, 0, 12,10,0,M_PI,50,0.138,0.158 );
   for(int ii = 0 ; ii<9; ++ii){
     //All TPC tracks
     mh3DstarDPhiInvMassPt_cent[ii] = new TH3F(Form("h3DstarDPhiInvMassPt_cent%i",ii), Form("dPhiVsInvMassVsPt_cent%i; d#phi ; m_{K#pi#pi - K#pi} (GeV/c^2); p_{T}(K#pi#pi) (GeV/c)",ii),120, 0.,12., 50, 0.138, 0.158, 10, 0., M_PI );
     mh3DstarCos2DPhiInvMassPt_cent[ii] = new TH3F(Form("h3DstarCos2DPhiInvMassPt_cent%i",ii), Form("Cos2dPhiVsInvMassVsPt_cent%i; cos(2d#phi); m_{K#pi#pi - K#pi} (GeV/c^2); p_{T}(K#pi#pi) (GeV/c)",ii),120, 0,12, 50, 0.138, 0.158, 200,-1.0,1.0 );
     mh2DstarInvMassPt_cent[ii] = new TH2F(Form("h2DstarInvMassPt_cent%i", ii), Form("InvMassVsPt_cent%i;  m_{K#pi#pi - K#pi} (GeV/c^2); p_{T}(K#pi#pi) (GeV/c)", ii), 120, 0, 12, 50, 0.138, 0.158);
   }
   //test
   nt = new TNtuple("nt","nt","flag:kaonPID:pionPID:softPionPID:"
		    "dcaKaon:dcaPion:dcaSoftPion:dcaDaugters_D0:kaonPt:pionPt:softPionPt:"
		    "theta_D0:decayL_D0:pt_D0:mass_D0:eta_D0:phi_D0:"
		    "pt_Dstar:mass_Dstar:eta_Dstar:phi_Dstar");
}
StPicoDstarAnaHists::~StPicoDstarAnaHists()
{
   delete mPrescales;
   // note that histograms are owned by mOutFile. They will be destructed
   // when the file is closed.
}
//-----------------------------------------------------------------------
void StPicoDstarAnaHists::addEvent(StPicoEvent const* const picoEvent)
{
   int runIndex = mPrescales->runIndex(picoEvent->runId());
   mh1TotalEventsInRun->Fill(runIndex);
}
//-----------------------------------------------------------------------
void StPicoDstarAnaHists::addEventBeforeCut(StPicoEvent const* const picoEvent)
{
   int runIndex = mPrescales->runIndex(picoEvent->runId());
   mh1TotalEventsInRunBeforeCut->Fill(runIndex);
}
//-----------------------------------------------------------------------
void StPicoDstarAnaHists::addCent(const double refmultCor, int centrality, const double reweight, const float vz)
{
   mh1gRefmultCor->Fill(refmultCor);
   mh1gRefmultCorWg->Fill(refmultCor, reweight);
   mh1Cent->Fill(centrality);
   mh1CentWg->Fill(centrality, reweight);
   mh2CentVz->Fill(centrality, vz);
   mh2CentVzWg->Fill(centrality, vz, reweight);
}
//-----------------------------------------------------------------------
void StPicoDstarAnaHists::addTpcDenom1(bool IsPion, bool IsKaon, float pt, int centrality, float Eta, float Phi, float Vz, float ZdcX)
{
   int EtaIndex = getEtaIndex(Eta);
   int PhiIndex = getPhiIndex(Phi);
   int VzIndex = getVzIndex(Vz);
   int ZdcxIndex = getZdcxIndex(ZdcX);
   if (IsPion)
   {
      mh2Tpc1PtCentPartEtaVz[0][EtaIndex][VzIndex]->Fill(pt, centrality);
      mh2Tpc1PtCentPartPhi[0][PhiIndex]->Fill(pt, centrality);
      mh2Tpc1PtCentPartEtaVzPhi[0][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
      mh2Tpc1PtCentPartZdcx[0][ZdcxIndex]->Fill(pt, centrality);
   }
   if (IsKaon)
   {
      mh2Tpc1PtCentPartEtaVz[1][EtaIndex][VzIndex]->Fill(pt, centrality);
      mh2Tpc1PtCentPartPhi[1][PhiIndex]->Fill(pt, centrality);
      mh2Tpc1PtCentPartEtaVzPhi[1][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
      mh2Tpc1PtCentPartZdcx[1][ZdcxIndex]->Fill(pt, centrality);
   }
   mh2Tpc1PtCent->Fill(pt, centrality);
   if (fabs(Eta) < 0.1 && pt > 3.0) mh2Tpc1PhiVz->Fill(Phi, Vz);

}
//-----------------------------------------------------------------------
void StPicoDstarAnaHists::addHFTNumer1(bool IsPion, bool IsKaon, float pt, int centrality, float Eta, float Phi, float Vz, float ZdcX)
{

   int EtaIndex = getEtaIndex(Eta);
   int PhiIndex = getPhiIndex(Phi);
   int VzIndex = getVzIndex(Vz);
   int ZdcxIndex = getZdcxIndex(ZdcX);
   if (IsPion)
   {
      mh2HFT1PtCentPartEtaVz[0][EtaIndex][VzIndex]->Fill(pt, centrality);
      mh2HFT1PtCentPartPhi[0][PhiIndex]->Fill(pt, centrality);
      mh2HFT1PtCentPartEtaVzPhi[0][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
      mh2HFT1PtCentPartZdcx[0][ZdcxIndex]->Fill(pt, centrality);
   }
   if (IsKaon)
   {
      mh2HFT1PtCentPartEtaVz[1][EtaIndex][VzIndex]->Fill(pt, centrality);
      mh2HFT1PtCentPartPhi[1][PhiIndex]->Fill(pt, centrality);
      mh2HFT1PtCentPartEtaVzPhi[1][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
      mh2HFT1PtCentPartZdcx[1][ZdcxIndex]->Fill(pt, centrality);
   }
   mh2HFT1PtCent->Fill(pt, centrality);
   if (fabs(Eta) < 0.1 && pt > 3.0) mh2HFT1PhiVz->Fill(Phi, Vz);
}
//-----------------------------------------------------------------------
void StPicoDstarAnaHists::addKaonPion(StKaonPion const* const kp, bool unlike, bool tpc, bool tof, int centrality, const double reweight)
{
   if (unlike)
   {
      if (tpc) mh2InvariantMassVsPt->Fill(kp->pt(), kp->m());
      if (tof) mh2InvariantMassVsPtTof->Fill(kp->pt(), kp->m());
      if (tpc) mh3InvariantMassVsPtVsCent->Fill(kp->pt(), centrality, kp->m(), reweight);
      if (tof) mh3InvariantMassVsPtVsCentTof->Fill(kp->pt(), centrality, kp->m(), reweight);
   }
   else
   {
      if (tpc) mh2InvariantMassVsPtLike->Fill(kp->pt(), kp->m());
      if (tof) mh2InvariantMassVsPtTofLike->Fill(kp->pt(), kp->m());
      if (tpc) mh3InvariantMassVsPtVsCentLike->Fill(kp->pt(), centrality, kp->m(), reweight);
      if (tof) mh3InvariantMassVsPtVsCentTofLike->Fill(kp->pt(), centrality, kp->m(), reweight);
   }
}
//---------------------------------------------------------------------
void StPicoDstarAnaHists::addDstar(const StLorentzVectorF &kpp, StKaonPion const & kp, float dPhi, int centrality, const double reweight, const bool isHft)
{
  mh3DstarInvariantMassVsPtVsCent->Fill(kpp.perp(), centrality, kpp.m()-kp.m(), reweight);
  mh3DstarInvMassPtdPhi->Fill(kpp.perp(), dPhi, kpp.m()-kp.m());
  if( isHft == true){
    mh3DstarInvariantMassVsPtVsCent_hft->Fill(kpp.perp(), centrality, kpp.m()-kp.m(), reweight);
    mh3DstarInvMassPtdPhi_hft->Fill(kpp.perp(), dPhi, kpp.m()-kp.m());
  }
  mh3DstarDPhiInvMassPt_cent[centrality]->Fill(kpp.perp(), kpp.m()-kp.m(), dPhi);;
  mh3DstarCos2DPhiInvMassPt_cent[centrality]->Fill(kpp.perp(), kpp.m()-kp.m(), cos(2.0*dPhi));
  mh2DstarInvMassPt_cent[centrality]->Fill(kpp.perp(), kpp.m()-kp.m());
}
//---------------------------------------------------------------------
void StPicoDstarAnaHists::addDzero(StKaonPion const & kp, float dPhi, int centrality, const double reweight)
{
  mh3D0InvariantMassVsPtVsCent->Fill(kp.pt(), centrality, kp.m(), reweight);
  mh3D0InvMassPtdPhi->Fill(kp.pt(), dPhi, kp.m());
}
//---------------------------------------------------------------------
void StPicoDstarAnaHists::addDstarNtuple(const StLorentzVectorF & kpp, StPicoTrack const & pion, StKaonPion const & kp)
{
  float ntVar[30];
  int flag=0;
  int ii = 0;
  if(pion.charge()>0) // Dplus
    flag = 1;
  ntVar[ii++] = flag;
  //PID
  //ntVar[ii++] = fabs(kaon->nSigmaKaon());
  //ntVar[ii++] = fabs(pion1->nSigmaPion());
  ntVar[ii++] = fabs(pion.nSigmaPion());
  //Track info
  //ntVar[ii++] = kp->kaonDca();
  //ntVar[ii++] = kp->pionDca();
  //ntVar[ii++] = (helix.origin()-primVtx).mag();
  //ntVar[ii++] = kp->dcaDaughters();
  //ntVar[ii++] = kaon->gPt();
  //ntVar[ii++] = pion1->gPt();
  ntVar[ii++] = pion.pMom().perp();
  //D0 information
  ntVar[ii++] = kp.pointingAngle();
  ntVar[ii++] = kp.decayLength();
  ntVar[ii++] = kp.pt();
  ntVar[ii++] = kp.m();
  ntVar[ii++] = kp.eta();
  ntVar[ii++] = kp.phi();
  //Dstar info
  ntVar[ii++] = kpp.perp();
  ntVar[ii++] = kpp.m();
  ntVar[ii++] = kpp.pseudoRapidity();
  ntVar[ii++] = kpp.phi();
  //Fill
  nt->Fill(ntVar);
}
//---------------------------------------------------------------------
void StPicoDstarAnaHists::addDcaPtCent(float dca, float dcaXy, float dcaZ, bool IsPion, bool IsKaon, float pt,  int centrality, float Eta, float Phi, float Vz, float ZdcX)
{
   int EtaIndex = getEtaIndex(Eta);
   int PhiIndex = getPhiIndex(Phi);
   int VzIndex = getVzIndex(Vz);
   int ZdcxIndex = getZdcxIndex(ZdcX);
   if (IsPion)
   {
      mh3DcaXyPtCentPartEtaVz[0][EtaIndex][VzIndex]->Fill(pt, centrality, dcaXy);
      mh3DcaXyPtCentPartPhi[0][PhiIndex]->Fill(pt, centrality, dcaXy);
      mh3DcaXyPtCentPartZdcx[0][ZdcxIndex]->Fill(pt, centrality, dcaXy);
      mh3DcaZPtCentPartEtaVz[0][EtaIndex][VzIndex]->Fill(pt, centrality, dcaZ);
      mh3DcaZPtCentPartPhi[0][PhiIndex]->Fill(pt, centrality, dcaZ);
      mh3DcaZPtCentPartZdcx[0][ZdcxIndex]->Fill(pt, centrality, dcaZ);
   }
   if (IsKaon)
   {
      mh3DcaXyPtCentPartEtaVz[1][EtaIndex][VzIndex]->Fill(pt, centrality, dcaXy);
      mh3DcaXyPtCentPartPhi[1][PhiIndex]->Fill(pt, centrality, dcaXy);
      mh3DcaXyPtCentPartZdcx[1][ZdcxIndex]->Fill(pt, centrality, dcaXy);
      mh3DcaZPtCentPartEtaVz[1][EtaIndex][VzIndex]->Fill(pt, centrality, dcaZ);
      mh3DcaZPtCentPartPhi[1][PhiIndex]->Fill(pt, centrality, dcaZ);
      mh3DcaZPtCentPartZdcx[1][ZdcxIndex]->Fill(pt, centrality, dcaZ);
   }
   mh3DcaXyPtCent->Fill(pt, centrality, dcaXy);
   mh3DcaZPtCent->Fill(pt, centrality, dcaZ);
}
//---------------------------------------------------------------------
int StPicoDstarAnaHists::getEtaIndex(float Eta)
{
   for (int i = 0; i < anaCuts::nEtas; i++)
   {
      if ((Eta >= anaCuts::EtaEdge[i]) && (Eta < anaCuts::EtaEdge[i + 1]))
         return i;
   }
   return anaCuts::nEtas - 1;
}
//---------------------------------------------------------------------
int StPicoDstarAnaHists::getPhiIndex(float Phi)
{
   for (int i = 0; i < anaCuts::nPhis; i++)
   {
      if ((Phi >= anaCuts::PhiEdge[i]) && (Phi < anaCuts::PhiEdge[i + 1]))
         return i;
   }
   return anaCuts::nPhis - 1;
}
//---------------------------------------------------------------------
int StPicoDstarAnaHists::getVzIndex(float Vz)
{
   for (int i = 0; i < anaCuts::nVzs; i++)
   {
      if ((Vz >= anaCuts::VzEdge[i]) && (Vz < anaCuts::VzEdge[i + 1]))
         return i;
   }
   return anaCuts::nVzs - 1;
}
//---------------------------------------------------------------------
int StPicoDstarAnaHists::getZdcxIndex(float ZdcX)
{
   for (int i = 0; i < anaCuts::nZdcxs; i++)
   {
      if ((ZdcX >= anaCuts::ZdcxEdge[i]) && (ZdcX < anaCuts::ZdcxEdge[i + 1]))
         return i;
   }
   return anaCuts::nZdcxs - 1;
}
//---------------------------------------------------------------------
void StPicoDstarAnaHists::addQaNtuple(int runnumber, float dca, float vz, float pt, float eta, float phi, int centrality, const double refmultCor, float zdcx, int tofMatchFlag, int hftMatchFlag)
{
//  nt->Fill(runnumber, dca, vz, pt, eta, phi, centrality, refmultCor, zdcx, tofMatchFlag, hftMatchFlag);
}
//---------------------------------------------------------------------
void StPicoDstarAnaHists::closeFile()
{
   mOutFile->cd();

   mh1TotalEventsInRun->Write();
   mh1TotalEventsInRunBeforeCut->Write();
   mh2InvariantMassVsPt->Write();
   mh2InvariantMassVsPtLike->Write();
   mh2InvariantMassVsPtTof->Write();
   mh2InvariantMassVsPtTofLike->Write();
   //centrality
   mh1Cent->Write();
   mh1CentWg->Write();
   mh1gRefmultCor->Write();
   mh1gRefmultCorWg->Write();
   mh2CentVz->Write();
   mh2CentVzWg->Write();
   mh3InvariantMassVsPtVsCent->Write();
   mh3InvariantMassVsPtVsCentLike->Write();
   mh3InvariantMassVsPtVsCentTof->Write();
   mh3InvariantMassVsPtVsCentTofLike->Write();
   //HFT ratio QA
   //mh2Tpc1PtCent->Write();
   // mh2Tpc1PhiVz->Write();
   //mh2HFT1PhiVz->Write();
   //mh2HFT1PtCent->Write();

   /*HFT DCA Ratio
   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < anaCuts::nEtas; iEta++)
      {
         for (int iVz = 0; iVz < anaCuts::nVzs; iVz++)
         {
            mh3DcaXyPtCentPartEtaVz[iParticle][iEta][iVz]->Write();
            mh3DcaZPtCentPartEtaVz[iParticle][iEta][iVz]->Write();
            mh2Tpc1PtCentPartEtaVz[iParticle][iEta][iVz]->Write();
            mh2HFT1PtCentPartEtaVz[iParticle][iEta][iVz]->Write();
         }
      }
   }
   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < anaCuts::nEtas; iEta++)
      {
         for (int iVz = 0; iVz < anaCuts::nVzs; iVz++)
         {
            for (int iPhi = 0; iPhi < anaCuts::nPhis; iPhi++)
            {
               mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]->Write();
               mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]->Write();
            }
         }
      }
   }
   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iPhi = 0; iPhi < anaCuts::nPhis; iPhi++)
      {
         mh3DcaXyPtCentPartPhi[iParticle][iPhi]->Write();
         mh3DcaZPtCentPartPhi[iParticle][iPhi]->Write();
         mh2Tpc1PtCentPartPhi[iParticle][iPhi]->Write();
         mh2HFT1PtCentPartPhi[iParticle][iPhi]->Write();
      }
   }
   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iZdcx = 0; iZdcx < anaCuts::nZdcxs; iZdcx++)
      {
         mh3DcaXyPtCentPartZdcx[iParticle][iZdcx]->Write();
         mh3DcaZPtCentPartZdcx[iParticle][iZdcx]->Write();
         mh2Tpc1PtCentPartZdcx[iParticle][iZdcx]->Write();
         mh2HFT1PtCentPartZdcx[iParticle][iZdcx]->Write();
      }
      }*/
   mh3DcaXyPtCent->Write();
   mh3DcaZPtCent->Write();
   //Dstar
   mh2EPresCentRandom->Write();
   mh2EPresCentEta->Write();
   mh2EventPlaneCent->Write();
   mh3D0InvariantMassVsPtVsCent->Write();
   mh3D0InvMassPtdPhi->Write();
   mh3DstarInvariantMassVsPtVsCent->Write();
   mh3DstarInvMassPtdPhi->Write();
   mh3DstarInvariantMassVsPtVsCent_hft->Write();
   mh3DstarInvMassPtdPhi_hft->Write();
   for(int ii=0; ii<9; ++ii){
     mh3DstarDPhiInvMassPt_cent[ii]->Write();
     mh3DstarCos2DPhiInvMassPt_cent[ii]->Write();
     mh2DstarInvMassPt_cent[ii]->Write();
   }
   //Test
   //nt->Write();
   
   mOutFile->Close();
   mOutFile->Delete();
}
