#ifndef StPicoDstarAnaHists__h
#define StPicoDstarAnaHists__h

/* **************************************************
 *  A class to create and save my D0 analysis histograms.
 *
 *  Authors: Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */

#include "StAnaCuts.h"
#include "THnSparse.h" 
#include "THn.h" 
#include "TH2.h"
#include "StLorentzVectorF.hh"

class TH1F;
class TH2F;
class TH3F;
class TFile;
class TString;
class StPicoPrescales;
class StPicoEvent;
class StPicoTrack;
class StKaonPion;
class TNtuple;


class StPicoDstarAnaHists
{
  public:
   StPicoDstarAnaHists(TString fileBaseName);
   virtual ~StPicoDstarAnaHists();
   void addEvent(StPicoEvent const *);
   void addEventBeforeCut(StPicoEvent const *);
   void addCent(const double refmultCor,int centrality, const double reweight, const float vz);
   void addKaonPion(StKaonPion const*, bool unlike, bool tpc, bool tof, int centrality, const double reweight);
   void addTpcDenom1(bool IsPion, bool IsKaon, float pt, int centrality, float Eta, float Phi, float Vz, float ZdcX);
   void addHFTNumer1(bool IsPion, bool IsKaon, float pt, int centrality, float Eta, float Phi, float Vz, float ZdcX);
   void addQaNtuple(int, float, float, float, float, float, int, const double, float, int, int);
   void addDcaPtCent(float dca, float dcaXy, float  dcaZ, bool IsPion, bool IsKaon, float pt,  int centrality, float Eta, float Phi, float Vz, float ZdcX);
   //
   void addEventPlane(int const centrality,  float const eventPlane, float const resolution_random, float const resolution_eta); 
   void addDstar(const StLorentzVectorF &kpp, StKaonPion const & kp,float dPhi, int centrality, const double reweight, const bool isHft);
   void addDzero(StKaonPion const & kp,float dPhi, int centrality, const double reweight);
   //test
   void addDstarNtuple(const StLorentzVectorF &kpp, StPicoTrack const & pion, StKaonPion const & kp);

   int getEtaIndex(float Eta) ;
   int getPhiIndex(float Phi) ;
   int getVzIndex(float Vz) ;
   int getZdcxIndex(float ZdcX) ;
   void closeFile();

  private:
   StPicoDstarAnaHists(){}

   StPicoPrescales* mPrescales;
   TFile* mOutFile;
   TH1F* mh1TotalEventsInRun;
   TH1F* mh1TotalEventsInRunBeforeCut;
   TH2F* mh2InvariantMassVsPt;
   TH2F* mh2InvariantMassVsPtLike;
   TH2F* mh2InvariantMassVsPtTof;
   TH2F* mh2InvariantMassVsPtTofLike;
   //centrality
   TH1F* mh1Cent;
   TH1F* mh1CentWg;
   TH1F* mh1gRefmultCor;
   TH1F* mh1gRefmultCorWg;
   TH2F* mh2CentVz;
   TH2F* mh2CentVzWg;
   TH3F* mh3InvariantMassVsPtVsCent;
   TH3F* mh3InvariantMassVsPtVsCentLike;
   TH3F* mh3InvariantMassVsPtVsCentTof;
   TH3F* mh3InvariantMassVsPtVsCentTofLike;
   //HFT ratio QA
   TH2F* mh2Tpc1PtCent;
   TH2F* mh2Tpc1PhiVz;
   TH2F* mh2HFT1PtCent;
   TH2F* mh2HFT1PhiVz;
   TH2F* mh2Tpc1PtCentPartEtaVz[anaCuts::nParticles][anaCuts::nEtas][anaCuts::nVzs];
   TH2F* mh2Tpc1PtCentPartEtaVzPhi[anaCuts::nParticles][anaCuts::nEtas][anaCuts::nVzs][anaCuts::nPhis];
   TH2F* mh2Tpc1PtCentPartPhi[anaCuts::nParticles][anaCuts::nPhis];
   TH2F* mh2Tpc1PtCentPartZdcx[anaCuts::nParticles][anaCuts::nZdcxs];
   TH2F* mh2HFT1PtCentPartEtaVz[anaCuts::nParticles][anaCuts::nEtas][anaCuts::nVzs];
   TH2F* mh2HFT1PtCentPartEtaVzPhi[anaCuts::nParticles][anaCuts::nEtas][anaCuts::nVzs][anaCuts::nPhis];
   TH2F* mh2HFT1PtCentPartPhi[anaCuts::nParticles][anaCuts::nPhis];
   TH2F* mh2HFT1PtCentPartZdcx[anaCuts::nParticles][anaCuts::nZdcxs];
   //HFT Dca 
   TH3F* mh3DcaXyPtCentPartEtaVz[anaCuts::nParticles][anaCuts::nEtas][anaCuts::nVzs];
   TH3F* mh3DcaXyPtCentPartPhi[anaCuts::nParticles][anaCuts::nPhis];
   TH3F* mh3DcaXyPtCentPartZdcx[anaCuts::nParticles][anaCuts::nZdcxs];
   TH3F* mh3DcaZPtCentPartEtaVz[anaCuts::nParticles][anaCuts::nEtas][anaCuts::nVzs];
   TH3F* mh3DcaZPtCentPartPhi[anaCuts::nParticles][anaCuts::nPhis];
   TH3F* mh3DcaZPtCentPartZdcx[anaCuts::nParticles][anaCuts::nZdcxs];
   TH3F* mh3DcaXyPtCent;
   TH3F* mh3DcaZPtCent;
   //Analysis plots
   TH2F* mh2EPresCentRandom;
   TH2F* mh2EPresCentEta;
   TH2F* mh2EventPlaneCent;
   TH1F* mh1EventPlaneEta;
   //all
   TH3F* mh3D0InvariantMassVsPtVsCent;
   TH3F* mh3D0InvMassPtdPhi;
   TH3F* mh3DstarInvariantMassVsPtVsCent;
   TH3F* mh3DstarInvMassPtdPhi;
   //with hft-test
   TH3F* mh3DstarInvariantMassVsPtVsCent_hft;
   TH3F* mh3DstarInvMassPtdPhi_hft;
   //centrality bins  
   TH3F* mh3DstarDPhiInvMassPt_cent[9];
   TH3F* mh3DstarCos2DPhiInvMassPt_cent[9];
   TH2F* mh2DstarInvMassPt_cent[9];
   //Tests
   TNtuple* nt;
};
inline void StPicoDstarAnaHists::addEventPlane(int const centrality,  float const eventPlane, float const resolution_random, float const resolution_eta)
{
  mh2EventPlaneCent->Fill(centrality, eventPlane);
  mh2EPresCentRandom->Fill(centrality, resolution_random);
  mh2EPresCentEta->Fill(centrality, resolution_eta);
}
#endif
