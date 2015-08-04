#ifndef StPicoDstarAnaMaker_h
#define StPicoDstarAnaMaker_h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoD0Event
 *  simultaneously and do analysis. 
 *
 *  Please write your analysis in the ::Make() function.
 *
 *  Authors:  Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */

#include "TChain.h"
#include "StMaker.h"
#include "StThreeVectorF.hh"
#include "StLorentzVectorF.hh"

class TString;
class TFile;
class TNtuple;
class StPicoEvent;
class StPicoD0Event;
class StKaonPion;
class StPicoTrack;
class StPicoDstMaker;
class StPicoDstarAnaHists;
class StRefMultCorr;
class StEventPlane;
class StPicoDstarAnaMaker : public StMaker
{
  public:
    StPicoDstarAnaMaker(char const * name, TString const inputFilesList, 
			TString const outBaseName,StPicoDstMaker* picoDstMaker, 
			StRefMultCorr* grefmultCorrUtil, StEventPlane *eventPlaneMaker);
    virtual ~StPicoDstarAnaMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    int getEntries() const;

  private:
    StPicoDstarAnaMaker() {}
    void readNextEvent();

    int  getD0PtIndex(StKaonPion const* ) const;
    bool isGoodEvent(StPicoEvent const*) const;
    bool isGoodQaTrack(StPicoTrack const * ,StThreeVectorF const momentum ,const double dca) const;
    bool isGoodTrack(StPicoTrack const*) const;
    bool isTpcPion(StPicoTrack const*) const;
    bool isTpcKaon(StPicoTrack const*) const;
    bool isTofPion(StPicoTrack const* const, float beta) const;
    bool isTofKaon(StPicoTrack const* const, float beta) const;
    bool isGoodPair(StKaonPion const*) const;
    bool isGoodDstar(StPicoTrack const*, StKaonPion const*, StThreeVectorF const &) const;
    float getTofBeta(StPicoTrack const*,StThreeVectorF const * pVtx) const;
    bool isSoftPion(StPicoTrack const*, StThreeVectorF const &) const;
    bool loadEventPlaneCorr(StEventPlane const * mEventPlane);

    StPicoDstMaker* mPicoDstMaker;
    StPicoD0Event* mPicoD0Event;
    StRefMultCorr* mGRefMultCorrUtil;
    StEventPlane* mEventPlane;

    TString mInputFilesList;
    TString mOutFileBaseName;
    TChain* mChain;
    int mEventCounter;

    // -------------- USER variables -------------------------
    // add your member variables here. 
    // Remember that ntuples size can be really big, use histograms where appropriate
    StPicoDstarAnaHists* mHists;

    ClassDef(StPicoDstarAnaMaker, 1)
};

inline int StPicoDstarAnaMaker::getEntries() const 
{
  return mChain? mChain->GetEntries() : 0;
}

inline void StPicoDstarAnaMaker::readNextEvent()
{
  mChain->GetEntry(mEventCounter++);
}
#endif
