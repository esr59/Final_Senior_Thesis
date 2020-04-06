// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StEshaMaker.h"
#include "StMemStat.h"

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include <THnSparse.h>
#include "TParameter.h"
#include <TProfile.h>
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"

// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoMtdTrigger.h"
#include "StRoot/StPicoEvent/StPicoBEmcPidTraits.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoMtdPidTraits.h"

// jet-framework includes
#include "StJetFrameworkPicoBase.h"
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StFemtoTrack.h"
#include "StEmcPosition2.h"
#include "StCentMaker.h"

// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StEshaMaker)

//________________________________________________________________________
StEshaMaker::StEshaMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", const char* jetMakerName = "", const char* rhoMakerName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
{
  fLeadingJet = 0x0; fSubLeadingJet = 0x0;
  fJets = 0x0 ;
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  JetMaker = 0;
  RhoMaker = 0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StEshaMaker::fRunFlagEnum
  doppAnalysis = kFALSE;
  fRequireCentSelection = kFALSE;
  fCentralitySelectionCut = -99;
  fDoEffCorr = kFALSE;
  doRejectBadRuns = kFALSE;
  fCorrJetPt = kFALSE;
  fMinPtJet = 0.0;
  fTrackBias = 0.0;
  fTowerBias = 0.0;
  fJetRad = 0.4;
  fEventZVtxMinCut = -40.0; fEventZVtxMaxCut = 40.0;
  fMaxEventTrackPt = 30.0;
  fMaxEventTowerEt = 1000.0; // 30.0
  fTrackPtMinCut = 0.2; fTrackPtMaxCut = 30.0;
  fTrackPhiMinCut = 0.0; fTrackPhiMaxCut = 2.0*TMath::Pi();
  fTrackEtaMinCut = -1.0; fTrackEtaMaxCut = 1.0;
  fTrackDCAcut = 3.0;
  fTracknHitsFit = 15; fTracknHitsRatio = 0.52;
  fTowerEMinCut = 0.2; fTowerEMaxCut = 100.0;
  fTowerEtaMinCut = -1.0; fTowerEtaMaxCut = 1.0;
  fTowerPhiMinCut = 0.0;  fTowerPhiMaxCut = 2.0*TMath::Pi();
  fCentralityScaled = 0.; ref16 = -99; ref9 = -99; // FIXME - maybe not make global
  Bfield = 0.0;
//  mVertex = 0x0;
  zVtx = 0.0;
  fEmcTriggerEventType = 0; fMBEventType = 2;
  fRho = 0x0;
  fRhoVal = 0;
  mEmcPosition = 0x0;
  mCentMaker = 0x0;
  mBaseMaker = 0x0;
  fAnalysisMakerName = name;
  fJetMakerName = jetMakerName;
  fRhoMakerName = rhoMakerName;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }

}

//
//________________________________________________________________________
StEshaMaker::~StEshaMaker()
{ /*  */
  // destructor
  if(hCentrality)  delete hCentrality;
  if(hMultiplicity)delete hMultiplicity;
  if(hJetPt)       delete hJetPt;
  if(hJetCorrPt)   delete hJetCorrPt;
  if(hJetMass)	   delete hJetMass;
  if(hTrackPt)     delete hTrackPt;
  if(hJetShape)    delete hJetShape;
  if(hJetShape1)   delete hJetShape1;
  if(hJetShape2)   delete hJetShape2;
  if(hJetShape3)   delete hJetShape3;
  if(hJetPt1)      delete hJetPt1;
  if(hJetPt2)      delete hJetPt2;
  if(hJetPt3)      delete hJetPt3;

  if(mEmcPosition) delete mEmcPosition;
}

//
//________________________________________________________________________
Int_t StEshaMaker::Init() {
  StJetFrameworkPicoBase::Init();

  // declare histograms
  DeclareHistograms();

  // position object for Emc
  mEmcPosition = new StEmcPosition2();

  // Jet TClonesArray
  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  //fJets->SetName(fJetsName);

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StEshaMaker::Finish() { 
  cout << "StEshaMaker::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StEshaMaker::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StEshaMaker::DeclareHistograms() {
  // binning for cent histograms
  int nHistCentBins = 20;

  // binning for mult histograms
  double kHistMultMax = 800.;
  int kHistMultBins = 400;

  // pp specific settings
  if(doppAnalysis) {
    kHistMultMax = 100.;
    kHistMultBins = 100.;
  }

  // histograms
  hCentrality = new TH1F("hCentrality", "No. events vs centrality", nHistCentBins, 0, 100);
  hMultiplicity = new TH1F("hMultiplicity", "No. events vs multiplicity", kHistMultBins, 0, kHistMultMax);

  // jet QA histos
  hJetPt = new TH1F("hJetPt", "Jet p_{T}", 100, 0, 100);
  hJetCorrPt = new TH1F("hJetCorrPt", "Corrected Jet p_{T}", 125, -25, 100);
  hJetMass = new TH1F("hJetMass", "Jet Mass", 20, 0, 20);

  // track pt histo
  hTrackPt = new TH1F("hTrackPt", "Track Pt", 50, 0, 50);

  // jet shape with 0.05 bins
  hJetShape = new TH1F("hJetShape", "Jet Shape", 8, 0, 0.4);

  //jet pt for different jet mass bins
  hJetPt1 = new TH1F("hJetPt1", "Jet p_{T} for < 2 GeV", 100, 0, 100);
  hJetPt2 = new TH1F("hJetPt2", "Jet p_{T} between 2 GeV and 4 GeV", 100, 0, 100);
  hJetPt3 = new TH1F("hJetPt3", "Jet p_{T} between 4 GeV and 6 GeV", 100, 0, 100);

  //jet shape for different jet mass bins
  hJetShape1 = new TH1F("hJetShape1", "Jet Shape for < 2 GeV", 8, 0, 0.4);
  hJetShape2 = new TH1F("hJetShape2", "Jet Shape between 2 GeV and 4 GeV", 8, 0, 0.4);
  hJetShape3 = new TH1F("hJetShape3", "Jet Shape between 4 GeV and 6 GeV", 8, 0, 0.4);

//trying to get the error bars
SetSumw2();

}
//
// write histograms
//_____________________________________________________________________________
void StEshaMaker::WriteHistograms() {
  // writing of histograms done here
  hCentrality->Write();
  hMultiplicity->Write();
  hJetPt->Write();
  hJetCorrPt->Write();
  hJetMass->Write();
  hTrackPt->Write();
  hJetShape->Write();
  hJetShape1->Write();
  hJetShape2->Write();
  hJetShape3->Write();
  hJetPt1->Write();
  hJetPt2->Write();
  hJetPt3->Write();
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StEshaMaker::Clear(Option_t *opt) {
  fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StEshaMaker::Make() {
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  //jets per event counter
  int nJetpE = 0;

  //initialize variables
  double _nJets = 0.0;
  double _totJets = 0.0;
  int _twojets = 0.0;
  int _fourjets = 0.0;
  int _sixjets = 0.0;

  //more global variables for jetshape
  double tracksum[8] = {0};
  double jetsum[8] = {0};

  // get PicoDstMaker 
  mPicoDstMaker = static_cast<StPicoDstMaker*>(GetMaker("picoDst"));
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  // construct PicoDst object from maker
  mPicoDst = static_cast<StPicoDst*>(mPicoDstMaker->picoDst());
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  // create pointer to PicoEvent 
  mPicoEvent = static_cast<StPicoEvent*>(mPicoDst->event());
  if(!mPicoEvent) {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }

  // get base class pointer
  mBaseMaker = static_cast<StJetFrameworkPicoBase*>(GetMaker("baseClassMaker"));
  if(!mBaseMaker) {
    LOG_WARN << " No baseMaker! Skip! " << endm;
    return kStWarn;
  }

  // get bad run, dead & bad tower lists
  badRuns = mBaseMaker->GetBadRuns();
  deadTowers = mBaseMaker->GetDeadTowers();
  badTowers = mBaseMaker->GetBadTowers();

  // get run number, check bad runs list if desired (kFALSE if bad)
  fRunNumber = mPicoEvent->runId();
  if(doRejectBadRuns) {
    if( !mBaseMaker->IsRunOK(fRunNumber) ) return kStOK;
  }

  //for(int i = 1; i<4801; i++) {  if(!mBaseMaker->IsTowerOK(i))  cout<<"tower: "<<i<<" is not good!!"<<endl;  }
  // can check for bad towers like this (if needed):
  // bool isTowOK = mBaseMaker->IsTowerOK(towerID);
  // ===========================================================================================

  // cut event on max track pt > 30.0 GeV
  if(GetMaxTrackPt() > fMaxEventTrackPt) return kStOK;

  // cut event on max tower Et > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField(); 

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();
  
  // Z-vertex cut: the Aj analysis cut on (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;


  // ============================ CENTRALITY ============================== //
  // get CentMaker pointer
  mCentMaker = static_cast<StCentMaker*>(GetMaker("CentMaker"));
  if(!mCentMaker) {
    LOG_WARN << " No CenttMaker! Skip! " << endm;
    return kStWarn;
  }

  // centrality variables
  int grefMult = mCentMaker->GetgrefMult(); // see StPicoEvent
  int refMult =  mCentMaker->GetrefMult();  // see StPicoEvent
  ref9 = mCentMaker->GetRef9();   // binning from central -> peripheral
  ref16 = mCentMaker->GetRef16(); // binning from central -> peripheral
  int cent16 = mCentMaker->GetCent16(); // centrality bin from StRefMultCorr (increasing bin corresponds to decreasing cent %) - Don't use except for cut below
  int centbin = mCentMaker->GetRef16();
  double refCorr2 = mCentMaker->GetRefCorr2();
  fCentralityScaled = mCentMaker->GetCentScaled();
  //double refCorr = mCentMaker->GetCorrectedMultiplicity(refMult, zVtx, zdcCoincidenceRate, 0); // example usage
  // for pp analyses:    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 

  // fill histograms
  hCentrality->Fill(fCentralityScaled);
  hMultiplicity->Fill(refCorr2);

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }

  // ============================ end of CENTRALITY ============================== //

  // ========================= Trigger Info =============================== //
  // looking at the EMCal triggers - used for QA and deciding on HT triggers
  FillEmcTriggers();

  // get trigger IDs from PicoEvent class and loop over them
  vector<unsigned int> mytriggers = mPicoEvent->triggerIds();
  if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<"EventTriggers: ";
  for(unsigned int i=0; i<mytriggers.size(); i++) {
    if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<"i = "<<i<<": "<<mytriggers[i] << ", ";
  }
  if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<endl;

  // check for MB/HT event
  bool fHaveMBevent = CheckForMB(fRunFlag, fMBEventType);
  bool fHaveMB5event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB5);
  bool fHaveMB30event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30);
  bool fHaveEmcTrigger = CheckForHT(fRunFlag, fEmcTriggerEventType);
  // ======================== end of Triggers ============================= //

  // =========================== JetMaker =============================== //
  // get JetMaker pointer
  JetMaker = static_cast<StJetMakerTask*>(GetMaker(fJetMakerName));
  const char *fJetMakerNameCh = fJetMakerName;
  if(!JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
    return kStWarn;
  }

  // get jet collection associated with JetMaker
  fJets = static_cast<TClonesArray*>(JetMaker->GetJets());
  if(!fJets) {
    LOG_WARN << Form(" No fJets object! Skip! ") << endm;
    return kStWarn;
  }

  // ============================= RhoMaker ============================== //
  // get RhoMaker pointer from event: old names "StRho_JetsBG"
  RhoMaker = static_cast<StRho*>(GetMaker(fRhoMakerName));
  const char *fRhoMakerNameCh = fRhoMakerName;
  if(!RhoMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fRhoMakerNameCh) << endm;
    return kStWarn;
  }

  // set rho object, alt fRho = GetRhoFromEvent(fRhoName);
  fRho = static_cast<StRhoParameter*>(RhoMaker->GetRho());
  if(!fRho) {
    LOG_WARN << Form("Couldn't get fRho object! ") << endm;
    return kStWarn;    
  } 
  
  // get rho/area value from rho object     fRho->ls("");
  fRhoVal = fRho->GetVal();
  // =======================================================================

  // get number of jets, tracks, and global tracks in events
  Int_t njets = fJets->GetEntries();
  const Int_t ntracks = mPicoDst->numberOfTracks();
  Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();

  // run Jets:
  RunJets();

  // run Tracks:
  RunTracks();

  // run Towers:
  RunTowers();

  return kStOK;
}

//
//
//_____________________________________________________________________________________________
void StEshaMaker::RunJets()
{
  // cache the leading + subleading jets within acceptance
  // first parameter is Jet Maker name, 2nd is Rho Parameter: fRho
  if(fCorrJetPt) {
    fLeadingJet = GetLeadingJet(fJetMakerName, fRho);
    fSubLeadingJet = GetSubLeadingJet(fJetMakerName, fRho);
  } else {
    fLeadingJet = GetLeadingJet(fJetMakerName);
    fSubLeadingJet = GetSubLeadingJet(fJetMakerName);
  }

  // ====================== Jet loop below ============================
  // loop over Jets in the event: initialize some parameter variables
  Int_t ijethi = -1;
  Double_t highestjetpt = 0.0;
  Int_t njets = fJets->GetEntries();

  // loop over jets
  for(int ijet = 0; ijet < njets; ijet++) {  // JET LOOP
    // get jet pointer
    StJet *jet = static_cast<StJet*>(fJets->At(ijet));
    if(!jet) continue;

    // get some jet parameters
    double jetarea = jet->Area();
    double jetpt = jet->Pt();
    double corrjetpt = jet->Pt() - jetarea*fRhoVal;
    double jetE = jet->E();
    double jetEta = jet->Eta();
    double jetPhi = jet->Phi();
    double jetNEF = jet->NEF();
    double jetmass = jet->M();

    //initial jet pt variables for each bin
    double jetpt1, jetpt2, jetpt3;
    
    //fill jet pt based on jet mass so we have a counter for each one
    if(jetmass<2){
    jetpt1 = jet->Pt();
    hJetPt1->Fill(jetpt1);
    }
    else if(jetmass>2 && jetmass<4){
    jetpt2 = jet->Pt();
    hJetPt2->Fill(jetpt2);
    }
    else if(jetmass>4 && jetmass<6){
    jetpt3 = jet->Pt();
    hJetPt3->Fill(jetpt3);
    }

 
    // initial tracksum
    double tracksum[8] = {0};
   
    //initialize jetsum for each of the separate jet mass bins (will do 0-2 Gev, 2-4 GeV, and 4-6 GeV for now)
    double jetsum[8] = {0};
    double jetsum1[8] = {0};
    double jetsum2[8] = {0};

    // get nTracks and maxTrackPt
    double maxtrackpt = jet->GetMaxTrackPt();
    double NtrackConstit = jet->GetNumberOfTracks();

    // get highest Pt jet in event (leading jet)
    if(highestjetpt < jetpt){
      ijethi = ijet;
      highestjetpt = jetpt;
    }

    // fill some basic histos
    hJetCorrPt->Fill(corrjetpt);
    hJetMass->Fill(jetmass);

/*
    // TEST - when using constituent subtractor
    vector<fastjet::PseudoJet> fConstituents = jet->GetJetConstituents();
    for(UInt_t ic = 0; ic < fConstituents.size(); ++ic) {
      // get user defined index
      Int_t uid = fConstituents[ic].user_index();
      double cpt = fConstituents[ic].perp();
      double ceta = fConstituents[ic].eta();
      double cphi = fConstituents[ic].phi();
      cout<<"ic = "<<ic<<", uid = "<<uid<<", cpt = "<<cpt<<", ceta = "<<ceta<<", cphi = "<<cphi<<endl;
    }
*/

    // get jet constituents: loop over constituent tracks
    for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {
      int trackid = jet->TrackAt(itrk);      
      
      //initialize variables
      double _nJets = 0.0;
      double _totJets = 0.0;
      _totJets++;
      int _twojets = 0.0;
      _twojets++;
      int _fourjets = 0.0;
      _fourjets++;
      int _sixjets = 0.0;
      _sixjets++;

      // get jet track pointer
      StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
      if(!trk){ continue; }

      // get momentum vector
      TVector3 mTrkMom;
      if(doUsePrimTracks) {
        if(!(trk->isPrimary())) continue; // check if primary
        // get primary track vector
        mTrkMom = trk->pMom();
      } else {
        // get global track vector
        mTrkMom = trk->gMom(mVertex, Bfield);
      }

      // track variables
      double pt = mTrkMom.Perp();
      double phi = mTrkMom.Phi();
      double eta = mTrkMom.PseudoRapidity();
      double px = mTrkMom.x();
      double py = mTrkMom.y();
      double pz = mTrkMom.z();
      short charge = trk->charge();
      double pi = TMath::Pi();

      if(phi < 0.0) phi += 2.0*pi;
      if(phi > 2.0*pi) phi -= 2.0*pi;

      //calculation for non reflected dphi
      double dphi = TMath::Abs(phi - jetPhi); // all  between 0 and 2pi
      double pi1 = TMath::Pi(); // pi value

      if(dphi > 1.5*pi1)                     dphi = 2*pi1 - dphi;  // quadrant 4 correction
      if(dphi < 1.5*pi1 && dphi > pi1)        dphi = 2*pi1 - dphi; // quandrant 3 correction
           

      //delta R calculation
      double distanceSig = TMath::Sqrt(TMath::Power(dphi,2) + TMath::Power(jetEta-eta,2));

      if(distanceSig > 0.4)     {continue;}

      double r = distanceSig;

      //fill track pt histo
      hTrackPt->Fill(pt);

      //separate the track pt fraction of total jet pt based on delta R
      if(r<0.05){
	tracksum[0]+=pt;}
      else if(r<0.1){
	tracksum[1]+=pt;}
      else if(r<0.15){
   	tracksum[2]+=pt;}
      else if(r<0.2){
	tracksum[3]+=pt;}
      else if(r<0.25){
        tracksum[4]+=pt;}
      else if(r<0.3){
        tracksum[5]+=pt;}
      else if(r<0.35){
        tracksum[6]+=pt;}
      else if(r<0.4){
        tracksum[7]+=pt;}
    
    } // track constituents loop

    // loop over constituents towers
    for(int itow = 0; itow < jet->GetNumberOfClusters(); itow++) {
      int towerid = jet->ClusterAt(itow);

      // get jet tower pointer
      StPicoBTowHit *tow = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(towerid));
      if(!tow){ continue; }

      // tower ID shifted by +1 from array index
      int towID = itow + 1;
    
    } // tower constituents loop

    //jets per event counter
    int nJetpE = 0;
          
    //number of jets per event
    nJetpE++;

    //initializing variables
    int _twojets = 0.0;
    int _fourjets = 0.0;
    int _sixjets = 0.0;
    
    //add up all the sums of the track pt in a specific range of jet mass to put into the jet sum
    	if(jetmass<2){
	  jetsum[0]+=tracksum[0];
	  jetsum[1]+=tracksum[1];
	  jetsum[2]+=tracksum[2];		
	  jetsum[3]+=tracksum[3];
	  jetsum[4]+=tracksum[4];
	  jetsum[5]+=tracksum[5];
	  jetsum[6]+=tracksum[6];
	  jetsum[7]+=tracksum[7];
	  _twojets++;}
	else if(jetmass>2 && jetmass<4){
	  jetsum1[0]+=tracksum[0];
	  jetsum1[1]+=tracksum[1];
	  jetsum1[2]+=tracksum[2];		
	  jetsum1[3]+=tracksum[3];
	  jetsum1[4]+=tracksum[4];
	  jetsum1[5]+=tracksum[5];
	  jetsum1[6]+=tracksum[6];
	  jetsum1[7]+=tracksum[7];
	  _fourjets++;}
	else if(jetmass>4 && jetmass<6){
	  jetsum2[0]+=tracksum[0];
	  jetsum2[1]+=tracksum[1];
	  jetsum2[2]+=tracksum[2];		
	  jetsum2[3]+=tracksum[3];
	  jetsum2[4]+=tracksum[4];
	  jetsum2[5]+=tracksum[5];
	  jetsum2[6]+=tracksum[6];
	  jetsum2[7]+=tracksum[7];
	  _sixjets++;}

    //fill jet pt plot instead of using counter to get entries as NJets for normalization later
    hJetPt->Fill(jetpt);

    //jet shape with bins of 0.05
    hJetShape->Fill(0.025, tracksum[0]/jetpt);
    hJetShape->Fill(0.075, tracksum[1]/jetpt);
    hJetShape->Fill(0.125, tracksum[2]/jetpt);
    hJetShape->Fill(0.175, tracksum[3]/jetpt);
    hJetShape->Fill(0.225, tracksum[4]/jetpt);
    hJetShape->Fill(0.275, tracksum[5]/jetpt);
    hJetShape->Fill(0.325, tracksum[6]/jetpt);
    hJetShape->Fill(0.375, tracksum[7]/jetpt);

    //jet shape with bins of 0.05 for a specific jet mass range
    hJetShape1->Fill(0.025, jetsum[0]/jetpt1);
    hJetShape1->Fill(0.075, jetsum[1]/jetpt1);
    hJetShape1->Fill(0.125, jetsum[2]/jetpt1);
    hJetShape1->Fill(0.175, jetsum[3]/jetpt1);
    hJetShape1->Fill(0.225, jetsum[4]/jetpt1);
    hJetShape1->Fill(0.275, jetsum[5]/jetpt1);
    hJetShape1->Fill(0.325, jetsum[6]/jetpt1);
    hJetShape1->Fill(0.375, jetsum[7]/jetpt1);
    
    hJetShape2->Fill(0.025, jetsum1[0]/jetpt2);
    hJetShape2->Fill(0.075, jetsum1[1]/jetpt2);
    hJetShape2->Fill(0.125, jetsum1[2]/jetpt2);
    hJetShape2->Fill(0.175, jetsum1[3]/jetpt2);
    hJetShape2->Fill(0.225, jetsum1[4]/jetpt2);
    hJetShape2->Fill(0.275, jetsum1[5]/jetpt2);
    hJetShape2->Fill(0.325, jetsum1[6]/jetpt2);
    hJetShape2->Fill(0.375, jetsum1[7]/jetpt2);

    hJetShape3->Fill(0.025, jetsum2[0]/jetpt3);
    hJetShape3->Fill(0.075, jetsum2[1]/jetpt3);
    hJetShape3->Fill(0.125, jetsum2[2]/jetpt3);
    hJetShape3->Fill(0.175, jetsum2[3]/jetpt3);
    hJetShape3->Fill(0.225, jetsum2[4]/jetpt3);
    hJetShape3->Fill(0.275, jetsum2[5]/jetpt3);
    hJetShape3->Fill(0.325, jetsum2[6]/jetpt3);
    hJetShape3->Fill(0.375, jetsum2[7]/jetpt3);
    
  } // jet loop

}

//
//
//________________________________________________________________________
void StEshaMaker::RunTracks()
{
  // constants: assume neutral pion mass
  double pi = 1.0*TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV
  unsigned int ntracks = mPicoDst->numberOfTracks();

  // loop over ALL tracks in PicoDst 
  for(unsigned short iTracks = 0; iTracks < ntracks; iTracks++){
    // get track pointer
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(iTracks));
    if(!trk){ continue; }

    // acceptance and kinematic quality cuts
    if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

    // primary track switch: get momentum vector of track - global or primary track
    TVector3 mTrkMom;
    if(doUsePrimTracks) {
      // get primary track vector
      mTrkMom = trk->pMom();
    } else {
      // get global track vector
      mTrkMom = trk->gMom(mVertex, Bfield);
    }

    // track variables
    double pt = mTrkMom.Perp();
    double phi = mTrkMom.Phi();
    double eta = mTrkMom.PseudoRapidity();
    double px = mTrkMom.x();
    double py = mTrkMom.y();
    double pz = mTrkMom.z();
    short charge = trk->charge();

    // do some things with tracks here
    // ......

    // fill track histograms here
    // .........

  } // track loop

}  // track function

//
//
//________________________________________________________________________
void StEshaMaker::RunTowers()
{
  // constants: assume neutral pion mass
  double pi = 1.0*TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV

  // looping over clusters - STAR: matching already done, get # of clusters and set variables
  unsigned int nBEmcPidTraits = mPicoDst->numberOfBEmcPidTraits();

  // loop over ALL clusters in PicoDst
  for(unsigned short iClus = 0; iClus < nBEmcPidTraits; iClus++){
    StPicoBEmcPidTraits *cluster = static_cast<StPicoBEmcPidTraits*>(mPicoDst->bemcPidTraits(iClus));
    if(!cluster){ cout<<"Cluster pointer does not exist.. iClus = "<<iClus<<endl; continue; }

    // cluster and tower ID
    int clusID = cluster->bemcId();  // index in bemc point array
    int towID = cluster->btowId();   // projected tower Id: 1 - 4800
    int towID2 = cluster->btowId2(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
    int towID3 = cluster->btowId3(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
    if(towID < 0) continue;

    // cluster and tower position - from vertex and ID
    TVector3 towPosition = mEmcPosition->getPosFromVertex(mVertex, towID);
    double towPhi = towPosition.Phi();
    double towEta = towPosition.PseudoRapidity();

    // matched track index
    int trackIndex = cluster->trackIndex();
    // get track pointer
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackIndex));
    if(!trk) { cout<<"No trk pointer...."<<endl; continue; }
    if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

  } // BEmc loop

  // loop over towers
  int nTowers = mPicoDst->numberOfBTowHits();
  for(int itow = 0; itow < nTowers; itow++) {
    // get tower pointer
    StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(itow));
    if(!tower) { cout<<"No tower pointer... iTow = "<<itow<<endl; continue; }

    // tower ID: tower ID shifted by +1 from array index
    int towerID = itow + 1;
    if(towerID < 0) continue; // double check these aren't still in the event list

    // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
    TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
    double towerPhi = towerPosition.Phi();
    double towerEta = towerPosition.PseudoRapidity();
    int towerADC = tower->adc();
    double towerE = tower->energy();
    double towerEt = towerE / (1.0*TMath::CosH(towerEta)); // THIS should be USED

    // do stuff with towers and fill histograms here
    // ........


  } // tower loop

}  // run towers function

//
//
// __________________________________________________________________________________
void StEshaMaker::SetSumw2() {
  hCentrality->Sumw2();
  hMultiplicity->Sumw2();
  hJetPt->Sumw2();
  hJetCorrPt->Sumw2();
  hJetMass->Sumw2();
  hTrackPt->Sumw2();
  hJetShape->Sumw2();
  hJetShape1->Sumw2();
  hJetShape2->Sumw2();
  hJetShape3->Sumw2();
  hJetPt1->Sumw2();
  hJetPt2->Sumw2();
  hJetPt3->Sumw2();
}

//
// Function: get relative phi of jet and track (-0.5pi, 1.5pi)
//________________________________________________________________________
Double_t StEshaMaker::RelativePhi(Double_t mphi,Double_t vphi) const
{ // function to calculate relative PHI
  double dphi = mphi-vphi;

  // set dphi to operate on adjusted scale
  if(dphi < -0.5*pi) dphi += 2.0*TMath::Pi();
  if(dphi >  1.5*pi) dphi -= 2.0*TMath::Pi();

  // test check
  if( dphi < -0.5*pi || dphi > 1.5*pi )
    Form("%s: dPHI not in range [-0.5*Pi, 1.5*Pi]!", GetName());

  return dphi; // dphi in [-0.5Pi, 1.5Pi]                                                                                   
}

//
// Function: calculate angle between jet and EP in the 1st quadrant (0,Pi/2)
//_________________________________________________________________________
Double_t StEshaMaker::RelativeEPJET(Double_t jetAng, Double_t EPAng) const
{
  Double_t pi = 1.0*TMath::Pi();
  Double_t dphi = 1.0*TMath::Abs(EPAng - jetAng);
  
  // ran into trouble with a few dEP<-Pi so trying this...
  if( dphi < -pi ){
    dphi = dphi + pi;
  } // this assumes we are doing full jets currently 
 
  if(dphi > 1.5*pi) dphi -= 2.0*pi;
  if((dphi > 1.0*pi) && (dphi < 1.5*pi)) dphi -= 1.0*pi;
  if((dphi > 0.5*pi) && (dphi < 1.0*pi)) dphi -= 1.0*pi;
  dphi = 1.0*TMath::Abs(dphi);

  // test check
  if( dphi < 0.0 || dphi > 0.5*pi ) {
    //Form("%s: dPHI not in range [0, 0.5*Pi]!", GetName());
    cout<<"dPhi not in range [0, 0.5*Pi]!"<<endl;
  }

  return dphi;   // dphi in [0, Pi/2]
}

//
//
//_________________________________________________________________________
void StEshaMaker::FillEmcTriggers() {
  // number of Emcal Triggers
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  Int_t nEmcTrigger = mPicoDst->numberOfEmcTriggers();

  // set kAny true to use of 'all' triggers
  fEmcTriggerArr[StJetFrameworkPicoBase::kAny] = 1;  // always TRUE, so can select on all event (when needed/wanted) 

  // loop over valid EmcalTriggers
  for(int i = 0; i < nEmcTrigger; i++) {
    // get trigger pointer
    StPicoEmcTrigger *emcTrig = static_cast<StPicoEmcTrigger*>(mPicoDst->emcTrigger(i));
    if(!emcTrig) continue;

    // fill for valid triggers
    if(emcTrig->isHT0()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT0] = 1; }
    if(emcTrig->isHT1()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT1] = 1; }
    if(emcTrig->isHT2()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT2] = 1; }
    if(emcTrig->isHT3()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT3] = 1; }
    if(emcTrig->isJP0()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP0] = 1; }
    if(emcTrig->isJP1()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP1] = 1; }
    if(emcTrig->isJP2()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP2] = 1; }
  }
}
//
// elems: sizeof(myarr)/sizeof(*myarr) prior to passing to function
// upon passing the array collapses to a pointer and can not get size anymore
//________________________________________________________________________
Bool_t StEshaMaker::DoComparison(int myarr[], int elems) {
  //std::cout << "Length of array = " << (sizeof(myarr)/sizeof(*myarr)) << std::endl;
  bool match = kFALSE;

  // loop over specific physics selection array and compare to specific event trigger
  for(int i=0; i<elems; i++) {
    if(mPicoEvent->isTrigger(myarr[i])) match = kTRUE;
    if(match) break;
  }

  return match;
}
