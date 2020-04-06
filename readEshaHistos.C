// ************************************************ //
// Author: Esha Rao
// Rutgers University
//
// March 19, 2020
// Updated to now make three canvases: one for jet mass of all pp data, one for jet shape of all pp data, and one for jet shape for bins of jet mass on the same canvas
//
// Macro for reading histograms and normalizing historgrams from readEshaDst
// March 6, 2020
// reading into Mar 4 Test 2 to grab the jetpt and jetshape histograms to get the number of events from jetpt to normalize the jetshape plot
//

#include <TSystem.h>

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include <THnSparse.h>
#include "TParameter.h"
#include <TProfile.h>
#include "TRandom.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TChain.h"

//C++ includes
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

// basic STAR classes
class StMemStat;
class StMaker;
class StChain;
class StPicoDstMaker;

// my added STAR classes
class StJetMakerTask;
class StRho;
class StRhoBase;
class StRhoSparse;
class StMyAnalysisMaker;

void readEshaHistos()
{
  
  //Get input
  TFile *finput = TFile::Open("jobhisto1.root");
  TDirectory* dir = gFile->GetDirectory("EshaMaker");

  //create output
  TFile *fout = new TFile("finalhisto1.root", "recreate");
  fout->cd();
  cout<<"Before declaring histograms"<<endl;

  //get any necessary histos
  TH1F*hJetPtALL = (TH1F*)dir->Get("hJetPt");
  hJetPtALL->SetName("hJetPtALL");
  TH1F*hJetMassALL = (TH1F*)dir->Get("hJetMass");
  hJetMassALL->SetName("hJetMassALL");
  TH1F*hJetShapeALL = (TH1F*)dir->Get("hJetShape");
  hJetShapeALL->SetName("hJetShapeNormalizedALL");

  //get necessary histos for each separate jet pt
  TH1F*hJetPtONE = (TH1F*)dir->Get("hJetPt1");
  hJetPtONE->SetName("hJetPtONE");
  TH1F*hJetPtTWO = (TH1F*)dir->Get("hJetPt2");
  hJetPtTWO->SetName("hJetPtTWO");
  TH1F*hJetPtTHREE = (TH1F*)dir->Get("hJetPt3");
  hJetPtTHREE->SetName("hJetPtTHREE");

  //get necessary histos for each separate jet shape
  TH1F*hJetShapeONE = (TH1F*)dir->Get("hJetShape1");
  hJetShapeONE->SetName("hJetShapeONE");
  TH1F*hJetShapeTWO = (TH1F*)dir->Get("hJetShape2");
  hJetShapeTWO->SetName("hJetShapeTWO");
  TH1F*hJetShapeTHREE = (TH1F*)dir->Get("hJetShape3");
  hJetShapeTHREE->SetName("hJetShapeTHREE");


  //jet pt count variable
  int jetptcount = hJetPtALL->GetEntries();

  //jet pt count for each separate histo variables
  int jetptcountONE = hJetPtONE->GetEntries();
  int jetptcountTWO = hJetPtTWO->GetEntries();
  int jetptcountTHREE = hJetPtTHREE->GetEntries();

  cout<<"Total Jet Pt Count: "<<jetptcount<<endl;
  cout<<"0-2 GeV Jet Pt Count: "<<jetptcountONE<<endl;
  cout<<"2-4 GeV Jet Pt Count: "<<jetptcountTWO<<endl;
  cout<<"4-6 GeV Jet Pt Count: "<<jetptcountTHREE<<endl;

  //scale jet shape with the bin width and jet pt count
  hJetShapeALL->Scale(1./jetptcount);
  hJetShapeALL->Scale(1./0.05);

  //scale jet shape for each bin of jet mass with the bin width and jet pt count for that bin
  hJetShapeONE->Scale(1./jetptcountONE);
  hJetShapeONE->Scale(1./0.05);
  hJetShapeTWO->Scale(1./jetptcountTWO);
  hJetShapeTWO->Scale(1./0.05);
  hJetShapeTHREE->Scale(1./jetptcountTHREE);
  hJetShapeTHREE->Scale(1./0.05);

  //make canvas for jet mass plot
  TCanvas * cjetmass = new TCanvas("cjetmass","Jet Mass STAR pp data 200 GeV", 600, 400);
  cjetmass->Divide(1, 1);
  cjetmass->cd(1);
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  hJetMassALL->SetXTitle("M_{J} (GeV/c^{2})");
  gPad->SetLogy();
  //hJetMassALL->GetYaxis()->SetRangeUser(0.0003, 35.);
  hJetMassALL->SetYTitle("(1/N_{jets})dN/dM_{J}");
  //hJetMassALL->SetMarkerColor(kBlack);
  hJetMassALL->SetTitle("Jet Mass STAR pp data 200 GeV");
  hJetMassALL->SetLineWidth(2);
  hJetMassALL->SetLineColor(kRed);
  hJetMassALL->Draw("E1");
  hJetMassALL->Draw("SAME");

  //legend for jet mass plot
  auto legend1 = new TLegend(0.45,0.6,0.8,0.8);
  legend1->SetBorderSize(0);
  //legend->SetHeader("Jet Pt > 10 GeV","C"); // option "C" allows to center the header
  legend1->AddEntry("hJetMassALL", "pp", "l");
  legend1->AddEntry("","anti-k_{T}, R = 0.4 Jets","");
  legend1->AddEntry("","Jet p_{T} > 10 GeV/c; #left|#eta#right| < 1","");
  legend1->Draw();


  //trying to make canvas to put jet shape normalized on
  TCanvas * cjetshape = new TCanvas("cjetshape","Jet Shape STAR pp data 200 GeV", 600, 400);
  cjetshape->Divide(1, 1);
  cjetshape->cd(1);
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky(); 
  gPad->SetGridx(0);
  gPad->SetGridy(0); 
  hJetShapeALL->SetXTitle("r from jet axis");
  gPad->SetLogy();
  //hJetShapeALL->GetYaxis()->SetRangeUser(0.003, 15.);
  hJetShapeALL->SetYTitle("#rho(r)");
  //hJetShapeALL->SetMarkerSize(100);
  hJetShapeALL->SetTitle("Jet Shape STAR pp data 200 GeV");
  hJetShapeALL->SetLineWidth(2);
  hJetShapeALL->SetLineColor(kRed);
  hJetShapeALL->Draw("");
  //hJetShapeALL->Draw("SAME");

  //jet shape ALL legend
  auto legend2 = new TLegend(0.15,0.2,0.5,0.4);
  //legend2->SetHeader("Jet Pt > 10 GeV","C"); // option "C" allows to center the header
  legend2->AddEntry("hJetShapeALL", "pp", "l");
  legend2->SetBorderSize(0);
  legend2->AddEntry("","anti-k_{T}, R = 0.4 Jets","");
  legend2->AddEntry("","Jet p_{T} > 10 GeV/c; #left|#eta#right| < 1","");
  legend2->Draw();

  //jet pt ALL canvas
  TCanvas * cjetpt = new TCanvas("cjetpt","Jet p_{T} STAR pp data 200 GeV", 600, 400);
  cjetpt->Divide(1, 1);
  cjetpt->cd(1);
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  hJetPtALL->SetXTitle("p_{T} (GeV)");
  gPad->SetLogy();
  //hJetPtALL->GetYaxis()->SetRangeUser(0.0003, 35.);
  hJetPtALL->SetYTitle("Counts");
  //hJetPtALL->SetMarkerColor(kBlack);
  hJetPtALL->SetTitle("Jet p_{T} STAR pp data 200 GeV");
  hJetPtALL->SetLineWidth(2);
  hJetPtALL->SetLineColor(kRed);
  hJetPtALL->Draw("");
  //hJetPtALL->Draw();
  
  //jet pt ALL legend
  auto legend3 = new TLegend(0.45,0.6,0.8,0.8);
  legend3->SetBorderSize(0);
  legend3->AddEntry("hJetPtALL", "pp", "l");
  legend3->AddEntry("","anti-k_{T}, R = 0.4 Jets","");
  legend3->AddEntry("","Jet p_{T} > 10 GeV/c; #left|#eta#right| < 1","");
  legend3->Draw();

  //jet shape for varying jet mass bins
  TCanvas * cjetshape1 = new TCanvas("cjetshape1","Jet Shape for Varying Jet Mass Bins", 600, 400);
  cjetshape1->Divide(1, 1);
  cjetshape1->cd(1);
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();
  hJetShapeONE->SetXTitle("r from jet axis");
  gPad->SetLogy();
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  //hJetShapeONE->GetYaxis()->SetRangeUser(0.005, 12.);
  hJetShapeONE->SetYTitle("#rho(r)");
  //hJetShapeONE->SetMarkerColor(kBlack);
  hJetShapeONE->SetLineColor(kBlack);
  hJetShapeTWO->SetLineColor(kRed);
  hJetShapeTHREE->SetLineColor(kBlue);
  hJetShapeONE->SetLineWidth(1);
  hJetShapeTWO->SetLineWidth(1);
  hJetShapeTHREE->SetLineWidth(1);
  //Style->SetEndErrorSize(8);
  hJetShapeONE->SetLineStyle(1);
  hJetShapeTWO->SetLineStyle(8);
  hJetShapeTHREE->SetLineStyle(3);
  hJetShapeONE->SetTitle("Jet Shape for Varying M_{J} STAR pp data");
  hJetShapeONE->Draw("");
  hJetShapeTWO->Draw("SAME");
  hJetShapeTHREE->Draw("SAME");
  //hJetShapeONE->Draw("SAME");
  //hJetShapeTWO->Draw("SAME");
  //hJetShapeTHREE->Draw("SAME");

  cout<<"Before writing histograms"<<endl;

  //write histos
  hJetPtALL->Write();
  hJetMassALL->Write();
  hJetShapeALL->Write();
  hJetPtONE->Write();
  hJetPtTWO->Write();
  hJetPtTHREE->Write();
  hJetShapeONE->Write();
  hJetShapeTWO->Write();
  hJetShapeTHREE->Write();

  //trying to make a legend
   auto legend = new TLegend(0.15,0.2,0.5,0.4);
   legend->SetBorderSize(0);
   legend->SetFillStyle(0);
   //legend->SetHeader("Jet Pt > 10 GeV","C"); // option "C" allows to center the header
   legend->AddEntry("hJetShapeONE","0 < M_{J} < 2","l");
   legend->AddEntry("hJetShapeTWO","2 < M_{J} < 4","l");
   legend->AddEntry("hJetShapeTHREE","4 < M_{J} < 6","l");
   legend->AddEntry("","anti-k_{T}, R = 0.4 Jets","");
   legend->AddEntry("","Jet p_{T} > 10 GeV/c; #left|#eta#right| < 1","");
   legend->Draw();

  cout<<"After writing histograms"<<endl;

  //close output file
  fout->Close();
}
