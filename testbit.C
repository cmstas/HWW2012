// ***************************************************************************
// This macro is made to reproduce Sync cut flow at smurf ntuple level 
// One can check following things
//  1) Bits are implemented correctly  
//  2) variables are stored correclty 
// 
// - Cut flow is made by combination of bits(variable : bit) and cuts(variable : addcut)
// - One can add bits and cuts successively to get yields 
// 
// Usage : 
// 		$ root 
// 	 	root[] .L testbit.C++
//  	root[] testbit("qqww.root")
//
// ***************************************************************************

#include "../Smurf/Core/SmurfTree.h"
#include "TChain.h"
#include "TTree.h"
#include "TCut.h"
#include "TH1F.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include <iostream> 

using namespace std;

// as a reference
enum Type {
	mm,
	me,
	em,
	ee
};

/*
unsigned int FullSelection	= 	SmurfTree::BaseLine |
								SmurfTree::ChargeMatch |
								SmurfTree::Lep1FullSelection |
								SmurfTree::Lep2FullSelection |
								SmurfTree::FullMET |
								SmurfTree::ZVeto |
								SmurfTree::TopVeto |
								SmurfTree::ExtraLeptonVeto;
*/

void printYield(TChain *ch, TString cutname, unsigned int bit, TString addcut) {

	TH1F *hmm 	= new TH1F("hmm", 	"hmm", 	10, -0.5, 19.5);
	TH1F *hme 	= new TH1F("hme", 	"hme", 	10, -0.5, 19.5);
	TH1F *hem 	= new TH1F("hem", 	"hem", 	10, -0.5, 19.5);
	TH1F *hee 	= new TH1F("hee", 	"hee", 	10, -0.5, 19.5);
	TH1F *hall 	= new TH1F("hall", 	"hall", 10, -0.5, 19.5);

	TString bitmm = Form("((cuts & %u)==%u) && type==0", bit, bit);
	TString bitme = Form("((cuts & %u)==%u) && type==1", bit, bit);
	TString bitem = Form("((cuts & %u)==%u) && type==2", bit, bit);
	TString bitee = Form("((cuts & %u)==%u) && type==3", bit, bit);
	TString bitall = Form("((cuts & %u)==%u)", bit, bit);

	ch->Draw("min(njets,19.4999)>>hmm", 	bitmm+"&&"+addcut, "goff");
	ch->Draw("min(njets,19.4999)>>hme", 	bitme+"&&"+addcut, "goff");
	ch->Draw("min(njets,19.4999)>>hem", 	bitem+"&&"+addcut, "goff");
	ch->Draw("min(njets,19.4999)>>hee", 	bitee+"&&"+addcut, "goff");
	ch->Draw("min(njets,19.4999)>>hall", bitall+"&&"+addcut, "goff");

	cout.width(40); 
    cout << cutname  << " \t"; 
	cout << hmm->Integral()  << " \t" 
		 << hee->Integral()  << " \t" 
		 << hem->Integral()  << " \t" 
		 << hme->Integral()  << " \t" 
		 << hall->Integral() << endl;

	delete hmm;
	delete hme;
	delete hem;
	delete hee;
	delete hall;
}

void testbit(TString inputsmurf="smurf/qqww_0_999999.root") {

  	TChain *ch = new TChain("tree");
    ch->Add(inputsmurf);

	// header
	cout.width(40); 
    cout << "Cuts \t\t";
	cout << "mm \tee \tem \tme \tall" << endl;

	// dummy additional cut
	TString addcut="1";

	//
	unsigned int bit = SmurfTree::BaseLine;
	printYield(ch, "Baseline", 	bit, addcut);

	//
	bit = bit | SmurfTree::ChargeMatch;
	printYield(ch, "Charge match", 	bit, addcut);

	//
	bit = bit | SmurfTree::Lep1FullSelection | SmurfTree::Lep2FullSelection;
	printYield(ch, "Full lepton selection", 	bit, addcut);

    //
	bit = bit | SmurfTree::ExtraLeptonVeto;
	printYield(ch, "Extralepton veto", 	bit, addcut);

	//
 	addcut = addcut+"&&met>20";
	printYield(ch, "met > 20 GeV", 	bit, addcut);
	
	//
 	addcut = addcut+"&&dilep.mass()>12";
	printYield(ch, "mll > 12 GeV", 	bit, addcut);
	
	//
	bit = bit | SmurfTree::ZVeto;
	printYield(ch, "|mll - mZ| > 15 GeV", 	bit, addcut);
	
	//
 	//addcut = addcut+"&&(type==0 || type==3 || min(pmet,pTrackMet)>20)";
	//printYield(ch, "minMET > 20 GeV for em/me", 	bit, addcut);
 	addcut = addcut+"&&min(pmet,pTrackMet)>20";
	printYield(ch, "minMET > 20 GeV", 	bit, addcut);
	
	//
 	addcut = addcut+"&&(type==1 || type==2  || min(pmet,pTrackMet)>40)";
	printYield(ch, "minMET > 40 GeV for ee/mm", 	bit, addcut);
	
	// needs some aliases for calculate dPhi(jet1+jet2, dilep)
	// because ROOT does not support (jet1+jet2).pt() in Draw
	ch->SetAlias("jpT1", 		"TMath::Sqrt(jet1.px()*jet1.px()+jet1.py()*jet1.py())");
 	addcut = addcut+"&&(type==1 || type==2 || njets>1 || jet1.pt()<15.|| (abs(TMath::ACos((jet1.px()*dilep.px()+jet1.py()*dilep.py())/jpT1/dilep.pt()))<(165.*TMath::Pi()/180.)))";
	//printYield(ch, "dPhiDiLepJet < 165 dg for ee/mm", 	bit, addcut);
	
    ch->SetAlias("jpx", 	"jet1.px()+jet2.px()");
	ch->SetAlias("jpy", 	"jet1.py()+jet2.py()");
	ch->SetAlias("dilpx", 	"dilep.px()");
	ch->SetAlias("dilpy", 	"dilep.py()");
	ch->SetAlias("jpT12", 		"TMath::Sqrt(jpx*jpx+jpy*jpy)");
	ch->SetAlias("dilpT", 	"TMath::Sqrt(dilpx*dilpx+dilpy*dilpy)");
 	addcut = addcut+"&&(type==1 || type==2 || njets<2 || (abs(TMath::ACos((jpx*dilpx+jpy*dilpy)/jpT12/dilpT))<(165.*TMath::Pi()/180.)))";
	printYield(ch, "dPhiDiLepJet < 165 dg for ee/mm", 	bit, addcut);

	//
 	addcut = addcut+"&&nSoftMuons==0";
	printYield(ch, "SoftMuons==0", 	bit, addcut);
	
    //
	bit = bit | SmurfTree::TopVeto;
	printYield(ch, "Top Veto", 	bit, addcut);

	//
 	addcut = addcut+"&&dilep.pt()>45";
	printYield(ch, "ptll > 45 GeV", 	bit, addcut);

	//
	//bit = bit | SmurfTree::FullMET;
	//printYield(ch, "cutname", 	bit, addcut);

	// 
	TString finalcut = addcut;

	// 0 jet bin
	cout.width(40); 
	cout << "---------- 0 jet ---------" << endl;
 	addcut = finalcut+"&&njets==0";
	printYield(ch, "njets == 0", 	bit, addcut);
 	
	//
	addcut = addcut+"&&max(lep1.pt(),lep2.pt())>30";
	printYield(ch, "max(lep1.pt(),lep2.pt())>30", 	bit, addcut);

	//
	addcut = addcut+"&&min(lep1.pt(),lep2.pt())>25";
	printYield(ch, "min(lep1.pt(),lep2.pt())>25", 	bit, addcut);
	
	// 1 jet bin
	cout.width(40); 
	cout << "---------- 1 jet ---------" << endl;
 	addcut = finalcut+"&&njets==1";
	printYield(ch, "njets==1", 	bit, addcut);

	// 2 jet bin
	cout.width(40); 
	cout << "---------- 2 jet ---------" << endl;
 	addcut = finalcut+"&&(njets==2 || njets==3)";
	printYield(ch, "njets==2 or 3", 	bit, addcut);
 	
	//
	addcut = addcut+"&&abs(jet1.eta())<4.7 && abs(jet2.eta())<4.7";
	printYield(ch, "abs(jet1.eta())<4.7 && abs(jet2.eta())<4.7", 	bit, addcut);
 	
	//
	addcut = addcut+"&& (njets==2 || jet3.pt()<30 || (((jet1.eta()-jet3.eta())<0 && (jet2.eta()-jet3.eta())<0 ) || ((jet2.eta()-jet3.eta())>0 && (jet1.eta()-jet3.eta())>0)))";
	printYield(ch, "no central jets", 	bit, addcut);

}
