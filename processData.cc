#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRegexp.h"
#include "TFile.h" 
#include "../HWW2012CORE/wwtypes.h"
#include "doAnalysis.h"

#include <stdlib.h>

int main(int argc, char *argv[])
{

  using namespace std;

  //
  // if the arguments are provided
  // then run one file
  //

  if (argc == 3) {
    SmurfTree::DataType dataType = (SmurfTree::DataType)atoi(argv[1]);
    std::string dataFile = argv[2];
    bool realData = false;
    if (dataType == 0) realData = true;
    const std::string cms2_json_file = "files/Cert_190456-195396_8TeV_PromptReco_Collisions12_JSON_v2.jmu";
    ProcessSample(dataFile, dataType, 3.6, -1, -1, false, realData, "");
	return 0;
  }
  
  if (argc == 4) {
  
    std::vector<string> dataSamples;

//  	// 2012A
// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/DoubleElectron_Run2012A-PromptReco-v1/V05-02-28/vvfilter/*.root");
// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/DoubleMu_Run2012A-PromptReco-v1/V05-02-28/vvfilter/*.root");
// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/MuEG_Run2012A-PromptReco-v1/V05-02-28/vvfilter/preprocessing/*.root");
// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/SingleElectron_Run2012A-PromptReco-v1/V05-02-28/vvfilter/*.root");
// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/SingleMu_Run2012A-PromptReco-v1/V05-02-28/vvfilter/*.root");

//  	// 2012B
// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/DoubleElectron_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
//  	dataSamples.push_back("/nfs-6/userdata/jaehyeok/DoubleMu_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
//  	dataSamples.push_back("/nfs-6/userdata/jaehyeok/MuEG_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
//  	dataSamples.push_back("/nfs-6/userdata/jaehyeok/SingleElectron_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
//  	dataSamples.push_back("/nfs-6/userdata/jaehyeok/SingleMu_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
 	
// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/jun11_from195538/DoubleElectron_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
//  	dataSamples.push_back("/nfs-6/userdata/jaehyeok/jun11_from195538/DoubleMu_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
//  	dataSamples.push_back("/nfs-6/userdata/jaehyeok/jun11_from195538/MuEG_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
//  	dataSamples.push_back("/nfs-6/userdata/jaehyeok/jun11_from195538/SingleElectron_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/jun11_from195538/SingleMu_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
 	
// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/jun17_from195904/DoubleElectron_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
//  	dataSamples.push_back("/nfs-6/userdata/jaehyeok/jun17_from195904/DoubleMu_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
//  	dataSamples.push_back("/nfs-6/userdata/jaehyeok/jun17_from195904/MuEG_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/jun17_from195904/SingleElectron_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/jun17_from195904/SingleMu_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
 	
// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/jun20_from196357/DoubleElectron_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
//  	dataSamples.push_back("/nfs-6/userdata/jaehyeok/jun20_from196357/DoubleMu_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
//  	dataSamples.push_back("/nfs-6/userdata/jaehyeok/jun20_from196357/MuEG_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
//  	dataSamples.push_back("/nfs-6/userdata/jaehyeok/jun20_from196357/SingleElectron_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/jun20_from196357/SingleMu_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");

// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/missing/DoubleMu_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/missing/DoubleElectron_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/missing/MuEG_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/missing/SingleElectron_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
// 	dataSamples.push_back("/nfs-6/userdata/jaehyeok/missing/SingleMu_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");

    SmurfTree::DataType dataType = (SmurfTree::DataType)atoi(argv[1]);
    int beginrun = atoi(argv[2]);
    int endrun = atoi(argv[3]);
    bool realData = false;
    if (dataType == 0) realData = true;
    const std::string cms2_json_file = "Cert_190456-202016_8TeV_PromptReco_Collisions12_JSON_cms2.txt";
	
    for(int run = beginrun; run<=endrun; run++) { 
//		dataSamples.push_back(Form("/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/DoubleElectron_Run2012C-PromptReco-v1_AOD/merged/*%d*.root",run));
//		dataSamples.push_back(Form("/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/DoubleMu_Run2012C-PromptReco-v1_AOD/merged/*%d*.root",run));
//		dataSamples.push_back(Form("/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/MuEG_Run2012C-PromptReco-v1_AOD/merged/*%d*.root",run));
//		dataSamples.push_back(Form("/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012C-PromptReco-v1_AOD/merged/*%d*.root",run));
//		dataSamples.push_back(Form("/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012C-PromptReco-v1_AOD/merged/*%d*.root",run));
		dataSamples.push_back(Form("/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/DoubleElectron_Run2012C-PromptReco-v2_AOD/merged/*%d*.root",run));
//		dataSamples.push_back(Form("/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/DoubleMu_Run2012C-PromptReco-v2_AOD/merged/*%d*.root",run));
//		dataSamples.push_back(Form("/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/MuEG_Run2012C-PromptReco-v2_AOD/merged/*%d*.root",run));
//		dataSamples.push_back(Form("/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012C-PromptReco-v2_AOD/merged/*%d*.root",run));
//		dataSamples.push_back(Form("/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012C-PromptReco-v2_AOD/merged/*%d*.root",run));
   	} 
	
	ProcessSample(dataSamples, dataType, 3.5, -1, -1, false, realData, "", beginrun, endrun);
    return 0;
  
  }


  //
  // Define how to handle complex Monte Carlo samples which may have 
  // more than one event type mixed in. Careful with this option. 
  // You don't need it for not mixed samples.
  //

  const bool identifyDYEvents = true;

  //
  // Flags for files to run over 
  // (0 and 1 are easier to modify)
  //
  bool runTest   = 1;
  bool runWW     = 0;
  bool runWWmcnlo= 0;
  bool runWWup   = 0;
  bool runWWdown = 0;
  bool runHWW    = 0;
  bool runGGWW   = 0;
  bool runWZ     = 0;
  bool runZZ     = 0;
  bool runZZpy   = 0;
  bool runWjets  = 0;
  bool runWgamma = 0;
  bool runWgstar = 0;
  bool runDYmgee = 0;
  bool runDYmgmm = 0;
  bool runDYee   = 0;
  bool runDYmm   = 0;
  bool runDYtt   = 0;
  bool runttbar  = 0;
  bool runttbarmg= 0;
  bool runtW     = 0;
  bool runtWds   = 0;
  bool runEmb    = 0;
  bool runData   = 0;

  // 
  // Ntuple version
  //
  string version = "wwfilter";
  if (gSystem->Getenv("VERSION")){
    version = gSystem->Getenv("VERSION");
    cout << "VERSION: " << version << endl;
  }
 
  // read dataset prefix
  string dataset = "data";
  if (gSystem->Getenv("DataDir")){
    dataset = gSystem->Getenv("DataDir");
    cout << "Dataset directory: " << dataset << endl;
  }
  
  const double integratedLumi = 1000.0; // pb^1  
  cout << "Integrated luminosity to scale to: " << integratedLumi << endl;

  if (runTest){
    std::vector<string> samples;
    //samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/merged_ntuple.root");
    //ProcessSample(samples, SmurfTree::hww160, integratedLumi, -1, -1);
    //samples.push_back("/nfs-6/userdata/jaehyeok/WW_TuneZ2star_8TeV_pythia6_tauola_Summer12-PU_S7_START50_V15-v1/V05-01-01/nofilter/postprocessing_forsync/*.root");
    //samples.push_back("/tas/jaehyeok/home/HWW/MakeCMS2ntuple/CMSSW_5_2_3_patch4/src/post_ntuple_A6DE4085-8191-E111-BF4E-001E67396D5B.root");
    //ProcessSample(samples, SmurfTree::qqww, integratedLumi, -1, -1);
    
	//samples.push_back("/home/users/jaehyeok/HWW/MakeCMS2ntuples/ICHEP2012/CMSSW_5_2_3_patch4_V05-02-28_dataJECfix/src/ntuple.root");
    //const std::string cms2_json = "files/Cert_190456-195396_8TeV_PromptReco_Collisions12_JSON_v2.jmu";
    //ProcessSample(samples, SmurfTree::data, 3.1, -1, -1, false, true, "");
    
	//ProcessSample("data/post_ntuple.root", SmurfTree::qqww, 1., -1, -1, false, false, "");
	ProcessSample("data/post_Slimntuple_START53_V15.root", SmurfTree::qqww, 1., -1, -1, false, false, "");
  }

  if (runWW)
    ProcessSample(dataset+"/WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1/V04-02-38/vvfilter/*.root", 
		  SmurfTree::qqww, integratedLumi, -1, -1, true);

  if (runWWmcnlo)
    ProcessSample(dataset+"/WWTo2L2Nu_CT10_7TeV-mcatnlo_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root", 
		  SmurfTree::qqww, integratedLumi, -1, -1, true);

  if (runWWup)
    ProcessSample(dataset+"/WWTo2L2Nu_scaleup_CT10_7TeV-mcatnlo_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root", 
		  SmurfTree::qqww, integratedLumi, -1, -1, true);

  if (runWWdown)
    ProcessSample(dataset+"/WWTo2L2Nu_scaledown_CT10_7TeV-mcatnlo_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root", 
		  SmurfTree::qqww, integratedLumi, -1, -1, true);

  if (runWgamma) {
    std::vector<string> samples;
    samples.push_back(dataset+"/WGToMuNuG_TuneZ2_7TeV-madgraph_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/WGToENuG_TuneZ2_7TeV-madgraph_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/WGToTauNuG_TuneZ2_7TeV-madgraph_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples,SmurfTree::wgamma, integratedLumi, -1, -1, true);
  }

  if (runWgstar) {
    std::vector<string> samples;
    samples.push_back(dataset+"/WGstarToLNu2E_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/WGstarToLNu2Mu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-36/DiLeptonFilter/*.root");
    ProcessSample(samples,SmurfTree::wgstar, integratedLumi, -1, -1, true);
  }

  if (runGGWW)
    ProcessSample(dataset+"/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root", 
		  SmurfTree::ggww, integratedLumi, -1, -1, true);

  if (runHWW){
    std::vector<string> samples;

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2LAndTau2Nu_M-110_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back("/nfs-7a/userdata/cerati/GluGluToHToWWTo2LAndTau2Nu_M-110_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww110, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2LAndTau2Nu_M-115_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back("/nfs-7a/userdata/cerati/GluGluToHToWWTo2LAndTau2Nu_M-115_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww115, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-120_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-120_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    //samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-120_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-120_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww120, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    //samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-130_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-130_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww130, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-140_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-140_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-140_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-140_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww140, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-150_7TeV-powheg-pythia6_Fall11-PU_S6-START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-150_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    //samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-150_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-150_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww150, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-6/userdata/cerati/VBF_HToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-160_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-160_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww160, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-170_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-170_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-170_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-170_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww170, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-180_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-180_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-180_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-180_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww180, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-190_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-190_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    //samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-190_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-190_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww190, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-200_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-200_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-200_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-200_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww200, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-210_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-210_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-210_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-210_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww210, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-220_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-220_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-220_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-220_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww220, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-230_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-230_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-230_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-230_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww230, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-250_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-250_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-250_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-250_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww250, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-300_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-300_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-300_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-300_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww300, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-350_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-350_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-350_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-350_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww350, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-400_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-400_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    //samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-400_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-400_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww400, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-450_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-450_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-450_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-450_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww450, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-500_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-500_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-500_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-500_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww500, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-550_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-550_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-550_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-550_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww550, integratedLumi, -1, -1); samples.clear();

    samples.push_back("/nfs-7a/userdata/cerati/VBF_HToWWTo2L2Nu_M-600_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2L2Nu_M-600_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/DiLeptonFilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWTo2Tau2Nu_M-600_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/GluGluToHToWWToLNuTauNu_M-600_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::hww600, integratedLumi, -1, -1); samples.clear();

  }

  if (runWZ)
    ProcessSample("/nfs-7a/userdata/cerati/WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root", SmurfTree::wz, integratedLumi, -1, -1);
  
  if (runZZ)
    ProcessSample(dataset+"/ZZJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root", SmurfTree::zz, integratedLumi, -1, -1);
  
  if (runZZpy)
    ProcessSample("/nfs-7a/userdata/cerati/ZZ_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/merged_ntuple*.root", SmurfTree::zz, integratedLumi, -1, -1);

  if (runWjets){
    std::vector<string> wSamples;
    wSamples.push_back(dataset+"/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(wSamples, SmurfTree::wjets,  integratedLumi , -1, -1);
  }

  if (runttbarmg)
    ProcessSample(dataset+"/TTJets_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/merged_ntuple*.root", 
		  SmurfTree::ttbar, integratedLumi, -1, -1);

  if (runttbar)
    ProcessSample(dataset+"/TTTo2L2Nu2B_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*root", 
		  SmurfTree::ttbar, integratedLumi, -1, -1);

  if (runtW) {
    std::vector<string> samples;
    samples.push_back(dataset+"/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back("/nfs-6/userdata/cerati/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Fall11-PU_S6_START42_V14B-v2/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::tw, integratedLumi, -1, -1);
  }

  if (runtWds) {
    std::vector<string> samples;
    samples.push_back("/nfs-7a/userdata/cerati/T_TuneZ2_tW-channel-DS_7TeV-powheg-tauola_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/Tbar_TuneZ2_tW-channel-DS_7TeV-powheg-tauola_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::tw, integratedLumi, -1, -1);
  }

  if (runDYmgee){
    std::vector<string> samples;
    samples.push_back("/nfs-7a/userdata/cerati/DYJetsToLL_M-10To50_TuneZ2_7TeV-madgraph_Fall11-PU_S6_START42_V14B-v1/V04-02-36/merged_ntuple*.root");
    samples.push_back("/nfs-7a/userdata/cerati/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1/V04-02-36/merged_ntuple*.root");
    ProcessSample(samples, SmurfTree::dyee, integratedLumi, -1, -1, identifyDYEvents);
  }

  if (runDYmgmm){
    std::vector<string> samples;
    samples.push_back("/nfs-7a/userdata/cerati/DYJetsToLL_M-10To50_TuneZ2_7TeV-madgraph_Fall11-PU_S6_START42_V14B-v1/V04-02-36/merged_ntuple*.root");
    samples.push_back("/nfs-7a/userdata/cerati/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1/V04-02-36/merged_ntuple*.root");
    ProcessSample(samples, SmurfTree::dymm, integratedLumi, -1, -1, identifyDYEvents);
  }

  if (runDYee){
    std::vector<string> samples;
    samples.push_back(dataset+"/DYToEE_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    samples.push_back(dataset+"/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::dyee, integratedLumi, -1, -1, identifyDYEvents);
  }

  if (runDYmm){
    std::vector<string> samples;
    samples.push_back(dataset+"/DYToMuMu_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/merged_ntuple*.root");
    samples.push_back(dataset+"/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/merged_ntuple*.root");
    ProcessSample(samples, SmurfTree::dymm, integratedLumi, -1, -1, identifyDYEvents);
  }

  if (runDYtt){
    std::vector<string> samples;
    samples.push_back(dataset+"/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Fall11-PU_S6_START42_V14B-v1/V04-02-36/vvfilter/*.root");
    ProcessSample(samples, SmurfTree::dytt, integratedLumi, -1, -1, identifyDYEvents);
  }
 
  // RealData
  TString cms2_json_file = "files/hww.Full2011.txt";
  cms2_json_file = "";

  std::vector<string> embSamples;
  embSamples.push_back("/nfs-7a/userdata/cerati/DoubleMu_StoreResults-DoubleMu_2011A_May10thRR_v1_embedded_trans1_tau123_pttau1_18tau2_8_v1-f456bdbb960236e5c696adfe9b04eaae/V04-02-36/vvfilter/*.root");
  embSamples.push_back("/nfs-7a/userdata/cerati/DoubleMu_StoreResults-DoubleMu_2011A_PR_v4_embedded_trans1_tau123_pttau1_18tau2_8_v1-f456bdbb960236e5c696adfe9b04eaae/V04-02-36/vvfilter/*.root");
  if (runEmb)
    ProcessSample(embSamples, SmurfTree::dyttDataDriven, 3.1, -1, -1, false, true, cms2_json_file);

  std::vector<string> dataSamples;
  
  dataSamples.push_back(dataset+"/DoubleElectron_Run2011A-May10ReReco-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/DoubleMu_Run2011A-May10ReReco-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/MuEG_Run2011A-May10ReReco-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/SingleElectron_Run2011A-May10ReReco-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/SingleMu_Run2011A-May10ReReco-v1/V04-02-36/vvfilter/*.root");

  dataSamples.push_back(dataset+"/DoubleElectron_Run2011A-PromptReco-v4/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/MuEG_Run2011A-PromptReco-v4/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/SingleElectron_Run2011A-PromptReco-v4/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/SingleMu_Run2011A-PromptReco-v4/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/DoubleMu_Run2011A-PromptReco-v4/V04-02-36/vvfilter/*.root");

  dataSamples.push_back(dataset+"/DoubleElectron_Run2011A-05Aug2011-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/MuEG_Run2011A-05Aug2011-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/SingleElectron_Run2011A-05Aug2011-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/SingleMu_Run2011A-05Aug2011-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/DoubleMu_Run2011A-05Aug2011-v1/V04-02-36/vvfilter/*.root");

  dataSamples.push_back(dataset+"/DoubleElectron_Run2011A-PromptReco-v6/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/MuEG_Run2011A-PromptReco-v6/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/SingleElectron_Run2011A-PromptReco-v6/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/SingleMu_Run2011A-PromptReco-v6/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/DoubleMu_Run2011A-PromptReco-v6/V04-02-36/vvfilter/*.root");

  dataSamples.push_back(dataset+"/DoubleElectron_Run2011B-PromptReco-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/MuEG_Run2011B-PromptReco-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/SingleElectron_Run2011B-PromptReco-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/SingleMu_Run2011B-PromptReco-v1/V04-02-36/vvfilter/*.root");
  dataSamples.push_back(dataset+"/DoubleMu_Run2011B-PromptReco-v1/V04-02-36/vvfilter/*.root");

  dataSamples.push_back("/nfs-6/userdata/jaehyeok/missing/DoubleMu_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
  dataSamples.push_back("/nfs-6/userdata/jaehyeok/missing/DoubleElectron_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
  dataSamples.push_back("/nfs-6/userdata/jaehyeok/missing/MuEG_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
  dataSamples.push_back("/nfs-6/userdata/jaehyeok/missing/SingleElectron_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");
  dataSamples.push_back("/nfs-6/userdata/jaehyeok/missing/SingleMu_Run2012B-PromptReco-v1/V05-02-28/vvfilter/*.root");

  if (runData)
    ProcessSample(dataSamples, SmurfTree::data, 5.1 -1, -1, false, true, "");
  
  return 0; 
}
