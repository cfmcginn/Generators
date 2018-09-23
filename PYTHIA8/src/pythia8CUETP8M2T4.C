//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TDatime.h"
#include "TFile.h"
#include "TTree.h"
#include "TNamed.h"
#include "TDirectoryFile.h"
#include "TRandom3.h"
#include "TMath.h"

//PYTHIA8 dependencies
#include "Pythia8/Pythia.h"

//FASTJET dependencies
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

//Local dependencies
#include "include/checkMakeDir.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

//TUNE LINKS:
//Common Block: https://github.com/cms-sw/cmssw/blob/master/Configuration/Generator/python/Pythia8CommonSettings_cfi.py
//CUETP8M1: https://github.com/cms-sw/cmssw/blob/master/Configuration/Generator/python/Pythia8CUEP8M1Settings_cfi.py
//CP5: https://github.com/cms-sw/cmssw/blob/master/Configuration/Generator/python/MCTunes2017/PythiaCP5Settings_cfi.py
//CUETP8M2T4: https://github.com/cms-sw/cmssw/blob/master/Configuration/Generator/python/Pythia8CUEP8M2T4Settings_cfi.py

int pythia8CUETP8M2T4(std::string outFileName, const std::string tuneStr, bool keepParticles, bool keepJets, bool keepWTA, const int nEvt = 10000, const bool doFlatPthat = false, const double flatPthatScale = 4.5, const double pthatMin = 15., const double pthatMax = 999999.)
{
  if(!keepParticles && !keepJets){
    std::cout << "Both keepParticles and keepJets are false. No point in running if nothing is saved! return 1" << std::endl;
    return 1;
  }

  const Int_t nTunes = 4;
  const std::string tunes[nTunes] = {"Vanilla", "CUETP8M1", "CP5", "CUETP8M2T4"};

  int tunePos = -1;
  
  for(Int_t tI = 0; tI < nTunes; ++tI){
    if(isStrSame(tuneStr, tunes[tI])){
      tunePos = tI;
      break;
    }
  }

  if(tunePos < 0){
    std::cout << "Given tuneStr, \'" << tuneStr << "\' is invalid. Please pick one of the following: " << std::endl;
    for(Int_t tI = 0; tI < nTunes; ++tI){
      std::cout << " " << tunes[tI] << ",";
    }
    std::cout << std::endl;
    std::cout << "Return 1" << std::endl;
    return 1;
  }

  TRandom3* randGen_p = new TRandom3(0);

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1,"");}
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");

  std::string descStr = tuneStr + "_NEvt" + std::to_string(nEvt) + "_DoFlatPthat";
  if(doFlatPthat) descStr = descStr + "True_FlatPower" + prettyString(flatPthatScale, 2, true);
  else descStr = descStr + "False";
  
  descStr = descStr + "_PthatMin" + prettyString(pthatMin, 1, true) + "_PthatMax" + prettyString(pthatMax, 1, true);
  
  outFileName = "output/" + dateStr + "/" + outFileName + "_" + descStr + "_" + dateStr + ".root";
  
  std::cout << "Creating " << outFileName << "..." << std::endl;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* genParticleTree_p = NULL; 
  TTree* ak4GenJetTree_ESchemeWTA_p = NULL;
  TTree* ak4GenJetTree_PureWTA_p = NULL;

  if(keepParticles) genParticleTree_p = new TTree("genParticleTree", "");
  if(keepJets){
    ak4GenJetTree_ESchemeWTA_p = new TTree("ak4GenJetTree_ESchemeWTA", "");
    if(keepWTA) ak4GenJetTree_PureWTA_p = new TTree("ak4GenJetTree_PureWTA", "");
  }
    
  TDirectoryFile* paramsDir_p = (TDirectoryFile*)outFile_p->mkdir("paramsDir");
  
  Float_t pthat_;
  Float_t weight_;

  //define particle params
  Int_t nPart_=0;
  std::vector<float>* genPt_p = new std::vector<float>;
  std::vector<float>* genPhi_p = new std::vector<float>;
  std::vector<float>* genEta_p = new std::vector<float>;
  std::vector<int>* genPDG_p = new std::vector<int>;
  std::vector<short>* genChg_p = new std::vector<short>;
  
  //define jet parameters
  const float minJtPtCut = 80.;
  const float maxJtAbsEtaCut = 5.1;

  const fastjet::JetAlgorithm jtAlgo = fastjet::antikt_algorithm;
  const double jtRVal = 0.4;
  const fastjet::RecombinationScheme jtRecombScheme = fastjet::E_scheme;
  const fastjet::JetDefinition jtDef(jtAlgo, jtRVal, jtRecombScheme);
  const fastjet::RecombinationScheme jtRecombScheme2 = fastjet::WTA_pt_scheme;
  const fastjet::JetDefinition jtDef2(jtAlgo, jtRVal+1., jtRecombScheme2);
  const fastjet::JetDefinition jtDef3(jtAlgo, jtRVal, jtRecombScheme2);
  
  Int_t nGenJt_ESchemeWTA_=0;
  std::vector<float>* genJtPt_ESchemeWTA_p = new std::vector<float>;
  std::vector<float>* toyRecoJtPt_ESchemeWTA_p = new std::vector<float>;
  std::vector<float>* genJtPhi_ESchemeWTA_p = new std::vector<float>;
  std::vector<float>* genJtEta_ESchemeWTA_p = new std::vector<float>;
  std::vector<float>* genJtPtWTA_ESchemeWTA_p = new std::vector<float>;
  std::vector<float>* genJtPhiWTA_ESchemeWTA_p = new std::vector<float>;
  std::vector<float>* genJtEtaWTA_ESchemeWTA_p = new std::vector<float>;

  Int_t nGenJt_PureWTA_=0;
  std::vector<float>* genJtPt_PureWTA_p = new std::vector<float>;
  std::vector<float>* genJtPhi_PureWTA_p = new std::vector<float>;
  std::vector<float>* genJtEta_PureWTA_p = new std::vector<float>;

  
  if(keepParticles){
    genParticleTree_p->Branch("pthat", &pthat_, "pthat/F");
    genParticleTree_p->Branch("weight", &weight_, "weight/F");  
    genParticleTree_p->Branch("nPart", &nPart_, "nPart/I");
    genParticleTree_p->Branch("pt", &genPt_p);
    genParticleTree_p->Branch("phi", &genPhi_p);
    genParticleTree_p->Branch("eta", &genEta_p);      
    genParticleTree_p->Branch("pdg", &genPDG_p);      
    genParticleTree_p->Branch("chg", &genChg_p);      
  }
  
  if(keepJets){
    if(!keepParticles){
      ak4GenJetTree_ESchemeWTA_p->Branch("pthat", &pthat_, "pthat/F");
      ak4GenJetTree_ESchemeWTA_p->Branch("weight", &weight_, "weight/F");  
    }
    
    ak4GenJetTree_ESchemeWTA_p->Branch("nGenJt", &nGenJt_ESchemeWTA_, "nGenJt/I");
    ak4GenJetTree_ESchemeWTA_p->Branch("genJtPt", &genJtPt_ESchemeWTA_p);
    ak4GenJetTree_ESchemeWTA_p->Branch("toyRecoJtPt", &toyRecoJtPt_ESchemeWTA_p);
    ak4GenJetTree_ESchemeWTA_p->Branch("genJtPhi", &genJtPhi_ESchemeWTA_p);
    ak4GenJetTree_ESchemeWTA_p->Branch("genJtEta", &genJtEta_ESchemeWTA_p);  
    if(keepWTA){
      ak4GenJetTree_ESchemeWTA_p->Branch("genJtPtWTA", &genJtPtWTA_ESchemeWTA_p);
      ak4GenJetTree_ESchemeWTA_p->Branch("genJtPhiWTA", &genJtPhiWTA_ESchemeWTA_p);
      ak4GenJetTree_ESchemeWTA_p->Branch("genJtEtaWTA", &genJtEtaWTA_ESchemeWTA_p);  
    
      ak4GenJetTree_PureWTA_p->Branch("nGenJt", &nGenJt_PureWTA_, "nGenJt/I");
      ak4GenJetTree_PureWTA_p->Branch("genJtPt", &genJtPt_PureWTA_p);
      ak4GenJetTree_PureWTA_p->Branch("genJtPhi", &genJtPhi_PureWTA_p);
      ak4GenJetTree_PureWTA_p->Branch("genJtEta", &genJtEta_PureWTA_p);
    }
  }
  
  Pythia8::Pythia pythia;
  pythia.readString("Beams:eCM = 5020.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");

  //Vanilla implement, shared

  pythia.readString("Tune:preferLHAPDF = 2");
  pythia.readString("Main:timesAllowErrors = 10000");
  pythia.readString("Check:epTolErr = 0.01");
  pythia.readString("Beams:setProductionScalesFromLHEF = off");
  pythia.readString("SLHA:keepSM = on");
  pythia.readString("SLHA:minMassSM = 1000.");
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10");
  pythia.readString("ParticleDecays:allowPhotonRadiation = on");

  if(tunePos == 1){//CUETP8M1
    pythia.readString("Tune:pp 14");
    pythia.readString("Tune:ee 7");
    pythia.readString("MultipartonInteractions:pT0Ref=2.4024");
    pythia.readString("MultipartonInteractions:ecmPow=0.25208");
    pythia.readString("MultipartonInteractions:expPow=1.6");
  }
  else if(tunePos == 2){//CP5
    pythia.readString("Tune:pp 14");
    pythia.readString("Tune:ee 7");
    pythia.readString("MultipartonInteractions:ecmPow=0.03344");
    pythia.readString("PDF:pSet=20");
    pythia.readString("MultipartonInteractions:bProfile=2");
    pythia.readString("MultipartonInteractions:pT0Ref=1.41");
    pythia.readString("MultipartonInteractions:coreRadius=0.7634");
    pythia.readString("MultipartonInteractions:coreFraction=0.63");
    pythia.readString("ColourReconnection:range=5.176");
    pythia.readString("SigmaTotal:zeroAXB=off");
    pythia.readString("SpaceShower:alphaSorder=2");
    pythia.readString("SpaceShower:alphaSvalue=0.118");
    pythia.readString("SigmaProcess:alphaSvalue=0.118");
    pythia.readString("SigmaProcess:alphaSorder=2");
    pythia.readString("MultipartonInteractions:alphaSvalue=0.118");
    pythia.readString("MultipartonInteractions:alphaSorder=2");
    pythia.readString("TimeShower:alphaSorder=2");
    pythia.readString("TimeShower:alphaSvalue=0.118");
  }
  else if(tunePos == 3){//CUETP8M2T4
    pythia.readString("Tune:pp 14");
    pythia.readString("Tune:ee 7");
    pythia.readString("MultipartonInteractions:ecmPow=0.25208");
    pythia.readString("SpaceShower:alphaSvalue=0.1108");
    pythia.readString("PDF:pSet=LHAPDF6:NNPDF30_lo_as_0130");
    pythia.readString("MultipartonInteractions:pT0Ref=2.20e+00");
    pythia.readString("MultipartonInteractions:expPow=1.60e+00");
    pythia.readString("ColourReconnection:range=6.59e+00");
  }
  
  if(doFlatPthat){
    pythia.readString("PhaseSpace:bias2Selection = on");
    pythia.readString("PhaseSpace:bias2SelectionPow = " + std::to_string(flatPthatScale));
    pythia.readString("PhaseSpace:bias2SelectionRef = " + std::to_string(pthatMin));
  }

  pythia.readString("PhaseSpace:pTHatMin = "  + std::to_string(pthatMin));
  pythia.readString("PhaseSpace:pTHatMax = "  + std::to_string(pthatMax));

  pythia.init();

  int totEvt = 0;
  
  while(totEvt < nEvt){
    if(!pythia.next()) continue;
    
    pthat_ = pythia.info.pTHat();
    if(doFlatPthat) weight_ = pythia.info.weight();
    else weight_ = 1.;

    //weight_ = 1.;

    std::vector<fastjet::PseudoJet> particles;
    
    for (int i = 0; i < pythia.event.size(); ++i){
      if(!pythia.event[i].isFinal()) continue;
      if(pythia.event[i].pT() < 0.5) continue; // assuming gev
      if(TMath::Abs(pythia.event[i].eta()) > maxJtAbsEtaCut) continue; // assuming rough detector geometry
      //continuing on neutrinos;
      if(TMath::Abs(pythia.event[i].id()) == 12) continue;
      if(TMath::Abs(pythia.event[i].id()) == 14) continue;
      if(TMath::Abs(pythia.event[i].id()) == 16) continue;

      if(keepParticles){
	genPt_p->push_back(pythia.event[i].pT());
	genPhi_p->push_back(pythia.event[i].phi());
	genEta_p->push_back(pythia.event[i].eta());
	genPDG_p->push_back(pythia.event[i].id());
	genChg_p->push_back(pythia.event[i].charge());
	++nPart_;
      }
      
      particles.push_back(fastjet::PseudoJet(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e()));
    }

    if(keepParticles){
      genParticleTree_p->Fill();

      nPart_ = 0;
      genPt_p->clear();
      genPhi_p->clear();
      genEta_p->clear();      
      genPDG_p->clear();      
      genChg_p->clear();      
    }
    
    if(keepJets){
      fastjet::ClusterSequence* cs = new fastjet::ClusterSequence(particles, jtDef);
      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs->inclusive_jets(minJtPtCut));
      
      for(unsigned int jI = 0; jI < jets.size(); ++jI){
	if(jets.at(jI).eta() > maxJtAbsEtaCut) continue;
	
	genJtPt_ESchemeWTA_p->push_back(jets.at(jI).pt());
	double sigma = TMath::Sqrt(0.06*0.06 + 0.9*0.9/jets.at(jI).pt());
	double tempToyRecoPt = jets.at(jI).pt()*randGen_p->Gaus(1.0, sigma);
	toyRecoJtPt_ESchemeWTA_p->push_back(tempToyRecoPt);
	
	
	genJtPhi_ESchemeWTA_p->push_back(jets.at(jI).phi_std());
	genJtEta_ESchemeWTA_p->push_back(jets.at(jI).eta());

	if(keepWTA){
	  std::vector<fastjet::PseudoJet> constituents = jets.at(jI).constituents();
	  
	  fastjet::ClusterSequence* cs2 = new fastjet::ClusterSequence(constituents, jtDef2);
	  std::vector<fastjet::PseudoJet> jets2 = fastjet::sorted_by_pt(cs->inclusive_jets());
	  
	  genJtPtWTA_ESchemeWTA_p->push_back(jets2.at(0).pt());
	  genJtPhiWTA_ESchemeWTA_p->push_back(jets2.at(0).phi_std());
	  genJtEtaWTA_ESchemeWTA_p->push_back(jets2.at(0).eta());
	
	  jets2.clear();
	  delete cs2;
	}
	
	++nGenJt_ESchemeWTA_;
      }

      jets.clear();
      delete cs;
    
      ak4GenJetTree_ESchemeWTA_p->Fill();

      nGenJt_ESchemeWTA_ = 0;
      genJtPt_ESchemeWTA_p->clear();
      toyRecoJtPt_ESchemeWTA_p->clear();
      genJtPhi_ESchemeWTA_p->clear();
      genJtEta_ESchemeWTA_p->clear();      
      genJtPtWTA_ESchemeWTA_p->clear();
      genJtPhiWTA_ESchemeWTA_p->clear();
      genJtEtaWTA_ESchemeWTA_p->clear();      

      if(keepWTA){
	cs = new fastjet::ClusterSequence(particles, jtDef3);
	jets = fastjet::sorted_by_pt(cs->inclusive_jets(minJtPtCut));
	
	for(unsigned int jI = 0; jI < jets.size(); ++jI){
	  if(jets.at(jI).eta() > maxJtAbsEtaCut) continue;
	  
	  genJtPt_PureWTA_p->push_back(jets.at(jI).pt());
	  genJtPhi_PureWTA_p->push_back(jets.at(jI).phi_std());
	  genJtEta_PureWTA_p->push_back(jets.at(jI).eta());	
	++nGenJt_PureWTA_;
	}

	jets.clear();
	delete cs;
      
	ak4GenJetTree_PureWTA_p->Fill();
      }
    }
    
    particles.clear();
    
    ++totEvt;
  }

  double sigmapb = pythia.info.sigmaGen() * 1.0E9;
  std::cout << "== Cross section for this run = " <<  sigmapb << " pb" << std::endl;
  
  outFile_p->cd();

  if(keepParticles) genParticleTree_p->Write("", TObject::kOverwrite);
  if(keepJets){
    ak4GenJetTree_ESchemeWTA_p->Write("", TObject::kOverwrite);
    if(keepWTA) ak4GenJetTree_PureWTA_p->Write("", TObject::kOverwrite);
  }

  paramsDir_p->cd();

  TNamed tuneName("tune", tuneStr.c_str());
  TNamed nEvtName("nEvt", std::to_string(nEvt).c_str());
  TNamed doFlatPthatName("doFlatPthat", std::to_string(doFlatPthat).c_str());
  TNamed flatPthatScaleName("flatPthatScale", prettyString(flatPthatScale, 2, false).c_str());
  TNamed pthatMinName("pthatMin", prettyString(pthatMin, 1, false).c_str());
  TNamed pthatMaxName("pthatMax", prettyString(pthatMax, 1, false).c_str());
  std::string crossSectionString = "crossSection_Pthat" + std::to_string((int)pthatMin);
  TNamed crossSectionName(crossSectionString.c_str(), prettyString(sigmapb, 3, false).c_str());
  
  tuneName.Write("", TObject::kOverwrite);
  nEvtName.Write("", TObject::kOverwrite);
  doFlatPthatName.Write("", TObject::kOverwrite);
  flatPthatScaleName.Write("", TObject::kOverwrite);
  pthatMinName.Write("", TObject::kOverwrite);
  pthatMaxName.Write("", TObject::kOverwrite);
  crossSectionName.Write("", TObject::kOverwrite);

  delete randGen_p;
  
  if(keepParticles) delete genParticleTree_p;
  if(keepJets){
    delete ak4GenJetTree_ESchemeWTA_p;
    if(keepWTA) delete ak4GenJetTree_PureWTA_p;
  }

  genPt_p->clear();
  genPhi_p->clear();
  genEta_p->clear();
  genPDG_p->clear();
  genChg_p->clear();

  delete genPt_p;
  delete genPhi_p;
  delete genEta_p;
  delete genPDG_p;
  delete genChg_p;

  genJtPt_ESchemeWTA_p->clear();
  toyRecoJtPt_ESchemeWTA_p->clear();
  genJtPhi_ESchemeWTA_p->clear();
  genJtEta_ESchemeWTA_p->clear();
  genJtPtWTA_ESchemeWTA_p->clear();
  genJtPhiWTA_ESchemeWTA_p->clear();
  genJtEtaWTA_ESchemeWTA_p->clear();

  delete genJtPt_ESchemeWTA_p;
  delete toyRecoJtPt_ESchemeWTA_p;
  delete genJtPhi_ESchemeWTA_p;
  delete genJtEta_ESchemeWTA_p;
  delete genJtPtWTA_ESchemeWTA_p;
  delete genJtPhiWTA_ESchemeWTA_p;
  delete genJtEtaWTA_ESchemeWTA_p;

  genJtPt_PureWTA_p->clear();
  genJtPhi_PureWTA_p->clear();
  genJtEta_PureWTA_p->clear();

  delete genJtPt_PureWTA_p;
  delete genJtPhi_PureWTA_p;
  delete genJtEta_PureWTA_p;

  outFile_p->Close();
  delete outFile_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 6 || argc > 11){
    std::cout << "Usage: ./bin/pythia8CUETP8M2T4.exe <outFileName> <tuneStr> <keepParticles> <keepJets> <keepWTA> <nEvt-opt> <doFlatPthat-opt> <flatPthatScale-opt> <pthatMin-opt> <pthatMax-opt>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 6) retVal += pythia8CUETP8M2T4(argv[1], argv[2], std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]));
  else if(argc == 7) retVal += pythia8CUETP8M2T4(argv[1], argv[2], std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]), std::stoi(argv[6]));
  else if(argc == 8) retVal += pythia8CUETP8M2T4(argv[1], argv[2], std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]), std::stoi(argv[6]), std::stoi(argv[7]));
  else if(argc == 9) retVal += pythia8CUETP8M2T4(argv[1], argv[2], std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]), std::stoi(argv[6]), std::stoi(argv[7]), std::stod(argv[8]));
  else if(argc == 10)retVal += pythia8CUETP8M2T4(argv[1], argv[2], std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]), std::stoi(argv[6]), std::stoi(argv[7]), std::stod(argv[8]), std::stod(argv[9]));
  else if(argc == 11) retVal += pythia8CUETP8M2T4(argv[1], argv[2], std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]), std::stoi(argv[6]), std::stoi(argv[7]), std::stod(argv[8]), std::stod(argv[9]), std::stod(argv[10]));
  return retVal;
}
