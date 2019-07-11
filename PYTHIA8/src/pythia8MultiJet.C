//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TDatime.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TNamed.h"
#include "TRandom3.h"
#include "TTree.h"

//PYTHIA8 dependencies
#include "Pythia8/Pythia.h"

//FASTJET dependencies
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

//Non-local (Utility) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/cppWatch.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/stringUtil.h"

//TUNE LINKS:
//Common Block: https://github.com/cms-sw/cmssw/blob/master/Configuration/Generator/python/Pythia8CommonSettings_cfi.py
//CUETP8M1: https://github.com/cms-sw/cmssw/blob/master/Configuration/Generator/python/Pythia8CUEP8M1Settings_cfi.py
//CP5: https://github.com/cms-sw/cmssw/blob/master/Configuration/Generator/python/MCTunes2017/PythiaCP5Settings_cfi.py
//CUETP8M2T4: https://github.com/cms-sw/cmssw/blob/master/Configuration/Generator/python/Pythia8CUEP8M2T4Settings_cfi.py

int pythia8MultiJet(const Int_t nEvt, const Int_t nSeed)
{  
  cppWatch timer;
  timer.start();

  TRandom3* randGen_p = new TRandom3(0);

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  
  std::string outFileName = "output/" + dateStr + "/multijet_Seed" + std::to_string(nSeed) + ".root";
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* ak4GenJetTree_p = new TTree("ak4GenJetTree", "");

  const Int_t nMaxPartons = 3;
  Int_t nPartons_ = nMaxPartons;
  Float_t partonPt_[nMaxPartons];
  Float_t partonPhi_[nMaxPartons];
  Float_t partonEta_[nMaxPartons];
  Int_t partonID_[nMaxPartons];
  
  Float_t pthat_;
  Float_t weight_;

  //define jet parameters
  const float minJtPtCut = 15.;
  const float maxJtAbsEtaCut = 5.1;

  const fastjet::JetAlgorithm jtAlgo = fastjet::antikt_algorithm;
  const double jtRVal = 0.4;
  const fastjet::RecombinationScheme jtRecombScheme = fastjet::E_scheme;
  const fastjet::JetDefinition jtDef(jtAlgo, jtRVal, jtRecombScheme);
  
  Int_t nGenJt_=0;
  std::vector<float>* genJtPt_p = new std::vector<float>;
  std::vector<float>* genJtPhi_p = new std::vector<float>;
  std::vector<float>* genJtEta_p = new std::vector<float>;

  ak4GenJetTree_p->Branch("nPartons", &nPartons_, "nPartons/I");
  ak4GenJetTree_p->Branch("partonPt", partonPt_, "partonPt[nPartons]/F");
  ak4GenJetTree_p->Branch("partonPhi", partonPhi_, "partonPhi[nPartons]/F");
  ak4GenJetTree_p->Branch("partonEta", partonEta_, "partonEta[nPartons]/F");
  ak4GenJetTree_p->Branch("partonID", partonID_, "partonID[nPartons]/I");
  
  ak4GenJetTree_p->Branch("pthat", &pthat_, "pthat/F");
  ak4GenJetTree_p->Branch("weight", &weight_, "weight/F");  


  ak4GenJetTree_p->Branch("nGenJt", &nGenJt_, "nGenJt/I");
  ak4GenJetTree_p->Branch("genJtPt", &genJtPt_p);
  ak4GenJetTree_p->Branch("genJtPhi", &genJtPhi_p);
  ak4GenJetTree_p->Branch("genJtEta", &genJtEta_p);  
  
  Pythia8::Pythia pythia;
  pythia.readString("Beams:eCM = 5020.");
  pythia.readString("HardQCD:all = off");
  pythia.readString("HardQCD:qq2qqgDiff = on");
  pythia.readString("HardQCD:qq2qqgSame = on");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(nSeed));
  pythia.readString("PhaseSpace:pTHatMin = 80.");
  //  pythia.readString("PhaseSpace:pTHatMax = 1000.");
  //  pythia.readString("PhaseSpace:bias2Selection = on");
  //  pythia.readString("PhaseSpace:bias2SelectionPow = 4.");
  //  pythia.readString("PhaseSpace:bias2SelectionRef = 80.");

  //Vanilla implement, shared

//  pythia.readString("Tune:preferLHAPDF = 2");
//  pythia.readString("Main:timesAllowErrors = 10000");
//  pythia.readString("Check:epTolErr = 0.01");
//  pythia.readString("Beams:setProductionScalesFromLHEF = off");
//  pythia.readString("SLHA:keepSM = on");
//  pythia.readString("SLHA:minMassSM = 1000.");
//  pythia.readString("ParticleDecays:limitTau0 = on");
//  pythia.readString("ParticleDecays:tau0Max = 10");
//  pythia.readString("ParticleDecays:allowPhotonRadiation = on");

  //CUETP8M1
//  pythia.readString("Tune:pp 14");
//  pythia.readString("Tune:ee 7");
//  pythia.readString("MultipartonInteractions:pT0Ref=2.4024");
//  pythia.readString("MultipartonInteractions:ecmPow=0.25208");
//  pythia.readString("MultipartonInteractions:expPow=1.6");
  
  pythia.init();

  int totEvt = 0;

  int nDiv = TMath::Max(nEvt/20, 1);
  
  while(totEvt < nEvt){
    if(!pythia.next()) continue;
    if(pythia.event[5].pT() < 60.) continue;
    if(pythia.event[5].id() == 21) continue;
    
    pthat_ = pythia.info.pTHat();
    weight_ = pythia.info.weight();
    //    weight_ = 1.;

    std::vector<fastjet::PseudoJet> particles;

    TLorentzVector sum123(0, 0, 0, 0);
    //    std::cout << "First 10 particles (" << pthat_ << ")" << std::endl;
    for (int i = 0; i < pythia.event.size(); ++i){
      if(i < 10){
	//	std::cout << " " << pythia.event[i].id() << ", " << pythia.event[i].pT() << ", " << pythia.event[i].phi() << ", " << pythia.event[i].eta() << std::endl;

	TLorentzVector tL(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
	if(i >= 5 && i <= 7){
	  sum123 += tL;

	  partonPt_[i-5] = pythia.event[i].pT();
	  partonPhi_[i-5] = pythia.event[i].phi();
	  partonEta_[i-5] = pythia.event[i].eta();
	  partonID_[i-5] = pythia.event[i].id();
	}
      
      }
      if(i == 10){
	//	std::cout << "SUM: " <<  sum123.Px() << ", " << sum123.Py() << ", " << sum123.Pz() << std::endl;
	//	std::cout << std::endl;
      }
      
      if(!pythia.event[i].isFinal()) continue;
      if(pythia.event[i].pT() < 0.5) continue; // assuming gev
      if(TMath::Abs(pythia.event[i].eta()) > maxJtAbsEtaCut) continue; // assuming rough detector geometry
      //continuing on neutrinos;
      if(TMath::Abs(pythia.event[i].id()) == 12) continue;
      if(TMath::Abs(pythia.event[i].id()) == 14) continue;
      if(TMath::Abs(pythia.event[i].id()) == 16) continue;
            
      particles.push_back(fastjet::PseudoJet(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e()));
    }

    fastjet::ClusterSequence* cs = new fastjet::ClusterSequence(particles, jtDef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs->inclusive_jets(minJtPtCut));
      
    for(unsigned int jI = 0; jI < jets.size(); ++jI){
      if(jets.at(jI).eta() > maxJtAbsEtaCut) continue;
	
      genJtPt_p->push_back(jets.at(jI).pt());
      genJtPhi_p->push_back(jets.at(jI).phi_std());
      genJtEta_p->push_back(jets.at(jI).eta());
      ++nGenJt_;
    }
    
    jets.clear();
    delete cs;
    
    ak4GenJetTree_p->Fill();

    nGenJt_ = 0;
    genJtPt_p->clear();
    genJtPhi_p->clear();
    genJtEta_p->clear();          
    particles.clear();
    
    ++totEvt;
    if(totEvt%nDiv == 0) std::cout << " Generated " << totEvt << "/" << nEvt << std::endl;
  }
  
  outFile_p->cd();

  ak4GenJetTree_p->Write("", TObject::kOverwrite);

  delete randGen_p;
  
  delete ak4GenJetTree_p;


  genJtPt_p->clear();
  genJtPhi_p->clear();
  genJtEta_p->clear();

  delete genJtPt_p;
  delete genJtPhi_p;
  delete genJtEta_p;

  outFile_p->Close();
  delete outFile_p;

  timer.stop();

  std::cout << "TIMER: " << timer.totalWall() << ", " << timer.totalCPU() << std::endl;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/pythia8MultiJet.exe <nEvt> <nSeed>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += pythia8MultiJet(std::stoi(argv[1]), std::stoi(argv[2]));
  return retVal;
}
