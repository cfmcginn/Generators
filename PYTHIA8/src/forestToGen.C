//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TDirectoryFile.h"

//Non-local (Utility) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/stringUtil.h"

const Int_t nPyt6PtHats = 11;
const Float_t pyt6PtHats[nPyt6PtHats] = {30, 50, 80, 100, 120, 170, 220, 280, 370, 460, 540};
//Taken from table here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiForest2015#PYTHIA_6_tune_Z2
//Note that 540 xsection based on PYTHIA8 value, renormalized to PYTHIA6 280, 370, 460, and averaged
const Double_t pyt6XSect[nPyt6PtHats] = {.03378,  .003778,  .0004412,  .0001511,  .00006147,  .00001018,  .000002477,  .0000006160,  .0000001088,  .00000002527, .000000007850};
const Int_t nPyt8PtHats = 12;
const Float_t pyt8PtHats[nPyt8PtHats] = {15, 30, 50, 80, 100, 120, 170, 220, 280, 370, 460, 540};
//Taken from table here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiForest2015#PYTHIA_8_tune_CUETP8M1
//Note that all values are taken from column 1 except for pthat100
const Double_t pyt8XSect[nPyt8PtHats] = {.5269, .03455, .004068, .0004959, .000173, .00007096, .00001223, .000003031, .0000007746, .0000001410, .00000003216, .00000001001};

int forestToGen(const std::string inFileName, bool isPyt6)
{
  if(!checkFile(inFileName)){
    std::cout << "Given inFileName \'" << inFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  std::vector<std::string> fileList;
  if(inFileName.find(".root") != std::string::npos) fileList.push_back(inFileName);
  else if(inFileName.find(".txt") != std::string::npos){

  }

  if(fileList.size() == 0){
    std::cout << "Given inFileName \'" << inFileName << "\' produces no valid root files. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;
  
  std::string jetTreeString = "";
  
  //Lets ID the jet tree and make sure it is present in all files
  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    const std::string fName = fileList.at(fI);
    TFile* inFile_p = new TFile(fName.c_str(), "READ");
    std::vector<std::string> listOfJetTrees = returnRootFileContentsList(inFile_p, "TTree", "JetAnalyzer");

    if(fI == 0){
      std::string tempJetTreeString = "";
      for(unsigned int jI = 0; jI < listOfJetTrees.size(); ++jI){
	if(listOfJetTrees.at(jI).find("ak4") != std::string::npos) tempJetTreeString = listOfJetTrees.at(jI);
	else if(listOfJetTrees.at(jI).find("akPu4") != std::string::npos) tempJetTreeString = listOfJetTrees.at(jI);
	else if(listOfJetTrees.at(jI).find("akCs4") != std::string::npos) tempJetTreeString = listOfJetTrees.at(jI);
	
	if(tempJetTreeString.size() != 0) break;
      }

      if(tempJetTreeString.size() == 0){
	std::cout << "First file \'" << fName << "\' contains no R=0.4 TTree. return 1" << std::endl;
	std::cout << " Jet trees are: ";
	for(unsigned int jI = 0; jI < listOfJetTrees.size(); ++jI){
	  std::cout << listOfJetTrees.at(jI) << ", ";
	}
	std::cout << std::endl;
	
	inFile_p->Close();
	delete inFile_p;
	return 1;
      }
      else jetTreeString = tempJetTreeString;
    }
    else{
      bool hasJetTreeString = false;
      for(unsigned int jI = 0; jI < listOfJetTrees.size(); ++jI){
	if(isStrSame(listOfJetTrees.at(jI), jetTreeString)){
	  hasJetTreeString = true;
	  break;
	}
      }

      if(!hasJetTreeString){
	std::cout << "File \'" << fName << "\' contains no R=0.4 TTree. return 1" << std::endl;
	std::cout << " Jet trees are: ";
	for(unsigned int jI = 0; jI < listOfJetTrees.size(); ++jI){
	  std::cout << listOfJetTrees.at(jI) << ", ";
	}
	std::cout << std::endl;
	
	inFile_p->Close();
	delete inFile_p;
	return 1;	
      }
    }
    
    inFile_p->Close();
    delete inFile_p;
  }

  std::cout << "Using jet tree \'" << jetTreeString << "\'..." << std::endl;

  std::string outFileName = inFileName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  if(outFileName.find(".") != std::string::npos) outFileName.replace(outFileName.find("."), outFileName.size(), "");

  checkMakeDir("output");
  outFileName = "output/" + outFileName + "_ForestToGen_" + dateStr + ".root";
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* genTree_p = new TTree("genTree", "");
  TDirectoryFile* paramsDir_p = (TDirectoryFile*)outFile_p->mkdir("paramsDir");
  Int_t nEvt = 0;
  Float_t pthatMin = 100000;

  Float_t pthat_;
  Float_t weight_;
  const Int_t nMaxJet = 500;
  Int_t nGenInJt_;
  Float_t genInJtPt_[nMaxJet];
  Float_t genInJtPhi_[nMaxJet];
  Float_t genInJtEta_[nMaxJet];

  Int_t nGenJt_;
  std::vector<float>* genJtPt_p = new std::vector<float>;
  std::vector<float>* genJtPhi_p = new std::vector<float>;
  std::vector<float>* genJtEta_p = new std::vector<float>;

  genTree_p->Branch("pthat", &pthat_, "pthat/F");
  genTree_p->Branch("weight", &weight_, "weight/F");
  genTree_p->Branch("nGenJt", &nGenJt_, "nGenJt/I");
  genTree_p->Branch("genJtPt", &genJtPt_p);
  genTree_p->Branch("genJtPhi", &genJtPhi_p);
  genTree_p->Branch("genJtEta", &genJtEta_p);
  
  weight_ = 1.;

  
  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    const std::string fName = fileList.at(fI);
    TFile* inFile_p = new TFile(fName.c_str(), "READ");
    TTree* inTree_p = (TTree*)inFile_p->Get(jetTreeString.c_str());
    inTree_p->SetBranchStatus("*", 0);
    inTree_p->SetBranchStatus("pthat", 1);
    inTree_p->SetBranchStatus("ngen", 1);
    inTree_p->SetBranchStatus("genpt", 1);
    inTree_p->SetBranchStatus("genphi", 1);
    inTree_p->SetBranchStatus("geneta", 1);

    inTree_p->SetBranchAddress("pthat", &pthat_);
    inTree_p->SetBranchAddress("ngen", &nGenInJt_);
    inTree_p->SetBranchAddress("genpt", genInJtPt_);
    inTree_p->SetBranchAddress("genphi", genInJtPhi_);
    inTree_p->SetBranchAddress("geneta", genInJtEta_);

    const Int_t nEntries = inTree_p->GetEntries();
    for(Int_t entry = 0; entry < nEntries; ++entry){
      inTree_p->GetEntry(entry);

      ++nEvt;
      if(pthat_ < pthatMin) pthatMin = pthat_;

      nGenJt_ = 0;
      genJtPt_p->clear();
      genJtPhi_p->clear();
      genJtEta_p->clear();

      for(Int_t jI = 0; jI < nGenInJt_; ++jI){
	genJtPt_p->push_back(genInJtPt_[jI]);
	genJtPhi_p->push_back(genInJtPhi_[jI]);
	genJtEta_p->push_back(genInJtEta_[jI]);
	++nGenJt_;
      }
      
      genTree_p->Fill();
    }
    
    inFile_p->Close();
    delete inFile_p;
  }

  std::cout << "pTHat Min " << pthatMin << std::endl;
  Float_t pthatOut = -1;
  Double_t xSecOut = -1;
  if(isPyt6){
    std::cout << "Searching PYTHIA6 pthats: " << std::endl;
    Int_t pos = -1;
    for(Int_t pI = 0; pI < nPyt6PtHats; ++pI){
      if(TMath::Abs(pthatMin - pyt6PtHats[pI]) < 1){
	std::cout << " Matched to " << pyt6PtHats[pI] << std::endl;
	pos = pI;
	break;
      }
    }
    
    if(pos >= 0){
      pthatOut = pyt6PtHats[pos];
      xSecOut = pyt6XSect[pos];
    }
  }
  else{
    std::cout << "Searching PYTHIA8 pthats: " << std::endl;
    Int_t pos = -1;
    for(Int_t pI = 0; pI < nPyt8PtHats; ++pI){
      if(TMath::Abs(pthatMin - pyt8PtHats[pI]) < 1){
	std::cout << " Matched to " << pyt8PtHats[pI] << std::endl;
	pos = pI;
	break;
      }
    }
    
    if(pos >= 0){
      pthatOut = pyt8PtHats[pos];
      xSecOut = pyt8XSect[pos];
    }
  }


  outFile_p->cd();

  genTree_p->Write("", TObject::kOverwrite);

  paramsDir_p->cd();
  std::string tuneStr = "";
  if(isPyt6) tuneStr = "PYT6Z2";
  else tuneStr = "CUETP8M1";
  //Cross sections above are given in mb, convert to pb
  Double_t multFact = 1;
  if(xSecOut >= 0) multFact = 1000000000;//Just so -1 doesnt get a crazy factor
  const Double_t sigmapb = xSecOut*multFact;
  TNamed tuneName("tune", tuneStr.c_str());
  TNamed nEvtName("nEvt", std::to_string(nEvt).c_str());
  //If forest, never flat, give 0 + default
  TNamed doFlatPthatName("doFlatPthat", std::to_string(0).c_str());
  TNamed flatPthatScaleName("flatPthatScale", prettyString(4.5, 2, false).c_str());
  TNamed pthatMinName("pthatMin", prettyString(pthatOut, 1, false).c_str());
  TNamed pthatMaxName("pthatMax", prettyString(999999., 1, false).c_str());
  std::string crossSectionString = "crossSection_Pthat" + std::to_string((int)pthatMin);
  TNamed crossSectionName(crossSectionString.c_str(), prettyString(sigmapb, 3, false).c_str());
  
  tuneName.Write("", TObject::kOverwrite);
  nEvtName.Write("", TObject::kOverwrite);
  doFlatPthatName.Write("", TObject::kOverwrite);
  flatPthatScaleName.Write("", TObject::kOverwrite);
  pthatMinName.Write("", TObject::kOverwrite);
  pthatMaxName.Write("", TObject::kOverwrite);
  crossSectionName.Write("", TObject::kOverwrite);

  delete genTree_p;

  delete genJtPt_p;
  delete genJtPhi_p;
  delete genJtEta_p;
  
  outFile_p->Close();
  delete outFile_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/forestToGen.exe <inFileName> <isPyt6>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += forestToGen(argv[1], std::stoi(argv[2]));
  return retVal;
}
