//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"

//Non-local (Utility) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/stringUtil.h"

int forestToGen(const std::string inFileName)
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

  outFile_p->cd();

  genTree_p->Write("", TObject::kOverwrite);
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
  if(argc != 2){
    std::cout << "Usage: ./bin/forestToGen.exe <inFileName>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += forestToGen(argv[1]);
  return retVal;
}
