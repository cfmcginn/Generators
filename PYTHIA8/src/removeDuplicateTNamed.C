//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TDirectoryFile.h"
#include "TH1D.h"

//Non-local (Utility) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/returnRootFileContentsList.h"

int recursiveRemove(TDirectory* inDir_p, TFile* outFile_p, std::vector<TDirectory*>* outDir_p, int outDirPos)
{
  TIter next(inDir_p->GetListOfKeys());
  TKey* key = NULL;

  inDir_p->cd();
  while( (key = (TKey*)next() ) ){
    const std::string name = key->GetName();
    const std::string className = key->GetClassName();

    if(className.find("TTree") != std::string::npos){
      inDir_p->cd();
      TTree* tree_p = (TTree*)key->ReadObj();

      outFile_p->cd();
      outDir_p->at(outDirPos)->cd();
      tree_p->CloneTree()->Write("", TObject::kOverwrite);
    } 
    else if(className.find("TNamed") != std::string::npos){
      inDir_p->cd();
      TNamed* tnamed_p = (TNamed*)key->ReadObj();

      outFile_p->cd();
      outDir_p->at(outDirPos)->cd();
      tnamed_p->Write("", TObject::kOverwrite);
    } 
    else if(className.find("TDirectory") != std::string::npos){
      inDir_p->cd();
      TDirectory* dir_p = (TDirectory*)key->ReadObj();
      outDir_p->push_back(outDir_p->at(outDirPos)->mkdir(name.c_str()));

      int pos = outDir_p->size() - 1;
      recursiveRemove(dir_p, outFile_p, outDir_p, pos);
    } 

    inDir_p->cd();   
  }  
  
  return 0;
}

int removeDuplicateTNamed(const std::string inFileName, std::string outFileName = "")
{
  if(!checkFile(inFileName) || inFileName.find(".root") == std::string::npos){
    std::cout << "Given inFileName \'" << inFileName << "\' is invaled. return 1" << std::endl;  
    return 1;
  }

  if(outFileName.size() == 0){
    outFileName = inFileName.substr(0, inFileName.find(".root"));
    outFileName = outFileName + "_RemovedDupTNamed.root";
  }

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  outFile_p->SetBit(TFile::kDevNull);
  TH1::AddDirectory(kFALSE);
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ"); 
  std::vector<std::string> dirList = returnRootFileContentsList(outFile_p, "TDirectory", "");
  std::vector<std::string> dirList2 = returnRootFileContentsList(outFile_p, "TDirectoryFile", "");

  dirList.insert(dirList.end(), dirList.begin(), dirList.end());
  const Int_t nDir = dirList.size();
  std::vector<TDirectory*> dirs_p;
  dirs_p.reserve(nDir);
  
  TIter next(inFile_p->GetListOfKeys());
  TKey* key = NULL;
  
  while( (key = (TKey*)next() ) ){
    const std::string name = key->GetName();
    const std::string className = key->GetClassName();

    if(className.find("TTree") != std::string::npos){
      inFile_p->cd();
      TTree* tree_p = (TTree*)key->ReadObj();

      outFile_p->cd();     
      tree_p->CloneTree()->Write("", TObject::kOverwrite);
    } 
    else if(className.find("TNamed") != std::string::npos){
      inFile_p->cd();
      TNamed* tnamed_p = (TNamed*)key->ReadObj();

      outFile_p->cd();
      tnamed_p->Write("", TObject::kOverwrite);
    } 
    else if(className.find("TDirectory") != std::string::npos){
      inFile_p->cd();
      TDirectory* dir_p = (TDirectory*)key->ReadObj();
      dirs_p.push_back(outFile_p->mkdir(name.c_str()));

      int pos = dirs_p.size() - 1;
      recursiveRemove(dir_p, outFile_p, &(dirs_p), pos);
    }

    inFile_p->cd();
  }  

  outFile_p->cd();
  
  outFile_p->Close();
  delete outFile_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 2 || argc > 3){
    std::cout << "Usage: ./bin/removeDuplicateTNamed.exe <inFileName> <outFileName-opt>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += removeDuplicateTNamed(argv[1]);
  else if(argc == 3) retVal += removeDuplicateTNamed(argv[1], argv[2]);
  return retVal;
}
