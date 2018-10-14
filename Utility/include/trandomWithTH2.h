#ifndef TRANDOMWITHTH2_H
#define TRANDOMWITHTH2_H

//cpp depdencies
#include <iostream>
#include <vector>
#include <string>

//ROOT dependencies
#include "TRandom3.h"
#include "TH2.h"

class trandomWithTH2
{
 public:
  TRandom3* randGen_p=NULL;
  TH2D* internalTH2_p=NULL;//X-axis is reco, y-axis is gen
  Int_t nBinsX = -1;
  Int_t nBinsY = -1;

  std::vector<double> xBinBounds;
  std::vector<double> yBinBounds;
  std::vector<std::vector<double> > binContentsNormed;
  std::vector<double> yBinIntegrals;
  
  bool isInit = false;
  
  trandomWithTH2(){isInit = false; return;};
  trandomWithTH2(TH2* inputTH2_p);
  void Init(TH2* inputTH2_p);
  void Clean();
  void Print();
  Double_t GetRecoPtFromMatrix(Double_t inGenPt); //Returns bin center -> this is tied to internal binning! can only be used w/ internally consistent bins
};

trandomWithTH2::trandomWithTH2(TH2* inputTH2_p)
{
  Init(inputTH2_p);
  return;
}

void trandomWithTH2::Init(TH2* inputTH2_p)
{
  if(isInit){
    std::cout << "WARNING: trandomWithTH2 is already INIT - if you want to INIT, CLEAN firest. return" << std::endl;
    return;
  }
  
  internalTH2_p = (TH2D*)inputTH2_p->Clone((std::string(inputTH2_p->GetName()) + "_internal").c_str());
  randGen_p = new TRandom3(0);

  nBinsX = internalTH2_p->GetXaxis()->GetNbins();
  nBinsY = internalTH2_p->GetYaxis()->GetNbins();

  for(Int_t bIX = 0; bIX < nBinsX+1; ++bIX){
    xBinBounds.push_back(internalTH2_p->GetXaxis()->GetBinLowEdge(bIX+1));
  }
  for(Int_t bIY = 0; bIY < nBinsY+1; ++bIY){
    yBinBounds.push_back(internalTH2_p->GetYaxis()->GetBinLowEdge(bIY+1));

    if(bIY < nBinsY+1){
      binContentsNormed.push_back({});
      yBinIntegrals.push_back(0.0);

      for(Int_t bIX = 0; bIX < nBinsX; ++bIX){
	binContentsNormed.at(bIY).push_back(internalTH2_p->GetBinContent(bIX+1, bIY+1));
	yBinIntegrals.at(bIY) += internalTH2_p->GetBinContent(bIX+1, bIY+1);
      }

      //Cumulative
      for(Int_t bIX = ((Int_t)binContentsNormed.at(bIY).size())-1; bIX >= 0; --bIX){
	for(Int_t bIX2 = 0; bIX2 < bIX; ++bIX2){
	  binContentsNormed.at(bIY).at(bIX) += binContentsNormed.at(bIY).at(bIX2);
	}
      }

      for(Int_t bIX = 0; bIX < nBinsX; ++bIX){
	binContentsNormed.at(bIY).at(bIX) /= yBinIntegrals.at(bIY);
      }
    }
  }

  isInit = true;
  
  return;
}

void trandomWithTH2::Clean()
{
  if(isInit){
    delete internalTH2_p;
    delete randGen_p;

    nBinsX = -1;
    nBinsY = -1;

    xBinBounds.clear();
    yBinBounds.clear();
    binContentsNormed.clear();
    yBinIntegrals.clear();
  }

  internalTH2_p = NULL;
  randGen_p = NULL;
  
  isInit = false;
  return;
}

void trandomWithTH2::Print()
{
  if(!isInit){
    std::cout << "WARNING: trandomWithTH2 PRINT called w/o INIT. return" << std::endl;
    return;
  }
  
  std::cout << "BinsX: ";
  for(Int_t bIX = 0; bIX < nBinsX; ++bIX){
    std::cout << xBinBounds.at(bIX) << ",";
  }
  std::cout << xBinBounds.at(nBinsX) << "." << std::endl;

  std::cout << "BinsY: ";
  for(Int_t bIY = 0; bIY < nBinsY; ++bIY){
    std::cout << yBinBounds.at(bIY) << ",";
  }
  std::cout << yBinBounds.at(nBinsY) << "." << std::endl;

  std::cout << "yBinIntegrals, normed contents: ";
  for(Int_t bIY = 0; bIY < nBinsY; ++bIY){
    std::cout << "Integral for " << yBinBounds.at(bIY) << "-" << yBinBounds.at(bIY+1) << ": " << yBinIntegrals.at(bIY) << std::endl;

    std::cout << " Normed contents: ";
    for(unsigned int bIX = 0; bIX < binContentsNormed.at(bIY).size()-1; ++bIX){
      std::cout << binContentsNormed.at(bIY).at(bIX) << ", ";
    }
    std::cout << binContentsNormed.at(bIY).at(binContentsNormed.at(bIY).size()-1) << "." << std::endl;
  }

  return;
}


Double_t trandomWithTH2::GetRecoPtFromMatrix(Double_t inGenPt)
{
  if(!isInit){
    std::cout << "WARNING: trandomWithTH2 is not INIT before GETRECOPTFROMMATRIX called. return -1" << std::endl;
    return -1;
  }

  if(inGenPt < internalTH2_p->GetYaxis()->GetBinLowEdge(1)) return -1;
  if(inGenPt >= internalTH2_p->GetYaxis()->GetBinLowEdge(nBinsY+1)) return -1;
  
  Int_t binPosY = -1;
  for(Int_t bIY = 0; bIY < nBinsY; ++bIY){
    if(inGenPt >= internalTH2_p->GetYaxis()->GetBinLowEdge(bIY+1) && inGenPt < internalTH2_p->GetYaxis()->GetBinLowEdge(bIY+2)){
      binPosY = bIY;
      break;
    }
  }

  std::vector<double> tempNorm = binContentsNormed.at(binPosY);
  Double_t randUni = randGen_p->Uniform(0.0, 1.0);

  if(randUni >= 0 && randUni < tempNorm.at(0)) return (xBinBounds.at(0) + xBinBounds.at(1))/2.;
  for(Int_t uI = 0; uI < tempNorm.size()-1; ++uI){
    if(randUni >= tempNorm.at(uI) && randUni < tempNorm.at(uI+1)){
      return (xBinBounds.at(uI+1) + xBinBounds.at(uI+2))/2.;
    }
  }
  if(randUni == 1.0) return (xBinBounds.at(nBinsX-1) + xBinBounds.at(nBinsX))/2.;
  
  return -1;
}

#endif
