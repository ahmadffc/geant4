#ifndef FakeModel_hh
#define FakeModel_hh

#include <vector>
#include <fstream>
#include <iostream>
#include <string>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class FakeModel {
public:
  FakeModel(const std::string file_name);
  ~FakeModel();

   size_t GetFirstFragmentIdInEvent();
   size_t GetLastFragmentIdInEvent();

  virtual int    GetEntry(Long64_t entry);
  // virtual int    NextEntry();
  
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  
  // Declaration of leaf types
  int           evt;
  float         b;
  int           A;
  int           Z;
  float         Eecc;
  float         spinx;
  float         spiny;
  float         spinz;
  float         px;
  float         py;
  float         pz;
  float         x;
  float         y;
  float         z;
  float         v;
  float         theta;
  float         phi;
  float         EK;

private:
  // List of branches
  TBranch        *b_evt;   //!
  TBranch        *b_b;   //!
  TBranch        *b_A;   //!
  TBranch        *b_Z;   //!
  TBranch        *b_Eecc;   //!
  TBranch        *b_spinx;   //!
  TBranch        *b_spiny;   //!
  TBranch        *b_spinz;   //!
  TBranch        *b_px;   //!
  TBranch        *b_py;   //!
  TBranch        *b_pz;   //!
  TBranch        *b_x;   //!
  TBranch        *b_y;   //!
  TBranch        *b_z;   //!
  TBranch        *b_v;   //!
  TBranch        *b_theta;   //!
  TBranch        *b_phi;   //!
  TBranch        *b_EK;   //!

  int           fCurrent; 
  virtual void     Init(TTree *tree);
  Long64_t fNentries;  
};
#endif
