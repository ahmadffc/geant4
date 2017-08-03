#include "FakeModel.hh"


FakeModel::FakeModel(const std::string file_name) :
  fChain(NULL)
{
  TFile *f = new TFile(file_name.c_str());
  TTree *tree = NULL;
  f->GetObject("blob", tree);
  Init(tree);
}

FakeModel::~FakeModel()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

void FakeModel::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("b", &b, &b_b);
   fChain->SetBranchAddress("A", &A, &b_A);
   fChain->SetBranchAddress("Z", &Z, &b_Z);
   fChain->SetBranchAddress("Eecc", &Eecc, &b_Eecc);
   fChain->SetBranchAddress("spinx", &spinx, &b_spinx);
   fChain->SetBranchAddress("spiny", &spiny, &b_spiny);
   fChain->SetBranchAddress("spinz", &spinz, &b_spinz);
   fChain->SetBranchAddress("px", &px, &b_px);
   fChain->SetBranchAddress("py", &py, &b_py);
   fChain->SetBranchAddress("pz", &pz, &b_pz);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("v", &v, &b_v);
   fChain->SetBranchAddress("theta", &theta, &b_theta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("EK", &EK, &b_EK);

}

Int_t FakeModel::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Int_t FakeModel::NextEntry()
{
  fCurrent++;
  return GetEntry(fCurrent);
}
