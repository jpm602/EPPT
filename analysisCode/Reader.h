//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb  1 19:24:20 2021 by ROOT version 6.14/04
// from TTree B5/Hits
// found on file: B5.root
//////////////////////////////////////////////////////////

#ifndef Reader_h
#define Reader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class Reader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Dc1Hits;
   Int_t           Dc2Hits;
   Double_t        ECEnergy;
   Double_t        HCEnergy;
   Double_t        Time1;
   Double_t        Time2;
   vector<double>  *ECEnergyVector;
   vector<double>  *HCEnergyVector;
   vector<double>  *Dc1HitsVector_x;
   vector<double>  *Dc1HitsVector_y;
   vector<double>  *Dc1HitsVector_z;
   vector<double>  *Dc2HitsVector_x;
   vector<double>  *Dc2HitsVector_y;
   vector<double>  *Dc2HitsVector_z;

   // List of branches
   TBranch        *b_Dc1Hits;   //!
   TBranch        *b_Dc2Hits;   //!
   TBranch        *b_ECEnergy;   //!
   TBranch        *b_HCEnergy;   //!
   TBranch        *b_Time1;   //!
   TBranch        *b_Time2;   //!
   TBranch        *b_ECEnergyVector;   //!
   TBranch        *b_HCEnergyVector;   //!
   TBranch        *b_Dc1HitsVector_x;   //!
   TBranch        *b_Dc1HitsVector_y;   //!
   TBranch        *b_Dc1HitsVector_z;   //!
   TBranch        *b_Dc2HitsVector_x;   //!
   TBranch        *b_Dc2HitsVector_y;   //!
   TBranch        *b_Dc2HitsVector_z;   //!

   Reader(TTree *tree=0);
   virtual ~Reader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Reader_cxx
Reader::Reader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("B5.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("B5.root");
      }
      f->GetObject("B5",tree);

   }
   Init(tree);
}

Reader::~Reader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Reader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Reader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Reader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ECEnergyVector = 0;
   HCEnergyVector = 0;
   Dc1HitsVector_x = 0;
   Dc1HitsVector_y = 0;
   Dc1HitsVector_z = 0;
   Dc2HitsVector_x = 0;
   Dc2HitsVector_y = 0;
   Dc2HitsVector_z = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Dc1Hits", &Dc1Hits, &b_Dc1Hits);
   fChain->SetBranchAddress("Dc2Hits", &Dc2Hits, &b_Dc2Hits);
   fChain->SetBranchAddress("ECEnergy", &ECEnergy, &b_ECEnergy);
   fChain->SetBranchAddress("HCEnergy", &HCEnergy, &b_HCEnergy);
   fChain->SetBranchAddress("Time1", &Time1, &b_Time1);
   fChain->SetBranchAddress("Time2", &Time2, &b_Time2);
   fChain->SetBranchAddress("ECEnergyVector", &ECEnergyVector, &b_ECEnergyVector);
   fChain->SetBranchAddress("HCEnergyVector", &HCEnergyVector, &b_HCEnergyVector);
   fChain->SetBranchAddress("Dc1HitsVector_x", &Dc1HitsVector_x, &b_Dc1HitsVector_x);
   fChain->SetBranchAddress("Dc1HitsVector_y", &Dc1HitsVector_y, &b_Dc1HitsVector_y);
   fChain->SetBranchAddress("Dc1HitsVector_z", &Dc1HitsVector_z, &b_Dc1HitsVector_z);
   fChain->SetBranchAddress("Dc2HitsVector_x", &Dc2HitsVector_x, &b_Dc2HitsVector_x);
   fChain->SetBranchAddress("Dc2HitsVector_y", &Dc2HitsVector_y, &b_Dc2HitsVector_y);
   fChain->SetBranchAddress("Dc2HitsVector_z", &Dc2HitsVector_z, &b_Dc2HitsVector_z);
   Notify();
}

Bool_t Reader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Reader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Reader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Reader_cxx
