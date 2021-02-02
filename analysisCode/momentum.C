#include "Reader.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TVectorD.h"
#include "TF1.h"
#include "TMath.h"

int momentum()
{
  // Open and read file
  TFile *file = new TFile("B5.root");
  TTree *tree = (TTree*)file->Get("B5");
  Reader reader(tree);
  std::cout << "file open and loaded" << std::endl;

  TCanvas *c1 = new TCanvas();
  TRandom3 *random = new TRandom3();
  std::cout << "canvas and RNG created" << std::endl;

  Long64_t nentries = 1;//reader.fChain->GetEntries(); // Temporarily set to 1 for testing momentum calculation
  // Loop over events
  for(Long64_t ientry=0; ientry<nentries; ++ientry)
    {
      reader.GetEntry(ientry);
      
      double overallHitX[10];
      double overallHitZ[10];

      double dc1HitX[5];
      double dc1HitZ[5];
      for(unsigned int i=0; i<reader.Dc1HitsVector_x->size(); ++i)
	{
	  if(i<5)
	    {
	      dc1HitX[i] = random->Gaus(reader.Dc1HitsVector_x->at(i), 0.1);
	      dc1HitZ[i] = ((reader.Dc2HitsVector_z->at(i)-5/2)*500)-5000; // Calculates z pos in mm for first arm
	  
	      overallHitX[i] = dc1HitX[i];
	      overallHitZ[i] = dc1HitZ[i];
	    }
	}
      
      double dc2HitX[reader.Dc2HitsVector_x->size()];
      double dc2HitZ[reader.Dc2HitsVector_x->size()];
      for(unsigned int i=0; i<reader.Dc2HitsVector_x->size(); ++i)
	{
	  if(i<5)
	    {
	      dc2HitX[i] = random->Gaus(reader.Dc2HitsVector_x->at(i), 0.1);
	      dc2HitZ[i] = ((reader.Dc2HitsVector_z->at(i)-5/2)*500 - 1500) + 5000; // Calculates z pos in mm for second arm
	      
	      overallHitX[i+reader.Dc1HitsVector_x->size()] = dc2HitX[i];
	      overallHitZ[i+reader.Dc1HitsVector_x->size()] = dc2HitZ[i];
	    }
	}
      // Create linear function for fitting later
      TF1 *f1 = new TF1("f1", "[0]*x + [1]", -6250,-1000);
      f1->SetParameter(0, 1);
      f1->SetParameter(1, 0);
      TF1 *f2 = new TF1("f2", "[0]*x + [1]", 1000, 6000);
      f2->SetParameter(0, 1);
      f2->SetParameter(1, 0);

      // Create graphs
      TGraph *gr2 = new TGraph(reader.Dc2HitsVector_x->size(), dc2HitZ, dc2HitX);
      TGraph *gr1 = new TGraph(reader.Dc1HitsVector_x->size(), dc1HitZ, dc1HitX);
      TGraph *gr = new TGraph(reader.Dc1HitsVector_x->size()+reader.Dc2HitsVector_x->size(), overallHitZ, overallHitX);

      // Fit graphs
      TFitResultPtr r1 = gr->Fit(f1, "RS");
      TFitResultPtr r2 = gr->Fit(f2, "RS+");

      // Axis labels, drawing, etc
      gr->SetTitle(";Z Position (mm);X Position (mm)");
      gr->SetMarkerStyle(8);
      gr->SetMarkerSize(1);
      gr->Draw("AP");

      // Momentum reconstruction
      // Retrieve fit parameter values
      double m1 = r1->Value(0);
      double c1 = r1->Value(1);
      double m2 = r2->Value(0);
      double c2 = r2->Value(1);
      // Finding deflection angle
      f1->SetParameter(0, m1);
      f1->SetParameter(1, c1);
      f2->SetParameter(0, m2);
      f2->SetParameter(1, c2);
      // Finding angle between the two lines using SOHCAHTOA and Pythagoras
      double theta = TMath::ATan((TMath::Sqrt(500*500 + (f1->Eval(-5000) - f1->Eval(-4500))*(f1->Eval(-5000) - f1->Eval(-4500)))) // Length of original track segment
				 /(TMath::Sqrt(500*500 + (f2->Eval(2500) - f2->Eval(3000))*(f2->Eval(2500) - f2->Eval(3000))))); // Length of deflected track segment;
      std::cout << "theta: " << theta << std::endl;
      double momentum = (0.3*0.5*2)/(2*TMath::Sin(theta/2)); // p = Bql/2sin(theta/2)
      std::cout << "momentum: " << momentum << std::endl;
    }
  
  c1->Show();

  return 0;
}
