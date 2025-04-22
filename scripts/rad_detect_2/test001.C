#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <vector>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TPad.h>
#include <TColor.h>


#include <fstream>
#include <filesystem>
#include <TSystem.h>

#include <iomanip>
#include <chrono>
#include <ctime> //時刻情報

#include <TLegend.h>

#include <fstream>

#include <TApplication.h>
#include <TChainElement.h>
#include <TObjArray.h>

#include <TFile.h>
#include <TRandom.h> // 乱数

Double_t test001(){
    Double_t beta_gamma = 0.0;
    Double_t step_beta_gamma = 0.01;
    Double_t beta_gamma_min = 0.1;
    Double_t beta_gamma_max = 10.0;

    Double_t m_proton = 0.938272;
    Double_t m_deuteron = 1.875612;
    Double_t m_alpha = 3.727379;
    //in MeV

    Int_t nsteps = 0;
    beta_gamma = beta_gamma_min;
    for(Int_t step=0; beta_gamma<beta_gamma_max; step++){
        beta_gamma += step_beta_gamma;
        // beta = sqrt(1 - 1/(beta_gamma*beta_gamma));
        nsteps++;
    }

    return nsteps;
}