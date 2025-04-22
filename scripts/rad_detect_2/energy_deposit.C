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

Double_t energy_deposit(){
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
        beta_gamma += step_beta_gamma*step;
        // beta = sqrt(1 - 1/(beta_gamma*beta_gamma));
        nsteps++;
    }

    Double_t E_thin_detector[3][nsteps];
    Double_t E_thick_detector[3][nsteps];

    std::vector<Double_t> mass_list = {m_proton, m_deuteron, m_alpha};

    Double_t beta, beta_buf;
    Double_t energy_deposit_buf;
    Double_t energy_deposit_thin_buf, energy_deposit_thick_buf;
    Double_t particle_energy_buf;
    Double_t mass_buf;
    for(Int_t particle_id=0; particle_id<3; particle_id++){
        mass_buf = mass_list[particle_id];
        beta_gamma = beta_gamma_min;
        // beta_gamma = 0.1;

        for(Int_t step=0; beta_gamma<beta_gamma_max; step++){
            beta_gamma += step_beta_gamma*step;
            beta_buf = sqrt(1 - 1/(beta_gamma*beta_gamma)); // various beta
            beta = beta_buf; // incident particle beta

            for(Int_t n_detector=0; n_detector<10; n_detector++){
                energy_deposit_buf = 0.18*pow(beta, -1.7); //in MeV
                if(n_detector == 0){
                    E_thin_detector[particle_id][step] = energy_deposit_buf;
                }
                energy_deposit_thick_buf += energy_deposit_buf;
                particle_energy_buf =+ - energy_deposit_buf;
                if(particle_energy_buf < 0){
                    particle_energy_buf = 0;
                    beta = 0;
                }
                else{
                    beta = sqrt(1-pow(mass_buf/particle_energy_buf, 2));
                }
            }
            E_thick_detector[particle_id][step] = energy_deposit_thick_buf;
        }
    }

    TGraph *graph1 = new TGraph(nsteps, E_thick_detector[0], E_thin_detector[0]);
    graph1->SetTitle("Energy deposit of proton");
    graph1->GetXaxis()->SetTitle("Energy deposit in thick detector [MeV]");
    graph1->GetYaxis()->SetTitle("Energy deposit in thin detector [MeV]");
    graph1->SetMarkerStyle(20);
    graph1->SetMarkerSize(0.5);
    graph1->SetLineColor(kRed);

    TGraph *graph2 = new TGraph(nsteps, E_thick_detector[1], E_thin_detector[1]);
    graph2->SetTitle("Energy deposit of deuteron");
    graph2->GetXaxis()->SetTitle("Energy deposit in thick detector [MeV]");
    graph2->GetYaxis()->SetTitle("Energy deposit in thin detector [MeV]");
    graph2->SetMarkerStyle(20);
    graph2->SetMarkerSize(0.5);
    graph2->SetLineColor(kBlue);

    TGraph *graph3 = new TGraph(nsteps, E_thick_detector[3], E_thin_detector[3]);
    graph3->SetTitle("Energy deposit of alpha");
    graph3->GetXaxis()->SetTitle("Energy deposit in thick detector [MeV]");
    graph3->GetYaxis()->SetTitle("Energy deposit in thin detector [MeV]");
    graph3->SetMarkerStyle(20);
    graph3->SetMarkerSize(0.5);
    graph3->SetLineColor(kGreen);

    TCanvas *c1 = new TCanvas("c1", "Energy deposit", 800, 600);
    c1->SetGrid();
    c1->SetLogx();
    c1->SetLogy();

    graph1->Draw("APL");
    graph2->Draw("PL SAME");
    graph3->Draw("PL SAME");
    c1->Update();

    TLegend *legend = new TLegend(0.1, 0.7, 0.3, 0.9);
    legend->AddEntry(graph1, "Proton", "p");
    legend->AddEntry(graph2, "Deuteron", "p");
    legend->AddEntry(graph3, "Alpha", "p");
    legend->Draw();
    
    
}