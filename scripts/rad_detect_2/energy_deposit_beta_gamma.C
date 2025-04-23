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

Double_t energy_deposit_beta_gamma(){
    Double_t beta_gamma = 0.0;
    Double_t step_beta_gamma = 0.001;
    Double_t beta_gamma_min = 0.26;
    Double_t beta_gamma_max = 10.0;

    Double_t m_proton = 0.938272*1000.0;
    Double_t m_deuteron = 1.875612*1000.0;
    Double_t m_alpha = 3.727379*1000.0;
    //in MeV

    Int_t nsteps = 0;
    beta_gamma = beta_gamma_min;
    for(Int_t step=0; beta_gamma<beta_gamma_max; step++){
        beta_gamma += step_beta_gamma;
        // beta = sqrt(1 - 1/(beta_gamma*beta_gamma));
        nsteps++;
    }

    std::vector<Double_t> beta_gamma_list(nsteps);
    std::vector<Double_t> beta_list(nsteps);
    std::vector<std::vector<Double_t>> E_thin_detector(3, std::vector<Double_t>(nsteps));
    std::vector<std::vector<Double_t>> E_thick_detector(3, std::vector<Double_t>(nsteps));

    std::vector<Double_t> mass_list = {m_proton, m_deuteron, m_alpha};

    Double_t beta = 0.0, beta_buf = 0.0;
    Double_t energy_deposit_buf = 0.0;
    Double_t energy_deposit_thin_buf = 0.0, energy_deposit_thick_buf = 0.0;
    Double_t particle_energy_buf = 0.0;
    Double_t mass_buf = 0.0;
    Double_t charge_buf = 0.0;
    Double_t energy_deposit_thin_max = 0.0, energy_deposit_thick_max = 0.0, energy_deposit_thin_min = 1e5, energy_deposit_thick_min = 1e5;
    for(Int_t particle_id=0; particle_id<3; particle_id++){
        mass_buf = mass_list[particle_id];
        charge_buf = pow(2.0, particle_id); // 1 for proton, 2 for deuteron, 4 for alpha
        beta_gamma = beta_gamma_min;
        // beta_gamma = 0.1;

        for(Int_t step=0; beta_gamma<beta_gamma_max; step++){
            energy_deposit_thick_buf = 0.0;

            beta_gamma += step_beta_gamma;
            beta_buf = sqrt(beta_gamma*beta_gamma/(1+beta_gamma*beta_gamma)); // various beta
            beta = beta_buf; // incident particle beta
            beta_list[step] = beta_buf;
            beta_gamma_list[step] = beta_gamma;

            particle_energy_buf = mass_buf*(beta_gamma/beta_buf); // in MeV

            for(Int_t n_detector=0; n_detector<10; n_detector++){
                energy_deposit_buf = charge_buf*0.18*pow(beta_buf, -1.7); //in MeV
                // energy_deposit_buf = 0.18*pow(beta_buf, -1.7); //in MeV
                if(n_detector == 0){
                    E_thin_detector[particle_id][step] = energy_deposit_buf;
                    if(step % 50 == 0){
                        // printf("thin detector: %d, %d, %f\n", particle_id, step, energy_deposit_buf);
                    }
                    if(energy_deposit_buf > energy_deposit_thin_max){
                        energy_deposit_thin_max = energy_deposit_buf;
                    }
                    if(energy_deposit_buf < energy_deposit_thin_min){
                        energy_deposit_thin_min = energy_deposit_buf;
                    }
                }
                energy_deposit_thick_buf += energy_deposit_buf;
                particle_energy_buf += - energy_deposit_buf;
                if(particle_energy_buf < mass_list[particle_id]){
                    particle_energy_buf = mass_list[particle_id];
                    beta_buf = 0.0;
                }
                else{
                    beta_buf = sqrt(1-pow(mass_buf/particle_energy_buf, 2));
                }
                if(step % 50 == 0){
                    printf("\tparticle: %f MeV/c2, incident beta gamma: %f, num of detector passed: %d, total energy: %f = %f * %f\n", mass_list[particle_id], beta_gamma, n_detector+1, particle_energy_buf, mass_buf, 1.0/sqrt(1-beta_buf*beta_buf));
                }
                
            }
            E_thick_detector[particle_id][step] = energy_deposit_thick_buf;
            if(energy_deposit_thick_buf > energy_deposit_thick_max){
                energy_deposit_thick_max = energy_deposit_thick_buf;
            }
            if(energy_deposit_thick_buf < energy_deposit_thick_min){
                energy_deposit_thick_min = energy_deposit_thick_buf;
            }
            if(step % 50 == 0){
                // printf("thick detector: %d, %d, %f\n", particle_id, step, energy_deposit_thick_buf);
            }
            // printf("thick detector: %d, %d, %f\n", particle_id, step, energy_deposit_thick_buf);
        }
    }

    // TGraph *graph1 = new TGraph(E_thin_detector[0].size(), &beta_list[0], &E_thin_detector[0][0]);
    // TGraph *graph2 = new TGraph(E_thin_detector[1].size(), &beta_list[0], &E_thin_detector[1][0]);
    // TGraph *graph3 = new TGraph(E_thin_detector[2].size(), &beta_list[0], &E_thin_detector[2][0]);

    // graph1->SetTitle("Energy deposit of proton");
    // graph1->GetXaxis()->SetTitle("beta");
    // graph1->GetYaxis()->SetTitle("Energy deposit in thin detector [MeV]");
    // graph1->SetMarkerStyle(20);
    // graph1->SetMarkerSize(0.5);
    // graph1->SetLineColor(kRed);

    // graph2->SetTitle("Energy deposit of deuteron");
    // graph2->GetXaxis()->SetTitle("beta");
    // graph2->GetYaxis()->SetTitle("Energy deposit in thin detector [MeV]");
    // graph2->SetMarkerStyle(20);
    // graph2->SetMarkerSize(0.5);
    // graph2->SetLineColor(kBlue);

    // graph3->SetTitle("Energy deposit of alpha");
    // graph3->GetXaxis()->SetTitle("beta");
    // graph3->GetYaxis()->SetTitle("Energy deposit in thin detector [MeV]");
    // graph3->SetMarkerStyle(20);
    // graph3->SetMarkerSize(0.5);
    // graph3->SetLineColor(kGreen);

    
    
    TGraph *graph1 = new TGraph(E_thin_detector[0].size(), &beta_gamma_list[0], &E_thick_detector[0][0]);
    TGraph *graph2 = new TGraph(E_thin_detector[1].size(), &beta_gamma_list[0], &E_thick_detector[1][0]);
    TGraph *graph3 = new TGraph(E_thin_detector[2].size(), &beta_gamma_list[0], &E_thick_detector[2][0]);
    
    graph1->SetTitle("Energy deposit of proton");
    graph1->GetXaxis()->SetTitle("beta gamma");
    graph1->GetYaxis()->SetTitle("Energy deposit in thick detector [MeV]");
    graph1->SetMarkerStyle(20);
    graph1->SetMarkerSize(0.5);
    graph1->SetLineColor(kRed);

    graph2->SetTitle("Energy deposit of deuteron");
    graph2->GetXaxis()->SetTitle("beta gamma");
    graph2->GetYaxis()->SetTitle("Energy deposit in thick detector [MeV]");
    graph2->SetMarkerStyle(20);
    graph2->SetMarkerSize(0.5);
    graph2->SetLineColor(kBlue);

    graph3->SetTitle("Energy deposit of alpha");
    graph3->GetXaxis()->SetTitle("beta gamma");
    graph3->GetYaxis()->SetTitle("Energy deposit in thick detector [MeV]");
    graph3->SetMarkerStyle(20);
    graph3->SetMarkerSize(0.5);
    graph3->SetLineColor(kGreen);

    TCanvas *c1 = new TCanvas("c1", "Energy deposit", 800, 600);
    c1->SetGrid();
    c1->SetLogx();
    // c1->SetLogy();

    graph1->Draw("AL");
    graph2->Draw("L SAME");
    graph3->Draw("L SAME");
    c1->Update();

    TLegend *legend = new TLegend(0.1, 0.7, 0.3, 0.9);
    legend->AddEntry(graph1, "Proton", "l");
    legend->AddEntry(graph2, "Deuteron", "l");
    legend->AddEntry(graph3, "Alpha", "l");
    legend->Draw();
    
    return nsteps;
}