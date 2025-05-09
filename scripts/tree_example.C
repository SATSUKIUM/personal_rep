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

Double_t tree_example(){
    Double_t nentries=0;
    // 保存ファイルと、その中に記録するツリー
    TFile *file = new TFile("data/tree.root", "RECREATE");
    TTree *tree = new TTree("tree", "example discription");

    // 各イベントごとの値を入れる変数
    Int_t time;
    Double_t energy_deposit[2];

    // ツリーにブランチを作る
    tree->Branch("time", &time, "time/I");
    tree->Branch("energy_deposit", energy_deposit, "energy_deposit[2]/D");

    // イベントごとに値を記録
    for(Int_t eventID=0; eventID<10000; eventID++){
        time = 2*eventID; // 例
        energy_deposit[0] = gRandom->Gaus(3,1);
        energy_deposit[1] = gRandom->Gaus(6,2);

        tree->Fill();
        nentries++;
    }

    // ツリーを書き出し
    tree->Write();
    file->Close();
    
    return nentries;
}