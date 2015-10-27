#include <iostream>
#include "Classes/hadtopBoost.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TEfficiency.h"
#include "TLegend.h"

#include <boost/progress.hpp>
using namespace std;

void ProduceWeightedHistos(TFile *f1, Double_t weight, TH1F &true_den, TH1F &true_num, TH1F &rela_den, TH1F &rela_num ){


    TDirectory * dir = (TDirectory*)f1->Get("hadtopBoost");
    TTree * tree;
    dir->GetObject("events",tree);
    hadtopBoost htBoostC;// = new hadtopBoost(tree);
    htBoostC.Init(tree);

    int nentries = htBoostC.fChain->GetEntries();
    boost::progress_display show_progress( nentries );

    for(int i=0; i<nentries; i++){


        htBoostC.GetEntry(i);
        if((*htBoostC.jetMassSoftDrop)[0]>50. && (*htBoostC.jetTau3)[0]/(*htBoostC.jetTau2)[0]<0.7)
            true_den.Fill((*htBoostC.jetPt)[0], weight);
        if((*htBoostC.jetMassSoftDrop)[0]>50. && (*htBoostC.triggerBit)[0] && (*htBoostC.jetTau3)[0]/(*htBoostC.jetTau2)[0]<0.7)
            true_num.Fill((*htBoostC.jetPt)[0], weight);

        if((*htBoostC.jetMassSoftDrop)[0]>50. && (*htBoostC.triggerBit)[2] && (*htBoostC.jetTau3)[0]/(*htBoostC.jetTau2)[0]<0.7)
            rela_den.Fill((*htBoostC.jetPt)[0], weight);
        if((*htBoostC.jetMassSoftDrop)[0]>50. && (*htBoostC.triggerBit)[0] && (*htBoostC.triggerBit)[2] && (*htBoostC.jetTau3)[0]/(*htBoostC.jetTau2)[0]<0.7)
            rela_num.Fill((*htBoostC.jetPt)[0], weight);

        ++show_progress;
    }


}

int BoostedTriggerEfficiency(){


    const int N = 6;
    const float LUMI = 1.;
    float XSEC[N] = {1.,1.,1.,1.,1.,1.};
    TString SAMPLE[N] = {
        "flatTree_QCD_HT200to300",
        "flatTree_QCD_HT300to500",
        "flatTree_QCD_HT500to700",
        "flatTree_QCD_HT700to1000",
        "flatTree_QCD_HT1000to1500",
        "flatTree_QCD_HT1500to2000"
    };

    TH1F * true_den_qcd = new TH1F("true_den_qcd","true_den_qcd",100,0,1000);
    TH1F * true_num_qcd = new TH1F("true_num_qcd","true_num_qcd",100,0,1000);
    TH1F * rela_den_qcd = new TH1F("rela_den_qcd","rela_den_qcd",100,0,1000);
    TH1F * rela_num_qcd = new TH1F("rela_num_qcd","rela_num_qcd",100,0,1000);
    true_den_qcd->Sumw2();
    true_num_qcd->Sumw2();
    rela_den_qcd->Sumw2();
    rela_num_qcd->Sumw2();

for(int i=0; i<N ;i++){
    cout<<"Analyzing File: "<<SAMPLE[i]<<".root"<<endl;
    TFile *f1 = TFile::Open("root://eoscms.cern.ch///store/cmst3/user/kkousour/ttbar/flat/"+SAMPLE[i]+".root");

    TH1F * true_den = new TH1F("true_den"+SAMPLE[i],"true_den"+SAMPLE[i],100,0,1000);
    TH1F * true_num = new TH1F("true_num"+SAMPLE[i],"true_num"+SAMPLE[i],100,0,1000);
    TH1F * rela_den = new TH1F("rela_den"+SAMPLE[i],"rela_den"+SAMPLE[i],100,0,1000);
    TH1F * rela_num = new TH1F("rela_num"+SAMPLE[i],"rela_num"+SAMPLE[i],100,0,1000);

    ProduceWeightedHistos(f1, XSEC[i]/LUMI, *true_den, *true_num, *rela_den, *rela_num );
    true_den_qcd->Add(true_den);
    true_num_qcd->Add(true_num);
    rela_den_qcd->Add(rela_den);
    rela_num_qcd->Add(rela_num);

}
    TEfficiency *true_eff = new TEfficiency("true_eff","true_eff",80,400,800);
    true_eff->SetPassedHistogram(*true_num_qcd,"f");
    true_eff->SetTotalHistogram(*true_den_qcd,"f");

    TEfficiency *rela_eff = new TEfficiency("rela_eff","rela_eff",80,400,800);
    rela_eff->SetPassedHistogram(*rela_num_qcd,"f");
    rela_eff->SetTotalHistogram(*rela_den_qcd,"f");

    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    c1->cd();
    c1->SetTitle("Trigger Efficiency");
    true_eff->Draw();
    rela_eff->SetLineColor(kRed);
    true_eff->SetTitle("Trigger Efficiency; p_{T}; Efficiency");
    //rela_eff->GetXaxis()->SetTitle("Leading Jet p_{T}");
    rela_eff->Draw("same");

    TLegend *leg;
    leg = new TLegend(0.64,0.17,0.88,0.31);
    leg->AddEntry(true_eff,"True", "l");
    leg->AddEntry(rela_eff,"Relative HLT_PFJet260", "l");
    leg->SetTextFont(42);
    leg->SetFillColor(kWhite);
    leg->SetLineColor(kWhite);
    leg->SetBorderSize(0);
    leg->Draw();

    c1->SaveAs("trigg_eff_jetmass50_tau3207_bit0_bit2.png");
    c1->SaveAs("trigg_eff_jetmass50_tau3207_bit0_bit2.pdf");

    TCanvas *c2 = new TCanvas("c2","c2",800,600);
    c2->cd();
    true_den_qcd->Draw();
    true_den_qcd->SetMarkerColor(kRed);
    true_den_qcd->SetLineColor(kRed);
    true_num_qcd->Draw("same");


return 0;
}
