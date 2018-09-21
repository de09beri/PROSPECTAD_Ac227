#include "PROSPECT_Style.cc"
#include "TROOT.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TPave.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TColor.h"
#include "TExec.h"
#include <sstream>
#include "TLatex.h"
#include "TMath.h"
#include "TGraphErrors.h"

TH1F *makeHist(TGraphErrors *gr){
	int numPt = gr->GetN();
	double grx, gry, grxErr, gryErr;

	TH1F *h = new TH1F("h","h",154,0,154);
	
	for(int i=0;i<numPt;i++){
		gr->GetPoint(i,grx,gry);
		gryErr = gr->GetErrorY(i);

		h->SetBinContent((int)grx,gry);
		h->SetBinError((int)grx,gryErr);
	}

	return h;
}

TH1F *makeEnHist(TGraphErrors *gr){
	int numPt = gr->GetN();
	double grx, gry, grxErr, gryErr;

	TH1F *h = new TH1F("h","h",154,0,154);
	
	for(int i=0;i<numPt;i++){
		gr->GetPoint(i,grx,gry);
		gryErr = gr->GetErrorY(i);

		gry = gry*1000.0;		//convert MeV to keV
		gryErr = gryErr*1000.0;

		h->SetBinContent((int)grx,gry);
		h->SetBinError((int)grx,gryErr);
	}

	return h;
}

TGraphErrors *makeRelGr(TGraphErrors *gr, double mean, double meanErr){
	TGraphErrors *grRel = (TGraphErrors*)gr->Clone();
	int numPt = gr->GetN();

	double rel, relErr; 
	double grx, gry, gryErr;

	for(int i=0;i<numPt;i++){
		gr->GetPoint(i,grx,gry);
		gryErr = gr->GetErrorY(i);

		rel = gry/mean;
		relErr = rel * sqrt(pow(gryErr/gry,2) + pow(meanErr/mean,2));	

		grRel->SetPoint(i,grx,rel);
		grRel->SetPointError(i,0,relErr);
	}

	return grRel;
}


int PlotRxOnVsRxOff(){
	int timeBin = 0;

	setup_PROSPECT_style();
  	gROOT->ForceStyle();

	//-------------------------------------------------------------------------------------------------------
	TFile *fOn_ESmear = new TFile("/g/g20/berish1/AD_Ac227Analysis/Data/2018A/AllCells/Results/RxOn/Ac227_HistsPerTime.root");
	if(!fOn_ESmear){
		printf("File not found. Exiting. \n");
		return -1;
	}

	TH1F *hRnPSD_On_ESmear 	= (TH1F*)fOn_ESmear->Get(Form("hRnPSD_%i",timeBin));
	TH1F *hPoPSD_On_ESmear 	= (TH1F*)fOn_ESmear->Get(Form("hPoPSD_%i",timeBin));
	TH1F *hRnEn_On_ESmear   = (TH1F*)fOn_ESmear->Get(Form("hRnEn_%i",timeBin));
	TH1F *hPoEn_On_ESmear 	= (TH1F*)fOn_ESmear->Get(Form("hPoEn_%i",timeBin));
	TH1F *hRnPos_On_ESmear 	= (TH1F*)fOn_ESmear->Get(Form("hRnPos_%i",timeBin));
	TH1F *hPoPos_On_ESmear 	= (TH1F*)fOn_ESmear->Get(Form("hPoPos_%i",timeBin));
	TH1F *hRnPoDz_On_ESmear = (TH1F*)fOn_ESmear->Get(Form("hRnPoDz_%i",timeBin));


	TFile *fOff_ESmear = new TFile("/g/g20/berish1/AD_Ac227Analysis/Data/2018A/AllCells/Results/RxOff/Ac227_HistsPerTime.root");
	if(!fOff_ESmear){
		printf("File not found. Exiting. \n");
		return -1;
	}

	TH1F *hRnPSD_Off_ESmear = (TH1F*)fOff_ESmear->Get(Form("hRnPSD_%i",timeBin));
	TH1F *hPoPSD_Off_ESmear = (TH1F*)fOff_ESmear->Get(Form("hPoPSD_%i",timeBin));
	TH1F *hRnEn_Off_ESmear  = (TH1F*)fOff_ESmear->Get(Form("hRnEn_%i",timeBin));
	TH1F *hPoEn_Off_ESmear 	= (TH1F*)fOff_ESmear->Get(Form("hPoEn_%i",timeBin));
	TH1F *hRnPos_Off_ESmear = (TH1F*)fOff_ESmear->Get(Form("hRnPos_%i",timeBin));
	TH1F *hPoPos_Off_ESmear = (TH1F*)fOff_ESmear->Get(Form("hPoPos_%i",timeBin));
	TH1F *hRnPoDz_Off_ESmear= (TH1F*)fOff_ESmear->Get(Form("hRnPoDz_%i",timeBin));

	//-------------------------------------------------------------------------------------------------------
	TFile *fOn_E = new TFile("/g/g20/berish1/AD_Ac227Analysis/Data/2018A/AllCells_E/Results/RxOn/Ac227_HistsPerTime.root");
	if(!fOn_E){
		printf("File not found. Exiting. \n");
		return -1;
	}

	TH1F *hRnPSD_On_E  = (TH1F*)fOn_E->Get(Form("hRnPSD_%i",timeBin));
	TH1F *hPoPSD_On_E  = (TH1F*)fOn_E->Get(Form("hPoPSD_%i",timeBin));
	TH1F *hRnEn_On_E   = (TH1F*)fOn_E->Get(Form("hRnEn_%i",timeBin));
	TH1F *hPoEn_On_E   = (TH1F*)fOn_E->Get(Form("hPoEn_%i",timeBin));
	TH1F *hRnPos_On_E  = (TH1F*)fOn_E->Get(Form("hRnPos_%i",timeBin));
	TH1F *hPoPos_On_E  = (TH1F*)fOn_E->Get(Form("hPoPos_%i",timeBin));
	TH1F *hRnPoDz_On_E = (TH1F*)fOn_E->Get(Form("hRnPoDz_%i",timeBin));


	TFile *fOff_E = new TFile("/g/g20/berish1/AD_Ac227Analysis/Data/2018A/AllCells_E/Results/RxOff/Ac227_HistsPerTime.root");
	if(!fOff_E){
		printf("File not found. Exiting. \n");
		return -1;
	}

	TH1F *hRnPSD_Off_E = (TH1F*)fOff_E->Get(Form("hRnPSD_%i",timeBin));
	TH1F *hPoPSD_Off_E = (TH1F*)fOff_E->Get(Form("hPoPSD_%i",timeBin));
	TH1F *hRnEn_Off_E  = (TH1F*)fOff_E->Get(Form("hRnEn_%i",timeBin));
	TH1F *hPoEn_Off_E  = (TH1F*)fOff_E->Get(Form("hPoEn_%i",timeBin));
	TH1F *hRnPos_Off_E = (TH1F*)fOff_E->Get(Form("hRnPos_%i",timeBin));
	TH1F *hPoPos_Off_E = (TH1F*)fOff_E->Get(Form("hPoPos_%i",timeBin));
	TH1F *hRnPoDz_Off_E= (TH1F*)fOff_E->Get(Form("hRnPoDz_%i",timeBin));

	double RxOnTime = 934.360379;	//hrs
	double RxOffTime = 864.126678;	//hrs

	//-------------------------------------------------------------------------------------------------------
	hRnEn_On_ESmear->Scale(1/RxOnTime);
	hPoEn_On_ESmear->Scale(1/RxOnTime);
	
	hRnEn_Off_ESmear->Scale(1/RxOffTime);
	hPoEn_Off_ESmear->Scale(1/RxOffTime);

	//~~~~~~~~~~~~~~~
	hRnEn_On_E->Scale(1/RxOnTime);
	hPoEn_On_E->Scale(1/RxOnTime);

	hRnEn_Off_E->Scale(1/RxOffTime);
	hPoEn_Off_E->Scale(1/RxOffTime);

	//-------------------------------------------------------------------------------------------------------
	TH1F *hRnEn_Sub_ESmear = (TH1F*)hRnEn_On_ESmear->Clone();
	hRnEn_Sub_ESmear->Sumw2();
	hRnEn_Sub_ESmear->Add(hRnEn_Off_ESmear,-1);
	hRnEn_Sub_ESmear->Divide(hRnEn_Off_ESmear);

	TH1F *hPoEn_Sub_ESmear = (TH1F*)hPoEn_On_ESmear->Clone();
	hPoEn_Sub_ESmear->Sumw2();
	hPoEn_Sub_ESmear->Add(hPoEn_Off_ESmear,-1);
	hPoEn_Sub_ESmear->Divide(hPoEn_Off_ESmear);

	TH1F *hRnEn_Sub_E = (TH1F*)hRnEn_On_E->Clone();
	hRnEn_Sub_E->Sumw2();
	hRnEn_Sub_E->Add(hRnEn_Off_E,-1);
	hRnEn_Sub_E->Divide(hRnEn_Off_E);

	TH1F *hPoEn_Sub_E = (TH1F*)hPoEn_On_E->Clone();
	hPoEn_Sub_E->Sumw2();
	hPoEn_Sub_E->Add(hPoEn_Off_E,-1);
	hPoEn_Sub_E->Divide(hPoEn_Off_E);


	//-------------------------------------------------------------------------------------------------------
	TPaveText *pt;
	int E_col = 1, ES_col = 2;
	int on_m = 21, off_m = 24;

	const char *yTitle = "(On-Off)/Off";

	TLegend *leg;

	TCanvas *cPoEn = new TCanvas("cPoEn","Po En",1);
	gPad->SetGrid();
	hPoEn_Sub_E->SetMarkerColor(E_col);
	hPoEn_Sub_E->SetLineColor(E_col);
	hPoEn_Sub_E->SetMarkerStyle(2);
	hPoEn_Sub_E->GetYaxis()->SetRangeUser(-0.4,0.4);
	hPoEn_Sub_E->GetXaxis()->SetRangeUser(0.6,1.2);
	hPoEn_Sub_E->Draw();
	hPoEn_Sub_ESmear->SetMarkerColor(ES_col);
	hPoEn_Sub_ESmear->SetLineColor(ES_col);
	hPoEn_Sub_ESmear->SetMarkerStyle(2);
	hPoEn_Sub_ESmear->Draw("same");
	leg = new TLegend(0.7,0.75,0.95,0.9);
	leg->AddEntry(hPoEn_Sub_E,"Po^{215} E");
	leg->AddEntry(hPoEn_Sub_ESmear,"Po^{215} ESmear");
	leg->Draw();
	cPoEn->SaveAs("/g/g20/berish1/AD_Ac227Analysis/Data/2018A/CompareE/PoEn.C");
	cPoEn->SaveAs("/g/g20/berish1/AD_Ac227Analysis/Data/2018A/CompareE/PoEn.png");
	

	//-------------------------------------------------------------------------------------------------------
	TFile *f_ESmear = new TFile("/g/g20/berish1/AD_Ac227Analysis/Data/2018A/AllCells/Results/Ac227_GraphsPerTime.root");
	if(!f_ESmear){
		printf("File not found. Exiting. \n");
		return -1;
	}

	TGraphErrors *grPoEnSigma_ESmear = (TGraphErrors*)f_ESmear->Get("grPoEnSigma");

	TFile *f_E = new TFile("/g/g20/berish1/AD_Ac227Analysis/Data/2018A/AllCells_E/Results/Ac227_GraphsPerTime.root");
	if(!f_E){
		printf("File not found. Exiting. \n");
		return -1;
	}

	TGraphErrors *grPoEnSigma_E = (TGraphErrors*)f_E->Get("grPoEnSigma");



	grPoEnSigma_ESmear->Fit("pol0","0");
	TGraphErrors *grRelPoEnSigma_ESmear = makeRelGr(grPoEnSigma_ESmear, grPoEnSigma_ESmear->GetFunction("pol0")->GetParameter(0), grPoEnSigma_ESmear->GetFunction("pol0")->GetParError(0));
	
	grPoEnSigma_E->Fit("pol0","0");
	TGraphErrors *grRelPoEnSigma_E = makeRelGr(grPoEnSigma_E, grPoEnSigma_E->GetFunction("pol0")->GetParameter(0), grPoEnSigma_E->GetFunction("pol0")->GetParError(0));


	TCanvas *cPoEnSigma = new TCanvas("cPoEnSigma","Po En Sigma",1);
	grPoEnSigma_E->SetMarkerColor(E_col);
	grPoEnSigma_E->SetLineColor(E_col);
	grPoEnSigma_E->SetMarkerStyle(20);
	grPoEnSigma_E->GetYaxis()->SetRangeUser(0.04,0.052);
	grPoEnSigma_E->GetXaxis()->SetTimeDisplay(1);
	grPoEnSigma_E->GetXaxis()->SetTimeFormat("%m/%d");
	grPoEnSigma_E->GetYaxis()->SetTitle("E_{#sigmaPo} [MeV]");
	grPoEnSigma_E->Draw("AP");	
	grPoEnSigma_ESmear->SetMarkerColor(ES_col);
	grPoEnSigma_ESmear->SetLineColor(ES_col);
	grPoEnSigma_ESmear->SetMarkerStyle(20);
	grPoEnSigma_ESmear->Draw("P");	
	leg = new TLegend(0.7,0.75,0.95,0.9);
	leg->AddEntry(grPoEnSigma_E,"Po^{215} E");
	leg->AddEntry(grPoEnSigma_ESmear,"Po^{215} ESmear");
	leg->Draw();
	cPoEnSigma->SaveAs("/g/g20/berish1/AD_Ac227Analysis/Data/2018A/CompareE/PoEnSigma.C");
	cPoEnSigma->SaveAs("/g/g20/berish1/AD_Ac227Analysis/Data/2018A/CompareE/PoEnSigma.png");
	

	TCanvas *cRelPoEnSigma = new TCanvas("cRelPoEnSigma","Relative Po En Sigma",1);
	grRelPoEnSigma_E->SetMarkerColor(E_col);
	grRelPoEnSigma_E->SetLineColor(E_col);
	grRelPoEnSigma_E->GetYaxis()->SetRangeUser(0.94,1.1);
	grRelPoEnSigma_E->GetXaxis()->SetTimeDisplay(1);
	grRelPoEnSigma_E->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoEnSigma_E->GetYaxis()->SetTitle("E_{#sigmaPo}/#LTE_{#sigmaPo}#GT");
	grRelPoEnSigma_E->Draw("AP");	
	grRelPoEnSigma_ESmear->SetMarkerColor(ES_col);
	grRelPoEnSigma_ESmear->SetLineColor(ES_col);
	grRelPoEnSigma_ESmear->Draw("P");	
	leg = new TLegend(0.7,0.75,0.95,0.9);
	leg->AddEntry(grRelPoEnSigma_E,"Po^{215} E");
	leg->AddEntry(grRelPoEnSigma_ESmear,"Po^{215} ESmear");
	leg->Draw();
	cRelPoEnSigma->SaveAs("/g/g20/berish1/AD_Ac227Analysis/Data/2018A/CompareE/RelPoEnSigma.C");
	cRelPoEnSigma->SaveAs("/g/g20/berish1/AD_Ac227Analysis/Data/2018A/CompareE/RelPoEnSigma.png");


	return 0;
} 	//end PlotDistributionsVsTime
