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


void PlotRnPoVsCell(){

	setup_PROSPECT_style();
    gROOT->ForceStyle();

	TFile *f = new TFile(Form("%s/Ac227_GraphsPerCell.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));

	TGraphErrors *grRate 		= (TGraphErrors*)f->Get("grRate");
	TGraphErrors *grRnPSDEff 	= (TGraphErrors*)f->Get("grRnPSDEff");
	TGraphErrors *grPoPSDEff 	= (TGraphErrors*)f->Get("grPoPSDEff");
	TGraphErrors *grRnEnEff 	= (TGraphErrors*)f->Get("grRnEnEff");
	TGraphErrors *grPoEnEff 	= (TGraphErrors*)f->Get("grPoEnEff");
	TGraphErrors *grRnPoDzEff 	= (TGraphErrors*)f->Get("grRnPoDzEff");
	TGraphErrors *grTotEff 		= (TGraphErrors*)f->Get("grTotEff");
	TGraphErrors *grLifetime 	= (TGraphErrors*)f->Get("grLifetime");
	TGraphErrors *grPoPSDMean 	= (TGraphErrors*)f->Get("grPoPSDMean");
	TGraphErrors *grPoPSDSigma 	= (TGraphErrors*)f->Get("grPoPSDSigma");
	TGraphErrors *grPoEnMean 	= (TGraphErrors*)f->Get("grPoEnMean");
	TGraphErrors *grPoEnSigma 	= (TGraphErrors*)f->Get("grPoEnSigma");
	TGraphErrors *grPoPosMean 	= (TGraphErrors*)f->Get("grPoPosMean");
	TGraphErrors *grPoPosSigma 	= (TGraphErrors*)f->Get("grPoPosSigma");
	TGraphErrors *grRnPoDzMean 	= (TGraphErrors*)f->Get("grRnPoDzMean");
	TGraphErrors *grRnPoDzSigma = (TGraphErrors*)f->Get("grRnPoDzSigma");

	f->Close();
	
	//-------------------------------------------------------------------------------------------------------
cout<<"\n RATE"<<endl;
	grRate->Fit("pol0","0");
	TGraphErrors *grRelRate = makeRelGr(grRate, grRate->GetFunction("pol0")->GetParameter(0), grRate->GetFunction("pol0")->GetParError(0));
cout<<"\n PSD MEAN"<<endl;
	grPoPSDMean->Fit("pol0","0");
	TGraphErrors *grRelPoPSDMean = makeRelGr(grPoPSDMean, grPoPSDMean->GetFunction("pol0")->GetParameter(0), grPoPSDMean->GetFunction("pol0")->GetParError(0));
cout<<"\n PSD SIGMA"<<endl;
	grPoPSDSigma->Fit("pol0","0");
	TGraphErrors *grRelPoPSDSigma = makeRelGr(grPoPSDSigma, grPoPSDSigma->GetFunction("pol0")->GetParameter(0), grPoPSDSigma->GetFunction("pol0")->GetParError(0));
cout<<"\n ENERGY MEAN"<<endl;
	grPoEnMean->Fit("pol0","0");
	TGraphErrors *grRelPoEnMean = makeRelGr(grPoEnMean, grPoEnMean->GetFunction("pol0")->GetParameter(0), grPoEnMean->GetFunction("pol0")->GetParError(0));
cout<<"\n ENERGY SIGMA"<<endl;
	grPoEnSigma->Fit("pol0","0");
	TGraphErrors *grRelPoEnSigma = makeRelGr(grPoEnSigma, grPoEnSigma->GetFunction("pol0")->GetParameter(0), grPoEnSigma->GetFunction("pol0")->GetParError(0));
cout<<"\n POSITION MEAN"<<endl;
	grPoPosMean->Fit("pol0","0");
	TGraphErrors *grRelPoPosMean = makeRelGr(grPoPosMean, grPoPosMean->GetFunction("pol0")->GetParameter(0), grPoPosMean->GetFunction("pol0")->GetParError(0));
cout<<"\n POSITION SIGMA"<<endl;
	grPoPosSigma->Fit("pol0","0");
	TGraphErrors *grRelPoPosSigma = makeRelGr(grPoPosSigma, grPoPosSigma->GetFunction("pol0")->GetParameter(0), grPoPosSigma->GetFunction("pol0")->GetParError(0));
cout<<"\n DZ MEAN"<<endl;
	grRnPoDzMean->Fit("pol0","0");
	TGraphErrors *grRelRnPoDzMean = makeRelGr(grRnPoDzMean, grRnPoDzMean->GetFunction("pol0")->GetParameter(0), grRnPoDzMean->GetFunction("pol0")->GetParError(0));
cout<<"\n DZ SIGMA"<<endl;
	grRnPoDzSigma->Fit("pol0","0");
	TGraphErrors *grRelRnPoDzSigma = makeRelGr(grRnPoDzSigma, grRnPoDzSigma->GetFunction("pol0")->GetParameter(0), grRnPoDzSigma->GetFunction("pol0")->GetParError(0));
	
	//-------------------------------------------------------------------------------------------------------
	TH2F *h2DRate = new TH2F("h2DRate","Rate per cell",14,0,14,11,0,11);
	TH2F *h2DRelRate = new TH2F("h2DRelRate","Rate per cell",14,0,14,11,0,11);
	int numPt = grRate->GetN();
	double grx,gry;
	int binx,biny;
	for(int i=0;i<numPt;i++){
		grRate->GetPoint(i,grx,gry);
		binx = (int)grx%14 + 1;
		biny = ((int)grx/14) + 1;
		if(grx!=111) h2DRate->SetBinContent(binx,biny,gry);

		grRelRate->GetPoint(i,grx,gry);
		if(grx!=111) h2DRelRate->SetBinContent(binx,biny,gry);
	}


	//-------------------------------------------------------------------------------------------------------
	const char *xLabel = "Cell";

	TCanvas *cRate = new TCanvas("cRate","Rate per cell",1000,400);
	grRate->GetXaxis()->SetTitle(xLabel);
	grRate->GetYaxis()->SetTitle("R_{RnPo} [mHz]");
	grRate->Draw("AP");
	cRate->SaveAs(Form("%s/RatePerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRate->SaveAs(Form("%s/RatePerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	gStyle->SetPaintTextFormat("2.2f");
	TCanvas *cRate2D = new TCanvas("cRate2D","Rate per cell",1);
	gPad->SetRightMargin(0.16);
	h2DRate->SetMarkerSize(1.25);
	h2DRate->SetMinimum(3.2);
	h2DRate->Draw("colz && text");	
	cRate2D->SaveAs(Form("%s/2DRatePerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRate2D->SaveAs(Form("%s/2DRatePerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelRate2D = new TCanvas("cRelRate2D","Relative rate per cell",1);
	gPad->SetRightMargin(0.16);
	h2DRelRate->SetMarkerSize(1.25);
	h2DRelRate->SetMinimum(0.93);
	h2DRelRate->Draw("colz && text");	
	cRelRate2D->SaveAs(Form("%s/2DRelRatePerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRate2D->SaveAs(Form("%s/2DRelRatePerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelRate = new TCanvas("cRelRate","Relative rate per cell",1000,400);
	grRelRate->GetXaxis()->SetTitle(xLabel);
	grRelRate->GetYaxis()->SetTitle("R_{RnPo}/#LTR_{RnPo}#GT");  
	grRelRate->GetYaxis()->SetRangeUser(0.93,1.07);
	grRelRate->Draw("AP");
	cRelRate->SaveAs(Form("%s/RelativeRatePerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRate->SaveAs(Form("%s/RelativeRatePerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoPSDMean = new TCanvas("cRelPoPSDMean","Relative Po PSD Mean",1000,400);
	grRelPoPSDMean->GetXaxis()->SetTitle(xLabel);
	grRelPoPSDMean->GetYaxis()->SetTitle("PSD_{#muPo}/#LTPSD_{#muPo}#GT");
	grRelPoPSDMean->Draw("AP");
	cRelPoPSDMean->SaveAs(Form("%s/RelativePoPSDMeanPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPSDMean->SaveAs(Form("%s/RelativePoPSDMeanPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoPSDSigma = new TCanvas("cRelPoPSDSigma","Relative Po PSD Sigma",1000,400);
	grRelPoPSDSigma->GetXaxis()->SetTitle(xLabel);
	grRelPoPSDSigma->GetYaxis()->SetTitle("PSD_{#sigmaPo}/#LTPSD_{#sigmaPo}#GT");
	grRelPoPSDSigma->Draw("AP");
	cRelPoPSDSigma->SaveAs(Form("%s/RelativePoPSDSigmaPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPSDSigma->SaveAs(Form("%s/RelativePoPSDSigmaPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoEnMean = new TCanvas("cRelPoEnMean","Relative Po Energy Mean",1000,400);
	grRelPoEnMean->GetXaxis()->SetTitle(xLabel);
	grRelPoEnMean->GetYaxis()->SetTitle("E_{#muPo}/#LTE_{#muPo}#GT");
	grRelPoEnMean->GetYaxis()->SetRangeUser(0.994,1.006);
	grRelPoEnMean->Draw("AP");
	cRelPoEnMean->SaveAs(Form("%s/RelativePoEnMeanPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoEnMean->SaveAs(Form("%s/RelativePoEnMeanPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoEnSigma = new TCanvas("cRelPoEnSigma","Relative Po Energy Sigma",1000,400);
	grRelPoEnSigma->GetXaxis()->SetTitle(xLabel);
	grRelPoEnSigma->GetYaxis()->SetTitle("E_{#sigmaPo}/#LTE_{#sigmaPo}#GT");
	grRelPoEnSigma->Draw("AP");
	cRelPoEnSigma->SaveAs(Form("%s/RelativePoEnSigmaPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoEnSigma->SaveAs(Form("%s/RelativePoEnSigmaPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoPosMean = new TCanvas("cRelPoPosMean","Relative Po Position Mean",1000,400);
	grRelPoPosMean->GetXaxis()->SetTitle(xLabel);
	grRelPoPosMean->GetYaxis()->SetTitle("z_{#muPo}/#LTz_{#muPo}#GT");
	grRelPoPosMean->Draw("AP");
	cRelPoPosMean->SaveAs(Form("%s/RelativePoPosMeanPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPosMean->SaveAs(Form("%s/RelativePoPosMeanPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoPosSigma = new TCanvas("cRelPoPosSigma","Relative Po Position Sigma",1000,400);
	grRelPoPosSigma->GetXaxis()->SetTitle(xLabel);
	grRelPoPosSigma->GetYaxis()->SetTitle("z_{RMS Po}/#LTz_{RMS Po}#GT");
	grRelPoPosSigma->Draw("AP");
	cRelPoPosSigma->SaveAs(Form("%s/RelativePoPosSigmaPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPosSigma->SaveAs(Form("%s/RelativePoPosSigmaPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelRnPoDzMean = new TCanvas("cRelRnPoDzMean","Relative RnPo Dz Mean",1000,400);
	grRelRnPoDzMean->GetXaxis()->SetTitle(xLabel);
	grRelRnPoDzMean->GetYaxis()->SetTitle("dz_{#mu}/#LTdz_{#mu}#GT");
	grRelRnPoDzMean->Draw("AP");
	cRelRnPoDzMean->SaveAs(Form("%s/RelativeRnPoDzMeanPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRnPoDzMean->SaveAs(Form("%s/RelativeRnPoDzMeanPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelRnPoDzSigma = new TCanvas("cRelRnPoDzSigma","Relative RnPo Dz Sigma",1000,400);
	grRelRnPoDzSigma->GetXaxis()->SetTitle(xLabel);
	grRelRnPoDzSigma->GetYaxis()->SetTitle("dz_{#sigma}/#LTdz_{#sigma}#GT");
	grRelRnPoDzSigma->Draw("AP");
	cRelRnPoDzSigma->SaveAs(Form("%s/RelativeRnPoDzSigmaPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRnPoDzSigma->SaveAs(Form("%s/RelativeRnPoDzSigmaPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cEff = new TCanvas("cEff","Total Efficiency",1000,400);
	grTotEff->GetXaxis()->SetTitle(xLabel);
	grTotEff->GetYaxis()->SetTitle("Efficiency");
	grTotEff->Draw("AP");
	cEff->SaveAs(Form("%s/EfficiencyPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cEff->SaveAs(Form("%s/EfficiencyPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

} 	//end PlotRnPoVsCell
