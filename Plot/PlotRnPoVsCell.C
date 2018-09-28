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
	TGraphErrors *grRnPoDzSigma 	= (TGraphErrors*)f->Get("grRnPoDzSigma");

	TGraph *grRnPSDChiSq = (TGraph*)f->Get("grRnPSDChiSq");
	TGraph *grPoPSDChiSq = (TGraph*)f->Get("grPoPSDChiSq");
	TGraph *grRnEnChiSq  = (TGraph*)f->Get("grRnEnChiSq");
	TGraph *grPoEnChiSq  = (TGraph*)f->Get("grPoEnChiSq");
	TGraph *grDzChiSq    = (TGraph*)f->Get("grDzChiSq");
	TGraph *grDtChiSq    = (TGraph*)f->Get("grDtChiSq");

	TGraph *grBGRate = (TGraph*)f->Get("grBGRate");

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
	TH1F *hRate = new TH1F("hRate","Rate per cell",50,2.7,3.1);
	TH2F *h2DRate = new TH2F("h2DRate","Rate per cell",14,0,14,11,0,11);
	TH2F *h2DRelRate = new TH2F("h2DRelRate","Rate per cell",14,0,14,11,0,11);
	int numPt = grRate->GetN();
	double grx,gry;
	int binx,biny;
	for(int i=0;i<numPt;i++){
		grRate->GetPoint(i,grx,gry);
		binx = (int)grx%14 + 1;
		biny = ((int)grx/14) + 1;

		h2DRate->SetBinContent(binx,biny,gry);
		hRate->Fill(gry);

		grRelRate->GetPoint(i,grx,gry);
		h2DRelRate->SetBinContent(binx,biny,gry);
	}


	//-------------------------------------------------------------------------------------------------------
	const char *xLabel = "Cell";

	TCanvas *cRate = new TCanvas("cRate","Rate per cell",1000,400);
	grRate->GetXaxis()->SetTitle(xLabel);
	grRate->GetYaxis()->SetTitle("R_{RnPo} [mHz]");
	grRate->Draw("AP");
	cRate->SaveAs(Form("%s/RatePerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRate->SaveAs(Form("%s/RatePerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *chRate = new TCanvas("chRate","Rate per cell",1);
	hRate->GetXaxis()->SetTitle("R_{RnPo} [mHz]");
	hRate->GetYaxis()->SetTitle("Counts");
	hRate->Draw("HIST");
	chRate->SaveAs(Form("%s/HistRatePerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	chRate->SaveAs(Form("%s/HistRatePerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	gStyle->SetPaintTextFormat("2.2f");
	TCanvas *cRate2D = new TCanvas("cRate2D","Rate per cell",1);
	gPad->SetRightMargin(0.16);
	h2DRate->SetMarkerSize(1.25);
	h2DRate->SetMinimum(grRate->GetHistogram()->GetMinimum());
	h2DRate->Draw("colz && text");	
	cRate2D->SaveAs(Form("%s/2DRatePerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRate2D->SaveAs(Form("%s/2DRatePerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cBGRate = new TCanvas("cBGRate","BGRate per cell",1000,400);
	gPad->SetLogy();
	gPad->SetGrid();
	grBGRate->GetXaxis()->SetTitle(xLabel);
	grBGRate->GetYaxis()->SetTitle("BG Rate [mHz]");
	grBGRate->Draw("AP");
	cBGRate->SaveAs(Form("%s/BGRatePerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cBGRate->SaveAs(Form("%s/BGRatePerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelRate = new TCanvas("cRelRate","Relative rate per cell",1000,400);
	gPad->SetGrid();
	grRelRate->GetXaxis()->SetTitle(xLabel);
	grRelRate->GetYaxis()->SetTitle("R_{RnPo}/#LTR_{RnPo}#GT");  
	grRelRate->Draw("AP");
	grRelRate->Fit("pol0");
	grRelRate->GetFunction("pol0")->SetLineStyle(2);
	TPaveText *pvR = new TPaveText(0.85,0.8,0.99,0.99,"NDC");
	pvR->AddText(Form("#Chi^{2}/NDF   %.1f/%d",grRelRate->GetFunction("pol0")->GetChisquare(),grRelRate->GetFunction("pol0")->GetNDF()));	
	pvR->AddText(Form("Prob   %f",grRelRate->GetFunction("pol0")->GetProb()));
	pvR->AddText(Form("p0   %.3f #pm %.3f",grRelRate->GetFunction("pol0")->GetParameter(0),grRelRate->GetFunction("pol0")->GetParError(0)));
	pvR->Draw();
	cRelRate->SaveAs(Form("%s/RelativeRatePerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRate->SaveAs(Form("%s/RelativeRatePerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelRate2D = new TCanvas("cRelRate2D","Relative rate per cell",1);
	gPad->SetRightMargin(0.16);
	h2DRelRate->SetMarkerSize(1.25);
	h2DRelRate->SetMinimum(grRelRate->GetHistogram()->GetMinimum());
	h2DRelRate->Draw("colz && text");	
	cRelRate2D->SaveAs(Form("%s/2DRelRatePerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRate2D->SaveAs(Form("%s/2DRelRatePerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));


	TCanvas *cRelPoPSDMean = new TCanvas("cRelPoPSDMean","Relative Po PSD Mean",1000,400);
	grRelPoPSDMean->GetXaxis()->SetTitle(xLabel);
	grRelPoPSDMean->GetYaxis()->SetTitle("PSD_{Po}/#LTPSD_{Po}#GT");
	grRelPoPSDMean->Draw("AP");
	cRelPoPSDMean->SaveAs(Form("%s/RelativePoPSDMeanPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPSDMean->SaveAs(Form("%s/RelativePoPSDMeanPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoPSDSigma = new TCanvas("cRelPoPSDSigma","Relative Po PSD Sigma",1000,400);
	grRelPoPSDSigma->GetXaxis()->SetTitle(xLabel);
	grRelPoPSDSigma->GetYaxis()->SetTitle("#sigma_{PSD}/#LT#sigma_{PSD}#GT");
	grRelPoPSDSigma->Draw("AP");
	cRelPoPSDSigma->SaveAs(Form("%s/RelativePoPSDSigmaPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPSDSigma->SaveAs(Form("%s/RelativePoPSDSigmaPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoEnMean = new TCanvas("cRelPoEnMean","Relative Po Energy Mean",1000,400);
	grRelPoEnMean->GetXaxis()->SetTitle(xLabel);
	grRelPoEnMean->GetYaxis()->SetTitle("E_{Po}/#LTE_{Po}#GT");
	grRelPoEnMean->GetYaxis()->SetRangeUser(0.994,1.006);
	grRelPoEnMean->Draw("AP");
	cRelPoEnMean->SaveAs(Form("%s/RelativePoEnMeanPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoEnMean->SaveAs(Form("%s/RelativePoEnMeanPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoEnSigma = new TCanvas("cRelPoEnSigma","Relative Po Energy Sigma",1000,400);
	grRelPoEnSigma->GetXaxis()->SetTitle(xLabel);
	grRelPoEnSigma->GetYaxis()->SetTitle("#sigma_{E}/#LT#sigma_{E}#GT");
	grRelPoEnSigma->Draw("AP");
	cRelPoEnSigma->SaveAs(Form("%s/RelativePoEnSigmaPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoEnSigma->SaveAs(Form("%s/RelativePoEnSigmaPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoPosMean = new TCanvas("cRelPoPosMean","Relative Po Position Mean",1000,400);
	grRelPoPosMean->GetXaxis()->SetTitle(xLabel);
	grRelPoPosMean->GetYaxis()->SetTitle("z_{Po}/#LTz_{Po}#GT");
	grRelPoPosMean->Draw("AP");
	cRelPoPosMean->SaveAs(Form("%s/RelativePoPosMeanPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPosMean->SaveAs(Form("%s/RelativePoPosMeanPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoPosSigma = new TCanvas("cRelPoPosSigma","Relative Po Position Sigma",1000,400);
	grRelPoPosSigma->GetXaxis()->SetTitle(xLabel);
	grRelPoPosSigma->GetYaxis()->SetTitle("RMS_{z}/#LTRMS_{z}#GT");
	grRelPoPosSigma->Draw("AP");
	cRelPoPosSigma->SaveAs(Form("%s/RelativePoPosSigmaPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPosSigma->SaveAs(Form("%s/RelativePoPosSigmaPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelRnPoDzMean = new TCanvas("cRelRnPoDzMean","Relative RnPo Dz Mean",1000,400);
	grRelRnPoDzMean->GetXaxis()->SetTitle(xLabel);
	grRelRnPoDzMean->GetYaxis()->SetTitle("dz_{RnPo}/#LTdz_{RnPo}#GT");
	grRelRnPoDzMean->Draw("AP");
	cRelRnPoDzMean->SaveAs(Form("%s/RelativeRnPoDzMeanPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRnPoDzMean->SaveAs(Form("%s/RelativeRnPoDzMeanPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelRnPoDzSigma = new TCanvas("cRelRnPoDzSigma","Relative RnPo Dz Sigma",1000,400);
	grRelRnPoDzSigma->GetXaxis()->SetTitle(xLabel);
	grRelRnPoDzSigma->GetYaxis()->SetTitle("#sigma_{dz}/#LT#sigma_{dz}#GT");
	grRelRnPoDzSigma->Draw("AP");
	cRelRnPoDzSigma->SaveAs(Form("%s/RelativeRnPoDzSigmaPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRnPoDzSigma->SaveAs(Form("%s/RelativeRnPoDzSigmaPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cEff = new TCanvas("cEff","Efficiency",1000,400);
	gPad->SetGrid();
	grTotEff->GetXaxis()->SetTitle(xLabel);
	grTotEff->GetYaxis()->SetTitle("Efficiency");

	grRnPSDEff->SetMarkerStyle(26);
	grRnPSDEff->SetMarkerColor(kBlack);
	grRnPSDEff->SetLineColor(kBlack);
	grPoPSDEff->SetMarkerStyle(26);
	grPoPSDEff->SetMarkerColor(kRed);
	grPoPSDEff->SetLineColor(kRed);
	grRnEnEff->SetMarkerStyle(25);
	grRnEnEff->SetMarkerColor(kBlack);
	grRnEnEff->SetLineColor(kBlack);
	grPoEnEff->SetMarkerStyle(25);
	grPoEnEff->SetMarkerColor(kRed);
	grPoEnEff->SetLineColor(kRed);
	grRnPoDzEff->SetMarkerStyle(27);
	grRnPoDzEff->SetMarkerColor(8);
	grRnPoDzEff->SetLineColor(8);

	grTotEff->Draw("AP");
	grRnPSDEff->Draw("P");
	grPoPSDEff->Draw("P");
	grRnEnEff->Draw("P");
	grPoEnEff->Draw("P");
	grRnPoDzEff->Draw("P");

	TLegend *leg = new TLegend(0.90,0.67,0.99,0.99);
	leg->AddEntry(grRnPSDEff,"Rn PSD","p");
	leg->AddEntry(grPoPSDEff,"Po PSD","p");
	leg->AddEntry(grRnEnEff,"Rn E","p");
	leg->AddEntry(grPoEnEff,"Po E","p");
	leg->AddEntry(grRnPoDzEff,"Dz","p");
	leg->AddEntry(grTotEff,"Total","p");
	leg->Draw();	

	cEff->SaveAs(Form("%s/EfficiencyPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cEff->SaveAs(Form("%s/EfficiencyPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	//--------------------------------------------------
	double grLifetimexStart, grLifetimexEnd, grLifetimeyStart, grLifetimeyEnd;
	grLifetime->GetPoint(0,grLifetimexStart,grLifetimeyStart);
	grLifetime->GetPoint(grLifetime->GetN()-1,grLifetimexEnd,grLifetimeyEnd);
	grLifetimexStart = grLifetimexStart - 3;
	grLifetimexEnd = grLifetimexEnd + 3;

	TLine *lLifetime = new TLine(grLifetimexStart,2.569,grLifetimexEnd,2.569);
	lLifetime->SetLineStyle(2);
	lLifetime->SetLineColor(30);

	TGraph *grlLifetime = new TGraph(4);
	grlLifetime->SetPoint(0,grLifetimexStart,2.562);
	grlLifetime->SetPoint(1,grLifetimexEnd,2.562);
	grlLifetime->SetPoint(2,grLifetimexEnd,2.576);
	grlLifetime->SetPoint(3,grLifetimexStart,2.576);
	grlLifetime->SetFillStyle(3001);
	grlLifetime->SetFillColor(30);
	
	TCanvas *cLifetime = new TCanvas("cLifetime","Lifetime",1000,400);
	gPad->SetGrid();
	grLifetime->GetXaxis()->SetTitle(xLabel);
	grLifetime->GetYaxis()->SetTitle("#tau_{Po} [ms]");	
	grLifetime->Draw("AP");
	lLifetime->Draw("same");
	grlLifetime->Draw("f");
	grLifetime->Draw("P");
	grLifetime->Fit("pol0");
	grLifetime->GetFunction("pol0")->SetLineStyle(2);
	TPaveText *pv = new TPaveText(0.85,0.8,0.99,0.99,"NDC");
	pv->AddText(Form("#Chi^{2}/NDF   %.1f/%d",grLifetime->GetFunction("pol0")->GetChisquare(),grLifetime->GetFunction("pol0")->GetNDF()));	
	pv->AddText(Form("Prob   %f",grLifetime->GetFunction("pol0")->GetProb()));
	pv->AddText(Form("p0   %.3f #pm %.3f",grLifetime->GetFunction("pol0")->GetParameter(0),grLifetime->GetFunction("pol0")->GetParError(0)));
	pv->Draw();
	cLifetime->SaveAs(Form("%s/LifetimePerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cLifetime->SaveAs(Form("%s/LifetimePerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	grRnPSDChiSq->SetMarkerStyle(22);
	grRnPSDChiSq->SetMarkerSize(1.1);
	grRnPSDChiSq->SetMarkerColor(kBlack);
	grRnPSDChiSq->SetLineColor(kBlack);
	grPoPSDChiSq->SetMarkerStyle(22);
	grPoPSDChiSq->SetMarkerSize(1.1);
	grPoPSDChiSq->SetMarkerColor(kRed);
	grPoPSDChiSq->SetLineColor(kRed);

	grRnEnChiSq->SetMarkerStyle(34);
	grRnEnChiSq->SetMarkerSize(1.1);
	grRnEnChiSq->SetMarkerColor(kBlack);
	grRnEnChiSq->SetLineColor(kBlack);
	grPoEnChiSq->SetMarkerStyle(34);
	grPoEnChiSq->SetMarkerSize(1.1);
	grPoEnChiSq->SetMarkerColor(kRed);
	grPoEnChiSq->SetLineColor(kRed);

	grDzChiSq->SetMarkerStyle(20);
	grDzChiSq->SetMarkerSize(1.1);
	grDzChiSq->SetMarkerColor(8);
	grDzChiSq->SetLineColor(8);

	grDtChiSq->SetMarkerStyle(21);
	grDtChiSq->SetMarkerSize(1.1);
	
	TCanvas *cFitChiSq = new TCanvas("cFitChiSq","Fit Chisquared",1000,400);
	gPad->SetGrid();
	grRnPSDChiSq->GetXaxis()->SetTitle(xLabel);
	grRnPSDChiSq->GetYaxis()->SetTitle("#Chi^{2}/NDF");
	grRnPSDChiSq->GetYaxis()->SetRangeUser(0,4);
	grRnPSDChiSq->Draw("AP");
	grPoPSDChiSq->Draw("P");
	grRnEnChiSq->Draw("P");
	grPoEnChiSq->Draw("P");
	grDzChiSq->Draw("P");
	grDtChiSq->Draw("P");	
	leg = new TLegend(0.90,0.67,0.99,0.99);
	leg->AddEntry(grRnPSDChiSq,"Prompt PSD","p");
	leg->AddEntry(grPoPSDChiSq,"Delay PSD","p");
	leg->AddEntry(grRnEnChiSq,"Prompt E","p");
	leg->AddEntry(grPoEnChiSq,"Delay E","p");
	leg->AddEntry(grDzChiSq,"Dz","p");
	leg->AddEntry(grDtChiSq,"Dt","p");
	leg->Draw();	
	cFitChiSq->SaveAs(Form("%s/FitChiSqPerCell.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cFitChiSq->SaveAs(Form("%s/FitChiSqPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));



} 	//end PlotRnPoVsCell
