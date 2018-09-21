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

TGraphErrors *makeCorrGr(TGraphErrors *gr, TF1 *fExp){
	TGraphErrors *grCorr = (TGraphErrors*)gr->Clone();
	int numPt = gr->GetN();

	double t0, N0;
	gr->GetPoint(0,t0,N0);

	double grx,gry,gryErr;
	double grCorry,grCorryErr;

	double lambda = log(2)/(21.772*365.0*24.0*60.0*60.0);
	

	for(int i=0;i<numPt;i++){
		gr->GetPoint(i,grx,gry);
		gryErr = gr->GetErrorY(i);	

		grCorry = gry/(N0*exp(-((grx-t0)*lambda)));
		grCorryErr = gryErr/(N0*exp(-((grx-t0)*lambda)));
	
		grCorr->SetPoint(i,grx,grCorry);
		grCorr->SetPointError(i,0,grCorryErr);	
	}

	return grCorr;
}


void PlotRnPoVsTime(){

	setup_PROSPECT_style();
    gROOT->ForceStyle();

	TFile *f = new TFile(Form("%s/Ac227_GraphsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));

	TGraphErrors *grRate 		= (TGraphErrors*)f->Get("grRate");
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

cout<<"\n POSITION RMS"<<endl;
	grPoPosSigma->Fit("pol0","0");
	TGraphErrors *grRelPoPosSigma = makeRelGr(grPoPosSigma, grPoPosSigma->GetFunction("pol0")->GetParameter(0), grPoPosSigma->GetFunction("pol0")->GetParError(0));

cout<<"\n DZ MEAN"<<endl;
	grRnPoDzMean->Fit("pol0","0");
	TGraphErrors *grRelRnPoDzMean = makeRelGr(grRnPoDzMean, grRnPoDzMean->GetFunction("pol0")->GetParameter(0), grRnPoDzMean->GetFunction("pol0")->GetParError(0));

cout<<"\n DZ SIGMA"<<endl;
	grRnPoDzSigma->Fit("pol0","0");
	TGraphErrors *grRelRnPoDzSigma = makeRelGr(grRnPoDzSigma, grRnPoDzSigma->GetFunction("pol0")->GetParameter(0), grRnPoDzSigma->GetFunction("pol0")->GetParError(0));
	
	//-------------------------------------------------------------------------------------------------------
	const char *xLabel = "Time";

	int numPt = grRate->GetN();
	double grxStart, grxEnd, gryStart, gryEnd;
	grRate->GetPoint(1,grxStart,gryStart);
	grRate->GetPoint(numPt-1,grxEnd,gryEnd);

	double May01_2018  = 1525168800;	//Rx On
	double May25_2018  = 1527285600;	//Rx Off
	double June12_2018 = 1528797600;	//Rx On	
	double July06_2018 = 1530914400;	//Rx Off
	double July24_2018 = 1532426400;	//Rx On
	double Aug17_2018  = 1534534200;	//Rx Off
	double Sept04_2018 = 1536055200;	//Rx On
	double Sept28_2018 = 1538172000; 	//Rx Off

	TGraph *grMayOn_2018  = new TGraph(4);
	TGraph *grJuneOn_2018 = new TGraph(4);
	TGraph *grJulyOn_2018 = new TGraph(4);
	TGraph *grSeptOn_2018 = new TGraph(4);

	grMayOn_2018->SetPoint(0,May01_2018,-50);
	grMayOn_2018->SetPoint(1,May01_2018,50);
	grMayOn_2018->SetPoint(2,May25_2018,50);
	grMayOn_2018->SetPoint(3,May25_2018,-50);
	
	grJuneOn_2018->SetPoint(0,June12_2018,-50);
	grJuneOn_2018->SetPoint(1,June12_2018,50);
	grJuneOn_2018->SetPoint(2,July06_2018,50);
	grJuneOn_2018->SetPoint(3,July06_2018,-50);

	grJulyOn_2018->SetPoint(0,July24_2018,-50);
	grJulyOn_2018->SetPoint(1,July24_2018,50);
	grJulyOn_2018->SetPoint(2,Aug17_2018,50);
	grJulyOn_2018->SetPoint(3,Aug17_2018,-50);

	grMayOn_2018->SetFillStyle(3002);
	grMayOn_2018->SetFillColor(16);
	grJuneOn_2018->SetFillStyle(3002);
	grJuneOn_2018->SetFillColor(16);
	grJulyOn_2018->SetFillStyle(3002);
	grJulyOn_2018->SetFillColor(16);

	TCanvas *cRate = new TCanvas("cRate","Rate vs time",1000,400);
	grRate->GetXaxis()->SetTitle(xLabel);
	grRate->GetYaxis()->SetTitle("R_{RnPo} [Hz]");
	grRate->GetXaxis()->SetTimeDisplay(1);
	grRate->GetXaxis()->SetTimeFormat("%m/%d");
	grRate->Draw("AP");
	grMayOn_2018->Draw("f");
	grJuneOn_2018->Draw("f");
	grJulyOn_2018->Draw("f");
	grRate->Draw("P");
	cRate->SaveAs(Form("%s/RateVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRate->SaveAs(Form("%s/RateVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));
/*
	TCanvas *cRelRate = new TCanvas("cRelRate","Relative rate vs time",1000,400);
	grRelRate->GetXaxis()->SetTitle(xLabel);
	grRelRate->GetYaxis()->SetTitle("R_{RnPo}/#LTR_{RnPo}#GT");  
	grRelRate->GetXaxis()->SetTimeDisplay(1);
	grRelRate->GetXaxis()->SetTimeFormat("%m/%d");
	grRelRate->Draw("AP");
	cRelRate->SaveAs(Form("%s/RelativeRateVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRate->SaveAs(Form("%s/RelativeRateVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));


	TF1 *fAcExp = new TF1("fAcExp","[0]*exp(-((x-[2])*log(2)/([1]*365.0*24.0*60.0*60.0)))",grxStart,grxEnd);
	fAcExp->SetParameters(2,21.772);
	fAcExp->FixParameter(2,grxStart);

	TF1 *fAcExpFixed = new TF1("fAcExp","[0]*exp(-((x-[2])*log(2)/([1]*365.0*24.0*60.0*60.0)))",grxStart,grxEnd);
	fAcExpFixed->FixParameter(1,21.772);
	fAcExpFixed->FixParameter(2,grxStart);


	TCanvas *cRateFit = new TCanvas("cRate","Rate vs time Fit",1000,400);
	grRate->GetXaxis()->SetTitle(xLabel);
	grRate->GetYaxis()->SetTitle("R_{RnPo} [Hz]");
	grRate->GetXaxis()->SetTimeDisplay(1);
	grRate->GetXaxis()->SetTimeFormat("%m/%d");
	grRate->Draw("AP");
	grRate->Fit(fAcExp,"0R");
	grRate->Fit(fAcExpFixed,"0R");
	fAcExp->Draw("same");
	fAcExpFixed->SetLineColor(8);
	fAcExpFixed->Draw("same");
	TPaveText *pv = new TPaveText(0.4,0.8,0.6,0.98,"NDC");
	pv->SetShadowColor(kRed);
	pv->AddText(Form("#Chi^{2}/NDF  %.1f/%d",fAcExp->GetChisquare(),fAcExp->GetNDF()));
	pv->AddText(Form("Prob  %f",fAcExp->GetProb()));
	pv->AddText(Form("R_{0}   %.2f #pm %.2f mHz",fAcExp->GetParameter(0)*1000.0,fAcExp->GetParError(0)*1000.0));
	pv->AddText(Form("t_{1/2}   %.2f #pm %.2f yrs",fAcExp->GetParameter(1),fAcExp->GetParError(1)));
	pv->Draw();
	TPaveText *pvFix = new TPaveText(0.65,0.8,0.85,0.98,"NDC");
	pvFix->SetShadowColor(8);
	pvFix->AddText(Form("#Chi^{2}/NDF  %.1f/%d",fAcExpFixed->GetChisquare(),fAcExpFixed->GetNDF()));
	pvFix->AddText(Form("Prob  %f",fAcExpFixed->GetProb()));
	pvFix->AddText(Form("R_{0}   %.2f #pm %.2f mHz",fAcExpFixed->GetParameter(0)*1000.0,fAcExpFixed->GetParError(0)*1000.0));
	pvFix->AddText(Form("t_{1/2}   %.2f #pm %.2f yrs",fAcExpFixed->GetParameter(1),fAcExpFixed->GetParError(1)));
	pvFix->Draw();
	cRateFit->SaveAs(Form("%s/RateVsTime_Fit.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRateFit->SaveAs(Form("%s/RateVsTime_Fit.png",gSystem->Getenv("AD_AC227_PLOTS")));


	TGraphErrors *grRateCorr = makeCorrGr(grRate,fAcExpFixed);
	TCanvas *cRateCorr = new TCanvas("cRateCorr","Corrected rate",1000,400);
	gPad->SetGrid();
	grRateCorr->GetXaxis()->SetTitle(xLabel);	
	grRateCorr->GetYaxis()->SetTitle("R_{RnPo}/R_{Ac}");
	grRateCorr->GetXaxis()->SetTimeDisplay(1);
	grRateCorr->GetXaxis()->SetTimeFormat("%m/%d");
	grRateCorr->Draw("AP");
	grRateCorr->Fit("pol0");
	TPaveText *pvPol = new TPaveText(0.8,0.8,0.99,0.99,"NDC");
	pvPol->AddText(Form("#Chi^{2}/NDF  %.1f/%d",grRateCorr->GetFunction("pol0")->GetChisquare(),grRateCorr->GetFunction("pol0")->GetNDF()));
	pvPol->AddText(Form("Prob  %f",grRateCorr->GetFunction("pol0")->GetProb()));
	pvPol->AddText(Form("p0  %.3f #pm %.3f",grRateCorr->GetFunction("pol0")->GetParameter(0),grRateCorr->GetFunction("pol0")->GetParError(0)));	
	pvPol->Draw();
	cRateCorr->SaveAs(Form("%s/RateVsTime_Corrected.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRateCorr->SaveAs(Form("%s/RateVsTime_Corrected.png",gSystem->Getenv("AD_AC227_PLOTS")));
	

	//======================================================================================

	TCanvas *cRelPoPSDMean = new TCanvas("cRelPoPSDMean","Relative Po PSD Mean",1000,400);
	grRelPoPSDMean->GetXaxis()->SetTitle(xLabel);
	grRelPoPSDMean->GetYaxis()->SetTitle("PSD_{#muPo}/#LTPSD_{#muPo}#GT");
	grRelPoPSDMean->GetXaxis()->SetTimeDisplay(1);
	grRelPoPSDMean->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoPSDMean->Draw("AP");
	cRelPoPSDMean->SaveAs(Form("%s/RelativePoPSDMeanVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPSDMean->SaveAs(Form("%s/RelativePoPSDMeanVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoPSDSigma = new TCanvas("cRelPoPSDSigma","Relative Po PSD Sigma",1000,400);
	grRelPoPSDSigma->GetXaxis()->SetTitle(xLabel);
	grRelPoPSDSigma->GetYaxis()->SetTitle("PSD_{#sigmaPo}/#LTPSD_{#sigmaPo}#GT");
	grRelPoPSDSigma->GetXaxis()->SetTimeDisplay(1);
	grRelPoPSDSigma->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoPSDSigma->Draw("AP");
	cRelPoPSDSigma->SaveAs(Form("%s/RelativePoPSDSigmaVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPSDSigma->SaveAs(Form("%s/RelativePoPSDSigmaVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoEnMean = new TCanvas("cRelPoEnMean","Relative Po Energy Mean",1000,400);
	grRelPoEnMean->GetXaxis()->SetTitle(xLabel);
	grRelPoEnMean->GetYaxis()->SetTitle("E_{#muPo}/#LTE_{#muPo}#GT");
	grRelPoEnMean->GetXaxis()->SetTimeDisplay(1);
	grRelPoEnMean->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoEnMean->Draw("AP");
	cRelPoEnMean->SaveAs(Form("%s/RelativePoEnMeanVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoEnMean->SaveAs(Form("%s/RelativePoEnMeanVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoEnSigma = new TCanvas("cRelPoEnSigma","Relative Po Energy Sigma",1000,400);
	grRelPoEnSigma->GetXaxis()->SetTitle(xLabel);
	grRelPoEnSigma->GetYaxis()->SetTitle("E_{#sigmaPo}/#LTE_{#sigmaPo}#GT");
	grRelPoEnSigma->GetXaxis()->SetTimeDisplay(1);
	grRelPoEnSigma->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoEnSigma->Draw("AP");
	cRelPoEnSigma->SaveAs(Form("%s/RelativePoEnSigmaVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoEnSigma->SaveAs(Form("%s/RelativePoEnSigmaVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoPosMean = new TCanvas("cRelPoPosMean","Relative Po Position Mean",1000,400);
	grRelPoPosMean->GetXaxis()->SetTitle(xLabel);
	grRelPoPosMean->GetYaxis()->SetTitle("z_{#muPo}/#LTz_{#muPo}#GT");
	grRelPoPosMean->GetXaxis()->SetTimeDisplay(1);
	grRelPoPosMean->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoPosMean->Draw("AP");
	cRelPoPosMean->SaveAs(Form("%s/RelativePoPosMeanVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPosMean->SaveAs(Form("%s/RelativePoPosMeanVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoPosSigma = new TCanvas("cRelPoPosSigma","Relative Po Position Sigma",1000,400);
	grRelPoPosSigma->GetXaxis()->SetTitle(xLabel);
	grRelPoPosSigma->GetYaxis()->SetTitle("z_{RMS Po}/#LTz_{RMS Po}#GT");
	grRelPoPosSigma->GetXaxis()->SetTimeDisplay(1);
	grRelPoPosSigma->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoPosSigma->Draw("AP");
	cRelPoPosSigma->SaveAs(Form("%s/RelativePoPosSigmaVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPosSigma->SaveAs(Form("%s/RelativePoPosSigmaVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelRnPoDzMean = new TCanvas("cRelRnPoDzMean","Relative RnPo Dz Mean",1000,400);
	grRelRnPoDzMean->GetXaxis()->SetTitle(xLabel);
	grRelRnPoDzMean->GetYaxis()->SetTitle("dz_{#mu}/#LTdz_{#mu}#GT");
	grRelRnPoDzMean->GetXaxis()->SetTimeDisplay(1);
	grRelRnPoDzMean->GetXaxis()->SetTimeFormat("%m/%d");
	grRelRnPoDzMean->Draw("AP");
	cRelRnPoDzMean->SaveAs(Form("%s/RelativeRnPoDzMeanVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRnPoDzMean->SaveAs(Form("%s/RelativeRnPoDzMeanVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelRnPoDzSigma = new TCanvas("cRelRnPoDzSigma","Relative RnPo Dz Sigma",1000,400);
	grRelRnPoDzSigma->GetXaxis()->SetTitle(xLabel);
	grRelRnPoDzSigma->GetYaxis()->SetTitle("dz_{#sigma}/#LTdz_{#sigma}#GT");
	grRelRnPoDzSigma->GetXaxis()->SetTimeDisplay(1);
	grRelRnPoDzSigma->GetXaxis()->SetTimeFormat("%m/%d");
	grRelRnPoDzSigma->Draw("AP");
	cRelRnPoDzSigma->SaveAs(Form("%s/RelativeRnPoDzSigmaVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRnPoDzSigma->SaveAs(Form("%s/RelativeRnPoDzSigmaVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cEff = new TCanvas("cEff","Total efficiency",1000,400);
	grTotEff->GetXaxis()->SetTitle(xLabel);
	grTotEff->GetYaxis()->SetTitle("Efficiency");
	grTotEff->GetXaxis()->SetTimeDisplay(1);
	grTotEff->GetXaxis()->SetTimeFormat("%m/%d");
	grTotEff->Draw("AP");
	cEff->SaveAs(Form("%s/EfficiencyVsTime.C",gSystem->Getenv("AD_AC227_PLOTS")));
	cEff->SaveAs(Form("%s/EfficiencyVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

*/
} 	//end PlotRnPoVsTime
