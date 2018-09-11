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
#include "Header.C"

using namespace std;

void PlotRnPoCellVsTime(){

	setup_PROSPECT_style();
    gROOT->ForceStyle();

	TGraphErrors *grRate;
	TF1 *fAcExp, *fAcExpFixed;
	bool exclude;
	int grRateNumPt;
	double grxStart, grxEnd, gryStart, gryEnd;
	
	int numPt = NUMCELLS - NUMEXCLUDECELLS;
	double x[numPt], yChiSq[numPt], yChiSqFixed[numPt];
	double yt[numPt], ytErr[numPt];
	double yR[numPt], yRErr[numPt], yRFixed[numPt], yRErrFixed[numPt];

	double ytChiSq[numPt], ytChiSqErr[numPt], ytChiSqFixed[numPt];

	//====================================================================
	TFile *f = new TFile(Form("%s/Ac227_GraphsCellvsTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));

	int pt = 0;
	for(int i=0;i<NUMCELLS;i++){
		exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), i) != end(ExcludeCellArr);
        if(exclude) continue;

		grRate = (TGraphErrors*)f->Get(Form("grRate_%i",i));	

		grRateNumPt = grRate->GetN();
		grRate->GetPoint(1,grxStart,gryStart);
		grRate->GetPoint(grRateNumPt-1,grxEnd,gryEnd);

		fAcExp = new TF1("fAcExp","[0]*exp(-((x-[2])*log(2)/([1]*365.0*24.0*60.0*60.0)))",grxStart,grxEnd);
		fAcExp->SetParameters(3.2,21.772);
		fAcExp->FixParameter(2,grxStart);
		fAcExp->SetParLimits(1,0,60);

		fAcExpFixed = new TF1("fAcExp","[0]*exp(-((x-[2])*log(2)/([1]*365.0*24.0*60.0*60.0)))",grxStart,grxEnd);
		fAcExpFixed->FixParameter(1,21.772);
		fAcExpFixed->FixParameter(2,grxStart);

		grRate->Fit(fAcExp,"RQ0");
		grRate->Fit(fAcExpFixed,"RQ0");

		x[pt] = i;

		yChiSq[pt] = (double)fAcExp->GetChisquare()/(double)fAcExp->GetNDF(); 
		yChiSqFixed[pt] = (double)fAcExpFixed->GetChisquare()/(double)fAcExpFixed->GetNDF(); 

		yt[pt] = fAcExp->GetParameter(1);
		ytErr[pt] = fAcExp->GetParError(1);
	
		yR[pt] = fAcExp->GetParameter(0);
		yRErr[pt] = fAcExp->GetParError(0);

		yRFixed[pt] = fAcExpFixed->GetParameter(0);
		yRErrFixed[pt] = fAcExpFixed->GetParError(0);

		ytChiSq[pt] = fAcExp->GetParameter(1)*((double)fAcExp->GetChisquare()/(double)fAcExp->GetNDF());
		ytChiSqErr[pt] = fAcExp->GetParError(1)*((double)fAcExp->GetChisquare()/(double)fAcExp->GetNDF());

		ytChiSqFixed[pt] = fAcExpFixed->GetParameter(1)*((double)fAcExpFixed->GetChisquare()/(double)fAcExpFixed->GetNDF());	
	
		pt++;

		TCanvas *c = new TCanvas("c","c",1000,400);
		grRate->GetYaxis()->SetTitle("R_{RnPo} [mHz]");
		grRate->GetXaxis()->SetTimeDisplay(1);
		grRate->GetXaxis()->SetTimeFormat("%m/%d");
		grRate->Draw("AP");
		fAcExp->Draw("same");
		fAcExpFixed->SetLineColor(8);
		fAcExpFixed->Draw("same");
		TPaveText *pv = new TPaveText(0.4,0.8,0.6,0.98,"NDC");
		pv->SetShadowColor(kRed);
		pv->AddText(Form("#Chi^{2}/NDF  %.1f/%d",fAcExp->GetChisquare(),fAcExp->GetNDF()));
		pv->AddText(Form("Prob  %f",fAcExp->GetProb()));
		pv->AddText(Form("R_{0}   %.2f #pm %.2f mHz",fAcExp->GetParameter(0),fAcExp->GetParError(0)));
		pv->AddText(Form("t_{1/2}   %.2f #pm %.2f yrs",fAcExp->GetParameter(1),fAcExp->GetParError(1)));
		pv->Draw();
		TPaveText *pvFix = new TPaveText(0.65,0.8,0.85,0.98,"NDC");
		pvFix->SetShadowColor(8);
		pvFix->AddText(Form("#Chi^{2}/NDF  %.1f/%d",fAcExpFixed->GetChisquare(),fAcExpFixed->GetNDF()));
		pvFix->AddText(Form("Prob  %f",fAcExpFixed->GetProb()));
		pvFix->AddText(Form("R_{0}   %.2f #pm %.2f mHz",fAcExpFixed->GetParameter(0),fAcExpFixed->GetParError(0)));
		pvFix->AddText(Form("t_{1/2}   %.2f #pm %.2f yrs",fAcExpFixed->GetParameter(1),fAcExpFixed->GetParError(1)));
		pvFix->Draw();
		c->SaveAs(Form("%s/CellVsTime/RateVsTime_Cell%03i.png",gSystem->Getenv("AD_AC227_PLOTS"),i));

		delete c;
		delete pv;
		delete pvFix;

	}	


	//====================================================================
	TGraph *grChiSq = new TGraphErrors(numPt,x,yChiSq);
	TGraph *grChiSqFixed = new TGraphErrors(numPt,x,yChiSqFixed);

	grChiSq->SetMarkerStyle(20);
	grChiSq->SetMarkerSize(0.9);
	grChiSq->SetMarkerColor(kBlack);

	grChiSqFixed->SetMarkerStyle(22);
	grChiSqFixed->SetMarkerSize(1.3);
	grChiSqFixed->SetMarkerColor(kRed);	

	TCanvas *cChiSq = new TCanvas("cChiSq","cChiSq",1000,400);
	grChiSq->Draw("AP");
	grChiSq->GetYaxis()->SetTitle("#Chi^{2}/NDF");
	grChiSq->GetXaxis()->SetTitle("Cell");
	cChiSq->Update();
	grChiSqFixed->Draw("P");
	TLegend *leg = new TLegend(0.8,0.85,0.95,0.95);
	leg->AddEntry(grChiSq,"Varying Lifetime");
	leg->AddEntry(grChiSqFixed,"Fixed Lifetime");
	leg->Draw();
	cChiSq->SaveAs(Form("%s/CellVsTime/ChiSqPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	TGraphErrors *grt = new TGraphErrors(numPt,x,yt,0,ytErr);
	grt->SetMarkerStyle(20);
	grt->SetMarkerSize(0.9);
	grt->SetMarkerColor(kBlack);
	
	TLine *lAc = new TLine(0,21.772,154,21.772);
	lAc->SetLineColor(8);
	lAc->SetLineStyle(2);
	lAc->SetLineWidth(2);

	TCanvas *cHalfLife = new TCanvas("cHalfLife","cHalfLife",1000,400);
	grt->Draw("AP");
	grt->GetYaxis()->SetTitle("t_{1/2} [yrs]");
	grt->GetXaxis()->SetTitle("Cell");
	grt->GetYaxis()->SetRangeUser(-5,50);
	lAc->Draw("same");
	cHalfLife->SaveAs(Form("%s/CellVsTime/HalfLifePerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	TGraphErrors *grR = new TGraphErrors(numPt,x,yR,0,yRErr);
	TGraphErrors *grRFixed = new TGraphErrors(numPt,x,yRFixed,0,yRErrFixed);

	grR->SetMarkerStyle(20);
	grR->SetMarkerSize(0.9);
	grR->SetMarkerColor(kBlack);
	grR->SetLineColor(kBlack);

	grRFixed->SetMarkerStyle(22);
	grRFixed->SetMarkerSize(1.3);
	grRFixed->SetMarkerColor(kRed);
	grRFixed->SetLineColor(kRed);

	TCanvas *cRate = new TCanvas("cRate","cRate",1000,400);
	grR->Draw("AP");
	grR->GetYaxis()->SetTitle("R_{0} [mHz]");
	grR->GetXaxis()->SetTitle("Cell");
	cRate->Update();
	grRFixed->Draw("P");
	leg = new TLegend(0.8,0.85,0.95,0.95);
	leg->AddEntry(grR,"Varying Lifetime");
	leg->AddEntry(grRFixed,"Fixed Lifetime");
	leg->Draw();
	cRate->SaveAs(Form("%s/CellVsTime/InitRatePerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));

		
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	TGraphErrors *grtChiSq = new TGraphErrors(numPt,x,ytChiSq,0,ytChiSqErr);
	TGraph *grtChiSqFixed = new TGraphErrors(numPt,x,ytChiSqFixed);
	
	grtChiSq->SetMarkerStyle(20);
	grtChiSq->SetMarkerSize(0.9);
	grtChiSq->SetMarkerColor(kBlack);
	grtChiSq->SetLineColor(kBlack);

	grtChiSqFixed->SetMarkerStyle(22);
	grtChiSqFixed->SetMarkerSize(1.3);
	grtChiSqFixed->SetMarkerColor(kRed);
	grtChiSqFixed->SetLineColor(kRed);

	TCanvas *ctChiSq = new TCanvas("ctChiSq","ctChiSq",1000,400);
	grtChiSq->Draw("AP");
	grtChiSq->GetYaxis()->SetTitle("t_{1/2}*(Chi^{2}/NDF) [yrs]");	
	grtChiSq->GetXaxis()->SetTitle("Cell");
	grtChiSq->GetYaxis()->SetRangeUser(-10,100);
	grtChiSqFixed->Draw("P");
	lAc->Draw("same");
	leg = new TLegend(0.8,0.85,0.95,0.95);
	leg->AddEntry(grtChiSq,"Varying Lifetime");
	leg->AddEntry(grtChiSqFixed,"Fixed Lifetime");
	leg->Draw();
	ctChiSq->SaveAs(Form("%s/CellVsTime/ChiSqTPerCell.png",gSystem->Getenv("AD_AC227_PLOTS")));



}	//end void

