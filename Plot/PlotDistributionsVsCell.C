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

TH1F *makeResidHist(TH1F *h, TF1 *fit){
	TH1F *hResid = (TH1F*)h->Clone();
	int numPt = h->GetSize() - 2;

	double hx, hy, hyErr;
	double fity;
	double resid, residErr;

	for(int i=1;i<numPt+1;i++){
		hx = h->GetXaxis()->GetBinCenter(i);
		hy = h->GetBinContent(i);	
		hyErr = h->GetBinError(i);

		fity = fit->Eval(hx);
		
		resid = hy-fity;
		residErr = hyErr;

		hResid->SetBinContent(i,resid);
		hResid->SetBinError(i,residErr);
	}

	hResid->GetYaxis()->SetTitle("Data - fit");
	hResid->GetYaxis()->SetLabelSize(0.07);
	hResid->GetYaxis()->SetTitleSize(0.09);
	hResid->GetYaxis()->SetTitleOffset(0.4);
	hResid->GetXaxis()->SetLabelSize(0.07);
	hResid->GetXaxis()->SetTitleSize(0.06);
	hResid->GetXaxis()->SetTitleOffset(0.7);


	return hResid;
}


int PlotDistributionsVsCell(int cellNum){

	setup_PROSPECT_style();
  	gROOT->ForceStyle();

	if(find(begin(ExcludeCellArr), end(ExcludeCellArr), cellNum) != end(ExcludeCellArr)){
		printf("No information for cell %i. Exiting. \n",cellNum);
		return -1;
	}

	int numLower = 0.0;
	for(int i=0;i<NUMEXCLUDECELLS;i++){
		int v = ExcludeCellArr[i];
		if(v<cellNum) numLower++;
		if(v>cellNum) break;
	}	

	int grIDX = cellNum - numLower;

	TFile *fg = new TFile(Form("%s/Ac227_GraphsPerCell.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS"))); 
	if(!fg){
		printf("Graph file not found. Exiting. \n");
		return -1;
	}
	TGraphErrors *grRnPSDEff    = (TGraphErrors*)fg->Get("grRnPSDEff");
    	TGraphErrors *grPoPSDEff    = (TGraphErrors*)fg->Get("grPoPSDEff");
    	TGraphErrors *grRnEnEff     = (TGraphErrors*)fg->Get("grRnEnEff");
    	TGraphErrors *grPoEnEff     = (TGraphErrors*)fg->Get("grPoEnEff");
    	TGraphErrors *grRnPoDzEff   = (TGraphErrors*)fg->Get("grRnPoDzEff");
	fg->Close();

	cout<<grIDX<<endl;
	cout<<grRnPSDEff->GetX()[grIDX]<<endl;
	//-------------------------------------------------------------------------------------------------------
	TFile *f = new TFile(Form("%s/Ac227_HistsPerCell.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));
	if(!f){
		printf("File not found. Exiting. \n");
		return -1;
	}


	TH1F *hSelectDt = (TH1F*)f->Get(Form("hSelectDt_%i",cellNum));
	TH1F *hBGDt 	= (TH1F*)f->Get(Form("hBGDt_%i",cellNum));
	TH1F *hRnPoDt 	= (TH1F*)f->Get(Form("hRnPoDt_%i",cellNum));

	TF1  *fRnPoDtExp = hRnPoDt->GetFunction("fRnPoDtExp");
	fRnPoDtExp->SetRange(0,2.57*5.0);

	TH1F *hSelectPromptPSD = (TH1F*)f->Get(Form("hSelectPromptPSD_%i",cellNum));
	TH1F *hBGPromptPSD 	   = (TH1F*)f->Get(Form("hBGPromptPSD_%i",cellNum));
	TH1F *hRnPSD 	 	   = (TH1F*)f->Get(Form("hRnPSD_%i",cellNum));

	TF1  *fRnPSDGaus = hRnPSD->GetFunction("fRnPSDGaus");
	fRnPSDGaus->SetRange(0.15,0.37);

	TH1F *hSelectDelayPSD = (TH1F*)f->Get(Form("hSelectDelayPSD_%i",cellNum));
	TH1F *hBGDelayPSD 	  = (TH1F*)f->Get(Form("hBGDelayPSD_%i",cellNum));
	TH1F *hPoPSD 		  = (TH1F*)f->Get(Form("hPoPSD_%i",cellNum));

	TF1  *fPoPSDGaus = hPoPSD->GetFunction("fPoPSDGaus");
	fPoPSDGaus->SetRange(0.15,0.37);

	TH1F *hSelectPromptEn = (TH1F*)f->Get(Form("hSelectPromptEn_%i",cellNum));
	TH1F *hBGPromptEn 	  = (TH1F*)f->Get(Form("hBGPromptEn_%i",cellNum));
	TH1F *hRnEn 		  = (TH1F*)f->Get(Form("hRnEn_%i",cellNum));

	TF1  *fRnEnGaus = hRnEn->GetFunction("fRnEnGaus");
	fRnEnGaus->SetRange(0.49,1.16);

	TH1F *hSelectDelayEn = (TH1F*)f->Get(Form("hSelectDelayEn_%i",cellNum));
	TH1F *hBGDelayEn 	 = (TH1F*)f->Get(Form("hBGDelayEn_%i",cellNum));
	TH1F *hPoEn 		 = (TH1F*)f->Get(Form("hPoEn_%i",cellNum));

	TF1  *fPoEnGaus = hPoEn->GetFunction("fPoEnGaus");
	fPoEnGaus->SetRange(0.49,1.16);

	TH1F *hSelectPromptPos = (TH1F*)f->Get(Form("hSelectPromptPos_%i",cellNum));
	TH1F *hBGPromptPos 	   = (TH1F*)f->Get(Form("hBGPromptPos_%i",cellNum));
	TH1F *hRnPos 		   = (TH1F*)f->Get(Form("hRnPos_%i",cellNum));

	TH1F *hSelectDelayPos = (TH1F*)f->Get(Form("hSelectDelayPos_%i",cellNum));
	TH1F *hBGDelayPos 	  = (TH1F*)f->Get(Form("hBGDelayPos_%i",cellNum));
	TH1F *hPoPos 		  = (TH1F*)f->Get(Form("hPoPos_%i",cellNum));

	TH1F *hSelectDz = (TH1F*)f->Get(Form("hSelectDz_%i",cellNum));
	TH1F *hBGDz 	= (TH1F*)f->Get(Form("hBGDz_%i",cellNum));
	TH1F *hRnPoDz 	= (TH1F*)f->Get(Form("hRnPoDz_%i",cellNum));

	TF1  *fRnPoDzGaus = hRnPoDz->GetFunction("fRnPoDzGaus");

	TH2F *hSelectPSDvsEn = (TH2F*)f->Get(Form("hSelectPSDvsEn_%i",cellNum));
	TH2F *hBGPSDvsEn 	 = (TH2F*)f->Get(Form("hBGPSDvsEn_%i",cellNum));
	TH2F *hRnPoPSDvsEn 	 = (TH2F*)f->Get(Form("hRnPoPSDvsEn_%i",cellNum));

	TH2F *hSelectPSDvsPos = (TH2F*)f->Get(Form("hSelectPSDvsPos_%i",cellNum));
	TH2F *hBGPSDvsPos 	  = (TH2F*)f->Get(Form("hBGPSDvsPos_%i",cellNum));
	TH2F *hRnPoPSDvsPos   = (TH2F*)f->Get(Form("hRnPoPSDvsPos_%i",cellNum));

	TH2F *hSelectEnvsPos = (TH2F*)f->Get(Form("hSelectEnvsPos_%i",cellNum));
	TH2F *hBGEnvsPos 	 = (TH2F*)f->Get(Form("hBGEnvsPos_%i",cellNum));
	TH2F *hRnPoEnvsPos   = (TH2F*)f->Get(Form("hRnPoEnvsPos_%i",cellNum));

	TH2F *hSelectDelayEnvsPromptEn = (TH2F*)f->Get(Form("hSelectDelayEnvsPromptEn_%i",cellNum));
	TH2F *hBGDelayEnvsPromptEn 	   = (TH2F*)f->Get(Form("hBGDelayEnvsPromptEn_%i",cellNum));
	TH2F *hPoEnvsRnEn 			   = (TH2F*)f->Get(Form("hPoEnvsRnEn_%i",cellNum));

	//-------------------------------------------------------------------------------------------------------
	TPaveText *pt;
	int p_col = 1, d_col = 4;
	int s_col = 8, bg_col = 1;

	TCanvas *cRnPoDt = new TCanvas("cRnPoDt","RnPo Dt",700,700);
	TPad *pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0,0);
	TPad *pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3,0);
	pad1->Draw();
	pad2->Draw();

	pad1->cd();
	gPad->SetGrid();
	gPad->SetLogy();
	hSelectDt->SetLineColor(8);
	hSelectDt->SetLineWidth(1);
	hSelectDt->SetMinimum(0.1);
	hSelectDt->Draw("HIST");
	hBGDt->SetLineColor(kBlack);
	hBGDt->SetLineWidth(1);
	hBGDt->Draw("same&HIST");
	hRnPoDt->SetLineWidth(1);
	hRnPoDt->Draw("same");
	hRnPoDt->GetXaxis()->SetTitle("t_{#alphaPo} - t_{#alphaRn} [ms]");
	fRnPoDtExp->SetLineStyle(2);
	fRnPoDtExp->SetRange(0,12.85);
	fRnPoDtExp->Draw("same");
	pt = new TPaveText(0.68,0.63,0.95,0.9,"NDCNB");
	pt->AddText(Form("Entries    %.0f",hRnPoDt->GetEntries()));
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fRnPoDtExp->GetChisquare(),fRnPoDtExp->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fRnPoDtExp->GetProb()));
	pt->AddText(Form("N    %.2f #pm %.2f",fRnPoDtExp->GetParameter(0),fRnPoDtExp->GetParError(0)));
	pt->AddText(Form("#tau_{Po^{215}}    %.2f #pm %.2f ms",fRnPoDtExp->GetParameter(1),fRnPoDtExp->GetParError(1)));
	pt->Draw();

	pad2->cd();
	gPad->SetGrid();
	TH1F *hRnPoDtResid = makeResidHist(hRnPoDt,fRnPoDtExp);
	hRnPoDtResid->GetXaxis()->SetTitle("dt [ms]");
	hRnPoDtResid->GetYaxis()->SetRangeUser(-50,50);
	hRnPoDtResid->Draw();
	hRnPoDtResid->Fit("pol0");
	hRnPoDtResid->GetFunction("pol0")->SetLineStyle(2);
	pt = new TPaveText(0.8,0.8,0.99,0.99,"NDCNB");
	pt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",hRnPoDtResid->GetFunction("pol0")->GetChisquare(),hRnPoDtResid->GetFunction("pol0")->GetNDF()));
	pt->AddText(Form("Prob   %.2f",hRnPoDtResid->GetFunction("pol0")->GetProb()));
	pt->AddText(Form("p0   %.2f #pm %.2f",hRnPoDtResid->GetFunction("pol0")->GetParameter(0),hRnPoDtResid->GetFunction("pol0")->GetParError(0)));
	pt->Draw();	

	cRnPoDt->SaveAs(Form("%s/RnPoDt_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cRnPoDt->SaveAs(Form("%s/RnPoDt_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));

	//----------------------------------------------------------
	TF1 *fRnPSDGaus0 = new TF1("fRnPSDGaus0","gaus",PSDMin,PSDMax);
	TF1 *fRnPSDGaus1 = new TF1("fRnPSDGaus1","gaus",PSDMin,PSDMax);

	fRnPSDGaus0->SetParameters(fRnPSDGaus->GetParameter(0),fRnPSDGaus->GetParameter(1),fRnPSDGaus->GetParameter(2));
	fRnPSDGaus1->SetParameters(fRnPSDGaus->GetParameter(3),fRnPSDGaus->GetParameter(4),fRnPSDGaus->GetParameter(5));

	fRnPSDGaus0->SetLineColor(8);
	fRnPSDGaus1->SetLineColor(kMagenta);

	TCanvas *cRnPSD = new TCanvas("cRnPSD","Rn PSD",700,700);
	pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0,0);
	pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3,0);
	pad1->Draw();
	pad2->Draw();

	pad1->cd();
	gPad->SetGrid();
	hSelectPromptPSD->SetLineColor(s_col);
	hSelectPromptPSD->Draw("HIST");
	hBGPromptPSD->SetLineColor(bg_col);
	hBGPromptPSD->Draw("HIST&SAME");
	hRnPSD->Draw("same");
	fRnPSDGaus->SetLineStyle(2);
	fRnPSDGaus->Draw("same");
	fRnPSDGaus0->Draw("same");
	fRnPSDGaus1->Draw("same");
	pt = new TPaveText(0.75,0.75,0.95,0.9,"NDCNB");
	TText *tRn = pt->AddText("Rn^{219}");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fRnPSDGaus->GetChisquare(),fRnPSDGaus->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fRnPSDGaus->GetProb()));
	pt->AddText(Form("Eff    %.3f #pm %.3f",grRnPSDEff->GetY()[grIDX],grRnPSDEff->GetEY()[grIDX]));
	pt->Draw();

	pad2->cd();
	gPad->SetGrid();
	TH1F *hRnPSDResid = makeResidHist(hRnPSD,fRnPSDGaus);
	hRnPSDResid->GetYaxis()->SetRangeUser(-70,70);
	hRnPSDResid->Draw();
	hRnPSDResid->Fit("pol0");
	hRnPSDResid->GetFunction("pol0")->SetLineStyle(2);
	pt = new TPaveText(0.8,0.8,0.99,0.99,"NDCNB");
	pt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",hRnPSDResid->GetFunction("pol0")->GetChisquare(),hRnPSDResid->GetFunction("pol0")->GetNDF()));
	pt->AddText(Form("Prob   %.2f",hRnPSDResid->GetFunction("pol0")->GetProb()));
	pt->AddText(Form("p0   %.2f #pm %.2f",hRnPSDResid->GetFunction("pol0")->GetParameter(0),hRnPSDResid->GetFunction("pol0")->GetParError(0)));
	pt->Draw();	

	cRnPSD->SaveAs(Form("%s/RnPSD_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cRnPSD->SaveAs(Form("%s/RnPSD_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));

	//----------------------------------------------------------
	TCanvas *cPoPSD = new TCanvas("cPoPSD","Po PSD",700,700);
	pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0,0);
	pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3,0);
	pad1->Draw();
	pad2->Draw();

	pad1->cd();
	gPad->SetGrid();
	hSelectDelayPSD->SetLineColor(s_col);
	hSelectDelayPSD->Draw("HIST");
	hBGDelayPSD->SetLineColor(bg_col);
	hBGDelayPSD->Draw("HIST&SAME");
	hPoPSD->Draw("same");
	fPoPSDGaus->SetLineStyle(2);
	fPoPSDGaus->Draw("same");
	pt = new TPaveText(0.75,0.75,0.95,0.9,"NDCNB");
	TText *tPo = pt->AddText("Po^{215}");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fPoPSDGaus->GetChisquare(),fPoPSDGaus->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fPoPSDGaus->GetProb()));
	pt->AddText(Form("Eff    %.3f #pm %.3f",grPoPSDEff->GetY()[grIDX],grPoPSDEff->GetEY()[grIDX]));
	pt->Draw();

	pad2->cd();
	gPad->SetGrid();
	TH1F *hPoPSDResid = makeResidHist(hPoPSD,fPoPSDGaus);
	hPoPSDResid->GetYaxis()->SetRangeUser(-70,70);
	hPoPSDResid->Draw();
	hPoPSDResid->Fit("pol0");
	hPoPSDResid->GetFunction("pol0")->SetLineStyle(2);
	pt = new TPaveText(0.8,0.8,0.99,0.99,"NDCNB");
	pt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",hPoPSDResid->GetFunction("pol0")->GetChisquare(),hPoPSDResid->GetFunction("pol0")->GetNDF()));
	pt->AddText(Form("Prob   %.2f",hPoPSDResid->GetFunction("pol0")->GetProb()));
	pt->AddText(Form("p0   %.2f #pm %.2f",hPoPSDResid->GetFunction("pol0")->GetParameter(0),hPoPSDResid->GetFunction("pol0")->GetParError(0)));
	pt->Draw();	

	cPoPSD->SaveAs(Form("%s/PoPSD_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cPoPSD->SaveAs(Form("%s/PoPSD_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));

	//----------------------------------------------------------
	TF1 *fRnEnGaus0 = new TF1("fRnEnGaus0","gaus",EnMin,EnMax);
	TF1 *fRnEnGaus1 = new TF1("fRnEnGaus1","gaus",EnMin,EnMax);
	TF1 *fRnEnGaus2 = new TF1("fRnEnGaus2","gaus",EnMin,EnMax);

	fRnEnGaus0->SetParameters(fRnEnGaus->GetParameter(0),fRnEnGaus->GetParameter(1),fRnEnGaus->GetParameter(2));
	fRnEnGaus1->SetParameters(fRnEnGaus->GetParameter(3),fRnEnGaus->GetParameter(4),fRnEnGaus->GetParameter(5));
	fRnEnGaus2->SetParameters(fRnEnGaus->GetParameter(6),fRnEnGaus->GetParameter(7),fRnEnGaus->GetParameter(8));

	fRnEnGaus0->SetLineColor(8);
	fRnEnGaus1->SetLineColor(kMagenta);
	fRnEnGaus2->SetLineColor(kViolet);

	TCanvas *cRnEn = new TCanvas("cRnEn","Rn En",700,700);
	pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0,0);
	pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3,0);
	pad1->Draw();
	pad2->Draw();

	pad1->cd();
	gPad->SetGrid();

	hSelectPromptEn->SetLineColor(s_col);
	hSelectPromptEn->Draw("HIST");
	hBGPromptEn->SetLineColor(bg_col);
	hBGPromptEn->Draw("HIST&SAME");
	hRnEn->Draw("same");
	fRnEnGaus->SetLineStyle(2);	
	fRnEnGaus->Draw("same");
	fRnEnGaus0->Draw("same");
	fRnEnGaus1->Draw("same");
	fRnEnGaus2->Draw("same");
	pt = new TPaveText(0.75,0.75,0.95,0.9,"NDCNB");
	tRn = pt->AddText("Rn^{219}");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fRnEnGaus->GetChisquare(),fRnEnGaus->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fRnEnGaus->GetProb()));
	pt->AddText(Form("Eff    %.3f #pm %.3f",grRnEnEff->GetY()[grIDX],grRnEnEff->GetEY()[grIDX]));
	pt->Draw();

	pad2->cd();
	gPad->SetGrid();
	TH1F *hRnEnResid = makeResidHist(hRnEn,fRnEnGaus);
	hRnEnResid->GetYaxis()->SetRangeUser(-70,70);
	hRnEnResid->Draw();
	hRnEnResid->Fit("pol0");
	hRnEnResid->GetFunction("pol0")->SetLineStyle(2);
	pt = new TPaveText(0.8,0.8,0.99,0.99,"NDCNB");
	pt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",hRnEnResid->GetFunction("pol0")->GetChisquare(),hRnEnResid->GetFunction("pol0")->GetNDF()));
	pt->AddText(Form("Prob   %.2f",hRnEnResid->GetFunction("pol0")->GetProb()));
	pt->AddText(Form("p0   %.2f #pm %.2f",hRnEnResid->GetFunction("pol0")->GetParameter(0),hRnEnResid->GetFunction("pol0")->GetParError(0)));
	pt->Draw();	

	cRnEn->SaveAs(Form("%s/RnEn_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cRnEn->SaveAs(Form("%s/RnEn_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));

	//----------------------------------------------------------
	TCanvas *cPoEn = new TCanvas("cPoEn","Po En",700,700);
	pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0,0);
	pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3,0);
	pad1->Draw();
	pad2->Draw();

	pad1->cd();
	gPad->SetGrid();
	hSelectDelayEn->SetLineColor(s_col);
	hSelectDelayEn->Draw("HIST");
	hBGDelayEn->SetLineColor(bg_col);
	hBGDelayEn->Draw("HIST&SAME");
	hPoEn->Draw("same");
	fPoEnGaus->SetLineStyle(2);	
	fPoEnGaus->Draw("same");
	pt = new TPaveText(0.75,0.75,0.95,0.9,"NDCNB");
	tPo = pt->AddText("Po^{215}");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fPoEnGaus->GetChisquare(),fPoEnGaus->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fPoEnGaus->GetProb()));
	pt->AddText(Form("Eff    %.3f #pm %.3f",grPoEnEff->GetY()[grIDX],grPoEnEff->GetEY()[grIDX]));
	pt->Draw();

	pad2->cd();
	gPad->SetGrid();
	TH1F *hPoEnResid = makeResidHist(hPoEn,fPoEnGaus);
	hPoEnResid->GetYaxis()->SetRangeUser(-70,70);
	hPoEnResid->Draw();
	hPoEnResid->Fit("pol0");
	hPoEnResid->GetFunction("pol0")->SetLineStyle(2);
	pt = new TPaveText(0.8,0.8,0.99,0.99,"NDCNB");
	pt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",hPoEnResid->GetFunction("pol0")->GetChisquare(),hPoEnResid->GetFunction("pol0")->GetNDF()));
	pt->AddText(Form("Prob   %.2f",hPoEnResid->GetFunction("pol0")->GetProb()));
	pt->AddText(Form("p0   %.2f #pm %.2f",hPoEnResid->GetFunction("pol0")->GetParameter(0),hPoEnResid->GetFunction("pol0")->GetParError(0)));
	pt->Draw();	

	cPoEn->SaveAs(Form("%s/PoEn_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cPoEn->SaveAs(Form("%s/PoEn_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));

	//----------------------------------------------------------
	TCanvas *cRnPoDz = new TCanvas("cRnPoDz","RnPo Dz",700,700);
	pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0,0);
	pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3,0);
	pad1->Draw();
	pad2->Draw();

	pad1->cd();
	gPad->SetGrid();
	hSelectDz->SetLineColor(s_col);
	hSelectDz->Draw("HIST");
	hBGDz->SetLineColor(bg_col);
	hBGDz->Draw("HIST&SAME");
	hRnPoDz->GetXaxis()->SetTitle("z_{#alphaPo} - z_{#alphaRn} [mm]");
	hRnPoDz->Draw("same");
	fRnPoDzGaus->SetLineStyle(2);
	fRnPoDzGaus->Draw("same");	
	pt = new TPaveText(0.75,0.75,0.95,0.9,"NDCNB");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fRnPoDzGaus->GetChisquare(),fRnPoDzGaus->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fRnPoDzGaus->GetProb()));
	pt->AddText(Form("Eff    %.2f #pm %.2f",grRnPoDzEff->GetY()[grIDX],grRnPoDzEff->GetEY()[grIDX]));
	pt->Draw();
	
	pad2->cd();
	gPad->SetGrid();
	TH1F *hRnPoDzResid = makeResidHist(hRnPoDz,fRnPoDzGaus);
	hRnPoDzResid->GetYaxis()->SetRangeUser(-70,70);
	hRnPoDzResid->Draw();
	hRnPoDzResid->Fit("pol0");
	hRnPoDzResid->GetFunction("pol0")->SetLineStyle(2);
	pt = new TPaveText(0.8,0.8,0.99,0.99,"NDCNB");
	pt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",hRnPoDzResid->GetFunction("pol0")->GetChisquare(),hRnPoDzResid->GetFunction("pol0")->GetNDF()));
	pt->AddText(Form("Prob   %.2f",hRnPoDzResid->GetFunction("pol0")->GetProb()));
	pt->AddText(Form("p0   %.2f #pm %.2f",hRnPoDzResid->GetFunction("pol0")->GetParameter(0),hRnPoDzResid->GetFunction("pol0")->GetParError(0)));
	pt->Draw();	
	cRnPoDz->SaveAs(Form("%s/RnPoDz_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cRnPoDz->SaveAs(Form("%s/RnPoDz_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));

	//----------------------------------------------------------
	TCanvas *cRnPos = new TCanvas("cRnPos","Rn Position",1);
	gPad->SetGrid();
	hSelectPromptPos->SetLineColor(s_col);
	hSelectPromptPos->Draw("HIST");
	hBGPromptPos->SetLineColor(bg_col);
	hBGPromptPos->Draw("HIST&SAME");
	hPoPos->Draw("same");
	pt = new TPaveText(0.75,0.75,0.95,0.9,"NDCNB");
    	tRn = pt->AddText("Rn^{219}");
  	pt->AddText(Form("#mu    %.2f #pm %.2f",hRnPos->GetMean(),hRnPos->GetMeanError()));
   	pt->AddText(Form("RMS    %.2f #pm %.2f",hRnPos->GetRMS(),hRnPos->GetRMSError()));
	pt->Draw();
	cRnPos->SaveAs(Form("%s/RnPos_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cRnPos->SaveAs(Form("%s/RnPos_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));

	TCanvas *cPoPos = new TCanvas("cPoPos","Po Position",1);
	gPad->SetGrid();
	hSelectDelayPos->SetLineColor(s_col);
	hSelectDelayPos->Draw("HIST");
	hBGDelayPos->SetLineColor(bg_col);
	hBGDelayPos->Draw("HIST&SAME");
	hPoPos->Draw("same");
	pt = new TPaveText(0.75,0.75,0.95,0.9,"NDCNB");
    	tPo = pt->AddText("Po^{215}");
  	pt->AddText(Form("#mu    %.2f #pm %.2f",hPoPos->GetMean(),hPoPos->GetMeanError()));
   	pt->AddText(Form("RMS    %.2f #pm %.2f",hPoPos->GetRMS(),hPoPos->GetRMSError()));
	pt->Draw();
	cPoPos->SaveAs(Form("%s/PoPos_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cPoPos->SaveAs(Form("%s/PoPos_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));

	//----------------------------------------------------------
	TCanvas *cRnPoPSDvsEn = new TCanvas("cRnPoPSDvsEn","RnPo PSD vs En",1);
	gPad->SetRightMargin(0.1);
	gPad->SetLeftMargin(0.12);
	hRnPoPSDvsEn->GetYaxis()->SetTitleOffset(1.1);
	hRnPoPSDvsEn->Draw("colz");
	cRnPoPSDvsEn->SaveAs(Form("%s/RnPoPSDvsEn_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cRnPoPSDvsEn->SaveAs(Form("%s/RnPoPSDvsEn_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	
	TCanvas *cRnPoPSDvsPos = new TCanvas("cRnPoPSDvsPos","RnPo PSD vs Pos",1);
	gPad->SetRightMargin(0.1);
	gPad->SetLeftMargin(0.12);
	hRnPoPSDvsPos->GetYaxis()->SetTitleOffset(1.1);
	hRnPoPSDvsPos->Draw("colz");
	cRnPoPSDvsPos->SaveAs(Form("%s/RnPoPSDvsPos_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cRnPoPSDvsPos->SaveAs(Form("%s/RnPoPSDvsPos_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));

	TCanvas *cRnPoEnvsPos = new TCanvas("cRnPoEnvsPos","RnPo En vs Pos",1);
	gPad->SetRightMargin(0.1);
	gPad->SetLeftMargin(0.12);
	hRnPoEnvsPos->GetYaxis()->SetTitleOffset(1.1);
	hRnPoEnvsPos->GetYaxis()->SetTitle("Energy [MeVee]");
	hRnPoEnvsPos->Draw("colz");
	cRnPoEnvsPos->SaveAs(Form("%s/RnPoEnvsPos_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cRnPoEnvsPos->SaveAs(Form("%s/RnPoEnvsPos_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));

	TCanvas *cPoEnvsRnEn = new TCanvas("cPoEnvsRnEn","Po En vs Rn En",650,600);
	gPad->SetRightMargin(0.12);
	gPad->SetLeftMargin(0.12);
	hPoEnvsRnEn->GetYaxis()->SetTitleOffset(1.1);
	hPoEnvsRnEn->Draw("colz");
	cPoEnvsRnEn->SaveAs(Form("%s/PoEnvsRnEn_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cPoEnvsRnEn->SaveAs(Form("%s/PoEnvsRnEn_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));

	return 0;
} 	//end PlotDistributionsVsCell
