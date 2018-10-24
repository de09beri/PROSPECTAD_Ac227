//This macro will create histograms for Ac227 coincidences
//according to cell number

#include "RNPO.C"

#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include <iostream>
#include <fstream>
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TMath.h"
#include "TLatex.h"
#include "TVectorD.h"

#include "Header.C"


void RnPoVsCell(double p_lowPSD, double d_lowPSD, double p_lowE, double d_lowE, double zLow, double zHigh, int dtFit){

	TFile *histFile = new TFile(Form("%s/Ac227_HistsPerCell.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")),"RECREATE");

	//---------------------------------------------------------------------------------
	//Initialize histograms for dt, PSD, E, dz, and position
	int N = NUMCELLS;

	TH1F *hSelectDt[N],		   *hBGDt[N],		 *hRnPoDt[N];
	TH1F *hSelectPromptPSD[N], *hBGPromptPSD[N], *hRnPSD[N];
	TH1F *hSelectDelayPSD[N],  *hBGDelayPSD[N],  *hPoPSD[N];
	TH1F *hSelectPromptEn[N],  *hBGPromptEn[N],  *hRnEn[N];
	TH1F *hSelectDelayEn[N],   *hBGDelayEn[N],   *hPoEn[N];
	TH1F *hSelectPromptPos[N], *hBGPromptPos[N], *hRnPos[N];
	TH1F *hSelectDelayPos[N],  *hBGDelayPos[N],  *hPoPos[N];
	TH1F *hSelectDz[N], 	   *hBGDz[N], 		 *hRnPoDz[N];

	TH1F *hSelectPromptTotEn[N],  *hBGPromptTotEn[N],  *hRnTotEn[N];

	TH2F *hSelectPSDvsEn[N],   *hBGPSDvsEn[N], 	 *hRnPoPSDvsEn[N];
	TH2F *hSelectPSDvsPos[N],  *hBGPSDvsPos[N],  *hRnPoPSDvsPos[N];
	TH2F *hSelectEnvsPos[N],   *hBGEnvsPos[N],	 *hRnPoEnvsPos[N];
	TH2F *hSelectDelayEnvsPromptEn[N], *hBGDelayEnvsPromptEn[N], *hPoEnvsRnEn[N];

	TH1F *hCell_tstamp[N];
	int timelow = 1519794000, timehigh = 1529952428 + (3600*10); 

	for(int i=0;i<NUMCELLS;i++){
		hSelectDt[i] 		= new TH1F(Form("hSelectDt_%i",i),";dt [ms];Counts/0.1 ms",numDtBins,dtMin,dtMax);	
		hBGDt[i] 	 		= new TH1F(Form("hBGDt_%i",i),";dt [ms];Counts/0.1 ms",numDtBins,dtMin,dtMax);	
		
		hSelectPromptPSD[i] = new TH1F(Form("hSelectPromptPSD_%i",i),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);
		hBGPromptPSD[i] 	= new TH1F(Form("hBGPromptPSD_%i",i),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);
	
		hSelectDelayPSD[i] 	= new TH1F(Form("hSelectDelayPSD_%i",i),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);
		hBGDelayPSD[i] 		= new TH1F(Form("hBGDelayPSD_%i",i),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);
			
		hSelectPromptEn[i] 	= new TH1F(Form("hSelectPromptEn_%i",i),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);
		hBGPromptEn[i] 		= new TH1F(Form("hBGPromptEn_%i",i),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);
	
		hSelectDelayEn[i] 	= new TH1F(Form("hSelectDelayEn_%i",i),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);
		hBGDelayEn[i] 		= new TH1F(Form("hBGDelayEn_%i",i),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);

		hSelectPromptPos[i] = new TH1F(Form("hSelectPromptPos_%i",i),";z [mm];Counts/cm",numPosBins,posMin,posMax);
		hBGPromptPos[i] 	= new TH1F(Form("hBGPromptPos_%i",i),";z [mm];Counts/cm",numPosBins,posMin,posMax);

		hSelectDelayPos[i] 	= new TH1F(Form("hSelectDelayPos_%i",i),";z [mm];Counts/cm",numPosBins,posMin,posMax);
		hBGDelayPos[i] 		= new TH1F(Form("hBGDelayPos_%i",i),";z [mm];Counts/cm",numPosBins,posMin,posMax);

		hSelectDz[i]		= new TH1F(Form("hSelectDz_%i",i),";z_{Po} - z_{Rn} [mm];Counts/0.25 cm",numDzBins,dzMin,dzMax);
		hBGDz[i]			= new TH1F(Form("hBGDz_%i",i),";z_{Po} - z_{Rn} [mm];Counts/0.25 cm",numDzBins,dzMin,dzMax);

		hSelectPromptTotEn[i] 	= new TH1F(Form("hSelectPromptTotEn_%i",i),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);
		hBGPromptTotEn[i] 		= new TH1F(Form("hBGPromptTotEn_%i",i),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);
		
		hSelectPSDvsEn[i] 	= new TH2F(Form("hSelectPSDvsEn_%i",i),";Energy [MeVee];PSD [arb]",numEnBins,EnMin,EnMax,numPSDBins,PSDMin,PSDMax);
		hBGPSDvsEn[i] 		= new TH2F(Form("hBGPSDvsEn_%i",i),";Energy [MeVee];PSD [arb]",numEnBins,EnMin,EnMax,numPSDBins,PSDMin,PSDMax);

		hSelectPSDvsPos[i] 	= new TH2F(Form("hSelectPSDvsPos_%i",i),";z [mm];PSD [arb]",numPosBins,posMin,posMax,numPSDBins,PSDMin,PSDMax);
		hBGPSDvsPos[i] 		= new TH2F(Form("hBGPSDvsPos_%i",i),";z [mm];PSD [arb]",numPosBins,posMin,posMax,numPSDBins,PSDMin,PSDMax);

		hSelectEnvsPos[i] 	= new TH2F(Form("hSelectEnvsPos_%i",i),";z [mm];En [MeVee]",numPosBins,posMin,posMax,numEnBins,EnMin,EnMax);
		hBGEnvsPos[i] 		= new TH2F(Form("hBGEnvsPos_%i",i),";z [mm];En [MeVee]",numPosBins,posMin,posMax,numEnBins,EnMin,EnMax);

		hSelectDelayEnvsPromptEn[i] = new TH2F(Form("hSelectDelayEnvsPromptEn_%i",i),";Rn Energy [MeVee];Po Energy [MeVee]",numEnBins,EnMin,EnMax,numEnBins,EnMin,EnMax);
		hBGDelayEnvsPromptEn[i] 	= new TH2F(Form("hBGDelayEnvsPromptEn_%i",i),";Rn Energy [MeVee];Po Energy [MeVee]",numEnBins,EnMin,EnMax,numEnBins,EnMin,EnMax);

		hCell_tstamp[i] = new TH1F(Form("hCell_tstamp_%i",i),Form("Cell %i Event Count;Time;Counts",i),400,timelow,timehigh);

	}	//end for loop creating histograms

	//---------------------------------------------------------------------------------
	//Initialize histograms for all cells
	
	TH1F *hSelectSeg = new TH1F("hSelectSeg",";Cell;Counts",NUMCELLS,0,NUMCELLS);
	TH1F *hBGSeg = new TH1F("hBGSeg",";Cell;Counts",NUMCELLS,0,NUMCELLS);

	//---------------------------------------------------------------------------------
	//Initialize TGraphs for results

	int numActiveCells = NUMCELLS - NUMEXCLUDECELLS;
	double x[numActiveCells], xErr[numActiveCells];
	double y[1], yErr[1];

	TGraphErrors *grRate 		= new TGraphErrors(numActiveCells,x,y,xErr,yErr);
	TGraphErrors *grRnPSDEff 	= new TGraphErrors(numActiveCells,x,y,xErr,yErr);
	TGraphErrors *grPoPSDEff 	= new TGraphErrors(numActiveCells,x,y,xErr,yErr);
	TGraphErrors *grRnEnEff 	= new TGraphErrors(numActiveCells,x,y,xErr,yErr);
	TGraphErrors *grPoEnEff 	= new TGraphErrors(numActiveCells,x,y,xErr,yErr);
	TGraphErrors *grRnPoDzEff 	= new TGraphErrors(numActiveCells,x,y,xErr,yErr);	
	TGraphErrors *grTotEff 		= new TGraphErrors(numActiveCells,x,y,xErr,yErr);
	TGraphErrors *grLifetime 	= new TGraphErrors(numActiveCells,x,y,xErr,yErr);
	TGraphErrors *grPoPSDMean 	= new TGraphErrors(numActiveCells,x,y,xErr,yErr);
	TGraphErrors *grPoPSDSigma 	= new TGraphErrors(numActiveCells,x,y,xErr,yErr);
	TGraphErrors *grPoEnMean   	= new TGraphErrors(numActiveCells,x,y,xErr,yErr);
	TGraphErrors *grPoEnSigma 	= new TGraphErrors(numActiveCells,x,y,xErr,yErr);
	TGraphErrors *grPoPosMean 	= new TGraphErrors(numActiveCells,x,y,xErr,yErr);
	TGraphErrors *grPoPosSigma 	= new TGraphErrors(numActiveCells,x,y,xErr,yErr);
	TGraphErrors *grRnPoDzMean 	= new TGraphErrors(numActiveCells,x,y,xErr,yErr);
	TGraphErrors *grRnPoDzSigma 	= new TGraphErrors(numActiveCells,x,y,xErr,yErr);

	TGraph *grBGRate = new TGraph(numActiveCells,x,y);

	TGraph *grRnPSDChiSq = new TGraph(numActiveCells,x,y);
	TGraph *grPoPSDChiSq = new TGraph(numActiveCells,x,y);
	TGraph *grRnEnChiSq  = new TGraph(numActiveCells,x,y);
	TGraph *grPoEnChiSq  = new TGraph(numActiveCells,x,y);
	TGraph *grDzChiSq    = new TGraph(numActiveCells,x,y);
	TGraph *grDtChiSq    = new TGraph(numActiveCells,x,y);

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	//---------------------------------------------------------------------------------
	//Fill histograms
	printf("=============== Filling Histograms =============== \n"); 

	RNPO *rnpo = new RNPO();
	
	bool exclude;
	Long64_t numEntries = Long64_t(rnpo->fChain->GetEntries());
	printf("Number of Ac-227 candidates: %lld \n",numEntries);

	double tstamp;
	double promptLowPSDCut, promptHighPSDCut, promptLowEnCut, promptHighEnCut;
	double delayLowPSDCut,  delayHighPSDCut,  delayLowEnCut,  delayHighEnCut;
	double dzCut;

	// Get cut values
	rnpo->GetEntry(0);
	tstamp 			 = rnpo->tstamp;
	promptLowPSDCut  = (p_lowPSD > rnpo->p_PSDCut[0]) ? p_lowPSD : rnpo->p_PSDCut[0]; 
	promptHighPSDCut = rnpo->p_PSDCut[1];
	delayLowPSDCut   = (d_lowPSD > rnpo->d_PSDCut[0]) ? d_lowPSD : rnpo->d_PSDCut[0];
	delayHighPSDCut  = rnpo->d_PSDCut[1];
	promptLowEnCut   = (p_lowE > rnpo->p_ECut[0]) ? p_lowE : rnpo->p_ECut[0];
	promptHighEnCut  = rnpo->p_ECut[1];
	delayLowEnCut    = (d_lowE > rnpo->d_ECut[0]) ? d_lowE : rnpo->d_ECut[0];
	delayHighEnCut   = rnpo->d_ECut[1];
	dzCut            = rnpo->dzCut;


	//Initialize variables
	int seg;
	double dt, dz;

	double livetime = 0.0;
	double lastTime = 0.0;

	double lastNumClusts = 0.0, numClusts = 0.0;
	double lastRuntime = 0.0, totRuntime = 0.0;

	for(Long64_t i=0;i<numEntries;i++){
		if(i%1000000==0) printf("Event: %lld  NumClusts: %f \n",i,numClusts);
		rnpo->GetEntry(i);

		if(rnpo->d_t*(1e-6) > ((double)((TVectorD*)rnpo->fChain->GetCurrentFile()->Get("runtime"))->Norm1()*1000.0 - (TIMEWINDOW+TIMEOFFSET))) continue;

		if(rnpo->d_t < lastTime){ 
			livetime += lastTime*(1e-6);		//livetime in ms	
			totRuntime += lastRuntime;
			numClusts += lastNumClusts;
		}
		lastTime = rnpo->d_t;
		lastNumClusts = rnpo->numClust;
		lastRuntime = ((TVectorD*)rnpo->fChain->GetCurrentFile()->Get("runtime"))->Norm1();	//[s]	

		seg = rnpo->d_seg;

/*
		double rnpo_p_E = rnpo->p_ESmear;	
		double rnpo_d_E = rnpo->d_ESmear;
		double rnpo_f_E = rnpo->f_ESmear;
*/

		double rnpo_p_E = rnpo->p_E;	
		double rnpo_d_E = rnpo->d_E;
		double rnpo_f_E = rnpo->f_E;

		if(rnpo->d_PSD < delayLowPSDCut || rnpo->d_E < delayLowEnCut) continue;
		if(rnpo->d_z < zLow || rnpo->d_z > zHigh) continue;

		if((rnpo->p_PSD>promptLowPSDCut && rnpo->p_E>promptLowEnCut) || (rnpo->f_PSD>promptLowPSDCut && rnpo->f_E>promptLowEnCut)) hCell_tstamp[seg]->Fill(rnpo->tstamp);

		exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), seg) != end(ExcludeCellArr);
		if(exclude) continue;


		//if prompt-delay pair
		dt = (rnpo->d_t - rnpo->p_t)*(1e-6);	//convert ns to ms	
		if(rnpo->p_seg > -1 && rnpo->p_PSD>promptLowPSDCut && rnpo->p_E>promptLowEnCut && rnpo->p_z>zLow && rnpo->p_z<zHigh && dt > 0.5){
			hSelectSeg->Fill(rnpo->d_seg);		
			
			dt = (rnpo->d_t - rnpo->p_t)*(1e-6);	//convert ns to ms	
			dz = rnpo->d_z - rnpo->p_z;


			hSelectDt[seg]->Fill(dt);
			hSelectPromptPSD[seg]->Fill(rnpo->p_PSD);
			hSelectDelayPSD[seg]->Fill(rnpo->d_PSD);
			hSelectPromptEn[seg]->Fill(rnpo_p_E);
			hSelectDelayEn[seg]->Fill(rnpo_d_E);
			hSelectPromptTotEn[seg]->Fill(rnpo->p_Etot);	
			hSelectPromptPos[seg]->Fill(rnpo->p_z);
			hSelectDelayPos[seg]->Fill(rnpo->d_z);
			hSelectDz[seg]->Fill(dz);

			hSelectPSDvsEn[seg]->Fill(rnpo_p_E,rnpo->p_PSD);
			hSelectPSDvsEn[seg]->Fill(rnpo_d_E,rnpo->d_PSD);
			hSelectPSDvsPos[seg]->Fill(rnpo->p_z,rnpo->p_PSD);
			hSelectPSDvsPos[seg]->Fill(rnpo->d_z,rnpo->d_PSD);
			hSelectEnvsPos[seg]->Fill(rnpo->p_z,rnpo_p_E);
			hSelectEnvsPos[seg]->Fill(rnpo->d_z,rnpo_d_E);
			hSelectDelayEnvsPromptEn[seg]->Fill(rnpo_p_E,rnpo_d_E);		
		}

		//if prompt-delay BG pair
		dt = (rnpo->f_t - rnpo->d_t)*(1e-6) - TIMEOFFSET;
		if(rnpo->f_seg > -1 && rnpo->f_PSD>promptLowPSDCut && rnpo->f_E>promptLowEnCut && rnpo->f_z>zLow && rnpo->f_z<zHigh && dt > 0.5){
			hBGSeg->Fill(rnpo->d_seg);
		
			dt = (rnpo->f_t - rnpo->d_t)*(1e-6) - TIMEOFFSET;
			dz = rnpo->d_z - rnpo->f_z;

			hBGDt[seg]->Fill(dt);
			hBGPromptPSD[seg]->Fill(rnpo->f_PSD);
			hBGDelayPSD[seg]->Fill(rnpo->d_PSD);
			hBGPromptEn[seg]->Fill(rnpo_f_E);
			hBGDelayEn[seg]->Fill(rnpo_d_E);
			hBGPromptTotEn[seg]->Fill(rnpo->f_Etot);
			hBGPromptPos[seg]->Fill(rnpo->f_z);
			hBGDelayPos[seg]->Fill(rnpo->d_z);
			hBGDz[seg]->Fill(dz);

			hBGPSDvsEn[seg]->Fill(rnpo_f_E,rnpo->f_PSD);
			hBGPSDvsEn[seg]->Fill(rnpo_d_E,rnpo->d_PSD);
			hBGPSDvsPos[seg]->Fill(rnpo->f_z,rnpo->f_PSD);
			hBGPSDvsPos[seg]->Fill(rnpo->d_z,rnpo->d_PSD);
			hBGEnvsPos[seg]->Fill(rnpo->f_z,rnpo_f_E);
			hBGEnvsPos[seg]->Fill(rnpo->d_z,rnpo_d_E);
			hBGDelayEnvsPromptEn[seg]->Fill(rnpo_f_E,rnpo_d_E);	
		}
	}	//end for loop over TChain


	livetime += lastTime*(1e-6);	//add time from last tree
	totRuntime += lastRuntime;
	numClusts += lastNumClusts;

	double pileupVetoTime = numClusts*pileupVetoT;	//[ms]

	printf("Total runtime: %f hours \n",totRuntime/(60.0*60.0));
	printf("Livetime: %f hours \n",livetime*(2.778e-7));
	printf("Pileup veto time: %f hours \n",pileupVetoTime*(2.778e-7));

	livetime = livetime - (2.0*pileupVetoTime);
	printf("Corrected Livetime: %f hours \n",livetime*(2.778e-7));

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	//---------------------------------------------------------------------------------
	//Subtract histograms
	printf("=============== Subtracting Histograms =============== \n"); 

	for(int i=0;i<NUMCELLS;i++){
		if(i%10==0) printf("Cell: %d \n",i);

		exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), i) != end(ExcludeCellArr);
		if(exclude) continue;

		hSelectDt[i]->Sumw2();
		hRnPoDt[i] = (TH1F*)hSelectDt[i]->Clone();
		hRnPoDt[i]->SetName(Form("hRnPoDt_%i",i));
		hBGDt[i]->Sumw2();
		hRnPoDt[i]->Add(hBGDt[i],-1);

		hRnPSD[i] = (TH1F*)hSelectPromptPSD[i]->Clone();
		hRnPSD[i]->SetName(Form("hRnPSD_%i",i));
		hRnPSD[i]->Sumw2();
		hRnPSD[i]->Add(hBGPromptPSD[i],-1);

		hPoPSD[i] = (TH1F*)hSelectDelayPSD[i]->Clone();
		hPoPSD[i]->SetName(Form("hPoPSD_%i",i));
		hPoPSD[i]->Sumw2();
		hPoPSD[i]->Add(hBGDelayPSD[i],-1);

		hRnEn[i] = (TH1F*)hSelectPromptEn[i]->Clone();
		hRnEn[i]->SetName(Form("hRnEn_%i",i));
		hRnEn[i]->Sumw2();
		hRnEn[i]->Add(hBGPromptEn[i],-1);

		hRnTotEn[i] = (TH1F*)hSelectPromptTotEn[i]->Clone();
		hRnTotEn[i]->SetName(Form("hRnTotEn_%i",i));
		hRnTotEn[i]->Sumw2();
		hRnTotEn[i]->Add(hBGPromptTotEn[i],-1);

		hPoEn[i] = (TH1F*)hSelectDelayEn[i]->Clone();
		hPoEn[i]->SetName(Form("hPoEn_%i",i));
		hPoEn[i]->Sumw2();
		hPoEn[i]->Add(hBGDelayEn[i],-1);
	
		hRnPos[i] = (TH1F*)hSelectPromptPos[i]->Clone();
		hRnPos[i]->SetName(Form("hRnPos_%i",i));
		hRnPos[i]->Sumw2();
		hRnPos[i]->Add(hBGPromptPos[i],-1);

		hPoPos[i] = (TH1F*)hSelectDelayPos[i]->Clone();
		hPoPos[i]->SetName(Form("hPoPos_%i",i));
		hPoPos[i]->Sumw2();
		hPoPos[i]->Add(hBGDelayPos[i],-1);
			
		hRnPoDz[i] = (TH1F*)hSelectDz[i]->Clone();
		hRnPoDz[i]->SetName(Form("hRnPoDz_%i",i));
		hRnPoDz[i]->Sumw2();
		hRnPoDz[i]->Add(hBGDz[i],-1);

		hRnPoPSDvsEn[i] = (TH2F*)hSelectPSDvsEn[i]->Clone();
		hRnPoPSDvsEn[i]->SetName(Form("hRnPoPSDvsEn_%i",i));
		hRnPoPSDvsEn[i]->Add(hBGPSDvsEn[i],-1);

		hRnPoPSDvsPos[i] = (TH2F*)hSelectPSDvsPos[i]->Clone();
		hRnPoPSDvsPos[i]->SetName(Form("hRnPoPSDvsPos_%i",i));
		hRnPoPSDvsPos[i]->Add(hBGPSDvsPos[i],-1);

		hRnPoEnvsPos[i] = (TH2F*)hSelectEnvsPos[i]->Clone();
		hRnPoEnvsPos[i]->SetName(Form("hRnPoEnvsPos_%i",i));
		hRnPoEnvsPos[i]->Add(hBGEnvsPos[i],-1);

		hPoEnvsRnEn[i] = (TH2F*)hSelectDelayEnvsPromptEn[i]->Clone();
		hPoEnvsRnEn[i]->SetName(Form("hPoEnvsRnEn_%i",i));
		hPoEnvsRnEn[i]->Add(hBGDelayEnvsPromptEn[i],-1);
	}	//end for loop to subtract hists


//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	//---------------------------------------------------------------------------------
	//Calculate results	
	printf("=============== Calculating Results =============== \n"); 

	double dtBinWidth = (dtMax - dtMin)/(double)numDtBins;

	TF1 *fRnPoDtExp;
	TF1 *fRnPSDGaus, *fPoPSDGaus;
	TF1 *fRnEnGaus,  *fPoEnGaus;
	TF1 *fRnPoDzGaus;

	double promptPSDEff, delayPSDEff, promptPSDEffErr, delayPSDEffErr;
	double promptEnEff,  delayEnEff,  promptEnEffErr,  delayEnEffErr;
	double dzEff, dzEffErr;
	double totEff, totEffErr;

	double NAlpha, NAlphaErr, lifetime, lifetimeErr;
	double rate, rateErr;

	int grPt = 0;

	for(int i=0;i<NUMCELLS;i++){
		if(i%10==0) printf("Cell: %d \n",i);

		exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), i) != end(ExcludeCellArr);
		if(exclude) continue;

		double BGRate = (hBGDt[i]->GetEntries()/livetime)*(1e6);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//Fit distributions
		fRnPoDtExp = new TF1("fRnPoDtExp",Form("[0]*exp(-x/[1])*(%f/[1])",dtBinWidth),0.0,dtMax);
		//fRnPoDtExp = new TF1("fRnPoDtExp",Form("[0]*exp(-x/[1])*(%f/[1])",dtBinWidth),0.5,11);
		fRnPoDtExp->SetParameter(1,POLIFETIME);
		hRnPoDt[i]->Fit(fRnPoDtExp,"R0");	

		double fitPSDMin = promptLowPSDCut;
		double maxValue = hRnPSD[i]->GetMaximum(), maxBin = hRnPSD[i]->GetBinCenter(hRnPSD[i]->GetMaximumBin());	
		fRnPSDGaus = new TF1("fRnPSDGaus","gaus(0)+gaus(3)",fitPSDMin,PSDMax);
		fRnPSDGaus->SetParameters(0.25*maxValue,maxBin+0.01,0.015,0.75*maxValue,maxBin-0.001,0.02);
		fRnPSDGaus->SetParLimits(0,0,maxValue);
		fRnPSDGaus->SetParLimits(3,0,maxValue);
		hRnPSD[i]->Fit(fRnPSDGaus,"RQ0B");
		fRnPSDGaus->SetRange(PSDMin,PSDMax);
	
		fitPSDMin = delayLowPSDCut;
		fPoPSDGaus = new TF1("fPoPSDGaus","gaus",fitPSDMin,PSDMax);
		hPoPSD[i]->Fit(fPoPSDGaus,"RQ0");
		fPoPSDGaus->SetRange(PSDMin,PSDMax);
		
		double fitRnEnMin = EnMin;	
		fRnEnGaus = new TF1("fRnEnGaus","gaus(0)+gaus(3)+gaus(6)",fitRnEnMin,EnMax);
		maxValue = hRnEn[i]->GetMaximum();
		maxBin = hRnEn[i]->GetBinCenter(hRnEn[i]->GetMaximumBin());	
		fRnEnGaus->SetParameters(0.8*maxValue,maxBin,0.035,0.05*maxValue,maxBin+0.14,0.08,0.15*maxValue,maxBin-0.03,0.05);
		hRnEn[i]->Fit(fRnEnGaus,"RQ0");
		fRnEnGaus->SetRange(EnMin,EnMax);

		double fitPoEnMin = delayLowEnCut;
		fPoEnGaus = new TF1("fPoEnGaus","gaus",fitPoEnMin,EnMax);
		hPoEn[i]->Fit(fPoEnGaus,"RQ0");
		fPoEnGaus->SetRange(EnMin,EnMax);

		fRnPoDzGaus = new TF1("fRnPoDzGaus","gaus",dzMin,dzMax);
		hRnPoDz[i]->Fit(fRnPoDzGaus,"RQ0");

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//Calculate efficiencies
		promptPSDEff = fRnPSDGaus->Integral(promptLowPSDCut,promptHighPSDCut)/fRnPSDGaus->Integral(PSDMin,PSDMax);
		promptPSDEffErr = sqrt((promptPSDEff*(1-promptPSDEff))/hRnPSD[i]->GetEntries()); 
			
		delayPSDEff = fPoPSDGaus->Integral(delayLowPSDCut,delayHighPSDCut)/fPoPSDGaus->Integral(PSDMin,PSDMax);
		delayPSDEffErr = sqrt((delayPSDEff*(1-delayPSDEff))/hPoPSD[i]->GetEntries()); 
		
		promptEnEff = fRnEnGaus->Integral(promptLowEnCut,promptHighEnCut)/fRnEnGaus->Integral(EnMin,EnMax);
		promptEnEffErr = sqrt((promptEnEff*(1-promptEnEff))/hRnEn[i]->GetEntries());

		delayEnEff = fPoEnGaus->Integral(delayLowEnCut,delayHighEnCut)/fPoEnGaus->Integral(EnMin,EnMax);
		delayEnEffErr = sqrt((delayEnEff*(1-delayEnEff))/hPoEn[i]->GetEntries()); 

		dzEff = fRnPoDzGaus->Integral(-dzCut,dzCut)/fRnPoDzGaus->Integral(dzMin,dzMax);
		dzEffErr = sqrt((dzEff*(1-dzEff))/hRnPoDz[i]->GetEntries());

		
		totEff = promptPSDEff * delayPSDEff * promptEnEff * delayEnEff * dzEff;	
		totEffErr = totEff * sqrt( pow(promptPSDEffErr/promptPSDEff,2) + pow(delayPSDEffErr/delayPSDEff,2) + pow(promptEnEffErr/promptEnEff,2) + pow(delayEnEffErr/delayEnEff,2) + pow(dzEffErr/dzEff,2) );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//Calculate rate
		NAlpha = fRnPoDtExp->GetParameter(0);
		NAlphaErr = fRnPoDtExp->GetParError(0);
		lifetime = fRnPoDtExp->GetParameter(1);
		lifetimeErr = fRnPoDtExp->GetParError(1);

		rate = (NAlpha/(livetime*totEff))*(1e6);		//mHz
		rateErr = rate * sqrt(pow(NAlphaErr/NAlpha,2) + pow(totEffErr/totEff,2));

		if(rate>5){
			printf("!=================================!\n");
			printf("ABNORMAL RATE: Cell %i \n",i);
			printf("Rate: %.2f +/- %.2f \n",rate,rateErr);
			printf("Prompt PSD Eff: %.6f +/- %.6f | En Eff: %.6f +/- %.6f \n", promptPSDEff,promptPSDEffErr,promptEnEff,promptEnEffErr);
			printf("Delay  PSD Eff: %.6f +/- %.6f | En Eff: %.6f +/- %.6f \n", delayPSDEff,delayPSDEffErr,delayEnEff,delayEnEffErr);
			printf("dz Eff: %.6f +/- %.6f \n",dzEff,dzEffErr);
			printf("!=================================!\n");
		}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//Populate graphs
		grRate->SetPoint(grPt,i,rate);
		grRate->SetPointError(grPt,0,rateErr);	

		grRnPSDEff->SetPoint(grPt,i,promptPSDEff);
		grRnPSDEff->SetPointError(grPt,0,promptPSDEffErr);
	
		grPoPSDEff->SetPoint(grPt,i,delayPSDEff);
		grPoPSDEff->SetPointError(grPt,0,delayPSDEffErr);

		grRnEnEff->SetPoint(grPt,i,promptEnEff);
		grRnEnEff->SetPointError(grPt,0,promptEnEffErr);
		
		grPoEnEff->SetPoint(grPt,i,delayEnEff);
		grPoEnEff->SetPointError(grPt,0,delayEnEffErr);

		grRnPoDzEff->SetPoint(grPt,i,dzEff);
		grRnPoDzEff->SetPointError(grPt,0,dzEffErr);
		
		grTotEff->SetPoint(grPt,i,totEff);
		grTotEff->SetPointError(grPt,0,totEffErr);

		grLifetime->SetPoint(grPt,i,lifetime);
		grLifetime->SetPointError(grPt,0,lifetimeErr);

		grPoPSDMean->SetPoint(grPt,i,fPoPSDGaus->GetParameter(1));	
		grPoPSDMean->SetPointError(grPt,0,fPoPSDGaus->GetParError(1));

		grPoPSDSigma->SetPoint(grPt,i,fPoPSDGaus->GetParameter(2));	
		grPoPSDSigma->SetPointError(grPt,0,fPoPSDGaus->GetParError(2));
	
		grPoEnMean->SetPoint(grPt,i,fPoEnGaus->GetParameter(1));
		grPoEnMean->SetPointError(grPt,0,fPoEnGaus->GetParError(1));
	
		grPoEnSigma->SetPoint(grPt,i,fPoEnGaus->GetParameter(2));
		grPoEnSigma->SetPointError(grPt,0,fPoEnGaus->GetParError(2));

		grPoPosMean->SetPoint(grPt,i,hPoPos[i]->GetMean());
		grPoPosMean->SetPointError(grPt,0,hPoPos[i]->GetMeanError());

		grPoPosSigma->SetPoint(grPt,i,hPoPos[i]->GetRMS());
		grPoPosSigma->SetPointError(grPt,0,hPoPos[i]->GetRMSError());

		grRnPoDzMean->SetPoint(grPt,i,fRnPoDzGaus->GetParameter(1));
		grRnPoDzMean->SetPointError(grPt,0,fRnPoDzGaus->GetParError(1));

		grRnPoDzSigma->SetPoint(grPt,i,fRnPoDzGaus->GetParameter(2));
		grRnPoDzSigma->SetPointError(grPt,0,fRnPoDzGaus->GetParError(2));

		grBGRate->SetPoint(grPt,i,BGRate);

		grRnPSDChiSq->SetPoint(grPt,i,fRnPSDGaus->GetChisquare()/(double)fRnPSDGaus->GetNDF());
		grPoPSDChiSq->SetPoint(grPt,i,fPoPSDGaus->GetChisquare()/(double)fPoPSDGaus->GetNDF());
	
		grRnEnChiSq->SetPoint(grPt,i,fRnEnGaus->GetChisquare()/(double)fRnEnGaus->GetNDF());
		grPoEnChiSq->SetPoint(grPt,i,fPoEnGaus->GetChisquare()/(double)fPoEnGaus->GetNDF());
		
		grDzChiSq->SetPoint(grPt,i,fRnPoDzGaus->GetChisquare()/(double)fRnPoDzGaus->GetNDF());	
	
		grDtChiSq->SetPoint(grPt,i,fRnPoDtExp->GetChisquare()/(double)fRnPoDtExp->GetNDF());

		grPt++;

	}	//end for loop to calculate results

	histFile->Write();
	histFile->Close();

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	//---------------------------------------------------------------------------------
	//Write TGraphs to file
	TFile *graphFile = new TFile(Form("%s/Ac227_GraphsPerCell.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")),"RECREATE");
	
	grRate->Write("grRate");
	grRnPSDEff->Write("grRnPSDEff");
	grPoPSDEff->Write("grPoPSDEff");
	grRnEnEff->Write("grRnEnEff");
	grPoEnEff->Write("grPoEnEff");
	grRnPoDzEff->Write("grRnPoDzEff");
	grTotEff->Write("grTotEff");	
	grLifetime->Write("grLifetime");
	grPoPSDMean->Write("grPoPSDMean");
	grPoPSDSigma->Write("grPoPSDSigma");
	grPoEnMean->Write("grPoEnMean");
	grPoEnSigma->Write("grPoEnSigma");
	grPoPosMean->Write("grPoPosMean");
	grPoPosSigma->Write("grPoPosSigma");
	grRnPoDzMean->Write("grRnPoDzMean");
	grRnPoDzSigma->Write("grRnPoDzSigma");	
	grBGRate->Write("grBGRate");
	grRnPSDChiSq->Write("grRnPSDChiSq");
	grPoPSDChiSq->Write("grPoPSDChiSq");
	grRnEnChiSq->Write("grRnEnChiSq");
	grPoEnChiSq->Write("grPoEnChiSq");
	grDzChiSq->Write("grDzChiSq");
	grDtChiSq->Write("grDtChiSq");

	graphFile->Close();

}	//end void RnPoVsCell

