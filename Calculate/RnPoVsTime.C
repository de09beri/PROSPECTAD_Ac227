//This macro will create histograms for Ac227 coincidences
//according to cell number

#include "RNPO.C"

#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include <iostream>
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TMath.h"
#include "TLatex.h"
#include "TVectorD.h"

#include "Header.C"


void RnPoVsTime(double p_lowPSD, double d_lowPSD, double p_lowE, double d_lowE, double zLow, double zHigh, double timeBin, int dtFit){

	const double TIMEBREAK = timeBin*(3.6e6);	//[ms]

	TFile *histFile = new TFile(Form("%s/Ac227_HistsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")),"RECREATE");

	//---------------------------------------------------------------------------------
	TH1F *hSelectDt, 		*hBGDt,        *hRnPoDt;
	TH1F *hSelectPromptPSD, *hBGPromptPSD, *hRnPSD,	 *hSelectDelayPSD, *hBGDelayPSD, *hPoPSD;
	TH1F *hSelectPromptEn,  *hBGPromptEn,  *hRnEn,   *hSelectDelayEn,  *hBGDelayEn,  *hPoEn;
	TH1F *hSelectPromptPos, *hBGPromptPos, *hRnPos,  *hSelectDelayPos, *hBGDelayPos, *hPoPos;
	TH1F *hSelectDz,        *hBGDz,        *hRnPoDz;

	TH1F *hSelectPromptTotEn, *hBGPromptTotEn, *hRnTotEn;

	TH2F *hSelectPSDvsEn, 			*hBGPSDvsEn, 		   *hRnPoPSDvsEn;
	TH2F *hSelectDelayEnvsPromptEn, *hBGDelayEnvsPromptEn, *hPoEnvsRnEn;

	//---------------------------------------------------------------------------------
	TF1 *fRnPoDtExp;
	TF1 *fRnPSDGaus, *fPoPSDGaus;
	TF1 *fRnEnGaus,  *fPoEnGaus;
	TF1 *fRnPoDzGaus;

	//---------------------------------------------------------------------------------
	vector<double> vRate,   vRateErr;
	vector<double> vTotEff, vTotEffErr;
	vector<double> vLifetime,   vLifetimeErr;
	vector<double> vLivetime,   vTimestamp;
	vector<double> vPoPSDMean,  vPoPSDMeanErr,  vPoPSDSigma,  vPoPSDSigmaErr;
	vector<double> vPoEnMean,   vPoEnMeanErr,   vPoEnSigma,   vPoEnSigmaErr;
	vector<double> vPoPosMean,  vPoPosMeanErr,  vPoPosSigma,  vPoPosSigmaErr;
	vector<double> vRnPoDzMean, vRnPoDzMeanErr, vRnPoDzSigma, vRnPoDzSigmaErr;

	vector<double> vTotLivetime,  vPileupVetoT;
	vector<double> vPileupVetoFrac;
	
	vector<double> vPromptEnEff,  vPromptEnEffErr,  vDelayEnEff,  vDelayEnEffErr;
	vector<double> vPromptPSDEff, vPromptPSDEffErr, vDelayPSDEff, vDelayPSDEffErr; 
	vector<double> vDzEff, 	      vDzEffErr;

	vector<double> vRnPSDChiSq, vPoPSDChiSq, vRnEnChiSq, vPoEnChiSq, vDzChiSq, vDtChiSq;

	vector<double> vBGRate;

	//---------------------------------------------------------------------------------
	RNPO *rnpo = new RNPO();
	
	bool exclude;
	Long64_t numEntries = Long64_t(rnpo->fChain->GetEntries());
	printf("Number of Ac-227 candidates: %lld \n",numEntries);

	double promptLowPSDCut, promptHighPSDCut, promptLowEnCut, promptHighEnCut;
	double delayLowPSDCut,  delayHighPSDCut,  delayLowEnCut,  delayHighEnCut;
	double dzCut;

	// Get cut values
	rnpo->GetEntry(0);
	promptLowPSDCut  = (p_lowPSD > rnpo->p_PSDCut[0]) ? p_lowPSD : rnpo->p_PSDCut[0];
	promptHighPSDCut = rnpo->p_PSDCut[1];
	delayLowPSDCut   = (d_lowPSD > rnpo->d_PSDCut[0]) ? d_lowPSD : rnpo->d_PSDCut[0];
	delayHighPSDCut  = rnpo->d_PSDCut[1];
	promptLowEnCut   = (p_lowE > rnpo->p_ECut[0]) ? p_lowE : rnpo->p_ECut[0];
    	promptHighEnCut  = rnpo->p_ECut[1];
    	delayLowEnCut    = (d_lowE > rnpo->d_ECut[0]) ? d_lowE : rnpo->d_ECut[0];
    	delayHighEnCut   = rnpo->d_ECut[1];
    	dzCut            = rnpo->dzCut;

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	//---------------------------------------------------------------------------------
	//Loop over all events, splitting into time periods of 24 hours

	Long64_t IDX = 0;
	int numTimeBin = 0;
	while(IDX<numEntries){

		//---------------------------------------------------------------------------------
		//Initialize histograms
		hSelectDt 			= new TH1F(Form("hSelectDt_%i",numTimeBin),";dt [ms];Counts/0.1 ms",numDtBins,dtMin,dtMax);
		hBGDt 				= new TH1F(Form("hBGDt_%i",numTimeBin),";dt [ms];Counts/0.1 ms",numDtBins,dtMin,dtMax);
	
		hSelectPromptPSD 	= new TH1F(Form("hSelectPromptPSD_%i",numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	
		hBGPromptPSD 		= new TH1F(Form("hBGPromptPSD_%i",numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	

		hSelectDelayPSD 	= new TH1F(Form("hSelectDelayPSD_%i",numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	
		hBGDelayPSD 		= new TH1F(Form("hBGDelayPSD_%i",numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	

		hSelectPromptEn 	= new TH1F(Form("hSelectPromptEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
		hBGPromptEn 		= new TH1F(Form("hBGPromptEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
	
		hSelectDelayEn 		= new TH1F(Form("hSelectDelayEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
		hBGDelayEn 			= new TH1F(Form("hBGDelayEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	

		hSelectPromptPos 	= new TH1F(Form("hSelectPromptPos_%i",numTimeBin),";z [mm];Counts/cm",numPosBins,posMin,posMax);
		hBGPromptPos 		= new TH1F(Form("hBGPromptPos_%i",numTimeBin),";z [mm];Counts/cm",numPosBins,posMin,posMax);

		hSelectDelayPos		= new TH1F(Form("hSelectDelayPos_%i",numTimeBin),";z [mm];Counts/cm",numPosBins,posMin,posMax);
		hBGDelayPos 		= new TH1F(Form("hBGDelayPos_%i",numTimeBin),";z [mm];Counts/cm",numPosBins,posMin,posMax);

		hSelectDz 			= new TH1F(Form("hSelectDz_%i",numTimeBin),";z_{Po} - z_{Rn} [mm];Counts/0.25 cm",numDzBins,dzMin,dzMax);
		hBGDz 				= new TH1F(Form("hBGDz_%i",numTimeBin),";z_{Po} - z_{Rn} [mm];Counts/0.25 cm",numDzBins,dzMin,dzMax);

		hSelectPromptTotEn 	= new TH1F(Form("hSelectPromptTotEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
		hBGPromptTotEn 		= new TH1F(Form("hBGPromptTotEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	

		hSelectPSDvsEn 		= new TH2F(Form("hSelectPSDvsEn_%i",numTimeBin),";Energy [MeVee];PSD [arb]",numEnBins,EnMin,EnMax,numPSDBins,PSDMin,PSDMax);
		hBGPSDvsEn 			= new TH2F(Form("hBGPSDvsEn_%i",numTimeBin),";Energy [MeVee];PSD [arb]",numEnBins,EnMin,EnMax,numPSDBins,PSDMin,PSDMax);

		hSelectDelayEnvsPromptEn 	= new TH2F(Form("hSelectDelayEnvsPromptEn_%i",numTimeBin),";Rn Energy [MeVee];Po Energy [MeVee]",numEnBins,EnMin,EnMax,numEnBins,EnMin,EnMax);
		hBGDelayEnvsPromptEn 		= new TH2F(Form("hBGDelayEnvsPromptEn_%i",numTimeBin),";Rn Energy [MeVee];Po Energy [MeVee]",numEnBins,EnMin,EnMax,numEnBins,EnMin,EnMax);

		//---------------------------------------------------------------------------------
		//Fill histograms
		printf("=============== Filling Histograms =============== \n"); 
	
		int seg;
		double dt, dz;
		
		double livetime = 0.0, tstamp;
		double lastTime = 0.0, lastRunTime = 0.0, lastTimestamp = 0.0;	
		double sumWeightedTimestamp = 0.0, sumRunTime = 0.0;

		double lastNumClusts = 0.0, numClusts = 0.0;

		for(Long64_t i=IDX;i<numEntries;i++){
			if(i%100000==0) printf("Event: %lld \n",i);
			rnpo->GetEntry(i);

			if(rnpo->d_t*(1e-6) > ((double)((TVectorD*)rnpo->fChain->GetCurrentFile()->Get("runtime"))->Norm1()*1000.0 - (TIMEWINDOW+TIMEOFFSET)) && i!=(numEntries-1)) continue;
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			//Calculate livetime and weighted timestamp
			if(rnpo->d_t < lastTime){ 
				livetime += lastTime*(1e-6);		//livetime in ms	

				sumWeightedTimestamp += lastRunTime * ((lastRunTime/2.0)+lastTimestamp);
				sumRunTime += lastRunTime;

				numClusts += lastNumClusts;
			}

			if(livetime>TIMEBREAK){

				tstamp = sumWeightedTimestamp/sumRunTime;
				vTimestamp.push_back(tstamp);	

				IDX = i;
				break;
			}

			lastTime = rnpo->d_t;
			lastRunTime = ((TVectorD*)rnpo->fChain->GetCurrentFile()->Get("runtime"))->Norm1();		//seconds
			lastTimestamp = rnpo->tstamp;			//epoch seconds	
			lastNumClusts = rnpo->numClust;

			//if we are at the last entry	
			if(i == (numEntries-1)){ 
				livetime += lastTime*(1e-6);

				//if livetime is less than 12 hours 
				if(livetime*(2.778e-7) < 12) goto endloop;

				sumWeightedTimestamp += lastRunTime * ((lastRunTime/2.0)+lastTimestamp);
				sumRunTime += lastRunTime;

				tstamp = sumWeightedTimestamp/sumRunTime;
				vTimestamp.push_back(tstamp);	
	
				numClusts += lastNumClusts;

				IDX = i+1;
			}

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			//Fill histograms

			seg = rnpo->d_seg;
/*
			double rnpo_p_E = rnpo->p_ESmear;	
			double rnpo_d_E = rnpo->d_ESmear;
			double rnpo_f_E = rnpo->f_ESmear;
*/
			double rnpo_p_E = rnpo->p_E;	
			double rnpo_d_E = rnpo->d_E;
			double rnpo_f_E = rnpo->f_E;

			exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), seg) != end(ExcludeCellArr);
			if(exclude) continue;

			if(rnpo->d_PSD < delayLowPSDCut || rnpo->d_E < delayLowEnCut) continue;	
			if(rnpo->d_z < zLow || rnpo->d_z > zHigh) continue;

			dt = (rnpo->d_t - rnpo->p_t)*(1e-6);	//convert ns to ms	
			if(rnpo->p_seg > -1 && rnpo->p_PSD>promptLowPSDCut && rnpo->p_E>promptLowEnCut && rnpo->p_z>zLow && rnpo->p_z<zHigh && dt>0.5){	
				dt = (rnpo->d_t - rnpo->p_t)*(1e-6);	//convert ns to ms	
				dz = rnpo->d_z - rnpo->p_z;

				hSelectDt->Fill(dt);
				hSelectPromptPSD->Fill(rnpo->p_PSD);
				hSelectDelayPSD->Fill(rnpo->d_PSD);
				hSelectPromptEn->Fill(rnpo_p_E);
				hSelectDelayEn->Fill(rnpo_d_E);
				hSelectPromptTotEn->Fill(rnpo->p_Etot);	
				hSelectPromptPos->Fill(rnpo->p_z);
				hSelectDelayPos->Fill(rnpo->d_z);
				hSelectDz->Fill(dz);

				hSelectPSDvsEn->Fill(rnpo_p_E,rnpo->p_PSD);
				hSelectPSDvsEn->Fill(rnpo_d_E,rnpo->d_PSD);
				hSelectDelayEnvsPromptEn->Fill(rnpo_p_E,rnpo_d_E);		
			}
			dt = (rnpo->f_t - rnpo->d_t)*(1e-6) - TIMEOFFSET;	
			if(rnpo->f_seg > -1 && rnpo->f_PSD>promptLowPSDCut && rnpo->f_E>promptLowEnCut && rnpo->f_z>zLow && rnpo->f_z<zHigh && dt>0.5){	
				dt = (rnpo->f_t - rnpo->d_t)*(1e-6) - TIMEOFFSET;	
				dz = rnpo->d_z - rnpo->f_z;

				hBGDt->Fill(dt);
				hBGPromptPSD->Fill(rnpo->f_PSD);
				hBGDelayPSD->Fill(rnpo->d_PSD);
				hBGPromptEn->Fill(rnpo_f_E);
				hBGDelayEn->Fill(rnpo_d_E);
				hBGPromptTotEn->Fill(rnpo->f_Etot);	
				hBGPromptPos->Fill(rnpo->f_z);
				hBGDelayPos->Fill(rnpo->d_z);
				hBGDz->Fill(dz);

				hBGPSDvsEn->Fill(rnpo_f_E,rnpo->f_PSD);
				hBGPSDvsEn->Fill(rnpo_d_E,rnpo->d_PSD);
				hBGDelayEnvsPromptEn->Fill(rnpo_f_E,rnpo_d_E);		
			}

		}	//end for loop over events

		double pileupVetoTime = numClusts*pileupVetoT;	//[ms]

		printf("Time bin: %i  |  Livetime: %f hrs \n",numTimeBin,livetime*(2.778e-7));
		printf("Pileup veto time: %f \n ms",pileupVetoTime);

		vTotLivetime.push_back(livetime/(1000.0*60.0));		//minutes
		vPileupVetoT.push_back(pileupVetoTime/(1000.0*60.0));	//minutes

		vPileupVetoFrac.push_back(pileupVetoTime/livetime);

		livetime = livetime - 2.0*pileupVetoTime;
		vLivetime.push_back(livetime);

		printf("Corrected Livetime: %f hours \n",livetime*(2.778e-7));

		double BGRate = (hBGDt->GetEntries()/livetime)*(1e3);	//Hz
		vBGRate.push_back(BGRate);

		//---------------------------------------------------------------------------------
		//Subtract histograms
		printf("=============== Subtracting Histograms =============== \n"); 

		hSelectDt->Sumw2();
		hRnPoDt = (TH1F*)hSelectDt->Clone();
		hRnPoDt->SetName(Form("hRnPoDt_%i",numTimeBin));
		hBGDt->Sumw2();
		hRnPoDt->Add(hBGDt,-1);

		hRnPSD = (TH1F*)hSelectPromptPSD->Clone();
		hRnPSD->SetName(Form("hRnPSD_%i",numTimeBin));
		hRnPSD->Sumw2();
		hRnPSD->Add(hBGPromptPSD,-1);

		hPoPSD = (TH1F*)hSelectDelayPSD->Clone();
		hPoPSD->SetName(Form("hPoPSD_%i",numTimeBin));
		hPoPSD->Sumw2();
		hPoPSD->Add(hBGDelayPSD,-1);

		hRnEn = (TH1F*)hSelectPromptEn->Clone();
		hRnEn->SetName(Form("hRnEn_%i",numTimeBin));
		hRnEn->Sumw2();
		hRnEn->Add(hBGPromptEn,-1);

		hRnTotEn = (TH1F*)hSelectPromptTotEn->Clone();
		hRnTotEn->SetName(Form("hRnTotEn_%i",numTimeBin));
		hRnTotEn->Sumw2();
		hRnTotEn->Add(hBGPromptTotEn,-1);

		hPoEn = (TH1F*)hSelectDelayEn->Clone();
		hPoEn->SetName(Form("hPoEn_%i",numTimeBin));
		hPoEn->Sumw2();
		hPoEn->Add(hBGDelayEn,-1);
	
		hRnPos = (TH1F*)hSelectPromptPos->Clone();
		hRnPos->SetName(Form("hRnPos_%i",numTimeBin));
		hRnPos->Sumw2();
		hRnPos->Add(hBGPromptPos,-1);

		hPoPos = (TH1F*)hSelectDelayPos->Clone();
		hPoPos->SetName(Form("hPoPos_%i",numTimeBin));
		hPoPos->Sumw2();
		hPoPos->Add(hBGDelayPos,-1);

		hRnPoDz = (TH1F*)hSelectDz->Clone();
		hRnPoDz->SetName(Form("hRnPoDz_%i",numTimeBin));
		hRnPoDz->Sumw2();
		hRnPoDz->Add(hBGDz,-1);

		hRnPoPSDvsEn = (TH2F*)hSelectPSDvsEn->Clone();
		hRnPoPSDvsEn->SetName(Form("hRnPoPSDvsEn_%i",numTimeBin));
		hRnPoPSDvsEn->Add(hBGPSDvsEn,-1);	

		hPoEnvsRnEn = (TH2F*)hSelectDelayEnvsPromptEn->Clone();
		hPoEnvsRnEn->SetName(Form("hPoEnvsRnEn_%i",numTimeBin));
		hPoEnvsRnEn->Add(hBGDelayEnvsPromptEn,-1);

		//---------------------------------------------------------------------------------
		//Calculate results
		printf("=============== Calculating Results =============== \n"); 

		double dtBinWidth = (dtMax - dtMin)/(double)numDtBins;

		double promptPSDEff, delayPSDEff, promptPSDEffErr, delayPSDEffErr;
		double promptEnEff,  delayEnEff,  promptEnEffErr,  delayEnEffErr;
		double dzEff, dzEffErr;
		double totEff, totEffErr;

		double NAlpha, NAlphaErr, lifetime, lifetimeErr;
		double rate, rateErr;		

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Fit distributions
		fRnPoDtExp = new TF1("fRnPoDtExp",Form("[0]*exp(-x/[1])*(%f/[1])",dtBinWidth),0.0,dtMax);
		//fRnPoDtExp = new TF1("fRnPoDtExp",Form("[0]*exp(-x/[1])*(%f/[1])",dtBinWidth),0.5,11);
		fRnPoDtExp->SetParameter(1,POLIFETIME);
		hRnPoDt->Fit(fRnPoDtExp,"R0");

		double fitPSDMin = promptLowPSDCut;
		double maxValue = hRnPSD->GetMaximum(), maxBin = hRnPSD->GetBinCenter(hRnPSD->GetMaximumBin());	
		fRnPSDGaus = new TF1("fRnPSDGaus","gaus(0)+gaus(3)",fitPSDMin,PSDMax);
		fRnPSDGaus->SetParameters(0.25*maxValue,maxBin+0.01,0.015,0.75*maxValue,maxBin-0.001,0.02);
		fRnPSDGaus->SetParLimits(0,0,maxValue);
		fRnPSDGaus->SetParLimits(3,0,maxValue);
		hRnPSD->Fit(fRnPSDGaus,"RQ0");
		fRnPSDGaus->SetRange(PSDMin,PSDMax);

		fitPSDMin = delayLowPSDCut;
		fPoPSDGaus = new TF1("fPoPSDGaus","gaus",fitPSDMin,PSDMax);
		hPoPSD->Fit(fPoPSDGaus,"RQ0");
		fPoPSDGaus->SetRange(PSDMin,PSDMax);
	
		double fitRnEnMin = EnMin;	
		fRnEnGaus = new TF1("fRnEnGaus","gaus(0)+gaus(3)+gaus(6)",fitRnEnMin,EnMax);
		maxValue = hRnEn->GetMaximum();
		maxBin = hRnEn->GetBinCenter(hRnEn->GetMaximumBin());	
		fRnEnGaus->SetParameters(0.8*maxValue,maxBin,0.035,0.05*maxValue,maxBin+0.14,0.08,0.15*maxValue,maxBin-0.03,0.05);
		hRnEn->Fit(fRnEnGaus,"RQ0");
		fRnEnGaus->SetRange(EnMin,EnMax);

		double fitPoEnMin = delayLowEnCut;
		fPoEnGaus = new TF1("fPoEnGaus","gaus",fitPoEnMin,EnMax);
		hPoEn->Fit(fPoEnGaus,"RQ0");
		fPoEnGaus->SetRange(EnMin,EnMax);

		fRnPoDzGaus = new TF1("fRnPoDzGaus","gaus",dzMin,dzMax);
		hRnPoDz->Fit(fRnPoDzGaus,"RQ0");

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Calculate efficiencies
		promptPSDEff = fRnPSDGaus->Integral(promptLowPSDCut,promptHighPSDCut)/fRnPSDGaus->Integral(PSDMin,PSDMax);
		promptPSDEffErr = sqrt((promptPSDEff*(1-promptPSDEff))/hRnPSD->GetEntries()); 
			
		delayPSDEff = fPoPSDGaus->Integral(delayLowPSDCut,delayHighPSDCut)/fPoPSDGaus->Integral(PSDMin,PSDMax);
		delayPSDEffErr = sqrt((delayPSDEff*(1-delayPSDEff))/hPoPSD->GetEntries()); 
		
		promptEnEff = fRnEnGaus->Integral(promptLowEnCut,promptHighEnCut)/fRnEnGaus->Integral(EnMin,EnMax);
		promptEnEffErr = sqrt((promptEnEff*(1-promptEnEff))/hRnEn->GetEntries());

		delayEnEff = fPoEnGaus->Integral(delayLowEnCut,delayHighEnCut)/fPoEnGaus->Integral(EnMin,EnMax);
		delayEnEffErr = sqrt((delayEnEff*(1-delayEnEff))/hPoEn->GetEntries()); 

		dzEff = fRnPoDzGaus->Integral(-dzCut,dzCut)/fRnPoDzGaus->Integral(dzMin,dzMax);
		dzEffErr = sqrt((dzEff*(1-dzEff))/hRnPoDz->GetEntries());

		
		totEff = promptPSDEff * delayPSDEff * promptEnEff * delayEnEff * dzEff;	
		totEffErr = totEff * sqrt( pow(promptPSDEffErr/promptPSDEff,2) + pow(delayPSDEffErr/delayPSDEff,2) + pow(promptEnEffErr/promptEnEff,2) + pow(delayEnEffErr/delayEnEff,2) + pow(dzEffErr/dzEff,2) );

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Calculate rate
		NAlpha = fRnPoDtExp->GetParameter(0);
		NAlphaErr = fRnPoDtExp->GetParError(0);
		lifetime = fRnPoDtExp->GetParameter(1);
		lifetimeErr = fRnPoDtExp->GetParError(1);
//		lifetimeErr = 0.007;
	
		rate = (NAlpha/(livetime*totEff))*(1e3);		//Hz
		rateErr = rate * sqrt(pow(NAlphaErr/NAlpha,2) + pow(totEffErr/totEff,2));

		printf("Rate: %.4f +/- %.4f \n",rate,rateErr);

		//---------------------------------------------------------------------------------
		//Populate vectors
		vRate.push_back(rate);
		vRateErr.push_back(rateErr);
		vTotEff.push_back(totEff);
		vTotEffErr.push_back(totEffErr);
		vLifetime.push_back(lifetime);
		vLifetimeErr.push_back(lifetimeErr);
		
		vPoPSDMean.push_back(fPoPSDGaus->GetParameter(1));	
		vPoPSDMeanErr.push_back(fPoPSDGaus->GetParError(1));
		vPoPSDSigma.push_back(fPoPSDGaus->GetParameter(2));
		vPoPSDSigmaErr.push_back(fPoPSDGaus->GetParError(2));

		vPoEnMean.push_back(fPoEnGaus->GetParameter(1));
		vPoEnMeanErr.push_back(fPoEnGaus->GetParError(1));
		vPoEnSigma.push_back(fPoEnGaus->GetParameter(2));
		vPoEnSigmaErr.push_back(fPoEnGaus->GetParError(2));

		vPoPosMean.push_back(hPoPos->GetMean());
		vPoPosMeanErr.push_back(hPoPos->GetMeanError());
		vPoPosSigma.push_back(hPoPos->GetRMS());
		vPoPosSigmaErr.push_back(hPoPos->GetRMSError());

		vRnPoDzMean.push_back(fRnPoDzGaus->GetParameter(1));
		vRnPoDzMeanErr.push_back(fRnPoDzGaus->GetParError(1));
		vRnPoDzSigma.push_back(fRnPoDzGaus->GetParameter(2));
		vRnPoDzSigmaErr.push_back(fRnPoDzGaus->GetParError(2));

		vPromptEnEff.push_back(promptEnEff);
		vPromptEnEffErr.push_back(promptEnEffErr);
		vDelayEnEff.push_back(delayEnEff);
		vDelayEnEffErr.push_back(delayEnEffErr);

		vPromptPSDEff.push_back(promptPSDEff);
		vPromptPSDEffErr.push_back(promptPSDEffErr);
		vDelayPSDEff.push_back(delayPSDEff);
		vDelayPSDEffErr.push_back(delayPSDEffErr);

		vDzEff.push_back(dzEff);
		vDzEffErr.push_back(dzEffErr);			

		vRnPSDChiSq.push_back(fRnPSDGaus->GetChisquare()/(double)fRnPSDGaus->GetNDF());
		vPoPSDChiSq.push_back(fPoPSDGaus->GetChisquare()/(double)fPoPSDGaus->GetNDF());
	
		vRnEnChiSq.push_back(fRnEnGaus->GetChisquare()/(double)fRnEnGaus->GetNDF());
		vPoEnChiSq.push_back(fPoEnGaus->GetChisquare()/(double)fPoEnGaus->GetNDF());
		
		vDzChiSq.push_back(fRnPoDzGaus->GetChisquare()/(double)fRnPoDzGaus->GetNDF());	
	
		vDtChiSq.push_back(fRnPoDtExp->GetChisquare()/(double)fRnPoDtExp->GetNDF());

		numTimeBin++;
	}	//end while loop IDX < numEntries

	endloop:
	histFile->Write();
	histFile->Close();

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	//---------------------------------------------------------------------------------
	//Initialize TGraphs for results

	int numPt = vRate.size();
	double x[numPt], xErr[numPt];
	double y[1], yErr[1];

	TGraphErrors *grRate 		= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grTotEff 		= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grLifetime 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoPSDMean 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoPSDSigma 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoEnMean   	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoEnSigma 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoPosMean 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoPosSigma 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grRnPoDzMean 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grRnPoDzSigma 	= new TGraphErrors(numPt,x,y,xErr,yErr);

	TGraph *grLivetime 	 = new TGraph(numPt,x,y);
	TGraph *grTotLivetime 	 = new TGraph(numPt,x,y);
	TGraph *grPileupVeto	 = new TGraph(numPt,x,y);
	TGraph *grPileupVetoFrac = new TGraph(numPt,x,y);

	TGraphErrors *grPromptEnEff  = new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grDelayEnEff   = new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPromptPSDEff = new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grDelayPSDEff  = new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grDzEff 	     = new TGraphErrors(numPt,x,y,xErr,yErr);

	TGraph *grRnPSDChiSq = new TGraph(numPt,x,y);
	TGraph *grPoPSDChiSq = new TGraph(numPt,x,y);
	TGraph *grRnEnChiSq  = new TGraph(numPt,x,y);
	TGraph *grPoEnChiSq  = new TGraph(numPt,x,y);
	TGraph *grDzChiSq    = new TGraph(numPt,x,y);
	TGraph *grDtChiSq    = new TGraph(numPt,x,y);

	TGraph *grBGRate = new TGraph(numPt,x,y);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//Fill TGraphs
	for(int i=0;i<numPt;i++){
		double time = vTimestamp[i];

		grRate->SetPoint(i,time,vRate[i]);
		grRate->SetPointError(i,0,vRateErr[i]);

		grTotEff->SetPoint(i,time,vTotEff[i]);
		grTotEff->SetPointError(i,0,vTotEffErr[i]);

		grLifetime->SetPoint(i,time,vLifetime[i]);
		grLifetime->SetPointError(i,0,vLifetimeErr[i]);

		grPoPSDMean->SetPoint(i,time,vPoPSDMean[i]);
		grPoPSDMean->SetPointError(i,0,vPoPSDMeanErr[i]);

		grPoPSDSigma->SetPoint(i,time,vPoPSDSigma[i]);
		grPoPSDSigma->SetPointError(i,0,vPoPSDSigmaErr[i]);

		grPoEnMean->SetPoint(i,time,vPoEnMean[i]);
		grPoEnMean->SetPointError(i,0,vPoEnMeanErr[i]);
	
		grPoEnSigma->SetPoint(i,time,vPoEnSigma[i]);
		grPoEnSigma->SetPointError(i,0,vPoEnSigmaErr[i]);

		grPoPosMean->SetPoint(i,time,vPoPosMean[i]);
		grPoPosMean->SetPointError(i,0,vPoPosMeanErr[i]);

		grPoPosSigma->SetPoint(i,time,vPoPosSigma[i]);
		grPoPosSigma->SetPointError(i,0,vPoPosSigmaErr[i]);		

		grRnPoDzMean->SetPoint(i,time,vRnPoDzMean[i]);
		grRnPoDzMean->SetPointError(i,0,vRnPoDzMeanErr[i]);

		grRnPoDzSigma->SetPoint(i,time,vRnPoDzSigma[i]);
		grRnPoDzSigma->SetPointError(i,0,vRnPoDzSigmaErr[i]);


		grLivetime->SetPoint(i,time,vLivetime[i]);
		grTotLivetime->SetPoint(i,time,vTotLivetime[i]);
		grPileupVeto->SetPoint(i,time,vPileupVetoT[i]);
		grPileupVetoFrac->SetPoint(i,time,vPileupVetoFrac[i]);

		grPromptEnEff->SetPoint(i,time,vPromptEnEff[i]);
		grPromptEnEff->SetPointError(i,0,vPromptEnEffErr[i]);

		grDelayEnEff->SetPoint(i,time,vDelayEnEff[i]);
		grDelayEnEff->SetPointError(i,0,vDelayEnEffErr[i]);

		grPromptPSDEff->SetPoint(i,time,vPromptPSDEff[i]);
		grPromptPSDEff->SetPointError(i,0,vPromptPSDEffErr[i]);

		grDelayPSDEff->SetPoint(i,time,vDelayPSDEff[i]);
		grDelayPSDEff->SetPointError(i,0,vDelayPSDEffErr[i]);

		grDzEff->SetPoint(i,time,vDzEff[i]);
		grDzEff->SetPointError(i,0,vDzEffErr[i]);	

		
		grRnPSDChiSq->SetPoint(i,time,vRnPSDChiSq[i]);
		grPoPSDChiSq->SetPoint(i,time,vPoPSDChiSq[i]);
		grRnEnChiSq->SetPoint(i,time,vRnEnChiSq[i]);
		grPoEnChiSq->SetPoint(i,time,vPoEnChiSq[i]);
		grDzChiSq->SetPoint(i,time,vDzChiSq[i]);
		grDtChiSq->SetPoint(i,time,vDtChiSq[i]);

		grBGRate->SetPoint(i,time,vBGRate[i]);

	}	//end for loop to populate TGraphs

	//---------------------------------------------------------------------------------
	//Write TGraphs to file
	TFile *graphFile = new TFile(Form("%s/Ac227_GraphsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")),"RECREATE");
	
	grRate->Write("grRate");
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
	grLivetime->Write("grLivetime");
	grTotLivetime->Write("grTotLivetime");
	grPileupVeto->Write("grPileupVeto");
	grPileupVetoFrac->Write("grPileupVetoFrac");
	grPromptEnEff->Write("grPromptEnEff");
	grDelayEnEff->Write("grDelayEnEff");
	grPromptPSDEff->Write("grPromptPSDEff");
	grDelayPSDEff->Write("grDelayPSDEff");
	grDzEff->Write("grDzEff");
	grRnPSDChiSq->Write("grRnPSDChiSq");
	grPoPSDChiSq->Write("grPoPSDChiSq");
	grRnEnChiSq->Write("grRnEnChiSq");
	grPoEnChiSq->Write("grPoEnChiSq");
	grDzChiSq->Write("grDzChiSq");
	grDtChiSq->Write("grDtChiSq");
	grBGRate->Write("grBGRate");

	graphFile->Close();

}	//end void RnPoVsTime

