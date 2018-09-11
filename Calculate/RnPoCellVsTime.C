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


void RnPoCellVsTime(double p_lowPSD, double d_lowPSD, double p_lowE, double d_lowE, double zLow, double zHigh, double timeBin, int dtFit){
	TH1::AddDirectory(kFALSE);

	const double TIMEBREAK = timeBin*(3.6e6);	//[ms]

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

	//---------------------------------------------------------------------------------
	TF1 *fRnPoDtExp;
	TF1 *fRnPSDGaus, *fPoPSDGaus;
	TF1 *fRnEnCB,    *fPoEnGaus;
	TF1 *fRnPoDzGaus;

	//---------------------------------------------------------------------------------
	vector<double> vLivetime,   vTimestamp;
	vector<vector<double>> vRate,   vRateErr;
	vector<vector<double>> vTotEff, vTotEffErr;
	vector<vector<double>> vLifetime,   vLifetimeErr;
	vector<vector<double>> vPoPSDMean,  vPoPSDMeanErr,  vPoPSDSigma,  vPoPSDSigmaErr;
	vector<vector<double>> vPoEnMean,   vPoEnMeanErr,   vPoEnSigma,   vPoEnSigmaErr;
	vector<vector<double>> vPoPosMean,  vPoPosMeanErr,  vPoPosSigma,  vPoPosSigmaErr;
	vector<vector<double>> vRnPoDzMean, vRnPoDzMeanErr, vRnPoDzSigma, vRnPoDzSigmaErr;

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	//---------------------------------------------------------------------------------
	//Loop over all events, splitting into time periods of 24 hours

	Long64_t IDX = 0;
	int numTimeBin = 0;
	while(IDX<numEntries){

		//---------------------------------------------------------------------------------
		//Initialize histograms
		int N = NUMCELLS;

		TH1F *hSelectDt[N],		   *hBGDt[N],		 *hRnPoDt[N];
		TH1F *hSelectPromptPSD[N], *hBGPromptPSD[N], *hRnPSD[N];
		TH1F *hSelectDelayPSD[N],  *hBGDelayPSD[N],  *hPoPSD[N];
		TH1F *hSelectPromptEn[N],  *hBGPromptEn[N],  *hRnEn[N];
		TH1F *hSelectDelayEn[N],   *hBGDelayEn[N],   *hPoEn[N];
		TH1F *hSelectPromptPos[N], *hBGPromptPos[N], *hRnPos[N];
		TH1F *hSelectDelayPos[N],  *hBGDelayPos[N],  *hPoPos[N];
		TH1F *hSelectDz[N], 	   *hBGDz[N], 		 *hRnPoDz[N];
	
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
		}	//end for loop creating histograms


		//---------------------------------------------------------------------------------
		//Fill histograms
		printf("=============== Filling Histograms =============== \n"); 
	
		int seg;
		double dt, dz;
		
		double livetime = 0.0, tstamp;
		double lastTime = 0.0, lastRunTime = 0.0, lastTimestamp = 0.0;	
		double sumWeightedTimestamp = 0.0, sumRunTime = 0.0;

		double lastOCSTime = 0.0, OCSTime = 0.0;
		double lastMuonVetoTime = 0.0, muonVetoTime = 0.0;

		double lastNumClusts = 0.0, numClusts = 0.0;

		for(Long64_t i=IDX;i<numEntries;i++){
			if(i%100000==0) printf("Event: %lld \n",i);
			rnpo->GetEntry(i);

			if(rnpo->d_t*(1e-6) > ((double)((TVectorD*)rnpo->fChain->GetCurrentFile()->Get("runtime"))->Norm1()*1000.0 - (TIMEWINDOW+TIMEOFFSET)) && i!=(numEntries-1)) continue;
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			//Calculate livetime and weighted timestamp
			if(rnpo->d_t < lastTime){ 
				livetime += lastTime*(1e-6);		//livetime in ms	

				OCSTime += lastOCSTime*(1e-6);
				muonVetoTime += lastMuonVetoTime*(1e-6);
	
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
			lastOCSTime = rnpo->OCSVeto_t;
			lastMuonVetoTime = rnpo->muonVeto_t;
			lastRunTime = ((TVectorD*)rnpo->fChain->GetCurrentFile()->Get("runtime"))->Norm1();		//seconds
			lastTimestamp = rnpo->tstamp;			//epoch seconds	
			lastNumClusts = rnpo->numClust;

			//if we are at the last entry	
			if(i == (numEntries-1)){ 
				livetime += lastTime*(1e-6);

				//if livetime is less than 12 hours 
				if(livetime*(2.778e-7) < 13) goto endloop;

				OCSTime += lastOCSTime*(1e-6);
				muonVetoTime += lastMuonVetoTime*(1e-6);
	
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
			double rnpo_p_E = rnpo->p_E;	
			double rnpo_d_E = rnpo->d_E;
			double rnpo_f_E = rnpo->f_E;
			exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), seg) != end(ExcludeCellArr);
			if(exclude) continue;

			if(rnpo->d_PSD < delayLowPSDCut || rnpo_d_E < delayLowEnCut) continue;	
			if(rnpo->d_z < zLow || rnpo->d_z > zHigh) continue;

			if(rnpo->p_seg > -1 && rnpo->p_PSD>promptLowPSDCut && rnpo_p_E>promptLowEnCut && rnpo->p_z>zLow && rnpo->p_z<zHigh){	
				dt = (rnpo->d_t - rnpo->p_t)*(1e-6);	//convert ns to ms	
				dz = rnpo->d_z - rnpo->p_z;

				hSelectDt[seg]->Fill(dt);
				hSelectPromptPSD[seg]->Fill(rnpo->p_PSD);
				hSelectDelayPSD[seg]->Fill(rnpo->d_PSD);
				hSelectPromptEn[seg]->Fill(rnpo_p_E);
				hSelectDelayEn[seg]->Fill(rnpo_d_E);
				hSelectPromptPos[seg]->Fill(rnpo->p_z);
				hSelectDelayPos[seg]->Fill(rnpo->d_z);
				hSelectDz[seg]->Fill(dz);
			}
			if(rnpo->f_seg > -1 && rnpo->f_PSD>promptLowPSDCut && rnpo_f_E>promptLowEnCut && rnpo->f_z>zLow && rnpo->f_z<zHigh){	
				dt = (rnpo->f_t - rnpo->d_t)*(1e-6) - TIMEOFFSET;	
				dz = rnpo->d_z - rnpo->f_z;

				hBGDt[seg]->Fill(dt);
				hBGPromptPSD[seg]->Fill(rnpo->f_PSD);
				hBGDelayPSD[seg]->Fill(rnpo->d_PSD);
				hBGPromptEn[seg]->Fill(rnpo_f_E);
				hBGDelayEn[seg]->Fill(rnpo_d_E);
				hBGPromptPos[seg]->Fill(rnpo->f_z);
				hBGDelayPos[seg]->Fill(rnpo->d_z);
				hBGDz[seg]->Fill(dz);
			}

		}	//end for loop over events

		double pileupVetoTime = numClusts*pileupVetoT;	//[ms]
		double pileupVetoCorr = (2.0*pileupVetoTime)/(livetime);

		printf("Time bin: %i  |  Livetime: %f hrs \n",numTimeBin,livetime*(2.778e-7));
		printf("OCS time: %f ms   |   Muon time: %f ms \n",OCSTime,muonVetoTime);
		printf("Pileup veto correction: %f \n",pileupVetoCorr);

		livetime = livetime - OCSTime;
		livetime = livetime - 2.0*muonVetoTime;
		livetime = livetime*(1-pileupVetoCorr);
		vLivetime.push_back(livetime);

		printf("Corrected Livetime: %f hours \n",livetime*(2.778e-7));

		//---------------------------------------------------------------------------------
		//Subtract histograms
		printf("=============== Subtracting Histograms =============== \n"); 

		for(int i=0;i<NUMCELLS;i++){
			exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), i) != end(ExcludeCellArr);
			if(exclude) continue;

			hRnPoDt[i] = (TH1F*)hSelectDt[i]->Clone();
			hRnPoDt[i]->SetName(Form("hRnPoDt_%i",numTimeBin));
			hRnPoDt[i]->Sumw2();
			hBGDt[i]->Sumw2();
			hRnPoDt[i]->Add(hBGDt[i],-1);
	
			hRnPSD[i] = (TH1F*)hSelectPromptPSD[i]->Clone();
			hRnPSD[i]->SetName(Form("hRnPSD_%i",numTimeBin));
			hRnPSD[i]->Sumw2();
			hBGPromptPSD[i]->Sumw2();
			hRnPSD[i]->Add(hBGPromptPSD[i],-1);

			hPoPSD[i] = (TH1F*)hSelectDelayPSD[i]->Clone();
			hPoPSD[i]->SetName(Form("hPoPSD_%i",numTimeBin));
			hPoPSD[i]->Sumw2();
			hBGDelayPSD[i]->Sumw2();
			hPoPSD[i]->Add(hBGDelayPSD[i],-1);
	
			hRnEn[i] = (TH1F*)hSelectPromptEn[i]->Clone();
			hRnEn[i]->SetName(Form("hRnEn_%i",numTimeBin));
			hRnEn[i]->Sumw2();
			hBGPromptEn[i]->Sumw2();
			hRnEn[i]->Add(hBGPromptEn[i],-1);
	
			hPoEn[i] = (TH1F*)hSelectDelayEn[i]->Clone();
			hPoEn[i]->SetName(Form("hPoEn_%i",numTimeBin));
			hPoEn[i]->Sumw2();
			hBGDelayEn[i]->Sumw2();
			hPoEn[i]->Add(hBGDelayEn[i],-1);
		
			hRnPos[i] = (TH1F*)hSelectPromptPos[i]->Clone();
			hRnPos[i]->SetName(Form("hRnPos_%i",numTimeBin));
			hRnPos[i]->Sumw2();
			hBGPromptPos[i]->Sumw2();
			hRnPos[i]->Add(hBGPromptPos[i],-1);
	
			hPoPos[i] = (TH1F*)hSelectDelayPos[i]->Clone();
			hPoPos[i]->SetName(Form("hPoPos_%i",numTimeBin));
			hPoPos[i]->Sumw2();
			hBGDelayPos[i]->Sumw2();
			hPoPos[i]->Add(hBGDelayPos[i],-1);
	
			hRnPoDz[i] = (TH1F*)hSelectDz[i]->Clone();
			hRnPoDz[i]->SetName(Form("hRnPoDz_%i",numTimeBin));
			hRnPoDz[i]->Sumw2();
			hBGDz[i]->Sumw2();
			hRnPoDz[i]->Add(hBGDz[i],-1);
		}	//end for loop subtracting histograms

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

		//---------------------------------------------------------------------------------
		vector<double> vRate_t, 	  vRateErr_t;
		vector<double> vTotEff_t, 	  vTotEffErr_t;
		vector<double> vLifetime_t,   vLifetimeErr_t;
		vector<double> vPoPSDMean_t,  vPoPSDMeanErr_t,  vPoPSDSigma_t,  vPoPSDSigmaErr_t;
		vector<double> vPoEnMean_t,   vPoEnMeanErr_t,   vPoEnSigma_t,   vPoEnSigmaErr_t;
		vector<double> vPoPosMean_t,  vPoPosMeanErr_t,  vPoPosSigma_t,  vPoPosSigmaErr_t;
		vector<double> vRnPoDzMean_t, vRnPoDzMeanErr_t, vRnPoDzSigma_t, vRnPoDzSigmaErr_t; 

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Fit distributions
		for(int i=0;i<NUMCELLS;i++){
			exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), i) != end(ExcludeCellArr);
			if(exclude){ 
				double nullPt = 0.0;
				vRate_t.push_back(nullPt);
				vRateErr_t.push_back(nullPt);
				vTotEff_t.push_back(nullPt);
				vTotEffErr_t.push_back(nullPt);
				vLifetime_t.push_back(nullPt);
				vLifetimeErr_t.push_back(nullPt);
			
				vPoPSDMean_t.push_back(nullPt);	
				vPoPSDMeanErr_t.push_back(nullPt);
				vPoPSDSigma_t.push_back(nullPt);
				vPoPSDSigmaErr_t.push_back(nullPt);
	
				vPoEnMean_t.push_back(nullPt);
				vPoEnMeanErr_t.push_back(nullPt);
				vPoEnSigma_t.push_back(nullPt);
				vPoEnSigmaErr_t.push_back(nullPt);
	
				vPoPosMean_t.push_back(nullPt);
				vPoPosMeanErr_t.push_back(nullPt);
				vPoPosSigma_t.push_back(nullPt);
				vPoPosSigmaErr_t.push_back(nullPt);
	
				vRnPoDzMean_t.push_back(nullPt);
				vRnPoDzMeanErr_t.push_back(nullPt);
				vRnPoDzSigma_t.push_back(nullPt);
				vRnPoDzSigmaErr_t.push_back(nullPt);

				continue;
			}
	
			fRnPoDtExp = new TF1("fRnPoDtExp",Form("[0]*exp(-x/[1])*(%f/[1])",dtBinWidth),dtBinWidth*dtFit,dtMax);
			fRnPoDtExp->SetParameter(1,POLIFETIME);
			hRnPoDt[i]->Fit(fRnPoDtExp,"RQ0");
	
			double fitPSDMin = 0.20;
			fRnPSDGaus = new TF1("fRnPSDGaus","gaus",fitPSDMin,PSDMax);
			hRnPSD[i]->Fit(fRnPSDGaus,"RQ0");
			fRnPSDGaus->SetRange(PSDMin,PSDMax);

			fPoPSDGaus = new TF1("fPoPSDGaus","gaus",fitPSDMin,PSDMax);
			hPoPSD[i]->Fit(fPoPSDGaus,"RQ0");
			fPoPSDGaus->SetRange(PSDMin,PSDMax);
	
			double fitRnEnMin = 0.57;	
			fRnEnCB = new TF1("fRnEnCB","crystalball",fitRnEnMin,EnMax);
			fRnEnCB->SetParameter(1,0.7);
			fRnEnCB->SetParameter(2,0.04);
			fRnEnCB->SetParameter(3,-1.7);
			fRnEnCB->SetParameter(4,1.2);
			hRnEn[i]->Fit(fRnEnCB,"RQ0");
			fRnEnCB->SetRange(EnMin,EnMax);
	
			double fitPoEnMin = 0.7;
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
		
			promptEnEff = fRnEnCB->Integral(promptLowEnCut,promptHighEnCut)/fRnEnCB->Integral(EnMin,EnMax);
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

			//---------------------------------------------------------------------------------
			//Populate vectors
			vRate_t.push_back(rate);
			vRateErr_t.push_back(rateErr);
			vTotEff_t.push_back(totEff);
			vTotEffErr_t.push_back(totEffErr);
			vLifetime_t.push_back(lifetime);
			vLifetimeErr_t.push_back(lifetimeErr);
		
			vPoPSDMean_t.push_back(fPoPSDGaus->GetParameter(1));	
			vPoPSDMeanErr_t.push_back(fPoPSDGaus->GetParError(1));
			vPoPSDSigma_t.push_back(fPoPSDGaus->GetParameter(2));
			vPoPSDSigmaErr_t.push_back(fPoPSDGaus->GetParError(2));

			vPoEnMean_t.push_back(fPoEnGaus->GetParameter(1));
			vPoEnMeanErr_t.push_back(fPoEnGaus->GetParError(1));
			vPoEnSigma_t.push_back(fPoEnGaus->GetParameter(2));
			vPoEnSigmaErr_t.push_back(fPoEnGaus->GetParError(2));

			vPoPosMean_t.push_back(hPoPos[i]->GetMean());
			vPoPosMeanErr_t.push_back(hPoPos[i]->GetMeanError());
			vPoPosSigma_t.push_back(hPoPos[i]->GetRMS());
			vPoPosSigmaErr_t.push_back(hPoPos[i]->GetRMSError());

			vRnPoDzMean_t.push_back(fRnPoDzGaus->GetParameter(1));
			vRnPoDzMeanErr_t.push_back(fRnPoDzGaus->GetParError(1));
			vRnPoDzSigma_t.push_back(fRnPoDzGaus->GetParameter(2));
			vRnPoDzSigmaErr_t.push_back(fRnPoDzGaus->GetParError(2));

		}	//end for loop fitting distributions

		//---------------------------------------------------------------------------------
		//Populate vectors
		vRate.push_back(vRate_t);
		vRateErr.push_back(vRateErr_t);
		vTotEff.push_back(vTotEff_t);
		vTotEffErr.push_back(vTotEffErr_t);
		vLifetime.push_back(vLifetime_t);
		vLifetimeErr.push_back(vLifetimeErr_t);

		vPoPSDMean.push_back(vPoPSDMean_t);
		vPoPSDMeanErr.push_back(vPoPSDMeanErr_t);
		vPoPSDSigma.push_back(vPoPSDSigma_t);
		vPoPSDSigmaErr.push_back(vPoPSDSigmaErr_t);

		vPoEnMean.push_back(vPoEnMean_t);
		vPoEnMeanErr.push_back(vPoEnMeanErr_t);
		vPoEnSigma.push_back(vPoEnSigma_t);
		vPoEnSigmaErr.push_back(vPoEnSigmaErr_t);
		
		vPoPosMean.push_back(vPoPosMean_t);
		vPoPosMeanErr.push_back(vPoPosMeanErr_t);
		vPoPosSigma.push_back(vPoPosSigma_t);
		vPoPosSigmaErr.push_back(vPoPosSigmaErr_t);
		
		vRnPoDzMean.push_back(vRnPoDzMean_t);
		vRnPoDzMeanErr.push_back(vRnPoDzMeanErr_t);
		vRnPoDzSigma.push_back(vRnPoDzSigma_t);
		vRnPoDzSigmaErr.push_back(vRnPoDzSigmaErr_t);

		numTimeBin++;

		vector<double>().swap(vRate_t);
		vector<double>().swap(vRateErr_t);
		vector<double>().swap(vTotEff_t);
		vector<double>().swap(vTotEffErr_t);
		vector<double>().swap(vLifetime_t);
		vector<double>().swap(vLifetimeErr_t);
	
		vector<double>().swap(vPoPSDMean_t);
		vector<double>().swap(vPoPSDMeanErr_t);
		vector<double>().swap(vPoPSDSigma_t);
		vector<double>().swap(vPoPSDSigmaErr_t);
	
		vector<double>().swap(vPoEnMean_t);
		vector<double>().swap(vPoEnMeanErr_t);
		vector<double>().swap(vPoEnSigma_t);
		vector<double>().swap(vPoEnSigmaErr_t);
	
		vector<double>().swap(vPoPosMean_t);
		vector<double>().swap(vPoPosMeanErr_t);
		vector<double>().swap(vPoPosSigma_t);
		vector<double>().swap(vPoPosSigmaErr_t);
	
		vector<double>().swap(vRnPoDzMean_t);
		vector<double>().swap(vRnPoDzMeanErr_t);
		vector<double>().swap(vRnPoDzSigma_t);
		vector<double>().swap(vRnPoDzSigmaErr_t);

		for(int i=0;i<NUMCELLS;i++){
			exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), i) != end(ExcludeCellArr);
			if(exclude) continue;

			delete hSelectDt[i];
			delete hBGDt[i];
			delete hRnPoDt[i];
	
			delete hSelectPromptPSD[i];
			delete hBGPromptPSD[i];
			delete hRnPSD[i];
			delete hSelectDelayPSD[i];
			delete hBGDelayPSD[i];
			delete hPoPSD[i];
	
			delete hSelectPromptEn[i];
			delete hBGPromptEn[i];
			delete hRnEn[i];
			delete hSelectDelayEn[i];
			delete hBGDelayEn[i];
			delete hPoEn[i];
	
			delete hSelectPromptPos[i];
			delete hBGPromptPos[i];
			delete hRnPos[i];
			delete hSelectDelayPos[i];
			delete hBGDelayPos[i];
			delete hPoPos[i];
	
			delete hSelectDz[i];
			delete hBGDz[i];
			delete hRnPoDz[i];
		}
	
	}	//end while loop IDX < numEntries
	endloop:

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	//---------------------------------------------------------------------------------
	//Initialize TGraphs for results

	int numPt = vLivetime.size();
	double x[numPt], xErr[numPt];
	double y[1], yErr[1];

	int N = NUMCELLS;
	TGraphErrors *grRate[N],	   *grTotEff[N],	 *grLifetime[N];
	TGraphErrors *grPoPSDMean[N],  *grPoPSDSigma[N];
	TGraphErrors *grPoEnMean[N],   *grPoEnSigma[N];
	TGraphErrors *grPoPosMean[N],  *grPoPosSigma[N];
	TGraphErrors *grRnPoDzMean[N], *grRnPoDzSigma[N];

	for(int i=0;i<NUMCELLS;i++){
		grRate[i] 		 = new TGraphErrors(numPt,x,y,xErr,yErr);
		grTotEff[i] 	 = new TGraphErrors(numPt,x,y,xErr,yErr);
		grLifetime[i] 	 = new TGraphErrors(numPt,x,y,xErr,yErr);
		grPoPSDMean[i] 	 = new TGraphErrors(numPt,x,y,xErr,yErr);
		grPoPSDSigma[i]  = new TGraphErrors(numPt,x,y,xErr,yErr);
		grPoEnMean[i]    = new TGraphErrors(numPt,x,y,xErr,yErr);
		grPoEnSigma[i] 	 = new TGraphErrors(numPt,x,y,xErr,yErr);
		grPoPosMean[i] 	 = new TGraphErrors(numPt,x,y,xErr,yErr);
		grPoPosSigma[i]  = new TGraphErrors(numPt,x,y,xErr,yErr);
		grRnPoDzMean[i]  = new TGraphErrors(numPt,x,y,xErr,yErr);
		grRnPoDzSigma[i] = new TGraphErrors(numPt,x,y,xErr,yErr);
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//Fill TGraphs
	for(int i=0;i<numPt;i++){
		double time = vTimestamp[i];
		for(int j=0;j<NUMCELLS;j++){
			grRate[j]->SetPoint(i,time,vRate[i][j]);
			grRate[j]->SetPointError(i,0,vRateErr[i][j]);

			grTotEff[j]->SetPoint(i,time,vTotEff[i][j]);
			grTotEff[j]->SetPointError(i,0,vTotEffErr[i][j]);
	
			grLifetime[j]->SetPoint(i,time,vLifetime[i][j]);
			grLifetime[j]->SetPointError(i,0,vLifetimeErr[i][j]);
	
			grPoPSDMean[j]->SetPoint(i,time,vPoPSDMean[i][j]);
			grPoPSDMean[j]->SetPointError(i,0,vPoPSDMeanErr[i][j]);
	
			grPoPSDSigma[j]->SetPoint(i,time,vPoPSDSigma[i][j]);
			grPoPSDSigma[j]->SetPointError(i,0,vPoPSDSigmaErr[i][j]);
	
			grPoEnMean[j]->SetPoint(i,time,vPoEnMean[i][j]);
			grPoEnMean[j]->SetPointError(i,0,vPoEnMeanErr[i][j]);
		
			grPoEnSigma[j]->SetPoint(i,time,vPoEnSigma[i][j]);
			grPoEnSigma[j]->SetPointError(i,0,vPoEnSigmaErr[i][j]);
	
			grPoPosMean[j]->SetPoint(i,time,vPoPosMean[i][j]);
			grPoPosMean[j]->SetPointError(i,0,vPoPosMeanErr[i][j]);
	
			grPoPosSigma[j]->SetPoint(i,time,vPoPosSigma[i][j]);
			grPoPosSigma[j]->SetPointError(i,0,vPoPosSigmaErr[i][j]);		
	
			grRnPoDzMean[j]->SetPoint(i,time,vRnPoDzMean[i][j]);
			grRnPoDzMean[j]->SetPointError(i,0,vRnPoDzMeanErr[i][j]);
	
			grRnPoDzSigma[j]->SetPoint(i,time,vRnPoDzSigma[i][j]);
			grRnPoDzSigma[j]->SetPointError(i,0,vRnPoDzSigmaErr[i][j]);
		}
	}

	//---------------------------------------------------------------------------------
	//Write TGraphs to file
	TFile *graphFile = new TFile(Form("%s/Ac227_GraphsCellvsTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")),"RECREATE");

	for(int i=0;i<NUMCELLS;i++){
		grRate[i]->Write(Form("grRate_%i",i));
	}
	
	graphFile->Close();

}	//end void RnPoVsTime

