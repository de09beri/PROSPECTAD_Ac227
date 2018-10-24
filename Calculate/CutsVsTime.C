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
#include <fstream>

#include "Header.C"

void CutsVsTime(double zLow, double zHigh, double timeBin){

	const double TIMEBREAK = timeBin*(3.6e6);	//[ms]
	
	ofstream cutFile("CutParameterVsTime.txt");

	//---------------------------------------------------------------------------------
	TH1F *hSelectPromptPSD, *hBGPromptPSD, *hRnPSD,	 *hSelectDelayPSD, *hBGDelayPSD, *hPoPSD;
	TH1F *hSelectPromptEn,  *hBGPromptEn,  *hRnEn,   *hSelectDelayEn,  *hBGDelayEn,  *hPoEn;

	//---------------------------------------------------------------------------------
	TF1 *fRnPSDGaus, *fPoPSDGaus;
	TF1 *fRnEnGaus,  *fPoEnGaus;

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
	promptLowPSDCut  = rnpo->p_PSDCut[0];
	promptHighPSDCut = rnpo->p_PSDCut[1];
	delayLowPSDCut   = rnpo->d_PSDCut[0];
	delayHighPSDCut  = rnpo->d_PSDCut[1];
	promptLowEnCut   = rnpo->p_ECut[0];
    	promptHighEnCut  = rnpo->p_ECut[1];
    	delayLowEnCut    = rnpo->d_ECut[0];
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
		hSelectPromptPSD 	= new TH1F(Form("hSelectPromptPSD_%i",numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	
		hBGPromptPSD 		= new TH1F(Form("hBGPromptPSD_%i",numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	

		hSelectDelayPSD 	= new TH1F(Form("hSelectDelayPSD_%i",numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	
		hBGDelayPSD 		= new TH1F(Form("hBGDelayPSD_%i",numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	

		hSelectPromptEn 	= new TH1F(Form("hSelectPromptEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
		hBGPromptEn 		= new TH1F(Form("hBGPromptEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
	
		hSelectDelayEn 		= new TH1F(Form("hSelectDelayEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
		hBGDelayEn 		= new TH1F(Form("hBGDelayEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	

		//---------------------------------------------------------------------------------
		//Fill histograms
		printf("=============== Filling Histograms =============== \n"); 
	
		int seg;
		double dt, dz;
		
		double livetime = 0.0, tstamp;
		double lastTime = 0.0, lastRunTime = 0.0, lastTimestamp = 0.0;	
		double sumWeightedTimestamp = 0.0, sumRunTime = 0.0;

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
			}

			if(livetime>TIMEBREAK){
				tstamp = sumWeightedTimestamp/sumRunTime;

				IDX = i;
				break;
			}

			lastTime = rnpo->d_t;
			lastRunTime = ((TVectorD*)rnpo->fChain->GetCurrentFile()->Get("runtime"))->Norm1();		//seconds
			lastTimestamp = rnpo->tstamp;			//epoch seconds	

			//if we are at the last entry	
			if(i == (numEntries-1)){ 
				livetime += lastTime*(1e-6);

				//if livetime is less than 12 hours 
				if(livetime*(2.778e-7) < 12) goto endloop;

				sumWeightedTimestamp += lastRunTime * ((lastRunTime/2.0)+lastTimestamp);
				sumRunTime += lastRunTime;

				tstamp = sumWeightedTimestamp/sumRunTime;

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
				hSelectPromptPSD->Fill(rnpo->p_PSD);
				hSelectDelayPSD->Fill(rnpo->d_PSD);
				hSelectPromptEn->Fill(rnpo_p_E);
				hSelectDelayEn->Fill(rnpo_d_E);
			}
			dt = (rnpo->f_t - rnpo->d_t)*(1e-6) - TIMEOFFSET;	
			if(rnpo->f_seg > -1 && rnpo->f_PSD>promptLowPSDCut && rnpo->f_E>promptLowEnCut && rnpo->f_z>zLow && rnpo->f_z<zHigh && dt>0.5){	
				hBGPromptPSD->Fill(rnpo->f_PSD);
				hBGDelayPSD->Fill(rnpo->d_PSD);
				hBGPromptEn->Fill(rnpo_f_E);
				hBGDelayEn->Fill(rnpo_d_E);
			}

		}	//end for loop over events

		printf("Time bin: %i  |  Livetime: %f hrs \n",numTimeBin,livetime*(2.778e-7));

		//---------------------------------------------------------------------------------
		//Subtract histograms
		printf("=============== Subtracting Histograms =============== \n"); 

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

		hPoEn = (TH1F*)hSelectDelayEn->Clone();
		hPoEn->SetName(Form("hPoEn_%i",numTimeBin));
		hPoEn->Sumw2();
		hPoEn->Add(hBGDelayEn,-1);
	
		//---------------------------------------------------------------------------------
		//Fit Histograms
		printf("=============== Fitting Histograms =============== \n"); 

		double RnPSD, RnPSDSigma, PoPSD, PoPSDSigma;
		double RnEn,  RnEnSigma,  PoEn,  PoEnSigma;

		double RnPSDMaxBin = hRnPSD->GetBinCenter(hRnPSD->GetMaximumBin());
		fRnPSDGaus = new TF1("fRnPSDGaus","gaus",RnPSDMaxBin-0.02,RnPSDMaxBin+0.02);
		hRnPSD->Fit(fRnPSDGaus,"RQ");
		RnPSD = fRnPSDGaus->GetParameter(1);
		RnPSDSigma = fRnPSDGaus->GetParameter(2);	

		double PoPSDMaxBin = hPoPSD->GetBinCenter(hPoPSD->GetMaximumBin());
		fPoPSDGaus = new TF1("fPoPSDGaus","gaus",PoPSDMaxBin-0.02,PoPSDMaxBin+0.02);
		hPoPSD->Fit(fPoPSDGaus,"RQ");
		PoPSD = fPoPSDGaus->GetParameter(1);
		PoPSDSigma = fPoPSDGaus->GetParameter(2);	

		double RnEnMaxBin = hRnEn->GetBinCenter(hRnEn->GetMaximumBin());
		fRnEnGaus = new TF1("fRnEnGaus","gaus",RnEnMaxBin-0.04,RnEnMaxBin+0.04);	
		hRnEn->Fit(fRnEnGaus,"RQ");
		RnEn = fRnEnGaus->GetParameter(1);
		RnEnSigma = fRnEnGaus->GetParameter(2);

		double PoEnMaxBin = hPoEn->GetBinCenter(hPoEn->GetMaximumBin());
		fPoEnGaus = new TF1("fPoEnGaus","gaus",PoEnMaxBin-0.04,PoEnMaxBin+0.04);	
		hPoEn->Fit(fPoEnGaus,"RQ");
		PoEn = fPoEnGaus->GetParameter(1);
		PoEnSigma = fPoEnGaus->GetParameter(2);
	

		cutFile<<std::fixed<<numTimeBin<<" "<<std::setprecision(0)<<tstamp<<" "<<std::setprecision(4)<<RnPSD<<" "<<std::setprecision(4)<<RnPSDSigma<<" "<<std::setprecision(4)<<PoPSD<<" "<<std::setprecision(4)<<PoPSDSigma<<" "<<std::setprecision(4)<<RnEn<<" "<<std::setprecision(4)<<RnEnSigma<<" "<<std::setprecision(4)<<PoEn<<" "<<std::setprecision(4)<<PoEnSigma<<"\n";	


		numTimeBin++;
	}	//end while loop IDX < numEntries

	endloop:

	cutFile.close();

}	//end void RnPoVsTime

