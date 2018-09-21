#!/bin/bash

p_lowPSD=0.20
d_lowPSD=0.20
p_lowE=0.52
d_lowE=0.61
zLow=-1000
zHigh=1000
timeBin=100000000;
dtFit=1

#==============================================
echo ======================================	
echo ============= REACTOR ON =============
echo ======================================

export AD_AC227ANALYSIS_RESULTS=$AD_ANALYSIS/Results/RxOn
export AD_AC227_PLOTS=$AD_ANALYSIS/Plots/RxOn

#==============================================
echo ========== Making TAc Class ==========

root -l -b <<EOF
int RxStatus = 1;
.L Calculate/MakeAcTreeClass_RxOnVsOff.C+
MakeAcTreeClass(RxStatus)
.q
EOF

mv RNPO* Calculate/

echo ========= Running RnPoVsCell =========	
root -l -b <<EOF 
.L Calculate/RnPoVsCell.C+
RnPoVsCell($p_lowPSD,$d_lowPSD,$p_lowE,$d_lowE,$zLow,$zHigh,$dtFit)
.q
EOF

echo ========= Running RnPoVsTime =========
root -l -b <<EOF 
.L Calculate/RnPoVsTime.C+
RnPoVsTime($p_lowPSD,$d_lowPSD,$p_lowE,$d_lowE,$zLow,$zHigh,$timeBin,$dtFit)
.q
EOF

#==============================================
echo =======================================
echo ============= REACTOR OFF =============
echo =======================================

export AD_AC227ANALYSIS_RESULTS=$AD_ANALYSIS/Results/RxOff
export AD_AC227_PLOTS=$AD_ANALYSIS/Plots/RxOff

#==============================================
echo ========== Making TAc Class ==========

root -l -b <<EOF
int RxStatus = 0;
.L Calculate/MakeAcTreeClass_RxOnVsOff.C+
MakeAcTreeClass(RxStatus)
.q
EOF

mv RNPO* Calculate/

echo ========= Running RnPoVsCell =========
root -l -b <<EOF 
.L Calculate/RnPoVsCell.C+
RnPoVsCell($p_lowPSD,$d_lowPSD,$p_lowE,$d_lowE,$zLow,$zHigh,$dtFit)
.q
EOF

echo ========= Running RnPoVsTime =========
root -l -b <<EOF 
.L Calculate/RnPoVsTime.C+
RnPoVsTime($p_lowPSD,$d_lowPSD,$p_lowE,$d_lowE,$zLow,$zHigh,$timeBin,$dtFit)
.q
EOF


#==============================================
echo ======================================
echo ============== PLOTTING ==============
echo ======================================

export AD_AC227ANALYSIS_RESULTS=$AD_ANALYSIS/Results
export AD_AC227_PLOTS=$AD_ANALYSIS/Plots

echo ======= Running PlotRxOnVsRxOff =======
root -l -b <<EOF 
.L Plot/PlotRxOnVsRxOff.C+
PlotRxOnVsRxOff()
.q
EOF

