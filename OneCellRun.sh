#!/bin/bash

# Run the complete analysis after the plugin stage of the Ac227 analysis

mkdir $AD_AC227_PLOTS/Cell104
export AD_AC227_PLOTS=$AD_AC227_PLOTS/Cell104

p_lowPSD=0.19
d_lowPSD=0.19
p_lowE=0.57
d_lowE=0.66
zLow=-1000
zHigh=1000
#timeBin=23.5
timeBin=47.5
dtFit=1


#==============================================
# Make Class
#echo ========== Making TAc Class ==========

#root -l -b <<EOF
#.L Calculate/MakeAcTreeClass.C+
#MakeAcTreeClass()
#.q
#EOF

#mv RNPO* Calculate/

#==============================================
# Run RnPoVsTime
echo ========= Running RnPoVsTime =========

root -l -b <<EOF 
.L Calculate/RnPoVsTime.C+
RnPoVsTime($p_lowPSD,$d_lowPSD,$p_lowE,$d_lowE,$zLow,$zHigh,$timeBin,$dtFit)
.q
EOF

#==============================================
# Run RnPoVsTime.C 
echo ======= Running PlotRnPoVsTime =======

root -l -b <<EOF
.L Plot/PlotRnPoVsTime.C+
PlotRnPoVsTime()
.q
EOF



