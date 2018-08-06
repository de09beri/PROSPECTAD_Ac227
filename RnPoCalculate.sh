#!/bin/bash

# ./RnPoCalculate.sh $1 $2 $3 $4
# $1: Run MakeClass
# $2: Run RnPoVsCell.C
# $3: Run RnPoVsTime.C
# $4: Run RnPoColVsTime.C and RnPoRowVsTime.C


p_lowPSD=0.19
d_lowPSD=0.19
p_lowE=0.57
d_lowE=0.66
zLow=-1000
zHigh=1000
timeBin=23.5
dtFit=1

#==============================================
# Make Class
echo ========== Making TAc Class ==========

if [ $1 -eq 1 ]
then

root -l -b <<EOF
.L Calculate/MakeAcTreeClass.C+
MakeAcTreeClass()
.q
EOF

fi

#==============================================
# Run RnPoVsCell
echo ========= Running RnPoVsCell =========

if [ $2 -eq 1 ]
then

root -l -b <<EOF 
.L Calculate/RnPoVsCell.C+
RnPoVsCell($p_lowPSD,$d_lowPSD,$p_lowE,$d_lowE,$zLow,$zHigh,$dtFit)
.q
EOF

fi

#==============================================
# Run RnPoVsTime
echo ========= Running RnPoVsTime =========

if [ $3 -eq 1]
then

root -l -b <<EOF 
.L Calculate/RnPoVsTime.C+
RnPoVsTime($p_lowPSD,$d_lowPSD,$p_lowE,$d_lowE,$zLow,$zHigh,$timeBin,$dtFit)
.q
EOF

fi

#==============================================
# Run RnPoColVsTime and RnPoRunVsTime
echo ========= Running RnPoVsTime =========

if [ $4 -eq 1]
then

echo ======= Running RnPoColVsTime =======

root -l -b <<EOF 
.L Calculate/RnPoColVsTime.C+
RnPoColVsTime($p_lowPSD,$d_lowPSD,$p_lowE,$d_lowE,$timeBin,$dtFit)
.q
EOF

echo ======= Running RnPoRowVsTime =======

root -l -b <<EOF 
.L Calculate/RnPoRowVsTime.C+
RnPoRowVsTime($p_lowPSD,$d_lowPSD,$p_lowE,$d_lowE,$timeBin,$dtFit)
.q
EOF

fi

