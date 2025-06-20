#!/bin/bash
sed -i "/d.ThresholdEnergy=.*/c\d.ThresholdEnergy=$1;" lib/SelectPhoton.cc
rivet-build TEST_ANALYSIS.cc lib/SelectPhoton.cc
rivet --pwd -a TEST_ANALYSIS Events/Events_120GeV/events_0__* -o Results/sherpa/Reweight/120GeV_Etr$1.yoda && echo "120GeV ended" &
rivet --pwd -a TEST_ANALYSIS Events/Events_1000GeV/events_0__1.* -o Results/sherpa/Reweight/1000GeVA_Etr$1.yoda && echo "1000GeV ended" &
rivet --pwd -a TEST_ANALYSIS Events/Events_1000GeV/events_0__10* -o Results/sherpa/Reweight/1000GeVB_Etr$1.yoda && echo "1000GeV ended" &
rivet --pwd -a TEST_ANALYSIS Events/Events_2000GeV/events_0.hepmc.gz -o Results/sherpa/Reweight/2000GeV_Etr$1.yoda && echo "2000GeV ended" &
wait
