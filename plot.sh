#!/bin/bash
yodastack -o Results/sherpa/Reweight/Etr_$1.yoda Results/sherpa/Reweight/120GeV_Etr$1.yoda Results/sherpa/Reweight/2000GeV_Etr$1.yoda Results/sherpa/Reweight/1000GeVA_Etr$1.yoda Results/sherpa/Reweight/1000GeVB_Etr$1.yoda
rivet-mkhtml --errs --no-weight Results/sherpa/Reweight/Etr_$1.yoda -o Plots/sherpa/Reweight/$1GeV
