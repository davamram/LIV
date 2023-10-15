#!/bin/bash
yodastack -o Results/MadGraph/Reweight/Etr_$1.yoda Results/MadGraph/Reweight/120GeV_Etr$1.yoda Results/MadGraph/Reweight/2000GeV_Etr$1.yoda Results/MadGraph/Reweight/1000GeVA_Etr$1.yoda Results/MadGraph/Reweight/1000GeVB_Etr$1.yoda
rivet-mkhtml --errs --no-weight Results/MadGraph/Reweight/Etr_$1.yoda -o Plots/MadGraph/Reweight/$1GeV
