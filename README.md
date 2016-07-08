# ISMRM-RF-Pulse-Design-Challenge
Scoring codes and examples for the 2016 ISMRM Challenge on RF Pulse Design. The challenge website is http://challenge.ismrm.org.

The files are:
```
pTx_phaseI/ - Phase I of the pTx sub-challenge
    |-- pTxMapProcessing_phaseI/ - Scripts to process the torso B1 maps for Phase I
            |-- interp_pTx_maps.m - Interpolates input B1 maps to the evaluation grid
            |-- loadAndProcMaps_phaseI.m - Loads input B1 maps and calls interp_pTx_maps.m
    |-- pTxChallengeCodeWalkthrough.pdf - A walkthrough of the structure of the Phase I pTx code
    |-- pTxExample_phaseI.m - Example script to design a 5-spoke pulse and evaluate it
```
