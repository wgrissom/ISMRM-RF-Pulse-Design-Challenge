# ISMRM-RF-Pulse-Design-Challenge
Scoring codes and examples for the 2016 ISMRM Challenge on RF Pulse Design. The challenge website is http://challenge.ismrm.org.

The files are:
```
pTx_phaseI/ - Phase I of the pTx sub-challenge
    +-- pTxMapProcessing_phaseI/ - Scripts to process the torso B1 maps for Phase I
    |   +-- interp_pTx_maps.m - Interpolates input B1 maps to the evaluation grid
    |   +-- loadAndProcMaps_phaseI.m - Loads and processes raw B1 maps for the challenge
    |   +-- pTxParams_torso_phaseI.m - Torso problem parameters
    +-- pTxChallengeCodeWalkthrough.pdf - A walkthrough of the structure of the Phase I pTx code
    +-- pTxExample_phaseI.m - Example script to design a 5-spoke pulse and evaluate it
pTx_phaseII/ - Phase II of the pTx sub-challenge
    +-- pTxMapProcessing_phaseII/ - Scripts to process the torso and brain B1 maps for Phase II
    |   +-- interp_pTx_maps_multislice.m - Interpolates input B1 maps to the evaluation grid, for multiple slices
    |   +-- interp_pTx_maps_sms.m - Interpolates input B1 maps to the evaluation grid, for SMS
    |   +-- loadAndProcMaps_brain_sms_axial_phaseII.m - Loads and processes raw B1 maps for the challenge
    |   +-- loadAndProcMaps_brain_sms_coronal_phaseII.m - Loads and processes raw B1 maps for the challenge
    |   +-- loadAndProcMaps_torso_phaseII.m - Loads and processes raw B1 maps for the challenge
    |   +-- pTxParams_brain_sms_axial_phaseII.m - Axial brain problem parameters
    |   +-- pTxParams_brain_sms_coronal_phaseII.m - Coronal brain problem parameters
    |   +-- pTxParams_torso_multislice_phaseII.m - Multislice torso problem parameters
    +-- pTxExamples_phaseII.m - Example script to design 5-spoke pulses for each problem and evaluate them
pTx_utils/ - Utility scripts for the pTx sub-challenge (both phases)
    +-- domintrap.m - Design a mininum duration trapezoid with 
    +-- doTrap.m - Design a minimum duration trapezoid
    +-- dzSpokes_sms.m - Design a spokes pulse for one or more slices
    +-- pTxEval.m - Evaluate and score a pTx pulse submission
    +-- subrf.mat - SLR slice-selective pulse shape, to make John Pauly's toolbox optional
sms_phaseI/ - Phase I of the SMS sub-challenge
    +-- MBChallengeCodeWalkthrough.pdf - A walkthrough of the structure of the Phase I SMS code
    +-- PINSRFandGrad.mat - PINS pulse shape, to make John Pauly's toolbox optional
    +-- diffParams_phaseI.m - Diffusion problem parameters
    +-- singleSliceRF.mat - SLR slice-selective pulse shape for conventional design, to make John Pauly's toolbox optional
    +-- smsExamples_phaseI.m - Design an example PINS pulse for TSE, and a conventional multiband pulse for diffusion, and evaluate them
    +-- tseParams_phaseI.m - TSE problem parameters
sms_phaseII/ - Phase II of the SMS sub-challenge
    +-- diffParams_phaseII.m - 
    +-- smsExamples_phaseII.m - 
    +-- tseParams_phaseII.m
sms_utils/ - Utility scripts for the SMS sub-challenge (both phases)
    +-- blochsim.m - 
    +-- dotrap.m - 
    +-- dz_pins.m - 
    +-- gen_sms_eval_prof.m - 
    +-- smsEval.m - Evaluate and score an SMS pulse submission
```
