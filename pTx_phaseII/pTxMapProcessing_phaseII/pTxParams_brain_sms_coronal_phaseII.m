function evalp = pTxParams_brain_sms_coronal_phaseII(rawMapsDir)

% Defines parameters for the parallel transmit pulse design challenge
% based on the 10.5T 16-channel brain array. 
% B1+ and SAR maps were provided by Nicolaus Boulant at Neurospin and
% collaborators at University of Minnesota.
% maps are 4 mm iso resolution

% scanner constants
evalp.maxg = 80; % mT/m
evalp.maxgslew = 200; % mT/m/ms
evalp.maxb1 = 30; % digital units, based on 1 kW per channel and B1 maps measured at 1 W
evalp.gamma = 2*pi*42.58; % (rad/s)/uT, proton gyro ratio

% source B1/df0 maps
evalp.rawmapsfname = 'brain_maps_orig.mat'; % Name of the .mat file containing the original source b1+SAR maps.
evalp.evalmapsfname = 'brain_maps_multiband_coronal.mat'; % name of the .mat file containing the processed maps
% File is assumed to contain a 'maps' structure with fields:
%   maps.b1: [Nx Ny Nz Nc] (uT), complex B1+ maps
%   maps.mask: [Nx Ny Nz] (logical), tissue mask
maps = load([rawMapsDir evalp.rawmapsfname]);
evalp.Nc = size(maps.B1Maps,4); % # Tx channels

% source VOPs and SAR specs
evalp.sourcevopsfname = 'brain_maps_orig.mat'; % Name of .mat file containing source VOPs
load([rawMapsDir evalp.sourcevopsfname]);
evalp.vop = permute(QVOPs_1,[2 3 1]);
evalp.vop(:,:,end+1) = QGlobal; 
evalp.maxSAR = 10*ones(size(evalp.vop,3)-1,1); % W/kg, max local SAR
evalp.maxSAR(end+1) = 3.2; % W/kg, max whole head SAR
evalp.TR = 0.06; % seconds, approximately equal to sequence TR divided by number of slices (from Moeller MRM 2009)

% global target pattern specs
evalp.flipAngle = 60; % target flip angle (degrees)
evalp.d1 = 0.01; % nominal fractional flip angle error in passband
evalp.d2 = 0.01; % nominal fractional flip angle error in stopband
evalp.tb = 4; % time-bandwidth product of slice profile
evalp.slThick = 1.5; % target slice thickness (mm)
evalp.dxyz = [0.4 0.4 (evalp.slThick/20)/10]; % resolution of target pattern for evaluation (cm)
  % Note: for NUFFT to maintain accuracy, in-plane grid size should be > 64 or so
% evaluation switches
evalp.maxMatrixDim = 2^10; % max eval system matrix dimension in any direction
evalp.useNUFFT = true; % use NUFFT (true), or build explicit system matrix (false)
evalp.maxPhsDev = 0.15; % radians, max amount of through-slice phase deviation from mean phase
evalp.maxOutOfSliceErr = 1; % max out-of-slice flip angle error (degrees, calculated as 1% of target flip angle)
evalp.maxOutOfSliceRMSE = 0.1; % max out-of-slice flip angle RMS error (degrees)
evalp.maxInSliceErr = 20; % max in-slice flip angle error (degrees)
evalp.maxInSliceRMSE = 4; % max in-slice flip angle RMS error (degrees)

% slice position-specific target pattern specs
evalp.slCent = 2.8*((0:5)-2.5);%3.5*((0:4)-2);%[-6 -2  2 6]; % slice positions in z (cm)

