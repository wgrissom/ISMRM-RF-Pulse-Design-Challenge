function evalp = paralleltxparams_torso_multislice

% Defines parameters for the parallel transmit pulse design challenge
% based on the 7T 8-channel torso array. 
% B1+ and SAR maps were provided by Arcan Erturk, University of Minnesota.

% scanner constants
evalp.maxg = 80; % mT/m
evalp.maxgslew = 200; % mT/m/ms
evalp.maxb1 = 30; % digital units, based on 1 kW per channel and B1 maps measured at 1 W
evalp.gamma = 2*pi*42.58; % (rad/s)/uT, proton gyro ratio

% source B1/df0 maps
evalp.rawmapsfname = 'torso_8ch_b1sar.mat'; % Name of the .mat file containing the original source b1+SAR maps.
evalp.evalmapsfname = 'torso_maps_phaseII.mat'; % name of the .mat file containing the processed maps
% File is assumed to contain a 'maps' structure with fields:
%   maps.b1: [Nx Ny Nz Nc] (uT), complex B1+ maps
%   maps.mask: [Nx Ny Nz] (logical), tissue mask
%   maps.fov: [3] (cm), fov in each dimension
load(evalp.rawmapsfname);
evalp.Nc = size(maps.b1,4); % # Tx channels

% source VOPs and SAR specs
evalp.sourcevopsfname = 'torso_VOPs_epspt5.mat'; % Name of .mat file containing source VOPs
load(evalp.sourcevopsfname);
evalp.vop = Sv; % the local SAR matrices
% get the global SAR matrix
evalp.vop(:,:,nVOP+1) = Sv(:,:,1)*nSv(1); % weight by number of spatial locs
for ii = 2:nVOP
    evalp.vop(:,:,nVOP+1) = evalp.vop(:,:,nVOP+1) + Sv(:,:,ii)*nSv(ii);
end
evalp.vop(:,:,nVOP+1) = evalp.vop(:,:,nVOP+1)/sum(nSv); % to get mean global SAR matrix
evalp.maxSAR = 20*ones(size(Sv,3),1); % W/kg, max local SAR
evalp.maxSAR(end+1) = 4; % W/kg, max global SAR
evalp.TR = 0.130; % seconds, sequence TR

% target pattern specs
evalp.flipAngle = 70; % target flip angle (degrees)
evalp.slCent = [-4 0 4]; % slice position in z (cm)
evalp.d1 = 0.01; % nominal fractional flip angle error in passband
evalp.d2 = 0.01; % nominal fractional flip angle error in stopband
evalp.maxInSliceErr = 20; % max in-slice flip angle error (degrees)
evalp.maxInSliceRMSE = 5; % max in-slice flip angle RMS error (degrees)
evalp.maxOutOfSliceErr = 1; % max out-of-slice flip angle error (degrees, calculated as 1% of target flip angle, with 2x fudge factor)
evalp.maxOutOfSliceRMSE = 0.1; % max out-of-slice flip angle RMS error (degrees)
evalp.maxPhsDev = 0.15; % radians, max amount of through-slice phase deviation from mean phase
evalp.tb = 4; % time-bandwidth product of slice profile
evalp.slThick = 2; % target slice thickness (mm)
evalp.dxyz = [maps.fov(1)/size(maps.b1,1)*2 maps.fov(2)/size(maps.b1,2)*2 ...
  (evalp.slThick/20)/10]; % resolution of target pattern for evaluation (cm)
  % Note: for NUFFT to maintain accuracy, in-plane grid size should be > 64 or so

% evaluation switches
evalp.maxMatrixDim = 2^10; % max eval system matrix dimension in any direction
evalp.useNUFFT = true; % use NUFFT (true), or build explicit system matrix (false)
