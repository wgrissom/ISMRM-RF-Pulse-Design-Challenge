%
% Case 1: 180 degree refocusing pulse for TSE
%

clear evalp

% set up evaluation parameters
evalp.maxg = 80; % mT/m
evalp.maxgslew = 200; % mT/m/ms
evalp.maxb1 = 18; % uT
evalp.gamma = 2*pi*42.58; % (rad/s)/uT, proton gyro ratio
evalp.maxSamples = 20000; % max # samples in the waveforms

% set up evaluation profile
de = 0.01; % target effective error in ref profile
evalp.maxErr = 2*de; % max amplitude deviation of ref profile from target
nb = 12; % number of bands
fov = 24; % cm, fov over which bands are evenly spaced
slthick = 1; % mm, slice thickness
dz = (slthick/40)/10; % deltaz of profile (cm)
tb = 3; % time-bandwidth product of slice profiles
evalp.maxPhsDev = 0.01; % radians, max amount of in-band phase deviation
[evalp.b2d,evalp.roi,evalp.z] = gen_mb_eval_prof(de,nb,fov,slthick,dz,tb);

% multibandeval function switches
evalp.fname = 'tse_eval';
evalp.genfig = true;