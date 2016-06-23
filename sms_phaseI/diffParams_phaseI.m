%
% Case: 180 degree refocusing pulse for twice-refocused/dPFG diffusion
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
evalp.maxErr = 2*de; % max deviation of refocusing profile from target
nb = 5; % number of bands
fov = 12; % cm, fov over which bands are evenly spaced
slthick = 1.25; % mm, slice thickness
dz = (slthick/40)/10; % deltaz of profile (cm)
tb = 4; % time-bandwidth product of slice profiles
evalp.maxPhsDev = Inf; % radians, max amount of in-band phase deviation
[evalp.b2d,evalp.roi,evalp.z] = gen_sms_eval_prof(de,nb,fov,slthick,dz,tb);

% multibandeval switches
evalp.fname = 'diff_eval';
evalp.genfig = true;
