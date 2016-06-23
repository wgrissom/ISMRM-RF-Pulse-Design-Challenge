%
% Case 2: 180 degree refocusing pulse for twice-refocused/dPFG diffusion
%

clear evalp

nb = 3:5; % number of bands (RANGE)
slthick = 1:0.25:2; % mm, slice thickness (RANGE)

% set up common evaluation parameters
evalpCommon.maxg = 80; % mT/m
evalpCommon.maxgslew = 200; % mT/m/ms
evalpCommon.maxb1 = 18; % uT
evalpCommon.gamma = 2*pi*42.58; % (rad/s)/uT, proton gyro ratio
evalpCommon.maxSamples = 20000; % max # samples in the waveforms
evalpCommon.maxSAR = 3.2; % W/kg, SAR constraint
evalpCommon.coilSARefficiency = 0.25; % W/kg/uT^2; representative for a 3T birdcage coil
evalpCommon.pulseFreq = 2/0.120; % pulses per second, for calculating SAR;
                           % two pulses per slice group, one slice group every 120 ms
                           % we are neglecting the 90 pulse's SAR here

% set up evaluation profile
de = 0.01; % target effective error in ref profile
evalpCommon.maxErr = 2*de; % max deviation of refocusing profile from target
fov = 12; % cm, fov over which bands are evenly spaced
tb = 4; % time-bandwidth product of slice profiles
evalpCommon.maxPhsDev = Inf; % radians, max amount of in-band phase deviation

% multibandeval switches
evalpCommon.fname = 'diff_eval';
%evalpCommon.genfig = true;

evalp = {};
% get the evaluation parameters for each # bands and each slice thickness
for ii = 1:length(nb)
    for jj = 1:length(slthick)
                
        evalp{ii,jj} = evalpCommon; 
        
        % set up evaluation profile
        dz = (slthick(jj)/40)/10; % deltaz of profile (cm)
        [evalp{ii,jj}.b2d,evalp{ii,jj}.roi,evalp{ii,jj}.z] = ...
            gen_mb_eval_prof(de,nb(ii),fov,slthick(jj),dz,tb);
        
    end
end
