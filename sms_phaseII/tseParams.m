%
% Case 1: 180 degree refocusing pulse for TSE
%

clear evalp

nb = 8:2:14; % number of bands (RANGE)
slthick = 0.5:0.5:2; % mm, slice thickness (RANGE)

% set up evaluation parameters
evalpCommon.maxg = 80; % mT/m
evalpCommon.maxgslew = 200; % mT/m/ms
evalpCommon.maxb1 = 18; % uT
evalpCommon.gamma = 2*pi*42.58; % (rad/s)/uT, proton gyro ratio
evalpCommon.maxSamples = 20000; % max # samples in the waveforms
evalpCommon.maxSAR = 3.2; % W/kg, SAR constraint
evalpCommon.coilSARefficiency = 0.25; % W/kg/uT^2; representative for a 3T birdcage coil
evalpCommon.pulseFreq = 12/0.220; % pulses per second, for calculating SAR; 
                                  % 12 pulses per slice group/ETL = 12, one slice group every 220 ms
                                  % we are neglecting the 90 pulse's SAR here
                            
% set up evaluation profile
de = 0.01; % target effective error in ref profile
evalpCommon.maxErr = 2*de; % max amplitude deviation of ref profile from target
fov = 24; % cm, fov over which bands are evenly spaced
tb = 3; % time-bandwidth product of slice profiles
evalpCommon.maxPhsDev = 0.01; % radians, max amount of in-band phase deviation

% multibandeval function switches
evalpCommon.fname = 'tse_eval';
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
