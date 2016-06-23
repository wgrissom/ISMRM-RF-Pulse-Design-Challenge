% This script will design and evaluate sms pulses for PHASE II of the
% ISMRM RF Pulse Design Challenge on multiband pulse design.
% It is intended to provide examples of pulses that meet the specifications
% for both cases. For more information on the challenge and the specific
% problems, check the website at challenge.ismrm.org
% 2016, Will Grissom and Kawin Setsompop
% Developed in MATLAB R2015a

% The main differences between this code and the Phase I code:
% 1) We now evaluate multiple numbers of bands and slice thicknesses for 
%    both TSE and diffusion pulses.
% 2) We now enforce SAR constraints (3.2 W/kg with a representative 3T head 
%    coil) for both TSE and diffusion pulses.
% Please note that we now store all the pulses in cell arrays of 
% structs, where the first dimension of each array is the number of bands
% index, and the second dimension is the slice thickness index. 

% Note that this code REQUIRES John Pauly's rf_tools SLR toolbox to run. 
% This can be downloaded at: http://rsl.stanford.edu/research/software.html

addpath ../sms_utils

maxScore = 999999; % max score that will be posted if pulses are not valid

%%
% Case 1: 180 degree refocusing pulse for TSE
%

% Get evaluation parameters
tseParams_phaseII; % defines evalp, for each number of bands and slice thickness

% design a pulse for each number of bands and slice thickness
tse = cell(length(nb),length(slthick));
for ii = 1:length(nb) % nb and slthick are defined by tseParams
    for jj = 1:length(slthick)

        % design and evaluate a PINS pulse
        tse{ii,jj}.dt = 2e-6; % dwell time
        mindurRF = 0; % switch to use min duration RF for all subpulses
        halfShift = true; % shift pattern by 1/2 slice gap to line up with target
        [tse{ii,jj}.rf,tse{ii,jj}.g] = dz_pins(tb,fov/nb(ii),slthick(jj),de,...
            0.9999*evalp{ii,jj}.maxb1/100,evalp{ii,jj}.maxg,evalp{ii,jj}.maxgslew,...
            tse{ii,jj}.dt,mindurRF,halfShift);
        tse{ii,jj}.rf = tse{ii,jj}.rf*100; % convert to uT

    end
end

% evaluate each TSE pulse
tseIsValid = zeros(length(nb),length(slthick));
tseDur = zeros(length(nb),length(slthick));
tseErrorCode = zeros(length(nb),length(slthick));
for ii = 1:length(nb) % nb and slthick are defined by tseParams
    for jj = 1:length(slthick)

        [tseIsValid(ii,jj),tseDur(ii,jj),tseErrorCode(ii,jj)] = smsEval(tse{ii,jj}.rf,tse{ii,jj}.g,tse{ii,jj}.dt,evalp{ii,jj});
        if tseIsValid(ii,jj) == true
            fprintf('%d-Band, %0.2f mm-thick TSE PINS pulse passed with duration %d us\n',...
                nb(ii),slthick(jj),tseDur(ii,jj));
        else
            fprintf('%d-Band, %0.2f mm-thick TSE PINS pulse FAILED with error code %d\n',...
                nb(ii),slthick(jj),tseErrorCode(ii,jj));
        end

    end
end

%%
% Case 2: 180 degree refocusing pulse for twice-refocused/dPFG diffusion
%

% Get evaluation parameters
diffParams_phaseII; % defines evalp, for each number of bands and slice thickness

% design and evaluate a conventional multiband pulse
n = 1024; % number of time points in pulse

% design and evaluate a pulse for each number of bands and slice thickness
diff = cell(length(nb),length(slthick));
for ii = 1:length(nb) % nb and slthick are defined by diffParams
    for jj = 1:length(slthick)
        
        % calculate the multiband modulation function
        bandsep = fov/nb(ii)/(slthick(jj)/10)*tb; % band separation (integer)
        
        % Direct multiband filter design, followed by inverse SLR
        % Compared to modulating the final RF waveforms, 
        % this approach mitigates between-band interference to some degree
        di = dinf(de/8,sqrt(sqrt(de)));
        w = di/tb; % fractional transition width
        centers = bandsep*((0:nb(ii)-1) - (nb(ii)-1)/2); % integer center of each band
        fNoShift = [-(1+w)*(tb/2) -(1-w)*(tb/2) (1-w)*(tb/2) (1+w)*(tb/2)];
        mNoShift = [0 1 1 0];
        f = -n/2;m = 0; % no passband at the most negative frequency
        % f: normalized band edges; m: filter response at each band edge
        for kk = 1:length(centers)
            f = [f centers(kk)+fNoShift];
            m = [m mNoShift];
        end
        f = [f (n/2)]/(n/2); % append nyquist and normalize to (0,1)
        m = [m 0]; % append final stopband
        wts = kron(ones(1,nb(ii)),[1 de/8/sqrt(sqrt(de))]); % band error weights
        wts = [de/8/sqrt(sqrt(de)) wts];
        f = f(length(f)/2:end);f(1) = 0;
        m = m(length(m)/2:end);if ~rem(nb(ii),2);m(1) = 0;end
        wts = wts(ceil(length(wts)/2):end);
        bmb = firls(n-1,f,m,wts); % get the multiband beta filter
        rfmb = b2rf(bmb); % call SLR to convert beta filter to RF
        
        % convert to uT
        rfmbut = rfmb(:)./max(abs(rfmb))*evalp{ii,jj}.maxb1*0.9999; % MB RF waveform (uT)
        
        % calculate the dwell time for the peak b1
        diff{ii,jj}.dt = max(abs(rfmb))/(0.9999*evalp{ii,jj}.maxb1*2*pi*42.58); % seconds
                
        % calculate the gradient waveform
        gmb = tb/(diff{ii,jj}.dt*n)/4258/(slthick(jj)/10)*ones(n,1)*10; % mT/m
        % put ramps on each end
        nramppts = ceil(gmb(1)/(diff{ii,jj}.dt*1000/2*evalp{ii,jj}.maxgslew));
        diff{ii,jj}.g = [gmb(1)*(0:nramppts-1)'/nramppts;gmb;...
            gmb(1)*(nramppts-1:-1:0)'/nramppts];
        % put zeros on the ends of the pulse to accommodate ramps
        diff{ii,jj}.rf = [zeros(nramppts,1);rfmbut;zeros(nramppts,1)];
        
    end
end

% evaluate each diffusion pulse
diffIsValid = zeros(length(nb),length(slthick));
diffDur = zeros(length(nb),length(slthick));
diffErrorCode = zeros(length(nb),length(slthick));
for ii = 1:length(nb) % nb and slthick are defined by diffParams
    for jj = 1:length(slthick)
       
        [diffIsValid(ii,jj),diffDur(ii,jj),diffErrorCode(ii,jj)] = ...
            smsEval(diff{ii,jj}.rf,diff{ii,jj}.g,diff{ii,jj}.dt,evalp{ii,jj});
        if diffIsValid(ii,jj) == true
            fprintf('%d-Band, %0.2f mm-thick diffusion MB pulse passed with duration %d us\n',...
                nb(ii),slthick(jj),diffDur(ii,jj));
        else
            fprintf('%d-Band, %0.2f mm-thick diffusion MB pulse FAILED with error code %d\n',...
                nb(ii),slthick(jj),diffErrorCode(ii,jj));
        end

    end
end

%%
% Check if all the pulses passed; if so, report the score
%
if ~any([tseIsValid(:);diffIsValid(:)] == false)
    % total score is sum of all pulse durations (us)
    totalScore = sum(tseDur(:)) + sum(diffDur(:));
    fprintf('Total score is %d us\n',totalScore);
else
    fprintf('One or more pulses failed; score is %d\n',maxScore);
end

