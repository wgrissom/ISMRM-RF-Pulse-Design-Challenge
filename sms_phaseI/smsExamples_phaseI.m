% This script will design and evaluate multiband pulses for the
% ISMRM RF Pulse Design Challenge on multiband pulse design for sms
% acquisitions.
% 
% This script provides examples of pulses that meet the specifications
% for both diffusion and TSE cases. It is best if you have John Pauly's
% rf_tools installed and in your path. They can be downloaded from:
% http://rsl.stanford.edu/research/software.html
%
% For more information on the challenge and the specific
% problems, visit the website at challenge.ismrm.org
%
% 2015, Will Grissom and Kawin Setsompop
% Developed in MATLAB R2015a

addpath ../sms_utils

%
% Case 1: 180 degree refocusing pulse for TSE
%

% Get evaluation parameters
tseParams_phaseI;

% design and evaluate a PINS pulse
tse.dt = 2e-6; % dwell time
mindurRF = 0; % switch to use min duration RF for all subpulses
if exist('dzrf')
  halfShift = true; % shift pattern by 1/2 slice gap to line up with target
  [tse.rf,tse.g] = dz_pins(tb,fov/nb,slthick,de,...
     0.9999*evalp.maxb1/100,evalp.maxg,evalp.maxgslew,...
     tse.dt,mindurRF,halfShift);
  tse.rf = tse.rf*100; % convert to uT
else
  load('PINSRFandGrad.mat');
  tse.rf = rfpins;
  tse.g = gpins;
end

% evaluate it
[tseIsValid,tseDur,tseErrorCode] = smsEval(tse.rf,tse.g,tse.dt,evalp);
if tseIsValid == true
    fprintf('TSE pulse passed with duration %d us\n',tseDur);
else
    fprintf('TSE pulse failed with error code %d\n',tseErrorCode);
end    

%
% Case 2: 180 degree refocusing pulse for twice-refocused/dPFG diffusion
%

% Get evaluation parameters
diffParams_phaseI;

% design and evaluate a conventional multiband pulse
n = 1024; % number of time points in pulse

% call JP's SLR code to get a single-band refocusing pulse
if exist('dzrf')
  % JP's SLR RF design tool; output is single-band RF in radians
  rf1b = real(dzrf(64,tb,'se','ls',de,sqrt(de)));
  % interpolate to the desired number of points
  rf1bi = interp1(0:1/63:1,rf1b,0:1/(n-1):1,'spline',0);
  rf1b = rf1bi./sum(rf1bi)*sum(rf1b);
else
  % load the waveform from a .mat file
  load('singleSliceRF.mat');
end

% calculate the multiband modulation function
bandsep = fov/nb/(slthick/10)*tb; % band separation (integer)
B = sum(exp(1i*2*pi/n*(-n/2:n/2-1)'*((0:nb-1)-(nb-1)/2)*bandsep),2);

% modulate the single band pulse to get the multiband pulse
rfmb = rf1b(:).*B;

% convert to uT
rfmbut = rfmb./max(abs(rfmb))*evalp.maxb1*0.9999; % MB RF waveform (uT)

% calculate the dwell time for the peak b1
diff.dt = max(abs(rfmb))/(0.9999*evalp.maxb1*2*pi*42.58); % seconds

% calculate the gradient waveform
gmb = tb/(diff.dt*n)/4258/(slthick/10)*ones(n,1)*10; % mT/m
% put ramps on each end
nramppts = ceil(gmb(1)/(diff.dt*1000/2*evalp.maxgslew));
diff.g = [gmb(1)*(0:nramppts-1)'/nramppts;gmb;...
    gmb(1)*(nramppts-1:-1:0)'/nramppts];
diff.rf = [zeros(nramppts,1);rfmbut;zeros(nramppts,1)];

% evaluate it
[diffIsValid,diffDur,diffErrorCode] = smsEval(diff.rf,diff.g,diff.dt,evalp);
if diffIsValid == true
    fprintf('Diffusion pulse passed with duration %d us\n',diffDur);
else
    fprintf('Diffusion pulse failed with error code %d\n',diffErrorCode);
end   

% total score is sum of two pulse durations (us)
totalScore = tseDur + diffDur;
fprintf('Total score is %d us\n',totalScore);

