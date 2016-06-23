% This script will design and evaluate spokes parallel transmit pulses for 
% the ISMRM RF Pulse Design Challenge on parallel transmit pulse design.
% It is intended to provide an example of a set of pulses that meet the 
% challenge specifications. To run efficiently, the scoring code
% will require Jeff Fessler's Gmri object, which is part of his Image
% Reconstruction Toolbox that can be downloaded at:
% http://web.eecs.umich.edu/~fessler/code/index.html
% 
% For more information on the challenge and the specific
% problems, check the website at challenge.ismrm.org
% 2015, Will Grissom and Kawin Setsompop
%
% Developed in MATLAB R2015a

%% load maps and evaluation parameters
disp 'Loading evaluation parameters and maps...'

addpath ../pTx_utils ../pTx_maps

% Load source b1 maps and tissue mask ('maps' structure)
% maps.b1: [Nx Ny Nz Nc] (uT), complex B1+ maps
% maps.mask: [Nx Ny Nz] (logical), tissue mask
% maps.fov: [3] (cm), fov in each dimension
% 
% This will also load the evalp structure defining the problem 
% constraints
load('torso_maps_phaseI.mat');


%% design a spokes pulse
disp 'Designing a spokes pulse...'

% Get scalar pulse and constraint parameters
probp.Nspokes = 5; % number of spokes in pulse
probp.flipAngle = evalp.flipAngle; % degrees
probp.beta = 10^-1; % RF power regularization parameter
probp.slthick = evalp.slThick; % slice thickness, mm
probp.tb = evalp.tb; % slice profile time-bandwidth product
probp.dt = 4e-6; % s, dwell time
probp.maxg = evalp.maxg; % max gradient amplitude, mT/m
probp.maxgslew = evalp.maxgslew; % max gradient slew, mT/m/ms
probp.deltaxmin = 10; % finest possible res of kx-ky spokes traj (cm)
probp.maxb1 = evalp.maxb1; % peak digital RF amplitude constraint
probp.TR = evalp.TR; % TR for SAR calculation
probp.slCent = 0;

% get the center slice b1+/df0 maps for the design
slMapInd = round((evalp.slCent + maps.fov(3)/2)/maps.fov(3)*size(maps.b1,3));
spokesmaps.b1 = squeeze(maps.b1(:,:,slMapInd,:));
spokesmaps.mask = squeeze(maps.mask(:,:,slMapInd));
spokesmaps.fov = maps.fov(1:2);
spokesmaps.xy = reshape(evalp.xyz(:,1:2),[size(maps.mask) 2]); % design with same xy grid as evaluation, to ensure agreement
spokesmaps.xy = squeeze(spokesmaps.xy(:,:,1,:));
spokesmaps.vop = evalp.vop;
spokesmaps.maxSAR = evalp.maxSAR;

% get a quadrature initial target phase pattern
tmp = 0;
for ii = 1:evalp.Nc
    tmp = tmp + exp(1i*(ii-1)/evalp.Nc*2*pi)*spokesmaps.b1(:,:,ii);
end
probp.phsinit = angle(tmp);

% run the design - we use the sms-enabled spokes function but with one slice
[rf,g,dt] = dzSpokes_sms(spokesmaps,probp);

%% evaluate the spokes pulse
disp 'Evaluating the spokes pulse...'

evalp.fName = 'spokes_eval';
evalp.genfig = true;
[isValid,dur,errorCode] = pTxEval(rf,g,dt,maps,evalp);
if isValid == true
    fprintf('Spokes pulses passed with duration %d us\n',dur);
else
    fprintf('Spokes pulses failed with error code %d\n',errorCode);
end   
