% This script will design and evaluate spokes parallel transmit pulses for
% the ISMRM RF Pulse Design Challenge on parallel transmit pulse design.
% It is intended to provide an example of a set of pulses that meet the
% challenge specifications.
%
% For more information on the challenge and the specific
% problems, check the website at challenge.ismrm.org
% 2016, Will Grissom and Kawin Setsompop
%
% Developed in MATLAB R2015a


%% New in Phase II:
%   1) Contestants will now design three sets of pulses for different slice
%      positions in the torso.
%   2) In addition, contestants will design two sets of multiband pulses 
%      for a brain volume at 10.5T with 16 tx channels, in axial and 
%      coronal orientations. 
%   3) The total score will be the sum of all these pulses' durations. 

maxScore = 999999; % max score that will be posted if pulses are not valid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Design and evaluate the multiband brain pulses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

orientations = {'axial','coronal'};

addpath ../pTx_utils ../pTx_maps

% cell array to store the solutions in
brain = cell(length(orientations),1);

for ii = 1:length(orientations)
    
    %% load maps and evaluation parameters
    fprintf('Loading brain evaluation parameters and maps for %s orientation\n',...
        orientations{ii});

    % Load source b1 maps and tissue mask ('maps' structure)
    % maps.b1: [Nx Ny Nz Nc] (uT), complex B1+ maps
    % maps.mask: [Nx Ny Nz] (logical), tissue mask
    % maps.fov: [3] (cm), fov in each dimension
    %
    % This will also load the evalp structure defining the problem
    % constraints
    load(['brain_maps_sms_' orientations{ii} '_phaseII']);
    
    %% design a spokes pulse for this slice
    fprintf('Designing %s multiband spokes pulses for the brain\n',orientations{ii});
    
    % Get scalar pulse and constraint parameters
    probp.Nspokes = 5; % number of spokes in pulse
    probp.flipAngle = evalp.flipAngle; % degrees
    probp.beta = 10^-1; % RF power regularization parameter
    probp.slthick = evalp.slThick; % slice thickness, mm
    probp.tb = evalp.tb; % slice profile time-bandwidth product
    probp.dt = 0.5e-6; % s, dwell time
    probp.maxg = evalp.maxg; % max gradient amplitude, mT/m
    probp.maxgslew = evalp.maxgslew; % max gradient slew, mT/m/ms
    probp.deltaxmin = 4; % finest possible res of kx-ky spokes traj (cm)
    probp.maxb1 = evalp.maxb1; % peak digital RF amplitude constraint
    probp.TR = evalp.TR; % TR for SAR calculation
    probp.slCent = evalp.slCent; % slice center position (cm)
    
    % get the center slice b1+/df0 maps for the design
    slMapInds = round((evalp.slCent + maps.fov(3)/2)/maps.fov(3)*size(maps.b1,3));
    spokesmaps.b1 = maps.b1(:,:,slMapInds,:);
    spokesmaps.mask = maps.mask(:,:,slMapInds);
    spokesmaps.fov = maps.fov(1:2);
    spokesmaps.xy = reshape(evalp.xyz(:,1:2),[size(maps.mask) 2]); % design with same xy grid as evaluation, to ensure agreement
    spokesmaps.xy = squeeze(spokesmaps.xy(:,:,1,:));
    spokesmaps.vop = evalp.vop;
    spokesmaps.maxSAR = evalp.maxSAR;
    
    % get a quadrature initial target phase pattern
    tmp = 0;
    for jj = 1:evalp.Nc
        tmp = tmp + exp(1i*(jj-1)/evalp.Nc*2*pi)*spokesmaps.b1(:,:,:,jj);
    end
    probp.phsinit = angle(tmp);
    
    % run the design
    [brain{ii}.rf,brain{ii}.g,brain{ii}.dt] = dzSpokes_sms(spokesmaps,probp);
    
end

brainIsValid = zeros(length(orientations),1);
brainDur = zeros(length(orientations),1);
brainErrorCode = zeros(length(orientations),1);
for ii = 1:length(orientations)
    
    %% evaluate the spokes pulses
    fprintf('Evaluating the brain pulses for %s orientation\n',orientations{ii});
    
    % Load source b1 maps and tissue mask ('maps' structure)
    % maps.b1: [Nx Ny Nz Nc] (uT), complex B1+ maps
    % maps.mask: [Nx Ny Nz] (logical), tissue mask
    % maps.fov: [3] (cm), fov in each dimension
    %
    % This will also load the evalp structure defining the problem
    % constraints
    load(['brain_maps_sms_' orientations{ii} '_phaseII']);
    
    %evalp.fName = 'spokes_eval';
    evalp.genfig = true;
    [brainIsValid(ii),brainDur(ii),brainErrorCode(ii)] = ...
        pTxEval(brain{ii}.rf,brain{ii}.g,brain{ii}.dt,maps,evalp);
    if brainIsValid(ii) == true
        fprintf('%s brain pulses passed with duration %d us\n',orientations{ii},brainDur(ii));
    else
        fprintf('%s brain pulses failed with error code %d\n',orientations{ii},brainErrorCode(ii));
    end
    
end % loop over orientations


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Design and evaluate the torso pulses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load source b1 maps and tissue mask ('maps' structure)
% maps.b1: [Nx Ny Nz Nc] (uT), complex B1+ maps
% maps.mask: [Nx Ny Nz] (logical), tissue mask
% maps.fov: [3] (cm), fov in each dimension
%
% This will also load the evalp structure defining the problem
% constraints, which is now a cell array in Phase II
disp 'Loading torso maps and evaluation parameters...'
load('torso_maps_phaseII.mat');

% cell array to store the solutions in
torso = cell(length(evalp),1);

for ii = 1:length(evalp)
    
    %% design a spokes pulse for this slice
    fprintf('Designing a spokes pulse for torso slice %d...\n',ii);
    
    % Get scalar pulse and constraint parameters
    probp.Nspokes = 5; % number of spokes in pulse
    probp.flipAngle = evalp{ii}.flipAngle; % degrees
    probp.beta = 10^-1; % RF power regularization parameter
    probp.slthick = evalp{ii}.slThick; % slice thickness, mm
    probp.tb = evalp{ii}.tb; % slice profile time-bandwidth product
    probp.dt = 1e-6; % s, dwell time
    probp.maxg = evalp{ii}.maxg; % max gradient amplitude, mT/m
    probp.maxgslew = evalp{ii}.maxgslew; % max gradient slew, mT/m/ms
    probp.deltaxmin = 10; % finest possible res of kx-ky spokes traj (cm)
    probp.maxb1 = evalp{ii}.maxb1; % peak digital RF amplitude constraint
    probp.TR = evalp{ii}.TR; % TR for SAR calculation
    probp.slCent = evalp{ii}.slCent; % slice center position (cm)
    
    % get the center slice b1+/df0 maps for the design
    slMapInd = round((evalp{ii}.slCent + maps.fov(3)/2)/maps.fov(3)*size(maps.b1,3));
    spokesmaps.b1 = squeeze(maps.b1(:,:,slMapInd,:));
    spokesmaps.mask = squeeze(maps.mask(:,:,slMapInd));
    spokesmaps.fov = maps.fov(1:2);
    spokesmaps.xy = reshape(evalp{ii}.xyz(:,1:2),[size(maps.mask) 2]); % design with same xy grid as evaluation, to ensure agreement
    spokesmaps.xy = squeeze(spokesmaps.xy(:,:,1,:));
    spokesmaps.vop = evalp{ii}.vop;
    spokesmaps.maxSAR = evalp{ii}.maxSAR;
    
    % get a quadrature initial target phase pattern
    tmp = 0;
    for jj = 1:evalp{ii}.Nc
        tmp = tmp + exp(1i*(jj-1)/evalp{ii}.Nc*2*pi)*spokesmaps.b1(:,:,jj);
    end
    probp.phsinit = angle(tmp);
    
    % run the design
    [torso{ii}.rf,torso{ii}.g,torso{ii}.dt] = dzSpokes_sms(spokesmaps,probp);
    
end

torsoIsValid = zeros(length(evalp),1);
torsoDur = zeros(length(evalp),1);
torsoErrorCode = zeros(length(evalp),1);
for ii = 1:length(evalp)
    
    %% evaluate the spokes pulses
    fprintf('Evaluating the pulses for torso slice %d...\n',ii);
    
    %evalp{ii}.fName = ['spokes_eval_slice' num2str(ii)];
    evalp{ii}.genfig = true;
    [torsoIsValid(ii),torsoDur(ii),torsoErrorCode(ii)] = ...
        pTxEval(torso{ii}.rf,torso{ii}.g,torso{ii}.dt,maps,evalp{ii});
    if torsoIsValid(ii) == true
        fprintf('Spokes pulses for torso slice %d passed with duration %d us\n',ii,torsoDur(ii));
    else
        fprintf('Spokes pulses for torso slice %d failed with error code %d\n',ii,torsoErrorCode(ii));
    end
    
end

%% Check if all the pulses passed; if so, report the score
if ~any([brainIsValid(:);torsoIsValid(:)] == false)
    % total score is sum of all pulse durations (us)
    totalScore = sum(brainDur(:)) + sum(torsoDur(:));
    fprintf('Total score is %d us\n',totalScore);
else
    fprintf('One or more pulses failed; score is %d\n',maxScore);
end

