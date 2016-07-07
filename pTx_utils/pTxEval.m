function [isValid,dur,errCode] = pTxEval(rf,g,dt,maps,evalp)

% function [isValid,dur,errCode] = pTxEval(rf,g,dt,evalp)
%
% Function to evaluate the validity and duration of a
% submitted small-tip-angle parallel transmit RF pulse
%
% In:
%   (contestant-supplied inputs)
%   rf      RF waveforms (digital units)
%   g       gradient waveforms (mT/m)
%   dt      dwell time (s)
%   (committee-supplied inputs)
%   maps    Structure containing B1+ maps and VOPs
%           maps Fields:
%               mask        Logical tissue mask for error calculation
%               b1          B1+ maps for the Tx coils
%               fov         Field-of-View of b1 map matrix in each dim
%   evalp   Parameter structure defining the constraints.
%
%
% Out:
%   isValid     Does the pulse meet all constraints (boolean)
%   dur         Pulse duration (ms)
%   errCode   If isValid == false, what is the error?
%               Values:     -1 =    pulse is valid, no error
%                           1 =     RF vector has wrong dimensions
%                           2 =     gradient vector has wrong dimensions
%                           3 =     gradient vector not real-valued
%                           4 =     RF and gradient vectors have ~= lengths
%                           5 =     dwell time is not a scalar value
%                           6 =     dwell time is not real-valued
%                           7 =     dwell time <= 0
%                           8 =     peak RF is too high
%                           9 =     peak gradient is too high
%                           10 =    peak gradient slew is too high
%                           11 =    gradient waveform doesn't start+end at 0
%                           12 =    SAR is too high
%                           13 =    In-slice profile does not meet max flip angle error constraint
%                           14 =    In-slice RMS flip angle error does not meet constraint
%                           15 =    Through-slice phase profile does not meet max deviation constraint
%                           16 =    Out-of-slice profile does not meet max flip angle error constraint
%                           17 =    Out-of-slice RMS flip angle error does not meet constraint
%

dur = Inf; % returned duration if the pulse is not valid
errCode = -1; % returned error code if pulse is valid
isValid = false;

try
    % check inputs
    if any(~isfinite(rf(:))) || any(~isfinite(g(:)))
        isValid = false;
        errCode = 0;
        errMsg = 'Error code 0: RF or gradient vectors contain NaNs or Infs';
        error(errMsg);
    end
    if size(rf,2) ~= evalp.Nc
        isValid = false;
        errCode = 1;
        error 'Error code 1: RF waveform vector must be Nt x Nc'
    end
    if size(g,2) ~= 3
        isValid = false;
        errCode = 2;
        error 'Error code 2: gradient waveform vector must be Nt x 3'
    end
    if ~isreal(g)
        isValid = false;
        errCode = 3;
        error 'Error code 3: gradient vector is not real-valued'
    end
    if size(rf,1) ~= size(g,1)
        isValid = false;
        errCode = 4;
        error 'Error code 4: RF and gradient vectors have ~= # time points'
    end
    if ~isscalar(dt)
        isValid = false;
        errCode = 5;
        error 'Error code 5: dwell time is not a scalar value'
    end
    if ~isreal(dt)
        isValid = false;
        errCode = 6;
        error 'Error code 6: dwell time is not real-valued'
    end
    if dt <= 0
        isValid = false;
        errCode = 7;
        error 'Error code 7: dwell time <= 0'
    end
    
    % check that peak RF meets constraint
    if any(max(abs(rf)) > evalp.maxb1)
        isValid = false;
        errCode = 8;
        error 'Error code 8: peak RF is too high'
    end
    
    % check gradient amplitude
    if any(max(abs(g)) > evalp.maxg)
        isValid = false;
        errCode = 9;
        error 'Error code 9: peak gradient is too high'
    end
    
    % check gradient slew
    gdiff = diff(g)/(dt*1000); % mT/m/ms
    if any(max(abs(gdiff)) > evalp.maxgslew)
        isValid = false;
        errCode = 10;
        error 'Error code 10: peak gradient slew is too high'
    end
    
    % check that gradient starts and ends at zero
    if any(g(1,:) ~= 0) || any(g(end,:) ~= 0)
        isValid = false;
        errCode = 11;
        error 'Error code 11: gradient waveform doesn''t start+end at 0'
    end
    
    % check that SAR is not too high
    if isfield(evalp,'vop')
        nVOP = size(evalp.vop,3); % # VOPs
        sarVOP = zeros(nVOP,size(rf,1));
        for ii = 1:size(rf,1)
            rft = rf(ii,:); % all coil's samples for this time point
            for jj = 1:nVOP
                sarVOP(jj,ii) = real(rft(:)'*(evalp.vop(:,:,jj)*rft(:)));
            end
        end
        sarVOP = sum(sarVOP,2) * dt./evalp.TR;
        sarVOPfrac = sarVOP./evalp.maxSAR(:);
        if any(sarVOPfrac > 1) % if the pulse violates any SAR constraint
            isValid = false;
            errCode = 12;
            error 'Error code 12: max SAR violated'
        end
    end
    
    % calculate the 3D excitation pattern
    k = -flipud(cumsum(flipud(g)))/10*dt*evalp.gamma*100; % exc k-space traj (rad/cm)
    %t = ((0:size(g,1)-1) - size(g,1))*dt; % time vector for off-resonance
    if evalp.useNUFFT && ~isempty(which('Gmri'))
        
        Kd = 2*size(maps.mask);
        nufft_args = {size(maps.mask), [7 7 7], Kd, ...
            size(maps.mask)/2, 'table', 2^10, 'minmax:kb'};
        mask = maps.mask;
        fov = maps.fov;
        G = Gmri(k/2/pi,mask,'fov',fov,'nufft',nufft_args)';
        b1 = maps.b1;
        b1 = permute(b1,[4 1 2 3]);
        b1 = b1(:,:).';
        theta = 0; % total excited flip angle pattern
        for ii = 1:evalp.Nc
            fprintf('Calculating excitation pattern for transmit coil %d of %d\n',ii,evalp.Nc);
            theta = theta + b1(mask,ii).*(G*rf(:,ii));
        end
        tmp = zeros(size(mask));
        tmp(mask) = theta;
        theta = tmp*1i*evalp.gamma*dt*180/pi;
        clear b1 mask
        
    else
        
        disp 'Warning: Using brute-force profile calculation; this could take a while...'
        
        b1 = maps.b1;b1 = permute(b1,[4 1 2 3]);b1 = b1(:,maps.mask(:)).';
        xyzm = evalp.xyz(maps.mask(:),:); % spatial locations, masked to tissue
        Ns = size(xyzm,1);Nt = size(rf,1);
        NsSeg = ceil(Ns/evalp.maxMatrixDim); % number of segments in space
        NtSeg = ceil(Nt/evalp.maxMatrixDim); % number of segments in time
        theta = {}; % excited flip angle pattern
        parfor ii = 1:NsSeg % loop over space segments
            
            % get spatial locations of this segment
            sInds = (ii-1)*evalp.maxMatrixDim+1:min(ii*evalp.maxMatrixDim,Ns);
            xyz = xyzm(sInds,:);
            %df0 = maps.df0(sInds);
            b1s = b1(sInds,:);
            
            thetat = 0;
            for jj = 1:NtSeg % loop over time segments
                
                % get temporal locations of this segment
                tInds = (jj-1)*evalp.maxMatrixDim+1:min(jj*evalp.maxMatrixDim,Nt);
                
                % calculate system matrix
                A = exp(1i*(xyz(:,1)*k(tInds,1)' + xyz(:,2)*k(tInds,2)' + xyz(:,3)*k(tInds,3)'));%.*...
                %exp(1i*2*pi*df0(:)*t(tInds)); % system matrix w/o B1 maps
                
                for kk = 1:evalp.Nc % sum excitation over coils
                    thetat = thetat + b1s(:,kk).*(A*rf(tInds,kk));
                end
                
            end
            
            theta{ii} = thetat;
            
            if rem(ii,100) == 0
                fprintf('Calculating excitation pattern for spatial segment %d of %d\n',ii,NsSeg);
            end
            
        end
        
        clear A;
        
        thetacol = [];
        for ii = 1:NsSeg
            thetacol = [thetacol;theta{ii}];
        end
        thetaall = double(maps.mask);
        thetaall(maps.mask) = thetacol;
        
        theta = thetaall*1i*evalp.gamma*dt*180/pi;
        
    end
    
    % show the profile in 3 axes
    if isfield(evalp,'genfig') 
        
        if evalp.genfig == true
            
            xyzpl = reshape(evalp.xyz,[size(theta) 3]);
            xpl = xyzpl(:,1,1,1);ypl = squeeze(xyzpl(1,:,1,2));zpl = squeeze(xyzpl(1,1,:,3));
            
            figure;
            %subplot(131)
            imagesc(ypl,zpl,abs(squeeze(theta(ceil(size(theta,1)/2),:,:))).');
            %axis image;
            colormap jet;title 'Coronal, x = 0'
            if isfield(evalp,'fName')
                savefig(gcf,[evalp.fName '_coronal']);
            end
            
            figure;
            %subplot(132)
            imagesc(xpl,zpl,abs(squeeze(theta(:,ceil(size(theta,2)/2),:))).');
            %axis image;
            colormap jet;
            title 'Sagittal, y = 0'
            if isfield(evalp,'fName')
                savefig(gcf,[evalp.fName '_sagittal']);
            end
            
            % show slice(s) in-plane
            figure;
            slMapInds = round((evalp.slCent + maps.fov(3)/2)/maps.fov(3)*size(maps.b1,3));
            tmp = theta(:,:,slMapInds);
            minTheta = min(abs(tmp(:)));maxTheta = max(abs(tmp(:)));
            for ii = 1:length(evalp.slCent)
                subplot(1,length(evalp.slCent),ii)
                imagesc(xpl,ypl,abs(squeeze(theta(:,:,slMapInds(ii))).'),[minTheta maxTheta]);
                axis image;
                colormap jet;title(sprintf('Axial, z = %d cm',evalp.slCent(ii)));
            end
            if isfield(evalp,'fName')
                savefig(gcf,[evalp.fName '_axial']);
            end
            
        end
        
    end
    
    % calculate the 3D excitation amplitude error
    err = abs(abs(theta) - evalp.thetad);
    
    % check the flip angle profile against the required profile
    % We check both max and RMS error in-slice to help ensure we don't get 'holes'
    if max(err(evalp.inSliceRoi)) > evalp.maxInSliceErr
        isValid = false;
        errCode = 13;
        error 'Error code 13: In-slice max flip angle error does not meet constraint'
    end
    
    if norm(err(evalp.inSliceRoi))/sqrt(sum(evalp.inSliceRoi(:))) > evalp.maxInSliceRMSE
        isValid = false;
        errCode = 14;
        error 'Error code 14: In-slice RMS flip angle error does not meet constraint'
    end
    
    for ii = 1:length(evalp.slCent) % check through-slice phase of each slice independently
        phs = unwrap(angle(theta(:,:,evalp.inSliceInds{ii})),[],3);
        phsd = phs - repmat(mean(phs,3),[1 1 size(phs,3)]);
        if any(abs(phsd(:)) > evalp.maxPhsDev)
            isValid = false;
            errCode = 15;
            errMsg = 'Error code 15: Through-slice phase profile does not meet constraints';
            error(errMsg);
        end
    end
    
    if max(err(evalp.outOfSliceRoi)) > evalp.maxOutOfSliceErr
        isValid = false;
        errCode = 16;
        error 'Error code 16: Out-of-slice profile does not meet amplitude constraints'
    end
    
    if norm(err(evalp.outOfSliceRoi))/sqrt(sum(evalp.outOfSliceRoi(:))) > evalp.maxOutOfSliceRMSE
        isValid = false;
        errCode = 17;
        error 'Error code 17: Out-of-slice RMS flip angle error does not meet constraint'
    end
    
    % calculate the pulse duration
    dur = ceil(dt*size(rf,1)*1000000); % milliseconds
    
    % if they made it this far, the pulse is valid
    isValid = true;
    
    % save everything in output file, exit
    if isfield(evalp,'fName')
        save(evalp.fName);
    end
    
    return
    
catch
    
    % save everything in output file, exit
    if isfield(evalp,'fName')
        fprintf('Error; Saving workspace to %s.mat and exiting...\n',evalp.fName);
        save(evalp.fName);
    else
        fprintf('Error; Exiting...\n');
    end
    
    return
    
end
