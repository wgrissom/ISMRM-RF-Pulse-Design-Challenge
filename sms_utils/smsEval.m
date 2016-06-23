function [isValid,dur,errCode] = smsEval(rf,g,dt,evalp)

% function [isValid,dur,errCode] = multibandeval(rf,g,dt,evalp)
%
% Function to evaluate the validity and duration of a
% submitted multiband RF pulse
%
% In:
%   (contestant-supplied inputs)
%   rf      RF waveform (uT)
%   g       gradient waveform (mT/m)
%   dt      dwell time (seconds)
%   (committee-supplied inputs)
%   evalp   Parameter structure defining the constraints.
%
%
% Out:
%   isValid     Does the pulse meet all constraints (boolean)
%   dur         Pulse duration (integer, microseconds)
%   errCode   If isValid == false, what is the error?
%               Values:     -1 =    pulse is valid, no error
%                           1 =     RF vector has > 1 dimension
%                           2 =     gradient vector has > 1 dimension
%                           3 =     too many samples in waveforms (> maxSamples)
%                           4 =     gradient vector not real-valued
%                           5 =     RF and gradient vectors have ~= lengths
%                           6 =     dwell time is not a scalar value
%                           7 =     dwell time is not real-valued
%                           8 =     dwell time <= 0
%                           9 =     peak RF is too high
%                           10 =     peak gradient is too high
%                           11 =    peak gradient slew is too high
%                           12 =    gradient waveform doesn't start+end at 0
%                           13 =    refocusing profile does not meet
%                                   amplitude constraints
%                           14 =    refocusing profile does not meet
%                                   phase constraints
%
%

dur = Inf; % returned duration if the pulse is not valid
errCode = -1; % returned error code if pulse is valid
isValid = false;
absb2 = [];
gdiff = [];

try

    % check inputs
    if any(~isfinite(rf(:))) || any(~isfinite(g(:)))
        isValid = false;
        errCode = 0;
        errMsg = 'Error code 0: RF or gradient vectors contain NaNs or Infs';
        error(errMsg);
    end
    if length(rf) ~= length(rf(:))
        isValid = false;
        errCode = 1;
        errMsg = 'Error code 1: RF vector has > 1 dimension';
        error(errMsg);
    end
    if length(g) ~= length(g(:))
        isValid = false;
        errCode = 2;
        errMsg = 'Error code 2: gradient vector has > 1 dimension';
        error(errMsg);
    end
    if max(length(rf),length(g)) > evalp.maxSamples
        isValid = false;
        errCode = 3;
        errMsg = 'Error code 3: too many samples in waveforms (> maxSamples)';
        error(errMsg);
    end
    if ~isreal(g)
        isValid = false;
        errCode = 4;
        errMsg = 'Error code 4: gradient vector is not real-valued';
        error(errMsg);
    end
    if length(rf(:)) ~= length(g(:))
        isValid = false;
        errCode = 5;
        errMsg = 'Error code 5: RF and gradient vectors have ~= lengths';
        error(errMsg);
    end
    if ~isscalar(dt)
        isValid = false;
        errCode = 6;
        errMsg = 'Error code 6: dwell time is not a scalar value';
        error(errMsg);
    end
    if ~isreal(dt)
        isValid = false;
        errCode = 7;
        errMsg = 'Error code 7: dwell time is not real-valued';
        error(errMsg);
    end
    if dt <= 0
        isValid = false;
        errCode = 8;
        errMsg = 'Error code 8: dwell time <= 0';
        error(errMsg);
    end

    % check that peak RF meets constraint
    if max(abs(rf)) > evalp.maxb1
        isValid = false;
        errCode = 9;
        errMsg = 'Error code 9: peak RF is too high';
        error(errMsg);
    end

    % check gradient amplitude
    if max(abs(g)) > evalp.maxg
        isValid = false;
        errCode = 10;
        errMsg = 'Error code 10: peak gradient is too high';
        error(errMsg);
    end

    % check gradient slew
    gdiff = diff(g)/(dt*1000); % mT/m/ms
    if max(abs(gdiff)) > evalp.maxgslew
        isValid = false;
        errCode = 11;
        errMsg = 'Error code 11: peak gradient slew is too high';
        error(errMsg);
    end

    % check that gradient starts and ends at zero
    if g(1) ~= 0 || g(end) ~= 0
        isValid = false;
        errCode = 12;
        errMsg = 'Error code 12: gradient waveform doesn''t start+end at 0';
        error(errMsg);
    end

    % simulate the pulse
    [~,b] = blochsim(rf(:)*dt*evalp.gamma,g(:)/10*dt*evalp.gamma*100,evalp.z);

    % calculate refocusing profile
    absb2 = abs(b.^2);

    % calculate the profile error
    err = abs(absb2 - evalp.b2d);

    % make a profile plot
    if isfield(evalp,'genfig') 
        
        if evalp.genfig == true
            figure
            plot(evalp.z,evalp.b2d);
            hold on
            plot(evalp.z,absb2);
            plot(evalp.z,evalp.roi);
            plot(evalp.z(evalp.roi),err(evalp.roi),'.');
            plot(evalp.z(evalp.roi),(0*evalp.roi(evalp.roi) + evalp.maxErr),'.');
            legend('Desired refocusing profile','Candidate refocusing profile',...
                'Error mask','Error','Max allowed error');
            xlabel 'z (cm)'
            ylabel 'Refocusing profile amplitude (\beta^2 = M_{xy}^+/M_{xy}^{-,*})'
            axis([evalp.z(1) evalp.z(end) -0.1 1.1]);
            title(sprintf('Submission: %s.',evalp.fname),'interpreter','none');
            
            if isfield(evalp,'fname')
                savefig(gcf,evalp.fname);
            end
        end
        
    end

    % check the refocusing profile against the required profile
    if max(err(evalp.roi)) > evalp.maxErr
        isValid = false;
        errCode = 13;
        errMsg = 'Error code 13: refocusing amplitude profile does not meet constraints';
        error(errMsg);
    end

    % if specified, evaluate in-band phase deviation
    if isfield(evalp,'maxPhsDev')
        % find the passband edges in the desired profile
        inds = find(abs(diff(evalp.b2d)) > 0.5);
        % fix right edge indices
        inds(1:2:end) = inds(1:2:end)+1;
        for ii = 1:length(inds)/2
            phs = unwrap(angle(b(inds(2*ii-1):inds(2*ii)).^2));
            if any(abs(phs - mean(phs)) > evalp.maxPhsDev)
                isValid = false;
                errCode = 14;
                errMsg = 'Error code 14: refocusing phase profile does not meet constraints';
                error(errMsg);
            end
        end
    end

    % calculate the pulse duration
    dur = ceil(dt*length(rf(:))*1000000); % microseconds

    % if they made it this far, the pulse is valid
    isValid = true;

    % save everything in output file, exit
    if isfield(evalp,'fname')
        save(evalp.fname);
    end

    return

catch

    % save everything in output file, exit
    if isfield(evalp,'fname')
        fprintf('Error; Saving workspace to %s.mat and exiting...\n',evalp.fname);
        save(evalp.fname);
    else
        fprintf('Error; Exiting...\n');
    end

    return

end
