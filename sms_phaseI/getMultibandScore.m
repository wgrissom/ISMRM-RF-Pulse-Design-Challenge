function totalScore = getMultibandScore(fname,maxScore)

% Input: tse and diff structs with fields:
%   rf: rf waveform (length-Nt vector)
%   g: gradient waveform (length-Nt vector)
%   dt: dwell time (scalar, seconds)

try
    
    printf('Attempting to load SMS pulse file %s\n',fname);
    pulses = load(fname);
    
    % get the tse evaluation parameters
    tseParams;
    
    tse = pulses.tse; % pull out submitted tse pulses
    
    % check that the fields we need are in the structure
    if ~isfield(tse,'rf') || ~isfield(tse,'g') || ~isfield(tse,'dt')
        error 'tse structure is missing one or more required fields'
    end
    
    % evaluate the tse pulse
    evalp = rmfield(evalp,{'genfig','fname'}); % suppress outputs
    [tseIsValid,tseDur,~] = multibandEval(tse.rf,tse.g,tse.dt,evalp);
    
    % get the diffusion evaluation parameters
    clear evalp;
    diffParams;
    
    diff = pulses.diff; % pull out submitted diff pulses
    
    % check that the fields we need are in the structure
    if ~isfield(diff,'rf') || ~isfield(diff,'g') || ~isfield(diff,'dt')
        error 'diff structure is missing one or more required fields'
    end
    
    % evaluate the diffusion pulse
    evalp = rmfield(evalp,{'genfig','fname'}); % suppress outputs
    [diffIsValid,diffDur,~] = multibandEval(diff.rf,diff.g,diff.dt,evalp);
    
    % calculate total score as sum of two pulse durations (us)
    if tseIsValid && diffIsValid
        totalScore = min(tseDur + diffDur,maxScore);
    else
        totalScore = maxScore;
    end
    
catch
    
    totalScore = maxScore;
    
end
