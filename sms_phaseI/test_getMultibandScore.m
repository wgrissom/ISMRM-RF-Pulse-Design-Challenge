clear

% Load the candidate waveform structures - 'tse' and 'diff'
load candidate.mat

% check that .mat file contains required structures
if ~exist('tse','var')
    error 'tse pulse structure not found'
end
if ~exist('diff','var')
    error 'diff pulse structure not found'
end
score = getMultibandScore(tse,diff);

printf('Total score of submission = %i us',score);