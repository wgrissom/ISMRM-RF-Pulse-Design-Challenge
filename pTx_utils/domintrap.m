function [trap,ramppts]=domintrap(area,gmax,dgdt,dt);
% dotrap(area,gmax,dgdt,dt);
%   area = pulse area in (g sec)/cm
%   gmax = max gradient in g/cm
%   dgdt = max slew in g/cm/sec
%   dt   = sample time in sec

if nargin < 5
  rampsamp = 1; % in case we are just making a rewinder
end

if abs(area) > 0

  % we get this solution for plateau amp by setting derivative of
  % duration as a function of amplitude to zero and solving.
  a = sqrt(dgdt*area/2);
  
  % now finish it with discretization
  % make a flat portion of magnitude a
  % and enough area for the entire swath
  flat = ones(1,floor(area/a/dt));
  flat = flat/sum(flat)*area/dt;
  if max(flat) > gmax
    flat = ones(1,ceil(area/gmax/dt));
    flat = flat/sum(flat)*area/dt;
  end
  
  % make attack and decay ramps
  ramppts = ceil(max(flat)/dgdt/dt);
  trap = [ (0:(ramppts-1))/ramppts*max(flat) flat ((ramppts-1):-1:0)/ ...
           ramppts*max(flat) ];
  
else
  trap = 0;
end

