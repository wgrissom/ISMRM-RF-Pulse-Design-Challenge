function [b2d,roi,z] = gen_mb_eval_prof(de,nb,fov,dthick,dz,tb)

% function [b2d,roi,z] = gen_mb_eval_prof(de,nb,fov,dthick,dz,tb)
% 
% Function to calculate a target refocusing profile for evaluation of 
% a submitted multiband RF pulse
% 
% In: 
%   de          target effective error in pass+stopbands of refocusing profile              
%   nb          number of bands
%   fov         cm, fov over which bands are evenly spaced
%   dthick      mm, slice thickness
%   dz          cm, spacing between profile samples
%   tb          time-bandwidth product of slice profiles
%
% Out:
%   b2d         target refocusing (beta^2 = Mxy+/conj(Mxy-)) profile,
%               assuming crushing
%   roi         region-of-interest for error calculation (excludes 
%               transition bands)
%   z           spatial coordinates of profile samples


d1 = de/4; % ripple in ref passband profile, converted to beta
d2 = sqrt(de); % ripple in ref stopband profile, converted to beta
z = -fov/2:dz:fov/2; % z-locations on which we define profile (cm)
w = dinf(d1,d2)/tb; % fractional transition width
d = zeros(size(z));s = d;
if rem(nb,2) % if Nb odd
    % start out the f, m and w vectors with the DC band
    f = [0 (1-w)*(tb/2) (1+w)*(tb/2)]*dthick/10/tb; % band edges (cm)
    d = abs(z) <= f(2); % target pattern
else
    f = 0;
end   

% add non-DC bands to the profiles
for ii = 1:floor(nb/2)
    cent = (ii - (rem(nb,2) == 0)/2)*fov/nb;
    f = [f (cent-(1+w)*(tb/2)*dthick/10/tb) (cent-(1-w)*(tb/2)*dthick/10/tb) ...
        (cent+(1-w)*(tb/2)*dthick/10/tb) (cent+(1+w)*(tb/2)*dthick/10/tb)];
    d = d | (z >= f(end-2) & z <= f(end-1));
    d = d | (-z >= f(end-2) & -z <= f(end-1));
    s = s | (z >= f(end-4) & z <= f(end-3));
    s = s | (-z >= f(end-4) & -z <= f(end-3));
end
% append the last stopband
s = s | (abs(z) >= f(end));
% erode the stopband mask once to give a little extra slop
s([abs(diff(double(s))) ~= 0 false]) = false;
s = fliplr(s);
s([abs(diff(double(s))) ~= 0 false]) = false;
b2d = d*(1-0.01); % target beta^2 profile
roi = d | s;

function di = dinf(d1,d2)

% code from John Pauly's function

a1 = 5.309e-3;
a2 = 7.114e-2;
a3 = -4.761e-1;
a4 = -2.66e-3;
a5 = -5.941e-1;
a6 = -4.278e-1;

l10d1 = log10(d1);
l10d2 = log10(d2);

di = (a1*l10d1.*l10d1+a2*l10d1+a3)*l10d2+(a4*l10d1.*l10d1+a5*l10d1+a6);

