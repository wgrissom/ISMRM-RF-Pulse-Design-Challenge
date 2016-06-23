function [a,b] = blochsim_simult(rf,gamgdt,xx)

[Ns,~] = size(xx); % Ns: # spatial locs. Nd: # of dimensions
Nt = size(gamgdt,1); % Nt: # time points
a = ones(Ns,1);
b = zeros(Ns,1);

for ii = 1:Nt
    
    % calculate rotation vector
    phi = -sqrt(abs(rf(ii))^2 + (gamgdt(ii)*xx).^2);
    nxy = rf(ii)./-phi;nxy(isinf(nxy)) = 0;nxy(isnan(nxy)) = 0;
    nz = gamgdt(ii)*xx./-phi;nz(isinf(nz)) = 0;nz(isnan(nz)) = 0;
    
    % calculate alpha and beta for this time point
    alpha = cos(phi/2) - 1i*nz.*sin(phi/2);
    beta = -1i*nxy.*sin(phi/2);
    
    % apply alpha and beta
    at = alpha.*a - conj(beta).*b;
    bt = beta.*a + conj(alpha).*b;
    a = at; b = bt;
  
end

