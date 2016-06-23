function [rf,g,dt,m,kx,ky] = dzSpokes_multiband(maps,probp)

% Spokes parallel transmit RF pulse design example
% Written for 2014 ISMRM Nuts and Bolts Course
% Adapted for ISMRM RF Pulse Design Challenge
% Will Grissom, Vanderbilt University, 2016
% will.grissom@vanderbilt.edu
%

if ndims(maps.b1) == 4 % we are designing a multiband spokes pulse
    [dim(1),dim(2),Nsl,Nc] = size(maps.b1); % get b1 map matrix size and # tx channels
else % we are designing a single band pulse
    [dim(1),dim(2),Nc] = size(maps.b1); % get b1 map matrix size and # tx channels
    Nsl = 1;
end
maps.b1 = reshape(maps.b1,[prod(dim) Nsl Nc]); % vectorize each slice's b1 maps
maps.maskv = reshape(maps.mask,[prod(dim) Nsl]); % vectorize each slice's mask

%% Step 1: set up design grid and (kx,ky) search grid
x = maps.xy(:,:,1);
y = maps.xy(:,:,2);
kmax = 1/probp.deltaxmin; % /cm, max spatial freq of traj
[kxs,kys] = meshgrid(-kmax/2:1/max(maps.fov):kmax/2-1/max(maps.fov)); % greedy kx-ky grid
kxs = kxs(:);kys = kys(:); % vectorize grid
dc = find((kxs == 0) & (kys == 0)); % find DC point
if ~isempty(dc)
    kxs = kxs([1:dc-1 dc+1:end]);kys = kys([1:dc-1 dc+1:end]); % remove DC point
end

%% Step 2: design the weights
kx = 0; % initial kx/ky location is DC
ky = 0;
if isfield(probp,'phsinit')
    for kk = 1:Nsl
        tmp = probp.phsinit(:,:,kk);
        phs{kk} = tmp(maps.maskv(:,kk));
    end
else
    for kk = 1:Nsl
        phs{kk} = zeros(sum(maps.maskv(:,kk)),1);
    end
end
m = {};

b1tmp = reshape(maps.b1,[prod(dim)*Nsl Nc]);
Afi = inv(b1tmp(maps.mask(:),:)'*b1tmp(maps.mask(:),:)); % precalc inverse of greedy matrix
for ii = 1:probp.Nspokes

    % build system matrix for each slice
    Afull = {};
    for kk = 1:Nsl
        A = exp(1i*2*pi*(x(maps.maskv(:,kk))*kx(:)' + y(maps.maskv(:,kk))*ky(:)')); % Fourier matrix
        Afull{kk} = zeros(sum(maps.maskv(:,kk)),ii*Nc);
        for jj = 1:Nc
            % multiply B1 maps and Fourier kernel - bsxfun does this quickly
            Afull{kk}(:,(jj-1)*ii+1:jj*ii) = bsxfun(@times,maps.b1(maps.maskv(:,kk),kk,jj),A);
        end
    end
    
    % calculate magnitude least-squares-optimal RF weights
    mAll = zeros(dim(1),dim(2),Nsl);
    wfull = zeros(ii*Nc,Nsl);
    favar = zeros(Nsl,1);
    for kk = 1:Nsl
        wfull(:,kk) = inv(Afull{kk}'*Afull{kk} + probp.beta*eye(ii*Nc))*(Afull{kk}'*exp(1i*phs{kk}));
        err = Afull{kk}*wfull(:,kk) - exp(1i*phs{kk});
        cost = real(err'*err + probp.beta*wfull'*wfull);
        costold = 10*cost; % to get loop going
        while abs(cost-costold) > 0.0001*costold
            costold = cost;
            phs{kk} = angle(Afull{kk}*wfull(:,kk)); % set phase to current excited phase
            wfull(:,kk) = (Afull{kk}'*Afull{kk} + probp.beta*eye(ii*Nc))\(Afull{kk}'*exp(1i*phs{kk}));
            err = Afull{kk}*wfull(:,kk) - exp(1i*phs{kk});
            cost = real(err'*err + probp.beta*wfull(:,kk)'*wfull(:,kk));
            m = zeros(dim);
            m(maps.maskv(:,kk)) = Afull{kk}*wfull(:,kk);
        end
        mAll(:,:,kk) = m;
    end
    
    fprintf('Number of spokes: %d. Normalized flip angle variance: %0.2d\n',ii,var(abs(mAll(maps.mask))));

    % add a spoke using greedy method
    if ii < probp.Nspokes
        r = [];
        for kk = 1:Nsl
            r = [r; exp(1i*phs{kk}) - Afull{kk}*wfull(:,kk)];
        end
        rfnorm = [];
        for jj = 1:length(kxs)
            Afull = [];
            for kk = 1:Nsl
                Afull = [Afull; bsxfun(@times,squeeze(maps.b1(maps.maskv(:,kk),kk,:)),exp(1i*2*pi*(x(maps.mask(:,:,kk))*kxs(jj)' + y(maps.mask(:,:,kk))*kys(jj))))];
            end
            rfnorm(jj) = norm(Afi*(Afull'*r));
        end
        [~,ind] = max(rfnorm);
        % alternate which end of the pulse we add the spoke to
        if rem(ii + 1,2)
            kx = [kx;kxs(ind)];ky = [ky;kys(ind)];
        else
            kx = [kxs(ind);kx];ky = [kys(ind);ky];
        end
        % remove selected point from search grid
        kxs = kxs([1:ind-1 ind+1:end]);
        kys = kys([1:ind-1 ind+1:end]);
    end

end

%% Step 3: build the full waveforms
area = probp.tb/(probp.slthick/10)/4258; % slthick*kwidth = tbw; kwidth = gamma*area
[subgz,nramp] = domintrap(area,probp.maxg/10,probp.maxgslew*100,probp.dt); % min-duration trapezoid
subgz = subgz*10; % mT/m
if ~isempty(which('dzrf')) % check if Pauly SLR tool is available
    subrf = dzrf(128,probp.tb,'st'); % get small-tip slice-selective subpulse
    save subrf subrf;
else
    load subrf;
end
Nplat = length(subgz) - 2*nramp; % # time points on trap plateau
subrf = interp1((0:127)/128,subrf,(0:Nplat-1)/Nplat,'spline',0);
subrf = [zeros(nramp,1);subrf(:);zeros(nramp,1)];
subrf = subrf./sum(subrf);

% Build gradient waveforms
gxarea = diff([kx; 0])/4258; % no (-) sign since we don't time-flip k
gyarea = diff([ky; 0])/4258;

gx = [];gy = [];gz = [];
for ii = 1:probp.Nspokes
    gz = [gz;((-1)^(ii-1))*subgz(:)];
    gx = [gx;zeros(length(subgz),1)];
    if abs(gxarea(ii)) > 0
        gxblip = sign(gxarea(ii))*dotrap(abs(gxarea(ii)),probp.maxg/10,probp.maxgslew*100,probp.dt)*10;
        gx(end-length(gxblip)+1:end) = gxblip;
    end
    gy = [gy;zeros(length(subgz),1)];
    if abs(gyarea(ii)) > 0
        gyblip = sign(gyarea(ii))*dotrap(abs(gyarea(ii)),probp.maxg/10,probp.maxgslew*100,probp.dt)*10;
        gy(end-length(gyblip)+1:end) = gyblip;
    end
end

gzref = ((-1)^probp.Nspokes)*dotrap(probp.dt*sum(subgz)/2/10,probp.maxg/10,probp.maxgslew*100,probp.dt)'*10;
gz = [gz;gzref];
gx = [gx;zeros(size(gzref))];
gy = [gy;zeros(size(gzref))];
g = [gx(:) gy(:) gz(:)];

% Build full RF waveforms
rf = zeros(length(gz),Nc);
for kk = 1:Nsl
    rft = reshape(kron(wfull(:,kk),subrf),[length(gz)-length(gzref) Nc]);
    rft = [rft;zeros(length(gzref),Nc)];
    
    % Scale to the target flip angle
    rft = rft*(probp.flipAngle/180*pi)/(probp.dt*42.57*2*pi);
    
    % modulate to slice position
    kz = -4258*flipud(cumsum(flipud(gz(:))))/10*probp.dt; % cycles/cm
    rft = rft.*repmat(exp(-1i*2*pi*kz*probp.slCent(kk)),[1 Nc]);
    
    % add on with the rest of the slices' waveforms
    rf = rf + rft;
end

% if RF violates peak B1 constraint, adjust dt
if max(abs(rf(:))) > probp.maxb1
    fprintf('Increasing spokes dwell time by factor of %f to meet peak RF constraints\n',max(abs(rf(:)))/(0.99*probp.maxb1));    
    dt = max(abs(rf(:)))/(0.99*probp.maxb1)*probp.dt;
    rf = rf*probp.dt/dt;
    g = g*probp.dt/dt;
else
    dt = probp.dt;
end

% if RF violates peak SAR, increase dt
if isfield(maps,'vop')
    nVOP = size(maps.vop,3); % # VOPs
    sarVOP = zeros(nVOP,size(rf,1));
    for ii = 1:size(rf,1)
        rft = rf(ii,:); % all coil's samples for this time point
        for jj = 1:nVOP
            sarVOP(jj,ii) = real(rft(:)'*(maps.vop(:,:,jj)*rft(:)));
        end
    end
    sarVOP = sum(sarVOP,2) * dt./probp.TR;
    sarVOPfrac = sarVOP./maps.maxSAR(:);
    if any(sarVOPfrac > 1) % if we are violating an SAR constraint
        fprintf('Increasing spokes dwell time by factor of %f to meet SAR constraints\n',max(sarVOPfrac)*0.99);
        dt = dt*max(sarVOPfrac)/0.99;
        rf = rf./max(sarVOPfrac)*0.99;
        g = g./max(sarVOPfrac)*0.99;
    end
end

%% Plot the results
figure
subplot(212)
plot((0:length(gx)-1)*dt*1000,[gx(:) gy(:) gz(:)])
axis([0 length(gx)*dt*1000 -probp.maxg probp.maxg]);
xlabel 'ms',ylabel 'mT/m'
title 'Gradient Waveforms'
subplot(211)
plot((0:length(gx)-1)*dt*1000,real(rf))
ysc = max(abs(min(real(rf(:)))),abs(max(real(rf(:)))));
axis([0 length(gx)*dt*1000 -1.1*ysc 1.1*ysc]);
xlabel 'ms',ylabel 'a.u.'
title 'RF Waveforms (real part)'

figure;
for ii = 1:Nsl
    subplot(1,Nsl,ii)
    imagesc(flipud(abs(mAll(:,:,ii).')*probp.flipAngle),...
        [min(abs(mAll(:))) max(abs(mAll(:)))]*probp.flipAngle);
    title(sprintf('Slice %d flip angle (degrees)',ii)),axis off,axis image;colorbar
    colormap jet
end