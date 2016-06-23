function [maps,evalpOut] = interp_ptx_maps_multislice(maps,evalp)

% function [maps,evalp] = proc_ptx_maps(maps,evalp)
%
% Given the source B1+/df0 maps, and the evaluation parameters,
% generate the target flip angle pattern, interpolate the maps,
% generate an error roi mask, and generate a spatial grid
%
% Will Grissom, 2015, for the ISMRM RF Pulse Design Challenge

%% get the evaluation grid
x = -maps.fov(1)/2:evalp.dxyz(1):maps.fov(1)/2;
y = -maps.fov(2)/2:evalp.dxyz(2):maps.fov(2)/2;
z = -maps.fov(3)/2:evalp.dxyz(3):maps.fov(3)/2;
[x3,y3,z3] = ndgrid(x,y,z);
evalp.xyz = [x3(:) y3(:) z3(:)]; % columnize for system matrix construction

%% interpolate B1+/df0 maps onto evaluation grid
[xin,yin,zin] = ndgrid(-maps.fov(1)/2:maps.fov(1)/size(maps.mask,1):maps.fov(1)/2-maps.fov(1)/size(maps.mask,1),...
        -maps.fov(2)/2:maps.fov(2)/size(maps.mask,2):maps.fov(2)/2-maps.fov(2)/size(maps.mask,2),...
        -maps.fov(3)/2:maps.fov(3)/size(maps.mask,3):maps.fov(3)/2-maps.fov(3)/size(maps.mask,3)); % input grid
%maps.df0 = interp3(xin,yin,zin,maps.df0,x,y,z,'linear',0); % interpolate
evalp.Nc = size(maps.b1,4); % # Tx channels
tmp = maps.b1;
maps = rmfield(maps,'b1');
for ii = 1:evalp.Nc
    % interpolate each tx coil's b1 map, after applying tissue mask
    maps.b1(:,:,:,ii) = interpn(xin,yin,zin,maps.mask.*tmp(:,:,:,ii),x3,y3,z3,'linear',0);
end
% recalculate tissue mask
maps.mask = interpn(xin,yin,zin,double(maps.mask),x3,y3,z3,'linear',0);
maps.mask = maps.mask > 0.9999; % exclude interpolated edges
%tmp = ssq(maps.b1,4);
%maps.mask = maps.mask & tmp > 0.02; % found some almost zero voxels in b1; this gets rid of those
%maps.mask = sum(abs(maps.b1),4) > 0;

%% get the target pattern and roi
w = dinf(evalp.d1,evalp.d2)/evalp.tb; % fractional transition width

%% get the f, m and w vectors with the DC band
f = [0 (1-w)*(evalp.tb/2) (1+w)*(evalp.tb/2)]*evalp.slThick/10/evalp.tb; % band edges (cm)

%% average the B1+ maps over the slice thickness so we can check for flat through-slice phase
evalp.inSliceInds = {}; % may have different # of in-slice inds for each band, depending on how grid falls
for ii = 1:length(evalp.slCent) % loop over slice positions
    inSliceIndst = find(abs(z - evalp.slCent(ii)) <= f(2)); % z-indices for this slice
    b1Ave = squeeze(mean(maps.b1(:,:,inSliceIndst,:),3)); % average b1 thru slice
    b1Ave = permute(repmat(b1Ave,[1 1 1 length(inSliceIndst)]),[1 2 4 3]); % duplicate to all in-slice locs
    maps.mask(:,:,inSliceIndst) = repmat(mean(double(maps.mask(:,:,inSliceIndst)),3) == 1,[1 1 length(inSliceIndst)]); % only include mask points that are in every sub-slice
    maps.mask(:,:,inSliceIndst) = repmat(mean(double(maps.mask(:,:,inSliceIndst)),3) == 1,[1 1 length(inSliceIndst)]); % only include mask points that are in every sub-slice
    maps.b1(:,:,inSliceIndst,:) = b1Ave;
    evalp.inSliceInds{ii} = inSliceIndst;
end
maps.b1 = maps.b1.*repmat(maps.mask,[1 1 1 evalp.Nc]); % re-apply mask to get rid of in-plane locs that may have had zeros in some sub-slices

%% remove first coil's phase from others - is best to do this after interpolation, since it can introduce rapid phase variations
maps.b1 = maps.b1.*repmat(exp(1i*angle(maps.b1(:,:,:,1))),[1 1 1 evalp.Nc]);

%% get passband and stopband masks
evalpOut = {};
for ii = 1:length(evalp.slCent)
    
    dt = abs(z3 - evalp.slCent(ii)) <= f(2); % target slice mask
    st = abs(z3 - evalp.slCent(ii)) >= f(3); % stopband mask
    
    % erode the stopband mask once to give a little extra slop
    st(abs(diff(double(st),1,3)) ~= 0) = false;
    st = flip(st,3);
    st(abs(diff(double(st),1,3)) ~= 0) = false;
    st = flip(st,3);
    
    evalpOut{ii} = evalp;
    evalpOut{ii}.slCent = evalp.slCent(ii);
    evalpOut{ii}.inSliceRoi = dt & maps.mask; % In-slice ROI for error calc
    evalpOut{ii}.outOfSliceRoi = st & maps.mask; % Out-of-Slice ROI for error calc
    evalpOut{ii}.thetad = (dt & maps.mask)*evalp.flipAngle; % target pattern
    evalpOut{ii}.inSliceInds{1} = evalp.inSliceInds{ii};
    
end

%% subtract off median frequency offset in slice (presumably the spectrometer would be set here)
%tmp = maps.df0(evalp.inSliceRoi);
%maps.df0 = maps.df0 - median(tmp);
