%% set up evaluation parameters
disp 'Getting evaluation parameters and processing maps...'

rawMapsDir = '../../pTx_rawMaps/';
targetMapsDir = '../../pTx_maps/';
evalp = pTxParams_torso_phaseI(rawMapsDir); % get evaluation parameter structure

% load source b1/df0 maps and tissue mask ('maps' structure)
% maps.b1: [Nx Ny Nz Nc] (uT), complex B1+ maps
% maps.df0: [Nx Ny Nz] (Hz), off-resonance map
% maps.mask: [Nx Ny Nz] (logical), tissue mask
% maps.fov: [3] (cm), fov in each dimension
% maps.vop: [Nc Nc Nvop], SAR VOPs
% maps.maxSAR: [Nvop], max SAR at each VOP
load([rawMapsDir evalp.rawmapsfname]);
maps = rmfield(maps,'mask_organs');

% get the target pattern, interpolated maps, error roi, and spatial grid
[maps,evalp] = interp_pTx_maps(maps,evalp);

save([targetMapsDir evalp.evalmapsfname]);
