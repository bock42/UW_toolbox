function [pRFs HDR opt misc seeds] = pRF_Estimate(vtc,stim,s,opt,seeds)
% [pRFs HDR opt misc seeds] = pRF_Estimate(vtc,stim,s,opt,seeds)
% pRF fitting function
% modified from Geoff and Ione's prf fitting & prf/hdr iterative fitting codes
% PB, 03/2013
%
% inputs:
%
% vtc -> fMRI data
%     vtc(1:nvoxels).tc (2D vector, rows=TRs, columns=runs)
%     vtc(1:nvoxels).id (Voxel indices, useful for eventually visualizing maps)
%
% stim ->
%       .secPerTR   - TR duration in seconds
%       .numTRs     - stim.durSec / stim.secPerTr
%
% s  -> stimulus movie with the following fields:
%     s.frames (the stimulus) [4D: x,y,timeinarun,differentruns] 
%     s.x and s.y (x and y basis where the stimulus is defined)
%     s.nx and s.ny (length of x and y)
% 	  s.dx and s.dy (steps in degrees)
%
% opt -> optional inputs -- pass as [] if you don't want to specify this
%     [opt.gof ('correl','mse') measure of goodness of fit ] default is 'correl'
%     [opt.corrthr (0..1) threshold corr with seeds for nonlinear fit, if Inf never do fit] default is 0.25
%     [opt.parallel (1/0) use parallelization toolbox] default is 1
%     [opt.hdr.function (hdr used to fit pRFs)] default is Boynton '96
%     [opt.vtcsize (4D size of original BV vtc, useful for vmp creation)] default is []
%     [opt.rate4seeds (sampling rate for seeds)] default is 0.5 
%     [opt.freelist (parameters to fit)] default is {'center','sig'}
%     [opt.normalize (logical: tc = (tc - mean)/mean)] default is 1
%
%
% seeds -> optional -- pass as [] if you don't want to specify this
%     seeds.rfs(1:nseeds).center : x,y position
%     seeds.rfs(1:nseeds).sig : sigma (of the pos,neg lobe in case pRF model is DoG; pos,NaN if model is Gauss)
%     seeds.predmat -> matrix all predicted responses for all seeds 
%
%
% outputs:
%
% pRFs -> 
%     pRFs(1:nvoxels).center : x,y position
%     pRFs(1:nvoxels).sig : sigma (of the pos,neg lobe in case pRF model is DoG; pos,NaN if model is Gauss)
%     pRFs(1:nvoxels).co : correlation
%           -- the params above are set even if only the initial best seed is selected --
%           -- the params below are NaNs if nonlinear fit is not performed --
%     pRFs(1:nvoxels).err4fit : quantity that is minimized for fit (-corr or mse)
%     pRFs(1:nvoxels).mse : mean square error
% 	  pRFs(1:nvoxels).amp : scaling factor (not a free param)
%     pRFs(1:nvoxels).offset : offset (not a free param, dc of fMRI timecourse)
%     pRFs(1:nvoxels).id : voxel id (copied from vtc(:).id)
%
% HDR -> used for estimating pRFs (as input, or Boynton '96)
%
% opt -> as in input, but updated with defaults
%
% misc ->
%     misc.pRFestimation_time (in s)
% 
% seeds -> seeds used to initialize non linear fit (as input, or computed prior to parfor loop)
%
%
% functions required:
% UW_toolbox/Optimization (fit interface)
% convUnwrap.m
% fitPRF.m
% Gauss.m
% mycorr.m (Paola's correlation)
% Gamma.m (Geoff's) -> rename the function capitalizing the initial to
%                      avoid conflict with Matlab built-in function
%
% discontinued: DoG fitting, iterative HDR fitting.



%% set defaults
if ~isfield(opt,'gof'), opt.gof = 'correl'; end
if ~isfield(opt,'corrthr'), opt.corrthr = 0.25; end
if ~isfield(opt,'parallel'), opt.parallel = 1; end
if isfield(opt,'hdr')
    HDR = opt.hdr;
else %  Boynton et al. '96
    HDR.tau = 1.5; %seconds
    HDR.delay = 2.25;  %seconds
    HDR.n = 3;
    HDR.dt = stim.secPerTR;
    HDR.th = 0:HDR.dt:30;
    HDR.function = shiftdim(Gamma(HDR.n,HDR.tau,HDR.th-HDR.delay),-1);  %use shiftdim to make it a 1x1xn vector
end
if ~isfield(opt,'vtcsize'), opt.vtcsize = []; end
if ~isfield(opt,'rate4seeds'), opt.rate4seeds = 0.5; end
if ~isfield(opt,'freelist'), opt.freelist = {'center','sig'}; end
if ~isfield(opt,'normalize'), opt.normalize = 1; end

%% Create stimulus timecourses convolved by HDR 
% siz = size(s.frames);
% siz(4) = number of runs
% siz(3) = time (1:stim.numTRs)
% siz(1,2) = space (1:s.nx,1:s.ny)
nruns = size(s.frames,4);

% unwrap stimulus timecourse [time,space1d,nruns]
UnwrappedConvStim = convUnwrap(s,HDR,stim);

% vtc(nVx).tc must be [time,nruns], UnwrappedConvStim must be [time,space,nruns]
if size(vtc(1).tc,1) ~= size(UnwrappedConvStim,1) || size(vtc(1).tc,2) ~= size(UnwrappedConvStim,3)
    error('stim dimension and vtc dimensions do not match')
end

%% Normalize fMRI timecourses (each run separately, each voxel separately)
if opt.normalize
    fprintf('normalizing timecourses ...\n')
    for voxNum = 1:length(vtc)
        ObsResp = vtc(voxNum).tc;
        
        dc = mean(ObsResp,1);
        dc = repmat(dc,size(ObsResp,1),1);
        ObsResp = ObsResp./dc;
        
        vtc(voxNum).tc = ObsResp;
    end
end
%% Seeds for fit initialization
tic
% define a series of potential pRFs that are used as starting
% points for the nonlinear optimization process (unless they were passed as input)
if ~isfield(seeds,'rfs')
    seedSig = linspace(1,10,5); %
    step = opt.rate4seeds;
    seedX = round(min(s.x(:))) : step : round(max(s.x(:)));
    seedY = round(min(s.y(:))) : step : round(max(s.y(:)));
    % create all possible combinations
    [xList,yList,sigList] = ndgrid(seedX,seedY,seedSig);
    % unwrap into 1D vectors length(seedX)*length(seedY)*length(seedSig) long
    xList = xList(:);
    yList = yList(:);
    sigList = sigList(:);
    for i = 1:length(xList)
        seeds.rfs(i).center = [xList(i),yList(i)];
        seeds.rfs(i).sig = [sigList(i) NaN];
    end
end
nseeds = length(seeds.rfs);
fprintf('number of seeds: %i x %i runs\n',nseeds,nruns)

if ~isfield(seeds,'predmat')
    disp('creating predicted timecourses ...')
    % columns of matSeeds will hold the unwrapped Gaussian pRFs
    % matrix size: [s.nx*s.ny, Number of runs, length(seeds.rfs)]
    seeds.predmat = zeros(size(UnwrappedConvStim,1),size(UnwrappedConvStim,3),length(seeds.rfs));
    createseedsmat = 1;
else
    createseedsmat = 0;
    disp('predicted timecourses passed as input')
end
for i = 1:nseeds
    % update fields
    seeds.rfs(i).shutup = 1;
    if createseedsmat
        [~,seeds.predmat(:,:,i)] = fitPRF(seeds.rfs(i),UnwrappedConvStim,[],s,opt);
        if ~mod(i,100), fprintf('..%i..',i); end
    end
end
fprintf('\n')

%% initialize pRFs struct
parfor voxNum = 1:length(vtc)
    
    pRFs(voxNum).center = [NaN NaN]; % x,y
    pRFs(voxNum).sig = [NaN NaN]; % pos and neg lobe (if pRF model is DoG)
    pRFs(voxNum).co = NaN;
    pRFs(voxNum).mse = NaN;
    pRFs(voxNum).err4fit = NaN;
    pRFs(voxNum).amp = [NaN NaN];
    pRFs(voxNum).offset = NaN;
    pRFs(voxNum).id = NaN;
    
end
toc

%% estimate pRFs
% initialize progress monitor bar -> dwnld free pkg "ParforProgMon"
if opt.parallel, ppm = ParforProgMon([ 'pRF fitting started at '  datestr(now,'HH:MM')], length(vtc) ); else ppm.increment = NaN; end
tic
disp('start pRFestimate, looping through voxels ...')

parfor voxNum = 1:length(vtc)
    
    % initializations required by parfor
    vN = voxNum;
    ObsResp = vtc(vN).tc;
    
    % show current status
    if ~opt.parallel, fprintf('voxNum %d of %d \n',vN,size(vtc,2)); end
    ppm.increment();
    
    % save voxel .id in pRF struct
    pRFs(voxNum).id = vtc(vN).id;
    
    % find best seed to initialize pRF fitting
    basiscorr = NaN(nruns,nseeds);
    for runNum = 1:nruns
        basiscorr(runNum,:) = mycorr(ObsResp(:,runNum),squeeze(seeds.predmat(:,runNum,:)));
    end
    basiscorr = mean(basiscorr,1); % average across runs
    [maxcorr,bestseed] = max(basiscorr); % find the best seed
    
    if maxcorr < opt.corrthr || isnan(maxcorr) 
        % only store best seed (NaNs happen for flat timecourses)
        pRFs(voxNum).co = maxcorr;
        pRFs(voxNum).center = seeds.rfs(bestseed).center;
        pRFs(voxNum).sig = seeds.rfs(bestseed).sig;
        
    else
        % perform non linear fitting
        disp(sprintf('..%.2f..',maxcorr))
%         disp(sprintf('..%i..',vN))
        
        [tmpRF err] = fit('fitPRF',seeds.rfs(bestseed),opt.freelist,UnwrappedConvStim,ObsResp,s,opt);
        pRFs(voxNum).err4fit = -err;
        
        pRFs(voxNum).center = tmpRF.center;
        pRFs(voxNum).sig = tmpRF.sig;

        [~,~,correl,mse,polyVals] = fitPRF(pRFs(voxNum),UnwrappedConvStim,ObsResp,s,opt);
        pRFs(voxNum).amp = [polyVals(1) NaN];
        pRFs(voxNum).offset = polyVals(2);
        pRFs(voxNum).co = -correl;
        pRFs(voxNum).mse = mse;
        
    end
end

misc.pRFestimation_time = toc;

disp(sprintf('\n prf estimation completed in %i min \n', round(misc.pRFestimation_time/60)))
fprintf('\n\n\n')


