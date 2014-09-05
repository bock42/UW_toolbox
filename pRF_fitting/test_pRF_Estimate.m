% create 'true pRFs' (generate responses, set options, etc) to recover with pRF_Estimate
% PB 03/2013

clear all; close all; clc;
% myseed = 1; rand('seed',myseed); randn('seed',myseed);

% set path to output file (where prfs are saved)
output_filename = 'myprfs.mat';

% set path to progress bar directory
progrbardir = '../../Parallelization/ParforProgMon/';

% % reload seeds for fit initialization (!!!)
% seeds_filename = output_filename;
% load(seeds_filename,'seeds')

% set number of runs
nruns = 1;

% set noiselevel
noiselevel = 0.05; %.1 corresponds to co about .6   .2 -> .4  .05 -> .85
howmanydims = 2; % dimensions of the stimulus (1,2 or 1/2=weight)

% set optional inputs to pRF_Estimate
% type help pRF_Estimate for the list of available options
opt.parallel = 1;
% opt.corrthr = Inf; % if Inf only finds best seed, no nonlinear fit
% opt.gof = 'mse';
% opt.rate4seeds = 1; % sampling rate of seeds matrix
ind = 1; % for Gauss (DoG discontinued)

%% define space & time
lims = [-8 8];
x = linspace(lims(1),lims(2),21);

if howmanydims == 2 % 2d
    [s.x,s.y] = meshgrid(x);
    s.nx = size( s.x,1 ) ;
    s.ny = size( s.x,2 ) ;
    s.dx = s.x(1,2) - s.x(1,1);
    s.dy = s.y(2,1) - s.y(1,1);
elseif howmanydims == 1 % 1D
    s.x = x'; s.x = Shuffle(s.x);
    s.y = 0;
    s.nx = length( s.x ) ;
    s.ny = 1 ;
    s.dx = s.x(2) - s.x(1);
    s.dy = 1;
elseif howmanydims == 1.5 % both x and y are defined, but both are 1D
    s.x = x'; s.x = Shuffle(s.x);
    s.y = x'; s.y = Shuffle(s.x);
    s.nx = length( s.x ) ;
    s.ny = 1;%length( s.y ) ;
    s.dx = 1; %s.x(2) - s.x(1);
    s.dy = 1;
else
    error('howmanydims out of range')
end

%% define stimulus

% spatiotemporal matrix (each position stimulated at one timepoint, alone)
stim.numTRs = numel(s.x);
s.frames = zeros(s.nx,s.ny,stim.numTRs,nruns);
for t = 1:stim.numTRs
    tmp = s.frames(:,:,t,1);
    tmp(t) = 1;
    for runNum = 1:nruns
        s.frames(:,:,t,runNum) = tmp;
    end
end
stim.secPerTR = 1;

%% stimulus movie
% visualize stimulus 's' with flag:
% 1 = show mean and var of the stimulus; 
% 2 = show (accelerated) stimulus movie; 
% 3 = 1D plot of stim present/absent vs time
% 4 = 1 screenshot
figure(10); clf
plotstim(s,1)

%% convolve with hdr & unwrap
dt = 1;
HDR.function = 1;
opt.hdr = HDR;
% unwrap stimulus timecourse [time,space1d,nruns]
UCStim = convUnwrap(s,HDR,stim);

%% define true pRFs
y = 0;
xList = [0.1 .2 .3 .4 .5 1 1.5 2 3 4];
sigList = [1 2 3 4];
ampList = 1 + 0*rand(size(xList));
offList = 1 + 0*rand(size(xList));
unwrap = 0;
vx = 0;
for cs = 1:length(sigList)
    for cx = 1:length(xList)
        vx = vx + 1;
        tRF(vx).center = [xList(cx) y];
        tRF(vx).sig = [sigList(cs) NaN];
        tRF(vx).co = NaN;
        tRF(vx).amp = [ampList(cx) NaN];
        tRF(vx).offset = offList(cx);
        
        myfun = Gauss( tRF(vx), s.x, s.y, ind, unwrap );
        
        % noise proportional to sigma and amp (its effect should be approx constant)
        tRF(vx).noisecoeff = noiselevel * sigList(cs)^.75 * max(myfun(:)) * tRF(vx).amp(ind);
%         tRF(vx).noisecoeff = noiselevel * sigList(cs).^.5 * max(myfun(:)) * tRF(vx).amp(ind);
    end
end
nVx = length(tRF);


%% create responses

for vx = 1:nVx
    vtc(vx).id = vx;
    % save the timecourse that is passed to pRF_Estimate
    [~,myresp] = fitPRF(tRF(vx),UCStim,[],s,opt);
    noise = randn(size(myresp)) * tRF(vx).noisecoeff;
    vtc(vx).tc = myresp * tRF(vx).amp(ind) + tRF(vx).offset + noise;
end

%% call pRF_estimate (and open parallel pool, if an option)
if opt.parallel
    if matlabpool('size') == 0 % checking to see if my pool is already open
        matlabpool open 6 % number of cores you have
    end
    currdir = cd(progrbardir); 
    pctRunOnAll javaaddpath(cd)
    cd(currdir)
    
    try  jm = findResource('scheduler', 'configuration', defaultParallelConfig);
    catch me
        parallel = 0;
    end
    
end
disp('calling pRF_estimate function...')
if ~exist('seeds')
    seeds = [];
end
[pRFs HDR opt misc seeds] = pRF_Estimate(vtc,stim,s,opt,seeds);

save(output_filename, 'pRFs', 'HDR', 'opt', 'misc', 'tRF', 'vtc','seeds');


%% plot results
tmp = cat(1,tRF.center); 
[pol, ecc] = cart2pol(tmp(:,1),tmp(:,2));
trueRFs = [tmp ecc pol cat(1,tRF.sig) cat(1,tRF.noisecoeff) cat(1,tRF.amp) cat(1,tRF.offset)];
tmp = cat(1,pRFs.center); 
[pol, ecc] = cart2pol(tmp(:,1),tmp(:,2));
estRFs = [tmp ecc pol cat(1,pRFs.sig) cat(1,pRFs.co) cat(1,pRFs.amp) cat(1,pRFs.offset)];
mycols = [3 4 5 7]; % ecc,polang,sig,corr
mycolors = hsv(length(sigList));

% all together
figure(1); clf; set(gcf,'position',[360,9,340,913;])
for i = 1:length(mycols)
    subplot(length(mycols),1,i); hold on
    myx = trueRFs(:,mycols(i)); % ecc
    myy = estRFs(:,mycols(i));
    plot(myx,myy,'.','color',mycolors(cs,:));
    plot(myx,myx,':k')
end

% group by sigma vals and plot against eccentricity
figure(2); clf; set(gcf,'position',[360+300,9,340,913;])
for cs = 1:length(sigList)
    thisSig = sigList(cs);
    myleg{cs} = num2str(thisSig);
    ind = trueRFs(:,5) == thisSig;
    
    mytrue = trueRFs(ind,:);
    myest = estRFs(ind,:);
    
    for i = 1:length(mycols)
        subplot(length(mycols),1,i); hold on
        myx = mytrue(:,mycols(1)); % ecc
        myy = myest(:,mycols(i));
        h(cs) = plot(myx,myy,'.-','color',mycolors(cs,:));
        if i == 1
            plot(myx,myx,':k')
            axis equal
%             ylim(xlim)
        else
            plot(myx,myx*0+myy(end),':','color',mycolors(cs,:))
        end
        
    end
    
end
legend(h,myleg)

%% plot correlation vs. noiselevel
figure(3); clf; hold on
myx = cat(1,tRF.sig); myx = myx(:,1);
myy = cat(1,pRFs.co); 
plot((myx),(myy),'.','color',mycolors(cs,:));
% plot(myx,myx,':k')
