% simulate runs and compute efficiency
clear all; close all; clc
stimList = {'DriftingBar'; 'RandomBar'};

stim.scotoma.flag = 1;
niter = 1;
%% default stimulus params
% may be overwritten for specific experiments

% carrier
stim.cyclePerDeg = .5; % checkerboard
stim.flickerHz = 8;

% stim
stim.radDeg = 8; % circular aperature
stim.widthDeg = 2; % bar width
stim.randAng = 1; % logical - if 0 bar orientation varies in steps of 45 deg, if 1 orient is random
stim.degPerSecRange = [ 1 2 ]; % range of speeds
if ~stim.randAng, stim.degPerSecRange = [1 1]*stim.widthDeg; end % non overlapping steps
stim.durBlanksSec = 15; % durantion of blanks between bar sweeps
stim.nSweepsBetweenBlanks = 4; % number of sweeps separating two blank periods
stim.begendBlankSec = 5; % duration of blanks at beginning and end of run

% time
stim.durSec = 240; %
stim.secPerTR = 1;
stim.numTRs = stim.durSec / stim.secPerTR;

% others
stim.fixMaskRadDeg = 0.25; % mask around fixation
%stim.seed = 8194; % always the same 'random' bar orient (Zach: 8194 best seed out of 1:1e5)
stim.scotoma.rad = [1 1] * 2 ;
stim.scotoma.center = [3 0];

% define display area
display.dist = 68;
display.screenAngle = [2 2]*stim.radDeg; % width,height in deg
% lower the resolution
stimres = 0.25*1; % deg
display.resolution = display.screenAngle / stimres;
display.width =  2 * display.dist * tan( pi * display.screenAngle(1) / (2*180) ) ; % reconvert (function uses width rather than angle)


me = NaN(niter,length(stimList));
for st = 1:length(stimList)
    stim.StimulusType = stimList{st};
    E{st} = NaN([display.resolution,niter]);
    for i = 1:niter
        switch stim.StimulusType
            
            case 'DriftingBar'
                
                % make bar seq
                stim = getBarSequence( stim );
                
            case 'RandomBar'
                
                % make bar seq
                stim = getBarSequence( stim );
                
                % shuffle bar positions (but preserve blanks)
                tbscrambled = find(stim.degPerSec); % blanks are marked as velocity = 0
                neworder = Shuffle(tbscrambled);
                stim.centersDeg(tbscrambled,:) = stim.centersDeg(neworder,:);
                stim.motDirsRad(tbscrambled) = stim.motDirsRad(neworder);
                stim.degPerSec(tbscrambled) = stim.degPerSec(neworder);
        end
        
        % make stim movie
        s = makeStimulusMovie( display, stim );
        
        
        % compute efficiency
        nhdr = 15; % try to estimate an hdr with 15 time points
        
        tmp = NaN(s.nx,s.ny);
        for ix = 1:s.nx
            for iy = 1:s.ny
                % compute a timecourse of this position (on/off)
                tc = squeeze(s.frames(ix,iy,:));
                if sum(tc) > 0
                    % compute efficiency of stimulation for this position
                    tmp(ix,iy) = efficiency(tc,nhdr);
                end
            end
        end
        E{st}(:,:,i) = tmp;
        me(i,st) = nanmean(tmp(:));
    end
    avgE{st} = nanmean(E{st},3);
    figure(10+st); clf
    plot(stim.motDirsRad/pi,'.-')
    xlim([1 stim.numTRs])

end


%% plot
figure(1); clf

clim = [.5 1.5];

for st = 1:length(stimList)
    mysubplot([2 length(stimList)],1,st)
    % map efficiency across space
    imagesc(s.x(1,:),s.y(:,1),avgE{st},clim); colorbar; set(gca,'ydir','normal')
    %         text(-5,9,sprintf('avg eff %.3f\n',nanmean(E{st,rd}(:))))
    title(sprintf('%s\nefficiency:\n(average across runs)\n%.2f',stimList{st},nanmean(avgE{st}(:))))
    axis tight square
end
colormap(jet)
if niter>1
    mysubplot([2 length(stimList)],2,1:2)
    plot(me,'.-')
    legend(stimList,'location','best')
    xlabel('run')
    ylabel('efficiency (average across space)')
end
