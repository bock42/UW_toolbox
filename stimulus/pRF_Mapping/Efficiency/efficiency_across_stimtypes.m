% compute efficiency of stimulus sequence
cls

stimList = {'DriftingBar'; 'RandomBar'; 'MultiFocal'}; 
randangList = {0;1};

stimList = stimList(1:2);
randangList = randangList(2);

homedir = cd;
targetdir = '/data/tmp/';

for st = 1:length(stimList)
    for rd = 1:length(randangList)
        
        % load stimulus file
        d = dir([homedir targetdir '*' stimList{st} '*scot0*.mat']);
        for cf = length(d):-1:1
            load([homedir targetdir d(cf).name])
            if stim.randAng == randangList{rd}
                disp(d(cf).name)
                break
            end
        end
        if stim.randAng
            shfstr = 'rndAng';
        else
            shfstr = '45steps';
        end
        titstr{st,rd} = [stim.StimulusType ' ' shfstr];
        
        % make stim movie
        % trim display area (so that the stimulus movie only represents the area of interest)
        display.screenAngle = [2 2]*stim.radDeg; % width,height in deg
        % lower the resolution
        stimres = 0.25*1; % deg
        display.resolution = display.screenAngle / stimres;
        display.width =  2 * display.dist * tan( pi * display.screenAngle(1) / (2*180) ) ; % reconvert (function uses width rather than angle)
        s = makeStimulusMovie( display, stim );
        
        % compute efficiency
        nhdr = 15; % try to estimate an hdr with 15 time points
        E{st,rd} = NaN(s.nx,s.ny);
        for ix = 1:s.nx
            for iy = 1:s.ny
                % compute a timecourse of this position (on/off)
                tc = squeeze(s.frames(ix,iy,:));
                if sum(tc) > 0
                    % compute efficiency of stimulation for this position
                    E{st,rd}(ix,iy) = efficiency(tc,nhdr);
                end
            end
        end
    end
end


%% plot
clim = [.5 1.5];

figure(3); clf
for st = 1:length(stimList)
    for rd = 1:length(randangList)
        % map efficiency across space
        mysubplot([length(stimList),length(randangList)],st,rd); hold on
        imagesc(s.x(1,:),s.y(:,1),E{st,rd},clim); colorbar; set(gca,'ydir','normal')
%         text(-5,9,sprintf('avg eff %.3f\n',nanmean(E{st,rd}(:))))
        title(sprintf('%s \navg eff %.2f',titstr{st,rd},nanmean(E{st,rd}(:))))
        axis off square
    end
end

colormap(jet)

