% compute efficiency of stimulus sequence
cls

% locate files
% homedir = cd;
% targetdir = '/data/tmp/';
homedir = '/media/VERBATIM/2013_03_19_SCOTOMA';
targetdir = '/20130319/axb/';


stimList = {'DriftingBar'; 'RandomBar'; 'MultiFocal'};
randangList = {0;1};
scotoma = {'0';'1'};

mystim = stimList{2};
randang = randangList{2};


% load stimulus file
avgE = cell(length(scotoma),1);
for sc = 1:length(scotoma)
    
    d = dir([homedir targetdir '*' mystim '*scot' scotoma{sc} '*.mat']);
    ct = 0;
    
    
    for cf = length(d):-1:1
        load([homedir targetdir d(cf).name])
        
        if stim.randAng == randang
            disp(d(cf).name)
            ct = ct+1;
            % make stim movie
            % trim display area (so that the stimulus movie only represents the area of interest)
            display.screenAngle = [2 2]*stim.radDeg; % width,height in deg
            % lower the resolution
            stimres = 0.25*1; % deg
            display.resolution = display.screenAngle / stimres;
            display.width =  2 * display.dist * tan( pi * display.screenAngle(1) / (2*180) ) ; % reconvert (function uses width rather than angle)
            s = makeStimulusMovie( display, stim );
            
            if ct == 1
                E = NaN(s.nx,s.ny,length(d));
            end
            
            titstr{ct} = num2str(cf);
            
            % compute efficiency
            nhdr = 15; % try to estimate an hdr with 15 time points
            
            for ix = 1:s.nx
                for iy = 1:s.ny
                    % compute a timecourse of this position (on/off)
                    tc = squeeze(s.frames(ix,iy,:));
                    if sum(tc) > 0
                        %                     return
                        % compute efficiency of stimulation for this position
                        E(ix,iy,ct) = efficiency(tc,nhdr);
                    end
                end
            end
            
            %         figure(10+ct); clf
            %         plot(stim.motDirsRad/pi,'.-')
            %         xlim([1 stim.numTRs])
            
        end
    end
    avgE{sc} = nanmean(E,3);
end

%% plot
figure(1); clf

clim = [.5 1.2];
ct = 1;
for sc = 1:length(scotoma)
    subplot(1,length(scotoma),sc)
    % map efficiency across space
    imagesc(s.x(1,:),s.y(:,1),avgE{sc},clim); colorbar; set(gca,'ydir','normal')
    %         text(-5,9,sprintf('avg eff %.3f\n',nanmean(E{st,rd}(:))))
    title(sprintf('%s\nefficiency:\n(average across runs)\n%.2f','',nanmean(avgE{sc}(:))))
    axis tight square
end
colormap(jet)
% if niter>1
%     mysubplot([2 length(stimList)],2,1:2)
%     plot(me,'.-')
%     legend(stimList,'location','best')
%     xlabel('run')
%     ylabel('efficiency (average across space)')
% end

% 
% 
% %% plot
% figure(3); clf
% for cf = 1:length(E)
%     % map efficiency across space
%     mysubplot([1,length(E)],1,cf); hold on
%     imagesc(s.x(1,:),s.y(:,1),E{cf},clim); colorbar; set(gca,'ydir','normal')
%     title(sprintf('%s \navg eff %.3f',titstr{cf},nanmean(E{cf}(:))))
%     axis off square
% end
% 
% colormap(jet)



