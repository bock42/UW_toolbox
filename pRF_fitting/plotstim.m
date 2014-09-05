function plotstim(s,flag,pausebetweenframes)

% visualize stimulus movie 's'
% flags -> 
% 0 = do nothing and return;
% 1 = show mean and var of the stimulus; 
% 2 = show stimulus movie; 
% 3 = 1D plot of stim present/absent vs time
% 4 = 1 screenshot
% 5 = efficiency over space

if nargin < 3
    pausebetweenframes = .1;
end

if 0
elseif flag==1 % statistics

    subplot(1,2,1); title('mean'); hold on
    imagesc(s.x(1,:),s.y(:,1),mean(mean(s.frames(:,:,:,:),3),4));
    axis equal %off
    colorbar
    subplot(1,2,2); title('variance'); hold on
    imagesc(s.x(1,:),s.y(:,1),mean(var(s.frames(:,:,:,:),[],3),4));
    axis equal %off
    colorbar
    drawnow
elseif flag==2 %  visualize the stimulus movie

    for fr = 1 : size(s.frames(:,:,:,1),3)
        title( fr )
        hold on
        image( s.x(1,:), s.y(:,1), s.frames(:,:,fr,1) * 255)
        drawnow
        pause(pausebetweenframes)
        axis equal
        
    end
elseif flag==3 %  space average over time

    avg = squeeze(mean(mean(s.frames(:,:,:,1),1),2));
    plot(avg,'.-')
    xlabel('time (fr)')
    ylabel('average stimulus contrast across space')
    ylim([0 .5])
%     sum(avg==0)/length(avg)
elseif flag==4 %  visualize one screenshot of the stimulus movie
    
    fr = 30;
    title( fr )
    hold on
    image( s.x(1,:), s.y(:,1), s.frames(:,:,fr,1) * 255)
    pause(.05)
    axis equal
    
elseif flag==5
    nhdr = 15; % TRs points where hdr is defined
    E = NaN(s.nx,s.ny,size(s.frames,4));
    for r = 1:size(s.frames,4)
    
    for ix = 1:s.nx
        for iy = 1:s.ny
            % compute a timecourse of this position (on/off)
            tc = squeeze(s.frames(ix,iy,:,r));
            if sum(tc) > 0
                % compute efficiency of stimulation for this position
                E(ix,iy,r) = efficiency(tc,nhdr);
            end
        end
    end
    
    end
    
    E = nanmean(E,3);
    
    imagesc(s.x(1,:),s.y(:,1),E); colorbar; set(gca,'ydir','normal')
    title(sprintf('efficiency: %.3f',nanmean(E(:))))
    
end

