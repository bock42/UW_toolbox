%function plotSimon(s)

dur = s.event(end).time;

figure(1)
clf
hold on

plotDur= 60;
nPlots = ceil(dur/plotDur);

%plot([0,dur],[0,0],'k-','LineWidth',2);

count =0;
curPlot = 0;
for i=1:s.numEvents
    newPlot = ceil((s.event(i).time+1e-5)/plotDur);
    if curPlot ~=newPlot
        curPlot = newPlot;
        subplot(nPlots,1,curPlot)
    end
    hold on
    switch(s.event(i).type);
        
        case 'switch to play';
            plot( [s.event(i).time,s.event(i).time],[-1,2.5],'r:','LineWidth',1);
        case 'switch to recall';
            plot( [s.event(i).time,s.event(i).time],[-1,2.5],'b:','LineWidth',1);
        case 'play: on'
            patch([s.event(i).time,s.event(i+1).time,s.event(i+1).time,s.event(i).time],[0,0,1,1],s.color{s.event(i).num}/255,'EdgeColor','none');
        case 'recall: on'
            patch([s.event(i).time,s.event(i+1).time,s.event(i+1).time,s.event(i).time],[0,0,1,1]+1,s.color{s.event(i).num}/255,'EdgeColor','none');
        case 'error: start';
            patch([s.event(i).time,s.event(i+1).time,s.event(i+1).time,s.event(i).time],[-1,-1,2.5,2.5],[.5,.5,.5]);
            count =count+1;
            text(s.event(i).time,2.5,num2str(s.maxLen(count)),'EdgeColor','k','VerticalAlignment','bottom');
    end
end
text(dur,2.5,sprintf('%d+',s.maxLen(end)),'HorizontalAlignment','center','EdgeColor','k','VerticalAlignment','bottom');

for i=1:nPlots
    subplot(nPlots,1,i)
    
    set(gca,'YTick',[-inf,inf]);
    set(gca,'XLim',[ (i-1)*plotDur,i*plotDur]);
    %set(gca,'XLim',[0,s.event(end).time]);
    pos = get(gca,'Position');
    pos(1) = .05; pos(3) = .9;
    set(gca,'Position',pos);
    %set(gca,'Position',[.05,.2,.9,.2]);
    set(gca,'YLim',[-1,2.5]);
    if i == nPlots 
        xlabel('Time (s)');
    end
    set(gca,'Color',[.7,.7,.7]);
end