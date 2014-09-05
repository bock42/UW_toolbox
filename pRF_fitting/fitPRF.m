function [err PredResp correl mse polyVals] = fitPRF(RF,UnwrappedConvStim,ObsResp,s,opt)
%
%  err = fitPRF(RF,UnwrappedConvStim,ObsResp,s,opt)
% 
% input:
%
% RF -> RF.center, RF.sig
% [ObsResp -> fMRI timecourse] optional
% s -> stimulus movie 
%        s.frames - x,y,time
%        s.x, x.y
% opt.gof -> corr or mse
% output:
% 
% err = error predicted and observed fMRI timecourse (average across runs)
% PredResp = timecourse predicted given RF and stimulus 
% correl,mse = -correlation and mse
% polyVals = amplitude and offset that best fit the observed fMRI timecourse

% stimulus([space,time]) * thisprf([time,1]) = predictedfMRItimecourse([time,1])

ind = 1; % used to be for DoG modeling (discontinued)
ntimes = size(UnwrappedConvStim,1);
nruns = size(UnwrappedConvStim,3);
PredResp = NaN(ntimes,nruns); %[time,nruns]
for runNum = 1:nruns
	PredResp(:,runNum) = UnwrappedConvStim(:,:,runNum)*Gauss(RF,s.x,s.y,ind,1);
end
if ~isempty(ObsResp)
    correl = NaN(1,nruns);
    mse = NaN(1,nruns);
    polyVals = NaN(2,nruns);
    for runNum = 1:nruns
        % compute correlation
        correl(:,runNum) = mycorr(ObsResp(:,runNum),PredResp(:,runNum));
%         % find linear coefficients
%         polyVals(:,runNum) = polyfit(PredResp(:,runNum),ObsResp(:,runNum),1);
%         % scale and shift predresp
%         PredResp(:,runNum) = polyval(polyVals(:,runNum),PredResp(:,runNum));
%         % compute mse
%         mse(:,runNum) = mean( (ObsResp(:,runNum) - PredResp(:,runNum)).^2 );
    end
    % average across runs
    correl = -mean(correl) ; 
    mse = mean(mse);
    polyVals = mean(polyVals,2)';
    eval(['err = ' opt.gof ';']) % 'gof' is either correl or mse
else
    correl = NaN;
    mse = NaN;
    polyVals = [NaN NaN];
    err = NaN;
end
return
%%
figure(100); clf
plot(PredResp,'r')
hold on
plot(ObsResp)
