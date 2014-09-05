cls

% load stimulus file
% homedir = cd;
% targetdir = '/data/bla/';
homedir = '/media/VERBATIM/fMRI/FreeSurfer/VisionLesion/';
targetdir = '2013_04_02_SCOTOMA_H32Coil_AB/stimulus/';
d = dir([homedir targetdir '*Ra*Bar*scot1*.mat']); 
stimulus_file = [homedir targetdir d(end).name];
disp(d(end).name)
load(stimulus_file)


%% compute the stimulus movie
% trim display area (so that the stimulus movie only represents the area of interest)
display.screenAngle = [2 2]*stim.radDeg; % width,height in deg
% lower the resolution
stimres = .5; % deg
display.resolution = display.screenAngle / stimres;
display.width =  2 * display.dist * tan( pi * display.screenAngle(1) / (2*180) ) ; % reconvert (function uses width rather than angle)
s = makeStimulusMovie( display, stim );


% %% plot bar orientation
% figure(10); clf
% plot(stim.motDirsRad/pi,'.-')
% xlim([1 stim.numTRs])
% 
% % return


figure(11); clf

% plotstim = visualize stimulus movie 's'
% flags -> 
% 0 = do nothing and return;
% 1 = show mean and var of the stimulus; 
% 2 = show stimulus movie; 
% 3 = 1D plot of stim present/absent vs time
% 4 = 1 screenshot
% 5 = efficiency over space

flag = 5;
plotstim(s,flag,.1);
