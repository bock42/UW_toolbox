function stim = getMeridiansSequence( stim )
% stim = getMeridiansSequence( stim );
%
% DESCRIPTION
% Computes orientation of the meridian over time. 
% INPUT
%   stim
%       .numTRs     - duration of sequence in TRs
%   	.widthDeg   - width of bar in degrees (Not degrees of visual angle, e.g. 30 deg)
%       .durBlanksSec - for blank periods: duration in sec
% OUTPUT
%   stim
%       .centersDeg - centers in deg, (0,0) is screen center (-y is up)
%   	.motDirsRad - direction of motion in radians (0 to the right, pi/2 is up)
%   	.degPerSec  - instantaneous speed at every TR

% VERSION HISTORY
%   1.0 (10.27.2011, ZRE)
%       - created
%   1.1 (02.22.2012, ZRE)
%       - the seed for rand() is now set here
%       - changed stim.ang -> stim.angsRad 
%       - revised comments
%   1.2 (02.24.2012, ZRE)
%       - change stim.angsRad -> stim.motDirsRad
%   0604.2012, PB
%       - blank periods added
%   0314.2013, PB
%       - switch to create a regular sequence added


% set vars
if stim.randAng
    % set seed to generate ''random'' sequence
    if isfield( stim, 'seed' )
        rand('seed',stim.seed);
    end
else
    myAngs = (0:1/4:7/4)*pi;
end
maxecc = stim.radDeg; % determines start and end positions of the bar within a sweep
nBlankFrames = round(stim.durBlanksSec/stim.secPerTR); % duration of blanks, in frames

% initialize matrices
stim.centersDeg = 999 * ones(stim.numTRs,2);
stim.motDirsRad = -1 * ones(stim.numTRs,1);
stim.degPerSec = 0 * ones(stim.numTRs,1);

curCenter = 999; dx = 0; dy = 0;
nSweeps = 0;
blankFr = nBlankFrames; % so that the first blank is skipped

for iTR = stim.begendBlankSec:(stim.numTRs-stim.begendBlankSec)
    % move stimulus every TR
    curCenter = curCenter + [dx,dy];
    
    % stimulus center has moved outside of range (this sweep is over || first TR)
    if norm(curCenter) >= maxecc % norm -> sqrt(c(1)^2+c(2)^2)
        
        if ~mod(nSweeps,stim.nSweepsBetweenBlanks) && blankFr < nBlankFrames
            % blank (set to initial val)
            curCenter = stim.centersDeg(iTR,:); 
            theta = stim.motDirsRad(iTR); 
            degPerSec = stim.degPerSec(iTR);
            blankFr = blankFr + 1; 
        else
            
            % pick new angle
            if stim.randAng
                theta = rand(1) * 2*pi; % random direction of motion
            else
                theta = myAngs(mod(nSweeps,length(myAngs))+1); % one of the set angles
            end
            
            % compute new starting bar position & dx,dy
            degPerSec = rand(1) * ( max(stim.degPerSecRange)-min(stim.degPerSecRange) ) + min(stim.degPerSecRange);
            curCenter = -maxecc * [cos(theta), -sin(theta)];
            dx = cos(theta) * degPerSec * stim.secPerTR;
            dy = -sin(theta) * degPerSec * stim.secPerTR;
            
            % stim
            blankFr = 0;
            nSweeps = nSweeps + 1;
        end
    end
    
    % store
    stim.centersDeg(iTR,:) = curCenter;
    stim.motDirsRad(iTR) = theta;
    stim.degPerSec(iTR) = degPerSec;
end

return

figure(10); clf
plot(stim.motDirsRad/pi,'.-')
% xlim([1 stim.numTRs])


sweepdur = maxecc*2/degPerSec
howmanysweeps = stim.numTRs / sweepdur

stim.numTRs / (sweepdur*8)