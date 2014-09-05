% stimulus to map pRFs
% Multifocal, DriftingBar, RandomBar, Meridians
% modified, from Zach's code
% task: Simon

% Multifocal: (Vanni et al., NeuroImage 2006) current
%   parameters courtesy of Omar Butt in Aguirre lab. 48 arcs
%   (6 rings, 8 wedges), a variable combination of which is presented every
%   second frame (TR). TR must be set to 3. Stimulus sequence is defined by
%   a file OB sent to GMB; it was originally created with m-sequences.
%   Ring eccentricities were modified to cover a smaller area.
%   IsoAreaArcs_createMFstimulus.m makes the necessary files and
%   illustrates the reasoning.
% 2013/03/21 changed output file (date first so that they are sorted by timestamp)
% 2013/04/05 changed TR for multifocal 
% 2013/04/15 added meridians

%PsychJavaTrouble
clear all
clear mex
KbName('UnifyKeyNames')
commandwindow

%% main params

% stimulus type
stim.StimulusType = 'DriftingBar'; %'MultiFocal'; %'Meridians'; % 'RandomBar'; %
% subject (names the directory where all runs are saved):
stim.scotoma.flag = 0; % comment this line for scotoma experiment
subj = 'foo'; % if *test* only 2 TRs
homedir = cd;
% save the variables even if the script didn't run until the end
always_save = 1;
display.gammaFileNm = 'luminanceGammaPowerFit_UWscanner.mat'; % this file must be in the path

%% screen variables

display.dist = 68; % distance from screen (cm)
display.width = 33; % width of screen (cm)
display.skipChecks = 2;
display.bkColor_PreGamma = [128 128 128];
display.screenNum = max(Screen('Screens'));

% Load luminance gamma table
if exist(display.gammaFileNm,'file')
    load( display.gammaFileNm, 'Linv' );
    display.gamma = Linv(:,4)-1; % pull out achromatic
else
    display.gamma = 0:255;
    fprintf( 'Warning! calibration file: %s\n not found\nusing linear look-up table instead\n', display.gammaFileNm )
end
Linearize = @(RGB) display.gamma( RGB+1 ); % pass values through gamma lookup table
rescale2RGB = @(img, pContrast) round( ( pContrast*img + 1)*127.5 ); % rescale img (-1 to 1) to RGB (0 to 255)

%% keyboard variables
a = cd;
if a(1)=='/' % mac or linux
    a = PsychHID('Devices');
    for i = 1:length(a), d(i) = strcmp(a(i).usageName, 'Keyboard'); end
    keybs = find(d);
else % windows
    keybs = [];
end


%% set seed manually
if strcmp(stim.StimulusType,'RandomBar') || strcmp(stim.StimulusType,'DriftingBar')
    stim.seed = input('seed for random: (Zach"s preferred: 8194)');
else
    stim.seed = NaN;
end
    

%% artificial scotoma
if ~isempty(strfind(subj,'test')) || ~isempty(strfind(subj,'tmp'))
    stim.scotoma.flag = 0;
end
if ~isfield(stim, 'scotoma')
    if ~isfield(stim.scotoma, 'flag')
        stim.scotoma.flag = (menu('Scotoma?','Yes','No') == 1);
    end
end
stim.scotoma.rad = [1 1] * 2 ;
stim.scotoma.center = [0 0];
try
    %% -- open window
    
    display.bkColor = Linearize( display.bkColor_PreGamma );
    display = OpenWindow(display);
    Screen( 'BlendFunction', display.windowPtr, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA' );
    rect = Screen('Rect', display.windowPtr );
    centerRect =  repmat( display.center, 1, 2);
    display.screenAngle = pix2angle( display, display.resolution );
    
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
    stim.durSec = 384; %
    stim.secPerTR = 2.4;
    stim.numTRs = stim.durSec / stim.secPerTR;
    
    % others
    stim.fixMaskRadDeg = 0.25; % mask around fixation (0.25 is taken up by the simon display)
    %stim.seed = 8194; % always the same 'random' bar orient (Zach: 8194 best seed out of 1:1e5)
    
    %% stimulus specific params
    switch stim.StimulusType
        
        case 'Meridians'
            
            % stimulus aperture
            blockdur = 5; 
            onecycle = [-ones(blockdur,1); zeros(blockdur,1); ones(blockdur,1); zeros(blockdur,1)];
            fulltc = repmat(onecycle,8,1); 
            stim.numTRs = length(fulltc); % 15*4*6 = 360
            stim.durSec = stim.numTRs * stim.secPerTR ;
            
            stim.meridianWidthRad = pi/6;
            
            stim.meridianPhaseRad = ones(stim.numTRs,1) * pi/4; 
            stim.meridianPhaseRad = stim.meridianPhaseRad + (pi/4 * fulltc) ./ abs(fulltc); % so that alternates between 0 and pi/2 (Inf for blanks)
            
            stim.centersDeg = zeros(stim.numTRs,2);
            stim.centersDeg(fulltc == 0) = 999; % so that blanks -> out of the screen == not shown
            
            % carrier
            stim.rotspeedDegPerSec = 30; % rotation speed of the carrier checkerboard
            
            % circular mask
            stim.radDeg = min(display.screenAngle)/2;

        case 'DriftingBar'
            
            % make bar seq
            stim = getBarSequence( stim );

        case 'RandomBar'
            
            % make bar seq
            stim = getBarSequence( stim );
            
            % shuffle bar positions (but preserve blanks)
            tbscrambled = find(stim.degPerSec); % blanks are marked as velocity = 0
            if isfield( stim, 'seed' )
                rand('seed',stim.seed);
            end
            neworder = Shuffle(tbscrambled);
            stim.centersDeg(tbscrambled,:) = stim.centersDeg(neworder,:);
            stim.motDirsRad(tbscrambled) = stim.motDirsRad(neworder);
            stim.degPerSec(tbscrambled) = stim.degPerSec(neworder);
            
            
        case 'MultiFocal'
            
            % carrier
            stim.rotspeedDegPerSec = 30; % rotation speed of the carrier checkerboard
            
            % stimulus
            stim.widthDeg = stim.radDeg*4; % bar width -> must cover whole screen
            load('MultifocalDefinitionFiles/MFStimLayoutAndSequence.mat', 'regDef', 'stimSeq', 'myTR')
            myTR = 2.4;
            stimSeq=stimSeq(:,1:160);
            if ~isempty(strfind(subj,'test'))
                stimSeq = stimSeq(:,1:2);
            end
            % time
            stim.activeRegions = stimSeq; % pre-set stimulus sequence
            stim.secPerTR = myTR;%3;
            stim.numTRs = size( stim.activeRegions,2 ) ;
            stim.durSec = stim.numTRs * stim.secPerTR ;
            % reset these as they are not used
            stim.durBlanksSec = NaN; % durantion of blanks between bar sweeps
            stim.nSweepsBetweenBlanks = NaN; % number of sweeps separating two blank periods
            stim.begendBlankSec = NaN; % duration of blanks at beginning and end of run

            % space
            stim.wedgeStartDeg = regDef(:,1);
            stim.wedgeSizeDeg = regDef(:,2);
            stim.ringOutDeg = regDef(:,3) + regDef(:,4); % radius of the outer circle where the arc is inscribed
            stim.ringSizeDeg = regDef(:,4); % arc width
            stim.radDeg = max( stim.ringOutDeg ); % (outer bound of) outermost of stimulated arcs
            stim.fixMaskRadDeg = min( stim.ringOutDeg - stim.ringSizeDeg ); % (inner bound of) innermost of stimulated arcs
            % add one more arc to delimit the fixation mask
            stim.wedgeStartDeg = [ stim.wedgeStartDeg ; 0 ] ;
            stim.wedgeSizeDeg = [ stim.wedgeSizeDeg ; 360 ] ;
            stim.ringOutDeg = [ stim.ringOutDeg ; regDef(1,3) ] ;
            stim.ringSizeDeg = [ stim.ringSizeDeg ; regDef(1,4) ] ;
            
            % modify a bar parameters and do deg2pix conversions
            % keep the bar fixed in the same pos and orient
            stim.motDirsRad = zeros(stim.numTRs,1);
            stim.centersDeg = zeros(stim.numTRs,2);
            stim.degPerSec = zeros(stim.numTRs,1);
            
            % define the rectangles where each arc is inscribed & and its
            % width
            ArcsRectPix = angle2pix( display, [ -repmat(stim.ringOutDeg,1,2)  +repmat(stim.ringOutDeg,1,2) ] ) + ...
                + repmat( display.center,size(stim.ringOutDeg,1),2 ) ; % [left, top, right, bottom]
            ArcsPenWidthPix = angle2pix( display,stim.ringSizeDeg );
            
            
        otherwise
            error('unknown stim.StimulusType')
    end
    if ~isempty(strfind(subj,'test'))
        stim.durSec = 2;
    end
    
    stim.nFrames = ceil( stim.durSec * display.frameRate );
    if any( 2*stim.radDeg > display.screenAngle )
        error( 'resize stimulus -- you are outside the bounds of your monitor' )
        Screen('CloseAll'); % ListenChar(1);
    end
        
    
    % carrier
    ndrfitpos = 10;
    stim_tmp = stim; stim_tmp.widthDeg = 100;%max(display.screenAngle)/2;
    
    if strcmp(stim.StimulusType,'Meridians')
        ph = 1;
        % 360 deg wheel
        carrier = makeMeridianImage( display, pi );
        img1 = rescale2RGB( carrier, 1 ); % rescale img to RGB vals (100% contrast)
        img2 = rescale2RGB( -carrier, 1 ); % phase-reversed
        Texture(1,ph) = Screen('MakeTexture', display.windowPtr, Linearize( img1 ) );
        Texture(2,ph) = Screen('MakeTexture', display.windowPtr, Linearize( img2 ) );
    else
        for ph = 1:ndrfitpos
            % bar covering the whole screen
            carrier = makeBarImage( display, stim_tmp, 0, [-ph/(ndrfitpos*stim.cyclePerDeg) 0] );
%             thisign = sign((~mod(ph,2))-0.5);
            img1 = rescale2RGB( carrier, 1 ); % rescale img to RGB vals (100% contrast)
            img2 = rescale2RGB( -carrier, 1 ); % phase-reversed
            Texture(1,ph) = Screen('MakeTexture', display.windowPtr, Linearize( img1 ) );
            Texture(2,ph) = Screen('MakeTexture', display.windowPtr, Linearize( img2 ) );
        end
    end
    clear stim_tmp
    
    % bar aperture
    if strcmp(stim.StimulusType,'Meridians')
        StimMask = abs(makeMeridianImage( display, stim.meridianWidthRad ));
    else
        StimMask = abs(makeBarImage( display, stim ));
    end
    stimImg = repmat( 128*ones( size(StimMask) ), [1 1 4] );
    stimImg(:,:,4) = 255*double(~StimMask); %alpha-channel - invert to make 1's into 0's (0=100% opaque)
    TextureStimMask = Screen('MakeTexture', display.windowPtr, Linearize( stimImg ) );
    
    % circular aperture
    Mask = makeMaskImage( display, stim );
    if stim.scotoma.flag 
        scotomaImg = makeScotomaImage( display, stim );
        Mask = ~Mask; % 0 = non stimulated (i.e. Masked)
        Mask = Mask .* scotomaImg; 
        Mask = ~Mask; % switch back -> 1 = non stimulated (i.e. Masked)
    end
    maskImg = repmat( 128*ones( size(Mask) ), [1 1 4] );
    maskImg(:,:,4) = 255*double(Mask); %alpha-channel - invert to make 1's into 0's (0=100% opaque)
    TextureMask = Screen('MakeTexture', display.windowPtr, Linearize( maskImg ) );
    
    
    % reformat stim.centersDeg in rect space
    stimCentersPix = angle2pix( display, stim.centersDeg );
    destRect = repmat( rect, stim.numTRs, 1 ) + [stimCentersPix stimCentersPix ] ;

    % deal with Simon task
    rand( 'seed', GetSecs ); % randomize seed for task
    task = doSimon(display);
    task.type = 'simon'; % task gets cleared inside of simon
    task.ISI = 1/3;  %seconds
    task.dur = .25;   %seconds
    task.pauseDur = .25; %seconds
    task.errDur = 1;  %seconds
    task.errFreq = 4;  %Hz
    task.keys = {'b','y','g','r'};%{'w','s','a','q'};
    
    % wait for ttl-pulse
    drawText( display, [0 1], '*waiting for ttl-pulse*' );
    Screen('Flip',display.windowPtr);
    wait4T(keybs);  %wait for 't' from scanner.

    SetMouse(0,0)

    %% -- play movie
    
    breakIt = 0;
    frameCnt = 1;
    timeStamp = 0;
    display.timeStamp = zeros( 1, stim.nFrames );
    startTime = GetSecs;  %read the clock
    while GetSecs-startTime < stim.durSec && ~breakIt  %loop until 'esc' pressed or time runs out
        
        % update timers
        elapsedTime = GetSecs-startTime;   
        
        curTR = ceil( elapsedTime / stim.secPerTR );

        flickState =  round( mod( elapsedTime * stim.flickerHz, 1) ) +1;

        driftState =  ceil( mod( elapsedTime * stim.flickerHz, ndrfitpos) );
        
        %% draw stim and mask
        % notes on DrawTexture:
        % first: DrawTexture rotates clockwise instead of
        % counter-clockwise (like I assume in radians) so I subtract
        % the direction from 2*pi first and everything is good
        % second: rotation is in deg instead of radians

        if strcmp(stim.StimulusType,'MultiFocal')
            % if Multifocal stimuli, draw arcs in all 'inactive' stimulus
            % regions => the underlying flickering checkerboard remains visible only in 'active' regions
            
            barAngDeg = mod(stim.rotspeedDegPerSec * elapsedTime,360);
            driftState = 1;

            % carrier
            Screen( 'DrawTexture', display.windowPtr, Texture(flickState,driftState), [], destRect(curTR,:), barAngDeg ); %
            % circular aperture
            Screen( 'DrawTexture', display.windowPtr, TextureMask )
            
            InActiveRegs = find( stim.activeRegions(:,curTR) == 0 ); % find inactive regions
            for ct = 1:length(InActiveRegs)
                i = InActiveRegs(ct);
                Screen('FrameArc',display.windowPtr,display.bkColor,ArcsRectPix(i,:) ,stim.wedgeStartDeg(i,:),...
                    stim.wedgeSizeDeg(i,:),ArcsPenWidthPix(i,:)); 
            end
            
            % add white lines to separate rings
            [bla,ind] = unique( stim.ringOutDeg );
            for ct = 1:length( ind)
                i = ind(ct);
                Screen('FrameArc',display.windowPtr, [1 1 1]*255, ArcsRectPix(i,:)+[-1 -1 +1 +1]*1, 0, 360, 2); % 2px width
            end
            
        elseif strcmp(stim.StimulusType,'Meridians')
            
            % rotating carrier
            barAngDeg = 0; %mod(stim.rotspeedDegPerSec * elapsedTime,360);
            driftState = 1;
            % orientation of meridian
            meridianAngDeg = stim.meridianPhaseRad(curTR) * 180/pi;
            
            % carrier
            Screen( 'DrawTexture', display.windowPtr, Texture(flickState,driftState), [], destRect(curTR,:), barAngDeg ); %
            % meridian aperture
            Screen( 'DrawTexture', display.windowPtr, TextureStimMask , [], destRect(curTR,:), meridianAngDeg )
            % circular aperture
            Screen( 'DrawTexture', display.windowPtr, TextureMask )
            
        else
            
            % Bar stimulus
            barAngDeg = (2*pi-stim.motDirsRad(curTR)) * 180/pi;
            
            % carrier
            Screen( 'DrawTexture', display.windowPtr, Texture(flickState,driftState), [], destRect(curTR,:), barAngDeg ); %
            % bar aperture
            Screen( 'DrawTexture', display.windowPtr, TextureStimMask , [], destRect(curTR,:), barAngDeg )
            % circular aperture
            Screen( 'DrawTexture', display.windowPtr, TextureMask )

        end
        
        % update simon
        task = doSimon(display, task, elapsedTime,keybs);
        Screen('Flip', display.windowPtr);
        display.timeStamp(frameCnt) = GetSecs-startTime;
        
        % check to see if the "esc" button was pressed
        breakIt = escPressed(keybs);
        frameCnt = frameCnt + 1;
    end
    frameCnt = frameCnt-1;
    % crop end off of timeStamp in case frames were dropped
    display.timeStamp = display.timeStamp( 1:frameCnt );
    fprintf( 'Elapsed time %.0f sec (expected: %.0f)\n', display.timeStamp(end), stim.durSec )
    fprintf( '# of frames displayed - expected = %.0f\n', frameCnt-stim.nFrames )
    
    % clean up the Simon process
    if strcmp( task.type, 'simon' )
        task.action = 'done';
        task = doSimon(display, task, elapsedTime,keybs);
        Screen('Flip', display.windowPtr);
    end
    
catch ME
    Screen('CloseAll');
    % ListenChar(0);
    rethrow(ME);
end
Screen('CloseAll');
% ListenChar(0);

%% save data
datadir = fullfile( homedir, 'data', subj );
if ~exist(datadir,'dir')
    mkdir(datadir)
end
fileName = sprintf( 'Stimulus_%s_%s_scot%i_seed%i_%s', datestr( now, 'YYYY_mm_dd_HH_MM_SS' ), stim.StimulusType, stim.scotoma.flag, stim.seed, subj );

if ~breakIt  || always_save
    save( fullfile( datadir, fileName ), 'display', 'stim', 'task');
    fprintf( 'DATA SAVED:\n\t%s\n\t%s\n', datadir, fileName )
end

%% plot performance & timing
if ~exist( 'breakIt', 'var' ); breakIt = false; end
if ~breakIt
    % simon results
    figure(1); clf
    set(gcf,'position', [ 1 31 1024 664]);
    s = task;
    plotSimon;
    % timings
    figure(2); clf
    plot( diff( display.timeStamp ), '.-' ); hold on
%     set( gca, 'YTick', 0:1/stim.flickerHz:1 ); grid on
    title( 'flip-to-flip time' )
    ylabel( 'time elapsed since previous frame (ms)' ); xlabel( 'frame' )
    ylim( [0, 1/stim.flickerHz ])
    legend(sprintf('frame duration: %.2f',1/stim.flickerHz))
end


