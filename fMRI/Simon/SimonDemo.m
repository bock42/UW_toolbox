%SimonDemo
%PsychJavaTrouble

clear all
clear MEX
%Load in the luminance gamma table for the scanner
if exist('scanner_luminance_gamma.mat','file')
    load scanner_luminance_gamma
    display.gamma = luminanceGamma;
else
    display.gamma = 0:255;
    disp('Warning! calibration file ''scanner_luminance_gamma.mat'' not found')
end

display.dist = 72;     %distance from screen (cm)
display.width = 27.5;  %width of screen (cm)

display.screenNum = 0;
display.skipChecks = 1;

p.dur = 60;  %duration of scan (seconds)

try
    display = OpenWindow(display);
    p.Hz = display.frameRate/8;
    ListenChar(2);

    %Put up blank image
    blankImg = 128*ones(512);
    blankTex =  makeTexture(display,blankImg);
    Screen('DrawTexture', display.windowPtr,blankTex );

    %set up Simon
    
    
    s = doSimon(display);
    
    s.ISI = 1/3;  %seconds
    s.dur = .25;   %seconds
    s.pauseDur = .25; %seconds
    s.errDur = 1;  %seconds
    s.errFreq = 4;  %Hz
    
    s.keys = {'w','s','a','q'};
    
    
    
    Screen('Flip',display.windowPtr);
    
    wait4T;  %wait for 't' from scanner.
    startTime = GetSecs;  %read the clock

    while GetSecs-startTime<p.dur & ~escPressed  %loop until 'esc' pressed or time runs out
        curTime = GetSecs-startTime;   %calculate the current time elapseed
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Do all of your regular stimulus presentation stuff here (except
        % for 'flip')

        Screen('DrawTexture', display.windowPtr,blankTex );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        s=doSimon(display,s,curTime);
        Screen('Flip',display.windowPtr);

    end
    %clean up the Simon process
    s.action = 'done';
    s= doSimon(display,s,curTime);
    Screen('Flip',display.windowPtr);
catch ME
    Screen('CloseAll');
    ListenChar(0);
    rethrow(ME);
end

Screen('CloseAll');
ListenChar(0);

%Show the results from the Simon task.
plotSimon;




