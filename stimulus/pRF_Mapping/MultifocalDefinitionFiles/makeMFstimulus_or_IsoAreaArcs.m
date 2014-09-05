% define arcs in the visual field that cover equal areas in V1
% set a desired area to be covered along the cortical surface, a number of
% radii and the eccentricity of the first on
% if 
% createMFstim = 1 -> load/save files required to run MultiFocal stimulus for
% pRF mapping (with parameters specific for the SCOTOMA experiment)
% 
% 2012/03/23, PB

cls
pwd
% set CMF main parameter (fovea expansion)
p.a = .5; % based on Duncan and Boynton -- Geoff gave me this value, not in the paper

createMFstim = 1; % 1 REQUIRES previously saved files and CREATES files for use in pRF_mapping.m

if ~createMFstim
    
    % pick a number for the cortical area that needs to be covered by each
    % arc (good values are around .1 .. .5, they need to be scaled by a
    % factor of about 20 because we are not modeling the full CMF equation
    % here--only reasoning in relative terms)
    desiredA = .5;
    
    % set the number of radii we want to define
    nRads = 6;
    
    % chose the eccentricity of the most foveal arc (inner radius)
    radStartValue = 0;
    
else
    
    % FOR SCOTOMA EXPERIMENT
    % set scaling factor for compatibility with the first run of this
    % experiment (see function scaling_for_scotomaexp for details)
    
    display.dist = 68;  % distance from screen (cm)
    display.width = 33; % screen width (cm)
    correctionFactor = 1; %scaling_for_scotomaexp( display ) ;

    % load the parameters of the arcs defined by Aguirre/Brainard, compute the
    % area and use it to set the desired areas of our own arcs 
    regDef = importdata('regiondefintionfile.txt'); regDef = regDef.data; % definition of regions
    stimSeq = importdata('48regionStimSeq.txt'); % definition of active regions per frame (not used for computations, only saved together with stimulus layout for convenience)
    
%     % DEBUG
%     regDef = unique(sort(regDef),'rows'); 
%     [~,ind] = unique(regDef(:,3));
%     stimSeq = ones(size(regDef,1),4); stimSeq(1:2:end,:) = 0;
%     % DEBUG
    
    % loop through all arcs and compute area
    desiredA = zeros(size(regDef,1),1);
    
    for i =1:size(regDef,1)
        
        % define arc parameters
        p.startRadDeg = regDef(i,3); % inner radius in deg
        p.widthRadDeg = regDef(i,4); % radius change in deg
        p.startAngRad  = 0; % starting angle in rad
        p.widthAngRad = pi/4; % angle change in rad
        p.shutup = 1;
        
        % find the area of the arc 
        [~,desiredA(i)] = fitArea(p);
        
    end
    % set desired area to their maximum (I picked this index instead of
    % mean or min)
    desiredA = max(desiredA); % about 0.16
    
    % set the number of radii we want to define
    nRads = length(unique(regDef(:,3)));  % 6
    
    % chose the eccentricity of the most foveal arc (inner radius)
    radStartValue = .99;% adjusted after repeatedly running this code to find a value that fits 2 arcs within the first 2 deg eccentricity
    
    % target scotoma radius
    scotoma = 2
    
end

% initialize the parameters of an arc
p.startAngRad  = 0; % starting angle in rad
p.widthAngRad = pi/4; % angle change in rad
p.startRadDeg = 0; % inner radius in deg -> will be updated
p.widthRadDeg = 0; % radius change in deg -> wil be updated
p.shutup = 1;

% iteratively define radii for all nRads arcs
radList = radStartValue; % radList holds the first inner radius now
widthList = []; % we'll find its width, which will also define the eccentricity of the second inner radius; we'll do the same with the second and the third and so on

for i = 1:nRads
    
    p.startRadDeg = radList(end); % last we've computed
    bestP{i} = fit('fitArea',p,{'widthRadDeg'},desiredA); % find the width that results in the arc covering the desired area
    
    widthList = [widthList, bestP{i}.widthRadDeg]; % save width for this arc
    if i <= (nRads-1)
        radList = [radList, p.startRadDeg + bestP{i}.widthRadDeg]; % define inner radius of next arc
    end
    
    % [~,bestA] = fitArea(bestP{i}); bestA % confirm the area is close to desired
end

disp('InnerRad OuterRad')
disp([radList' (radList+widthList)'])
myecc = [radList' widthList'];

% save values & make stimulus layout/sequence matrix for MultiFocal
% stimuli experiment
if createMFstim 
    
    ecc = unique(regDef(:,3));
    for i = 1:length(ecc)
        regDef(regDef(:,3)==ecc(i),4) = myecc(i,2);
        regDef(regDef(:,3)==ecc(i),3) = myecc(i,1);
    end
    save MFStimLayoutAndSequence.mat regDef stimSeq correctionFactor 
end

%% check that the arcs indeed cover equal areas on the cortex
colors = 'gcbymkr';
figure(11); clf; 

% loop through arcs
for i = 1:length(bestP)
    
    p = bestP{i};
    
    % define patch in visual space
    angList = linspace(0,p.widthAngRad,21);
    radList = linspace(0,p.widthRadDeg,5);
        
    % 1: isoOr line, centrifugal direction, starting Angle
    % 2: isoEcc line, CCW direction, Outer radius
    % 3: isoOr line, centripetal direction, starting Angle + Angle width
    % 4: isoEcc line, CW direction, Inner radius
    z1 = (p.startRadDeg+radList)*exp(sqrt(-1)*p.startAngRad);
    z2 = (p.startRadDeg+p.widthRadDeg)*exp(sqrt(-1)*angList);
    z3 = (p.startRadDeg+fliplr(radList))*exp(sqrt(-1)*(p.startAngRad+p.widthAngRad));
    z4 = p.startRadDeg*exp(sqrt(-1)*fliplr(angList));
    
    z = [z1,z2,z3,z4];
    
    % plot it
    subplot(1,2,1); hold on
    patch(real(z),imag(z),colors(1+mod(i,length(colors))))
   
    % transform into cortical space
    w = log(z + p.a);
    
    % plot it
    subplot(1,2,2); hold on
    patch(real(w),imag(w),colors(1+mod(i,length(colors))))
    
end

% now plot isoEccentriciy and isoOrientation lines in both spaces (for
% reference)
radList = myecc(:,1)+myecc(:,2);
angList = pi * ( -.5 : .25 : .5 ) ;

isoOr = linspace( 0 , max(radList) , 201 )' * exp( angList * sqrt(-1) );
isoEcc = radList * exp( linspace( min( angList ) , max( angList ) , 201 ) * sqrt(-1) ); isoEcc = isoEcc';

w_isoOr = log( isoOr + p.a );
w_isoEcc = log (isoEcc + p.a );

subplot(1,2,1)
plot( real( isoOr ), imag( isoOr ), 'k' )
plot( real( isoEcc ), imag( isoEcc ), 'k' )
axis equal
xlabel('deg');
ylabel('deg');
title('Visual space');

subplot(1,2,2)
plot( real( w_isoOr ), imag( w_isoOr ), 'k' )
plot( real( w_isoEcc ), imag( w_isoEcc ), 'k' )
axis equal
xlabel('mm');
ylabel('mm');
title('V1');

