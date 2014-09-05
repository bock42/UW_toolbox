function img = makeBarImage( display, stim, motDirRad, center )
% img = makeBarImage( display, stim [,ang] [,center] )
%
% DESCRIPTION
% returns the basic image of the stimulus to be used
% INPUT
%   display
%       .screenAngle - resolution of screen in degrees
%       .resolution  - resolution of screen in pixels
%   stim
%       .widthDeg       - width of bar stim in degrees
%       [.cyclePerDeg]  - number of black and white chex per degree, ELSE white bar
%   [motDirRad]      - direction of motion in radians (default theta = 0, rightward motion)
%   [center]            - bar center - OTHERWISE generate standard image
% OUTPUT
%   img - stimulus image [display.resolution] ranging between -1 and 1

% VERSION HISTORY
%   1.0 (10.27.2011, ZRE)
%       - created
%   1.1 (02.22.2012, ZRE)
%       - made checkerboard pattern optional
%       - revised comments
%   1.2 (02.24.2012, ZRE)
%       - changed ang -> motDirRad to reflect that stim.angRad is the
%       direction of motion NOT the angle of the bar. To adjust bar angle
%       we must add pi/2 to the motion angle!
%       - removed negative from sin component (thought to fix the fact that
%       positive is down along the y-axis but really just confusing things)
% 03/22/2012 PB
% changed definition of x and y -> see XYdefinition_explain.m 

if ~exist( 'motDirRad', 'var' )
    motDirRad = 0; % default = rightward motion 
end
if ~exist( 'center', 'var' )
    center = [0 0];
end

% rotate bar angle 90 deg so that it is perpendicular to motion direction
barOrienationRad = motDirRad + pi/2;

[xpx ypx] = meshgrid(   ( 1 : display.resolution(1) ) - display.resolution(1)/2 ,...
                    ( 1 : display.resolution(2) )  - display.resolution(2)/2 );
x = pix2angle( display, xpx );
y = pix2angle( display, ypx ) ;

a = cos(barOrienationRad)*(x-center(1)) + sin(barOrienationRad)*(y-center(2));
b = sin(barOrienationRad)*(x-center(1)) + cos(barOrienationRad)*(y-center(2));

mask = abs(b) < stim.widthDeg/2;

% add checkerboard pattern
if isfield( stim, 'cyclePerDeg' )
    % checkerboard bar
    chex = sign( cos(2*pi*a*stim.cyclePerDeg) .* cos(2*pi*b*stim.cyclePerDeg) );
else
    % solid white bar
    chex = ones( size( mask ) );
end
img = chex .* mask;