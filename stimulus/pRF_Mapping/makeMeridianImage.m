function img = makeMeridianImage( display, meridianWidthRad, meridianPhaseRad, center )
% img = makeMeridianImage( display, meridianWidthRad, meridianPhaseRad, center )
%
% DESCRIPTION
% returns the basic image of the stimulus to be used
% INPUT
%   display
%       .screenAngle - resolution of screen in degrees
%       .resolution  - resolution of screen in pixels
%   meridianWidthRad - width of merdian stim in radians
%   [meridianPhaseRad] - orientation in radians , default is 0
%   [center]         - bar center - default is 0
% OUTPUT
%   img - stimulus image [display.resolution] ranging between 0 and 1
% 
% 04/2013 PB

if ~exist( 'center', 'var' )
    center = [0 0];
end

if ~exist( 'meridianPhaseRad', 'var' )
    meridianPhaseRad = 0;
end

% make basis functions
[xpx ypx] = meshgrid(   ( 1 : display.resolution(1) ) - display.resolution(1)/2 ,...
                    ( 1 : display.resolution(2) )  - display.resolution(2)/2 );
x = pix2angle( display, xpx );
y = pix2angle( display, ypx ) ;

[ang,rad] = cart2pol(x,y);
basis = mod(ang-meridianPhaseRad,pi);

mask = basis < (meridianWidthRad/2) | basis > (pi - meridianWidthRad/2);

% add checkerboard pattern
chex = sign( sin(24*ang) .* sin(4*rad) );

img = chex .* mask;
