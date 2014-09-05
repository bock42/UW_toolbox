function mask = makeMaskImage( display, stim )
% mask = makeMaskImage( display, stim )
%
% DESCRIPTION
% Creates logical mask which will be used to make a texture and mask off
% an annulus around the stimulus and a disc around fixation.
% INPUT
%   display
%       .screenAngle - resolution of screen in degrees
%   	.resolution  - resolution of screen in pixels
%   stim
%       .radDeg         - radius of stimulus aperture in degrees
%   	[.fixMaskRadDeg]- radius of fixation mask in degrees
% OUTPUT
%   mask - binary image [display.resolution]

% VERSION HISTORY
%   1.0 (10.27.2011, ZRE)
%       - created
%   1.1 (02.22.2012, ZRE)
%       - allowed the stim field "fixMaskRadDeg" to be optional
%       - revised comments
% 03/22/2012 PB
% changed definition of x and y -> see XYdefinition_explain.m 

% if there is no fixation mask field set to zero
if ~isfield( stim, 'fixMaskRadDeg' )
    stim.fixMaskRadDeg = 0;
end

[xpx ypx] = meshgrid(   ( 1 : display.resolution(1) ) - display.resolution(1)/2 ,...
                    ( 1 : display.resolution(2) )  - display.resolution(2)/2 );
x = pix2angle( display, xpx );
y = pix2angle( display, ypx ) ;

r = sqrt( x.^2 + y.^2 );

mask = r > stim.radDeg | r < stim.fixMaskRadDeg;
% figure(99); clf; imagesc(mask); colormap(gray)
% unique(mask)