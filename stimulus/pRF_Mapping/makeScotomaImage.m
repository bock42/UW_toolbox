function scotomaImg = makeScotomaImage( display, stim )
% mask = makeMaskImage( display, stim )
%
% DESCRIPTION
% Creates logical mask which will be used to make a texture and mask off
% scotoma area.
% INPUT
%   display
%       .screenAngle - resolution of screen in degrees
%   	.resolution  - resolution of screen in pixels
%   stim
%       .scotoma
%           .flag   - [bool] toggle scotoma on/off
%           .center - [x y] center in degrees
%           .rad    - radius in degrees
%   
% OUTPUT
%   scotomaImg - double size = resolution

% VERSION HISTORY
%   1.0 (02.22.2012, ZRE)
%       - created as a function (previously internal code)
%       - created a field, stim.scotoma.flag to toggle scotoma on/off
%       - [host functions should also check for the field stim.scotoma as a
%       second alternative method for toggling the scotoma on/off]
%       - WARNING: this function is untested
% 03/22/2012 PB
% changed definition of x and y -> see XYdefinition_explain.m 

% added for backwards compatability
if ~isfield( stim.scotoma, 'flag' )
    stim.scotoma.flag = true;
end

scotomaImg = ones(display.resolution(2),display.resolution(1));
if stim.scotoma.flag

    [xpx ypx] = meshgrid(   ( 1 : display.resolution(1) ) - display.resolution(1)/2 ,...
        ( 1 : display.resolution(2) )  - display.resolution(2)/2 );
    x = pix2angle( display, xpx );
    y = pix2angle( display, ypx ) ;

    r = sqrt((x-stim.scotoma.center(1)).^2 + (y-stim.scotoma.center(2)).^2);
    scotomaImg(r<stim.scotoma.rad(1)) = 0;
    id = r>stim.scotoma.rad(1) & r<stim.scotoma.rad(2);
    scotomaImg(id) = .5*(1-cos( pi*(r(id)-stim.scotoma.rad(1))/1.5));
end

% figure(99); clf; imagesc(scotomaImg); colormap(gray)

