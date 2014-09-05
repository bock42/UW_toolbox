function [err,A] = fitArea(p,desiredA)

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

% project to cortical space
w = log(z+p.a); % w = p.k * log(z+p.a) would take the appropriate scaling factor p.k into account

% compute area of the polygon
A = polyarea(real(w),imag(w));

% compute squared distance from desired area, if defined
if exist('desiredA','var')
    err = (A-desiredA)^2;
else
    err = NaN;
end


