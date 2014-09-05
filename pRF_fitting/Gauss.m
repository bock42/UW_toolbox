function out = Gauss(RF,x,y,ind,unwrap)

% out = Gauss(RF,x,y,ind,unwrap)
% ind is the index to the RF structure

out = exp( -( (x-RF.center(1)).^2 + (y-RF.center(2)).^2 ) / (2*RF.sig(ind)^2) );
out = (out / sum(out(:))) * 100; % area constant
% out = (out / max(out(:))) * 1; % peak constant

% unwrap 2D into 1D vector
if exist('unwrap','var')
    if unwrap
        out = out(:);
    end
end

