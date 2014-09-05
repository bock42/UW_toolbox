function UnwrappedConvStim = convUnwrap(s,HDR,stim)

% convolve stimulus timecourse with HDR, collapse the 2 spatial dimensions into one
% inputs:
% s.frames [x,y,time]
% s.nx s.ny (size of spatial dimensions)
% stim.numTRs (size of time)
% HDR.function [1,1,time]
%
% output: UnwrappedConvStim [time,space,runs]
% PB 03/2013

UnwrappedConvStim = [];
for nR = 1:size(s.frames,4)
    % convolve the [spacex,spacey,time] by the HDR [1,1,time]
    ConvStim = stim.secPerTR*convn(s.frames(:,:,:,nR),HDR.function); 
    % unwrap 2D space into 1D -> unwrapped [time,space] (if space is 1D, this is simply a shiftdim)
    unwrapped = s.dx*s.dy*reshape(ConvStim(:,:,1:stim.numTRs),[s.nx*s.ny,stim.numTRs])';
    % concatenate across runs -> UnwrappedConvStim [time,space,runs]
    UnwrappedConvStim = cat( 3,UnwrappedConvStim,unwrapped);
end
