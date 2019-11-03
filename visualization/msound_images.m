function msound_images(RF, axial, mgrid)
% MSOUND_IMAGES( RF, axial, mgrid )
% 
% Generate images of the raw RF data and a log-compressed B-mode
% based on the beamformed output of an mSOUND simulation.

% 2019-10-28 - Keita Yokoyama (UNC/NCSU)
%              initial version

    nD=msound_nDim(mgrid);
    switch nD
        case 3, RF=squeeze(RF(:,:,round(size(RF,3)/2)));
    end

    % envelope-detect and log-compress RF data
    bimg=20*log10(abs(hilbert(RF)));
    
    % convert mgrid label vectors into more useful scales
    x=1000*mgrid.x;  y=1000*axial; %1e6*mgrid.t;
    
    % generate RF figure
    figure
    subplot(1,2,1)
        imagesc(x, y, RF); title('RF')
        axis image; colormap(gca,'gray'); xlabel('Lateral (mm)'); ylabel('Axial (mm)')
        caxis([-1e4 1e4])

    % generate B-mode
    subplot(1,2,2)
        imagesc(x, y, bimg); title('B-mode')
        axis image; colormap(gca,'gray'); xlabel('Lateral (mm)'); ylabel('Axial (mm)')
        
end