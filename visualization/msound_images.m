function msound_images(RF, mgrid)
% MSOUND_IMAGES( RF, mgrid )
% 
% Using an mSOUND set_grid object "mgrid" and RF data as processed using
% "msound_beamform", generate the raw RF data and envelope-detected,
% log-compressed B-mode image.

% 2019-10-28 - Keita Yokoyama (UNC/NCSU)
%              initial version

    % envelope-detect and log-compress RF data
    bimg=20*log10(abs(hilbert(RF)));
    
    % convert mgrid label vectors into more useful scales
    x=1000*mgrid.x;  y=1e6*mgrid.t;
    
    % generate RF figure
    figure
    subplot(1,2,1)
        imagesc(x, y, RF)
        title('RF')
        axis image; colormap(gca,'gray'); xlabel('Lateral (mm)');
        caxis([-1e4 1e4])

    % generate B-mode
    subplot(1,2,2)
        imagesc(x, y, bimg)
        title('B-mode')
        axis image; colormap(gca,'gray'); xlabel('Lateral (mm)');
        
end