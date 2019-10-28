function xdcr=msound_xdcr(mgrid, fc, focus)
% xdcr = MSOUND_XDCR( mgrid, fc, focus )
% 
% Using an mSOUND set_grid object "mgrid", and transducer center frequency
% "fc", create a structure that encodes information about a planar
% transducer at a scalar focal depth "focus".

% 2019-10-26 - Keita Yokoyama (UNC/NCSU)
%              initial version
    
    % characterize center frequency + focal depth + aperture radius
    xdcr.fc=fc;
    xdcr.focus=focus;
    xdcr.dim=[10, 10]*1e-3; % transducer size (lateral, elevation)
    xdcr.plane_depth=0.005; % depth of "transducer" plane
    
    % initialize transducer mask
    xdcr.mask=zeros(mgrid.num_x, mgrid.num_y+1);
    
    % create new axial vector with respect to transducer position
    depth=mgrid.y' - mgrid.y(1) - xdcr.plane_depth;
    
    % position "transducer" plane
    x_xdcr_loc=knnsearch(mgrid.x', -xdcr.dim(1)/2) : knnsearch(mgrid.x', xdcr.dim(1)/2);
    y_xdcr_loc=knnsearch(depth,0);
    xdcr.mask(x_xdcr_loc,y_xdcr_loc)=1;