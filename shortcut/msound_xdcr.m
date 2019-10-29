function xdcr=msound_xdcr(mgrid, fc, varargin)
% xdcr = MSOUND_XDCR( mgrid, fc )
% xdcr = MSOUND_XDCR( mgrid, fc, focus, dim, plane_depth )
% 
% Create a structure that encodes information about a planar transducer
% for an mSOUND simulation.
%
% REQUIRED INPUT:
%          mgrid = mSOUND set_grid object
%             fc = center frequency of transducer
%                  [Hz] (scalar)
%
% OPTIONAL INPUT:
%          focus = focal depth
%                  [m] (scalar)
%                  DEFAULT: where axial depth = 0 in "mgrid"
%            dim = transducer size
%                  [m] (1x2 vector - [lat,ele])
%                  DEFAULT: 10mm x 10mm
%    plane_depth = transducer position, axial direction
%                  [m] (scalar)
%                  DEFAULT: 0.5mm from top edge of simulation grid

% 2019-10-26 - Keita Yokoyama (UNC/NCSU)
%              initial version
% 2019-10-29 - Keita Yokoyama (UNC/NCSU)
%              generalized function for 1/2/3-D environments;
%              added more helpful help texts

    % define number of dimensions in simulation
    nD=msound_nDim(mgrid);
    
    % characterize center frequency
    xdcr.fc=fc;
    
    % ditto, aperture radius (lateral x elevational size)...
    xdcr.dim=[10, 10]*1e-3;
    
    % ...transducer plane depth...
    depthRaw=0.0005;
    switch nD
        case 1, xdcr.plane_depth=mgrid.x( knnsearch(mgrid.x', mgrid.x(1)+depthRaw) );
        case 2, xdcr.plane_depth=mgrid.y( knnsearch(mgrid.y', mgrid.y(1)+depthRaw) );
        case 3, xdcr.plane_depth=mgrid.z( knnsearch(mgrid.z', mgrid.z(1)+depthRaw) );
    end
    
    % ...and focal depth...
    xdcr.focus=(mgrid.y_length/2);
    
    % ...but change default values according to optional inputs
    for nOpt=0:length(varargin)
        switch nOpt
            case 1, xdcr.focus=varargin{1};
            case 2, xdcr.dim=varargin{2};
            case 3, xdcr.plane_depth=varargin{3};
        end
    end
    
    switch nD
        case 1
            % initialize transducer mask
            xdcr.mask=zeros(mgrid.num_x+1, 1);
            
            % create new axial vector with respect to transducer position
            depth=mgrid.x' - mgrid.x(1) - xdcr.plane_depth;
            
            % define position of "transducer" point
            axi_xdcr_loc=knnsearch(depth,0);
            
            % mark position of "transducer"
            xdcr.mask(axi_xdcr_loc)=1;
            
        case 2
            % initialize transducer mask
            xdcr.mask=zeros(mgrid.num_x, mgrid.num_y+1);
            
            % create new axial vector with respect to transducer position
            depth=mgrid.y' - mgrid.y(1) - xdcr.plane_depth;
            
            % define position of "transducer" line
            lat_xdcr_loc=knnsearch(mgrid.x',-xdcr.dim(1)/2) : knnsearch(mgrid.x',xdcr.dim(1)/2);
            axi_xdcr_loc=knnsearch(depth,0);
            
            % mark position of "transducer"
            xdcr.mask(lat_xdcr_loc,axi_xdcr_loc)=1;
            
        case 3
            % initialize transducer mask
            xdcr.mask=zeros(mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
            
            % create new axial vector with respect to transducer position
            depth=mgrid.z' - mgrid.z(1) - xdcr.plane_depth;
            
            % define position of "transducer" plane
            lat_xdcr_loc=knnsearch(mgrid.x',-xdcr.dim(1)/2) : knnsearch(mgrid.x',xdcr.dim(1)/2);
            ele_xdcr_loc=knnsearch(mgrid.y',-xdcr.dim(2)/2) : knnsearch(mgrid.y',xdcr.dim(2)/2);
            axi_xdcr_loc=knnsearch(depth,0);
            
            % mark position of "transducer"
            xdcr.mask(lat_xdcr_loc,ele_xdcr_loc,axi_xdcr_loc)=1;
    end
end