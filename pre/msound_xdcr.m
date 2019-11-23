function xdcr=msound_xdcr(mgrid, fc, varargin)
% xdcr = MSOUND_XDCR( mgrid, fc )
% xdcr = MSOUND_XDCR( mgrid, fc, ... )
% 
% Create a structure that encodes information about a transducer
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
%
%            dim = transducer size
%                  [m] (1x2 vector - [lat,ele])
%                  DEFAULT: 10mm x 5mm
%
%    plane_depth = transducer position with respect to the
%                  axial vector in "mgrid"
%                  [m] (scalar)
%                  DEFAULT: top edge of simulation grid
%
%       elemsize = size (width/height) of a single element
%                  [m] (1x2 vector - [lat,ele])
%                  DEFAULT: 0.2mm x 1mm
%
%           kerf = spacing between elements
%                  [m] (1x2 vector - [lat,ele])
%                  DEFAULT: single spatial "period" in each dim.
%
%       dynrange = dynamic range
%                  [dB] (scalar)
%                  DEFAULT: -60dB
%                  TRANSMIT ONLY
%
%        angsens = angular sensitivity of element, in degrees
%                  (scalar)
%                  DEFAULT: 45 degrees
%                  NOT YET ACTIVE

% 2019-10-26 - Keita Yokoyama (UNC/NCSU)
%              initial version
% 2019-10-29 - Keita Yokoyama (UNC/NCSU)
%              generalized function for 1/2/3-D environments;
%              added more helpful help texts
% 2019-11-14 - Keita Yokoyama (UNC/NCSU)
%              added channel maps etc. for multiple transmits +
%              easier beamforming

% define number of dimensions in simulation
    nD=msound_nDim(mgrid);
    
% default parameters for whole transducer
    xdcr.dim=[10, 5]*1e-3;         % [lat x ele] size of smallest transducer cross-section
    xdcr.focus=(mgrid.y_length/2); % focal depth
    xdcr.kerf=[mgrid.dx, mgrid.dy];% [lat x ele] kerf (spacing between elements)
    
% default parameters for each element
    xdcr.fc=fc;                    % center frequency
    xdcr.elemsize=[0.2, 1]*1e-3;   % [lat x ele] element size
    xdcr.dynrange=60;              % dynamic range of element
    xdcr.angsens=45;               % angular sensitivity of element, in degrees
    
% default parameters to be converted
    numLayers=1;                   % axial thickness of pressure recording region
    depthRaw=0;                    % offset depth of transducer from grid edge
    switch nD
        case 1, xdcr.plane_depth=mgrid.x( knnsearch(mgrid.x', mgrid.x(1)+depthRaw) );
        case 2, xdcr.plane_depth=mgrid.y( knnsearch(mgrid.y', mgrid.y(1)+depthRaw) );
        case 3, xdcr.plane_depth=mgrid.z( knnsearch(mgrid.z', mgrid.z(1)+depthRaw) );
    end
    
% add optional inputs, if defined
    for nOpt=0:length(varargin)
        switch nOpt
            case 1, if ~isempty(varargin{1}),       xdcr.focus=varargin{1};  end
            case 2, if ~isempty(varargin{2}),         xdcr.dim=varargin{2};  end
            case 3, if ~isempty(varargin{3}), xdcr.plane_depth=varargin{3};  end
            case 4, if ~isempty(varargin{4}),    xdcr.elemsize=varargin{4};  end
            case 5, if ~isempty(varargin{5}),        xdcr.kerf=varargin{5};  end
            case 6, if ~isempty(varargin{6}),    xdcr.dynrange=varargin{6};  end
            case 7, if ~isempty(varargin{7}),     xdcr.angsens=varargin{7};  end
        end
    end
    
% override manual inputs, if dimensions don't make sense
    switch nD
        case 1
            xdcr.dim=[0, 0];
            xdcr.elemsize=[0, 0];
            xdcr.kerf=[0, 0]; 
        case 2
            xdcr.dim=[xdcr.dim(1), 0];
            xdcr.elemsize=[xdcr.elemsize(1), 0];
            xdcr.kerf=[xdcr.kerf(1), 0];
    end
    
% derive remaining parameters about the input
% (for values that do not depend on the simulation dimension)
    xdcr.pitch=xdcr.kerf+xdcr.elemsize;
    
    if nD==1
            % initialize transducer mask
            xdcr.mask=zeros(mgrid.num_x+1, 1);
            
            % create new axial vector with respect to transducer position
            depth=mgrid.x' - mgrid.x(1) - xdcr.plane_depth;
            
            % define position of "transducer" point
            axi_xdcr_loc=knnsearch(depth,0);
            
            % mark position of "transducer"
            xdcr.mask(axi_xdcr_loc)=1;
            xdcr.thickness=mgrid.dx; % a point is only one unit thick
            
            % create map of transducer
            % (...except there's only one channel to begin with)
            xdcr.chanmap=1;
            xdcr.Nelem=[1, 1];
            xdcr.chanID=1;
            xdcr.chanLoc={0,0};
    else
% identify positions and indices of helpful features
        if nD==2
            % initialize transducer mask
            xdcr.mask=zeros(mgrid.num_x, mgrid.num_y+1);
            
            % create new axial vector with respect to transducer position
            depth=mgrid.y' - mgrid.y(1) - xdcr.plane_depth;
            
            % define thickness of "transducer" layer
            xdcr.thickness=mgrid.dy * (knnsearch(depth,0)+numLayers-1);
            
            % define position of "transducer" line
            lat_xdcr_loc=knnsearch(mgrid.x',-xdcr.dim(1)/2) : ...
                         knnsearch(mgrid.x', xdcr.dim(1)/2);
            axi_xdcr_loc=knnsearch(depth,0) + (1:numLayers);
            
            % mark position of "transducer"
            xdcr.mask(lat_xdcr_loc,axi_xdcr_loc)=1;
            
            % get number of indices for element size/kerf/pitch
            id_eL=xdcr.elemsize(1)/mgrid.dx; % size
            id_kL=xdcr.kerf(1)/mgrid.dx;     % kerf
            id_pL=id_eL+id_kL;               % pitch
            
            % identify grid indices for edges of transducer
            id_edgeL=[knnsearch(mgrid.x',-xdcr.dim(1)/2),...
                      knnsearch(mgrid.x', xdcr.dim(1)/2)];
            id_edgeE=[1, 1];  id_eE=1;  id_pE=1;
            
        else
            % initialize transducer mask
            xdcr.mask=zeros(mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
            
            % create new axial vector with respect to transducer position
            depth=mgrid.z' - mgrid.z(1) - xdcr.plane_depth;
            
            % define thickness of "transducer" layer
            xdcr.thickness=mgrid.dz * (knnsearch(depth,0)+numLayers-1);
            
            % define position of "transducer" plane
            lat_xdcr_loc=knnsearch(mgrid.x',-xdcr.dim(1)/2) :...
                         knnsearch(mgrid.x', xdcr.dim(1)/2);
            ele_xdcr_loc=knnsearch(mgrid.y',-xdcr.dim(2)/2) :...
                         knnsearch(mgrid.y', xdcr.dim(2)/2);
            axi_xdcr_loc=knnsearch(depth,0) + (1:numLayers);
            
            % mark position of "transducer"
            xdcr.mask(lat_xdcr_loc,ele_xdcr_loc,axi_xdcr_loc)=1;
            
            % get number of indices for element size/kerf/pitch
            id_eL=xdcr.elemsize(1)/mgrid.dx;   id_eE=xdcr.elemsize(2)/mgrid.dy; % size
            id_kL=xdcr.kerf(1)/mgrid.dx;       id_kE=xdcr.kerf(2)/mgrid.dy;     % kerf
            id_pL=id_eL+id_kL;                 id_pE=id_eE+id_kE;               % pitch
            
            % identify grid indices for edges of transducer
            id_edgeL=[knnsearch(mgrid.x',-xdcr.dim(1)/2),...
                      knnsearch(mgrid.x', xdcr.dim(1)/2)];
            id_edgeE=[knnsearch(mgrid.y',-xdcr.dim(2)/2),...
                      knnsearch(mgrid.y', xdcr.dim(2)/2)];
        end
        
        % count number of elements in each direction
        xdcr.Nelem=[length( id_edgeL(1):id_pL:id_edgeL(2) ),...
                    length( id_edgeE(1):id_pE:id_edgeE(2) )];

% map elements of the transducer to the "transducer mask" 
% (which is just a spatial range to record the mSOUND simulation)
        chanID=0;
        idxLat=id_edgeL(1):id_pL:id_edgeL(2);
        idxEle=id_edgeE(1):id_pE:id_edgeE(2);

% create map of transducer, looking along axial direction
        if nD==2, xdcr.chanmap=zeros(size( mean(xdcr.mask(:,  axi_xdcr_loc), 2) ));
        else,     xdcr.chanmap=zeros(size( mean(xdcr.mask(:,:,axi_xdcr_loc), 3) ));
        end
        xdcr.chanID=nan(length(idxLat), length(idxEle));
        xdcr.chanLoc=cell(length(idxLat), length(idxEle));

% map out position of elements
        for i=1:length(idxLat)
        for j=1:length(idxEle)
            % designate ID of current channel
            chanID=chanID+1;  la=idxLat(i);  el=idxEle(j);

            % designate indices in "chanmap" for this channel
            lat_inElem=          la                  :min([ la+id_eL-1, id_edgeL(2) ]);
            latspacing=min([ la+id_eL, id_edgeL(2) ]):min([ la+id_pL,   id_edgeL(2) ]);
            ele_inElem=          el                  :min([ el+id_eE-1, id_edgeE(2) ]);
            elespacing=min([ el+id_eE, id_edgeE(2) ]):min([ el+id_pE,   id_edgeE(2) ]);

            % apply channel ID + mask for kerf in current element
            xdcr.chanmap(lat_inElem, elespacing)=0;
            xdcr.chanmap(latspacing, ele_inElem)=0;
            xdcr.chanmap(lat_inElem, ele_inElem)=chanID;%assignment overrules zeros

            % add current channel to map of channel IDs
            xdcr.chanID(i,j)=chanID;

            % identify position of the center of current element (lateral)
            if mod(length(lat_inElem),2)~=0
                latID=lat_inElem( ceil(length(lat_inElem)/2) );
                locLat=mgrid.x(latID);
            else
                latID=lat_inElem( floor(length(lat_inElem)/2) );
                locLat=mean(mgrid.x([latID latID+1]));
            end
            
            % identify position of the center of current element (elev.)
            if nD==2
                locEle=0;
            elseif mod(length(ele_inElem),2)==0
                locEle=mgrid.y(ele_inElem( length(ele_inElem)/2 ));
            else
                ele1=ele_inElem( floor(length(ele_inElem)/2) );
                ele2=ele_inElem(  ceil(length(ele_inElem)/2) );
                locEle=mean(mgrid.x([ele1 ele2]));
            end
            
            xdcr.chanLoc{i,j}=[locLat, locEle];
        end
        end

        % FOR DEBUG: create map of channels
        %figure;
        %imagesc(1000*mgrid.x(lat_xdcr_loc), 1000*mgrid.y(ele_xdcr_loc), xdcr.chanmap');
        %axis image;  ylabel('Elev');  xlabel('Lat')
    end
end