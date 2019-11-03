function [RF,axial] = msound_beamform(chanData, mgrid, medium, xdcr, exci)
% [RF,axial] = MSOUND_BEAMFORM(chanData, mgrid, medium, xdcr, exci)
% 
% Beamform RF data using channel data acquired from an mSOUND simulation.
%
% This script only works for delay-sum beamforming with dynamic-focus
% receive.
%
% CAUTION: this code is intended only for simulations that are done
%          ENTIRELY using mSOUND. Simulations that need to be coupled with
%          other software, such as Field-II or k-Wave, should be
%          beamformed on there instead.
%
% REQUIRED INPUT:
%       chanData = mSOUND result, re-organized as channel data for the
%                  transducer "xdcr" using "msound_get_chan.m"
%          mgrid = mSOUND set_grid object
%         medium = mSOUND medium structure
%           xdcr = transducer settings, created using "msound_xdcr.m"
%           exci = excitation settings, created using "msound_excite.m"

% 2019-10-26 - Keita Yokoyama (UNC/NCSU)
%              initial version
% 2019-10-29 - Keita Yokoyama (UNC/NCSU)
%              laid out groundworks to eventually generalize this function
%              for 3-D simulations; added more helpful help texts
% 2019-10-30 - Keita Yokoyama (UNC/NCSU)
%              added Hanning apodization; re-formatted to deal with
%              matrix probe-like datasets (via msound_reform.m) by default

    % define number of dimensions in simulation
    nD=msound_nDim(mgrid);
    
    % reconstruct axial label, based on time vector
    axial=medium.c0*(mgrid.t'./2);%
    
    % correct for position of transducer plane
    switch nD
        case 1, axial=axial-(xdcr.plane_depth-mgrid.x(1));
        case 2, axial=axial-(xdcr.plane_depth-mgrid.y(1));
        case 3, axial=axial-(xdcr.plane_depth-mgrid.z(1));
    end
    
    % define scan lines
    numLineL=xdcr.Nelem(1);  dx=xdcr.pitch(1);
    numLineE=xdcr.Nelem(2);  dy=xdcr.pitch(2);
    RF=zeros(size(chanData,1), numLineL, numLineE);
    
    % create reference apodization matrix, which contains all possible
    % apodization weights over channels AND lines in lat/elev dimensions
    apodLat=repmat( hann(size(chanData,2)+numLineL-1), [1, 2*numLineE-1]);
    apodEle=repmat( hann(size(chanData,3)+numLineE-1)',[2*numLineL-1, 1]);
    apodAll=apodEle.*apodLat;
    
    for eleline=1:numLineE
    for latline=1:numLineL
        % get apodization weights for current line
        if length(0:numLineL-1)==1
            apodNow=1;
        elseif length(0:numLineE-1)==1
            apodNow=apodAll( latline+(0:numLineL-1), 1 );
        else
            apodNow=apodAll( latline+(0:numLineL-1), eleline+(0:numLineE-1) );
        end
        
        % get label vectors for all channels of current line
        lat1D=dx*( (1:numLineL)-latline );  if isempty(lat1D), lat1D=0; end
        ele1D=dy*( (1:numLineE)-eleline );  if isempty(ele1D), ele1D=0; end

        % create matrix version of label vectors
        axi3D=repmat(axial, [1, size(chanData,2), size(chanData,3)]);
        lat3D=repmat(lat1D, [size(chanData,1), 1, size(chanData,3)]);
        ele3D=repmat(ele1D, [size(chanData,1), size(chanData,2), 1]);
        
        % dynamic-focus receive: find indices to shift pressure data
        delay=sqrt( axi3D.^2 + lat3D.^2 + ele3D.^2 )./medium.c0;
        delay_shift=repmat( min(min(delay,[],3),[],2),...
                           [1, size(chanData,2), size(chanData,3)]);
        delay=delay-delay_shift;
        delay_idx=round( delay./mgrid.dt );

        % dynamic-focus receive: shift channel data for all positions
        idx3D=repmat( (1:size(chanData,1))', [1, size(chanData,2), size(chanData,3)]);
        idx3D=idx3D-delay_idx;

        % dynamic-focus receive: correct for over-shifted values
        idx3D( idx3D>size(chanData,1) )=size(chanData,1);
        idx3D( idx3D<1 )=1;
        
        % angular sensitivity
        
        % sum together the delayed channel data
        chanSum=zeros(size(idx3D,1),1);
        for chanX=1:size(chanData,2)
        for chanY=1:size(chanData,3)
            idx_Now=squeeze( idx3D(:,chanX,chanY) );            % ID this channeel
            % angular shift of ID
            chanNow=squeeze( chanData(idx_Now, chanX, chanY) ); % get delayed channel data
            chanNow=squeeze( chanNow.*apodNow(chanX, chanY) );  % apodize delayed channels
            chanSum=chanSum+chanNow;                            % sum channels
        end
        end
        RF(:,latline,eleline)=chanSum;
    end
    end
    
    % truncate RF data to remove excitation
    idx_offset=floor( size(exci.Pi,1)/2 );           % define duration of excitation       
    axial=axial(idx_offset:end)-idx_offset*mean(diff(axial));
    RF=RF(idx_offset:end, :, :);
    
    % truncate RF data to remove "echoes" from impossible axial depths
    switch nD
        case 1, axlim=knnsearch(axial, min(axial(axial>mgrid.x_length)));
        case 2, axlim=knnsearch(axial, min(axial(axial>mgrid.y_length)));
        case 3, axlim=knnsearch(axial, min(axial(axial>mgrid.z_length)));
    end
    axial=axial(1:axlim);  RF=RF(1:axlim,:,:);
end