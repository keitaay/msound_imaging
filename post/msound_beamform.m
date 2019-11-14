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
%          beamformed on there, instead.
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
% 2019-11-13 - Keita Yokoyama (UNC/NCSU)
%              reformatted delay-sum scheme to only output individual
%              A-lines (will support B-mode outputs after reformatting
%              transmit sequence)

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
    
% extract lateral/elevational positions of each scan line
    numLines=size(chanData,4);  lineL=zeros(numLines,1);  lineE=lineL;
    for lineID=1:numLines
        lineL(lineID)=exci(lineID).linepos(1);
        lineE(lineID)=exci(lineID).linepos(2);
    end
    lineL=unique(lineL);  numLineL=length(lineL);
    lineE=unique(lineE);  numLineE=length(lineE);
    
% extract lateral/elevational positions of each channel
    chanL=zeros(xdcr.Nelem(1)*xdcr.Nelem(2), 1); chanE=chanL;
    for chanY=1:xdcr.Nelem(2)
    for chanX=1:xdcr.Nelem(1)
        chanL( chanX + xdcr.Nelem(1)*(chanY-1) )=xdcr.chanLoc{chanX,chanY}(1);
        chanE( chanX + xdcr.Nelem(1)*(chanY-1) )=xdcr.chanLoc{chanX,chanY}(2);
    end
    end
    chanL=unique(chanL);  chanE=unique(chanE);
    
% create matrix of apodizations for each channel
    apodLat=repmat( hann(xdcr.Nelem(1)), [1, xdcr.Nelem(2)]);
    apodEle=repmat( hann(xdcr.Nelem(2))',[xdcr.Nelem(1), 1]);
    apodAll=apodEle.*apodLat;
    
% initialize storage variable for RF data
    RF=zeros(size(chanData,1), numLineL, numLineE);
    
    for eleline=1:numLineE
    for latline=1:numLineL
% convert lateral/elevational line indices to line ID for "chanData"
        lineID=latline+(eleline-1)*numLineL;
        
% get label vectors for all channels with respect to line position
        lat1D=chanL'+lineL(latline); if isempty(lat1D), lat1D=0; end
        ele1D=chanE+lineE(eleline);  if isempty(ele1D), ele1D=0; end

% create matrix version of label vectors
        axi3D=repmat(axial, [1, size(chanData,2), size(chanData,3)]);
        lat3D=repmat(lat1D, [size(chanData,1), 1, size(chanData,3)]);
        ele3D=repmat(ele1D, [size(chanData,1), size(chanData,2), 1]);
        
% dynamic-focus receive: calculate time delays for each depth/line position
        delay=sqrt( axi3D.^2 + lat3D.^2 + ele3D.^2 )./medium.c0;
        
% dynamic-focus receive: find indices to shift pressure data
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
        
% apodize, then sum channel data
        chanSum=zeros(size(chanData,1),1);
        for chanY=1:xdcr.Nelem(2)
        for chanX=1:xdcr.Nelem(1)
            idxNow=squeeze( idx3D(:,chanX,chanY) );
            chanNow=squeeze( chanData(idxNow, chanX, chanY, lineID) );
            chanNow=squeeze( chanNow.*apodAll(chanX, chanY) );
            chanSum=chanSum+chanNow;
        end
        end
        
% append beamformed line to RF output
        RF(:,latline,eleline)=chanSum;
    end
    end
    
% truncate outputs to remove "echoes" from impossible depths
    switch nD
        case 1, axlim=knnsearch(axial, min(axial(axial>mgrid.x_length)));
        case 2, axlim=knnsearch(axial, min(axial(axial>mgrid.y_length)));
        case 3, axlim=knnsearch(axial, min(axial(axial>mgrid.z_length)));
    end
    axial=axial(1:axlim);  RF=RF(1:axlim,:,:);
end