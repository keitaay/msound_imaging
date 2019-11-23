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
%              (1) laid out groundworks to generalize for 3-D simulations
%              (2) added more helpful help texts
% 2019-10-30 - Keita Yokoyama (UNC/NCSU)
%              added Hanning apodization; re-formatted to deal with
%              matrix probe-like datasets (via msound_reform.m) by default
% 2019-11-13 - Keita Yokoyama (UNC/NCSU)
%              (1) reformatted delay-sum scheme to accept channel data for
%                  multiple transmit-receive sequences
%              (2) added capability for plane-wave and phased-array
%                  beam steering
%              (3) separated time-delay index finder into separate func.

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
    
% initialize storage variable for RF data
    RF=zeros(size(chanData,1), length(exci));
    
    for lineID=1:length(exci)
% get positions, steering angles, & apodization of current line
        posNow=exci(lineID).linepos;
        apodRx=exci(lineID).apodRx;
        
% calculate time delays
        delay=getTimeDelays( xdcr, posNow, axial, medium.c0, 'focused' );

% shift channel data for all positions
        idx3D=repmat( (1:size(chanData,1))',...
                      [1, size(chanData,2), size(chanData,3)]...
                    )-round( delay./mgrid.dt );

% correct time delay indices for over-shifted values
        idx3D( idx3D>size(chanData,1) )=size(chanData,1);
        idx3D( idx3D<1 )=1;
        
% apodize, then sum channel data
        chanSum=zeros(size(chanData,1),1);
        for chanY=1:xdcr.Nelem(2)
        for chanX=1:xdcr.Nelem(1)
            idxNow=squeeze( idx3D(:,chanX,chanY) );
            chanNow=squeeze( chanData(idxNow, chanX, chanY, lineID) );
            chanNow=squeeze( chanNow.*apodRx(chanX, chanY) );
            chanSum=chanSum+chanNow;
        end
        end
        
% append beamformed line to RF output
        RF(:,lineID)=chanSum;
    end
    
% truncate outputs to remove "echoes" from impossible depths
    switch nD
        case 1, axlim=knnsearch(axial, min(axial(axial>mgrid.x_length)));
        case 2, axlim=knnsearch(axial, min(axial(axial>mgrid.y_length)));
        case 3, axlim=knnsearch(axial, min(axial(axial>mgrid.z_length)));
    end
    axial=axial(1:axlim);  RF=RF(1:axlim,:);
end