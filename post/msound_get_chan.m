function chanData = msound_get_chan(P, mgrid, xdcr)
% chanData = MSOUND_GET_CHAN(P, mgrid, xdcr)
% 
% Consolidate reshaped mSOUND pressure output so that there is only one
% waveform per channel.
%
% Each position in an mSOUND simulation output encodes a point in the
% mSOUND grid; in Field-II terminology, this amounts to raw simulation
% outputs accounting for each individual subelement. This function converts
% them into element-wise readouts.
%
% CAUTION: this code is intended only for simulations that are done
%          ENTIRELY using mSOUND. Simulations that need to be coupled with
%          other software, such as Field-II or k-Wave, should be
%          beamformed on there instead.
%
% REQUIRED INPUT:
%              P = mSOUND result, reshaped with "msound_reform.m"
%          mgrid = mSOUND set_grid object
%           xdcr = transducer settings, created using "msound_xdcr.m"

% 2019-10-31 - Keita Yokoyama (UNC/NCSU)
%              initial version
% 2019-11-14 - Keita Yokoyama (UNC/NCSU)
%              added support for pressure inputs of multiple transmits
%              (element in cell array = individual transmit sequence)
    
% define number of dimensions in simulation
    nD=msound_nDim(mgrid);
    
% concatenate waveforms from all lines into single matrix
%ACCOUNT FOR ANGULAR SENSITIVITY HERE?
    if iscell(P) % presumes all pressure output matrices have same dimensions
        Nlines=length(P);
        switch nD
            case 1, P_red=nan([size(P{1}),1,1,Nlines]);
            case 2, P_red=nan([size(P{1}), 1, Nlines]);
            case 3, P_red=nan([size(P{1}),    Nlines]);
        end
        for line=1:Nlines, P_red(:,:,:,line)=P{line}; end
    else, Nlines=1; P_red=P;
    end
    
% average pressure waveform readings across transducer's
% axial thickness (AKA axial range of region of simulation recording)
% - dimensions: 1=time, 2=axial, 3=lateral, 4=elev., 5=line
    P_red=permute( mean(P_red,2), [1,3,4,5,2]);
    
% get channel ID
    channels=xdcr.chanID;
    
% for each line + channel...
    chanData=zeros(size(P_red,1), xdcr.Nelem(1), xdcr.Nelem(2), Nlines);
    for line=1:Nlines
    for chL=1:size(channels,1)
    for chE=1:size(channels,2)
        
        % from "chanmap", isolate out a mask for the current channel
        chIDnow=channels(chL,chE);
        chanNow=false(size(xdcr.chanmap));     % all false as default
        chanNow(xdcr.chanmap==chIDnow)=true;   % current channel = true
        
        % truncate parts of channel map that are not included in "P"
        switch nD
            case 2
                idL=knnsearch(mgrid.x',-xdcr.dim(1)/2):knnsearch(mgrid.x', xdcr.dim(1)/2);
                chanNow=chanNow(idL);
            case 3
                idL=knnsearch(mgrid.x',-xdcr.dim(1)/2):knnsearch(mgrid.x', xdcr.dim(1)/2);
                idE=knnsearch(mgrid.y',-xdcr.dim(2)/2):knnsearch(mgrid.y', xdcr.dim(2)/2);
                chanNow=chanNow(idL,idE);
        end
        
        % copy "chanNow" to match dimensions of "P_red"
        chan4D=false(size(P_red));
        chan4D(:,:,:,line)=repmat( permute(chanNow,[3,1,2]), [size(P_red,1),1,1] );
        
        % get pressure readings, only for the current channel, and
        % reformat into a 2D matrix (time x subelements)
        P_now=reshape(P_red(chan4D), size(P_red,1), []);
        
        % combine all pressure readings in "P_now" into a single vector,
        % and save that sum as the current channel
        chanData(:,chL,chE,line)=mean(P_now,2);
    end
    end
    end
    
end