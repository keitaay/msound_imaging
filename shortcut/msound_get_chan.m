function chanData = msound_get_chan(P, mgrid, xdcr)
% chanData = MSOUND_BEAMFORM(P, mgrid, xdcr)
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
%              P = mSOUND result, reshaped with "msound_reform.m"
%          mgrid = mSOUND set_grid object
%           xdcr = transducer settings, created using "msound_xdcr.m"

% 2019-10-31 - Keita Yokoyama (UNC/NCSU)
%              initial version
    
    % define number of dimensions in simulation
    nD=msound_nDim(mgrid);
    
    % average pressure waveform readings across the axial
    % thickness of the transducer (region of pressure recording)
    P_red=permute( mean(P,2), [1,3,4,2]);
    
    % get channel ID from transducer map
    channels=unique(xdcr.chanmap(:));  channels=channels( channels>0 );
    
    % reshape list of channels into shape of transducer
    channels=reshape(channels, [xdcr.Nelem(1),xdcr.Nelem(2)]);
    
    % for each channel in the lateral/elevation directions...
    chanData=zeros(size(P_red,1), xdcr.Nelem(1), xdcr.Nelem(2));
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
        
        % copy "chanNow" over time
        chanNow=repmat( permute(chanNow,[3,1,2]), [size(P_red,1),1,1] );
        
        % get pressure readings, only for the current channel
        P_now=reshape(P_red(chanNow), size(P_red,1), []);
        
        % combine all pressure readings in "P_now" into a single vector,
        % and save that sum as the current channel
        chanData(:,chL,chE)=mean(P_now,2);
    end
    end
    
end