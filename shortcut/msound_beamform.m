function RF = msound_beamform(P, mgrid, medium, xdcr, exci)
% RF = MSOUND_BEAMFORM(P, mgrid, medium, xdcr, exci)
% 
% Beamform RF data using channel data acquired from an mSOUND simulation.
%
% REQUIRED INPUT:
%              P = pressure waveform from mSOUND forward simulation
%          mgrid = mSOUND set_grid object
%         medium = mSOUND medium structure
%           xdcr = transducer settings, created using "msound_xdcr.m"
%           exci = excitation settings, created using "msound_excite.m"

% 2019-10-26 - Keita Yokoyama (UNC/NCSU)
%              initial version
% 2019-10-29 - Keita Yokoyama (UNC/NCSU)
%              laid out groundworks to eventually generalize this function
%              for 3-D simulations; added more helpful help texts

    % define number of dimensions in simulation
    nD=msound_nDim(mgrid);
    if nD~=2, error(['msound_beamform is not yet ready for ',num2str(nD),'-D simulation results!']); end
    
    % define time shift from transducer position
    % (in case transducer plane is not placed at zero depth)
    idx_offset=size(exci.t,1)+knnsearch(mgrid.y', mgrid.y(1)+xdcr.plane_depth);
    
    % define expected axial label vector
    axial=medium.c0 * mgrid.t(idx_offset+1:end)';
    
    % isolate channel data from input
    chanData=P( idx_offset+knnsearch(axial,0) : end,   : );
    
    % initialize storage variables
    numLineL=size(chanData,2);
    RF=zeros(size(chanData,1), numLineL);
    
    
    switch nD
        case 1
            % beamforming is trivial when there are no spatial
            % dimensions to beamform through...
            RF=chanData;
            
        case 2
            for latline=1:numLineL
                % spacing between A-lines = lateral spatial freq
                pitch=mgrid.dx;
                
                % get index for current line position
                lineIDnow=floor( latline - size(chanData,2)/2 );
                
                % define lateral vector
                x_orig=pitch*lineIDnow - (-xdcr.dim(1)/2 : pitch : xdcr.dim(1)/2);

                % create 2-D matrix version of label vectors
                axi2D=repmat(axial,  [1, size(chanData,2)]);
                lat2D=repmat(x_orig, [size(chanData,1), 1]);

                % dynamic-focus receive
                delay=sqrt( axi2D.^2 + lat2D.^2 )./medium.c0;
                delay_idx=round( (delay-min(delay,[],2))./mgrid.dt );

                % shift channel data by above index, for all positions
                idx2D=repmat( (1:size(chanData,1))', [1, size(chanData,2)]);
                idx2D=idx2D-delay_idx;

                % correct for over-shifted values
                idx2D( idx2D>size(chanData,1) )=size(chanData,1);
                idx2D( idx2D<1 )=1;

                % sum together the delayed channel data
                chanSum=zeros(size(idx2D,1),1);
                for chan=1:size(chanData,2)
                    chanNow=chanData( idx2D(:,chan), chan );
                    chanSum=chanSum+chanNow;
                end
                RF(:,latline)=chanSum;
            end
            
        case 3
            error('How did you even get here???')
            %for eleline=1:numLineE
            %for latline=1:numLineL
            %end
            %end
            
    end
end