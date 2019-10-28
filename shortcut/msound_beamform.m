function RF = msound_beamform(P, mgrid, medium, xdcr, exci)
% RF = MSOUND_BEAMFORM(P, mgrid, medium, xdcr, exci)
% 
% Beamform RF data using channel data acquired from an mSOUND simulation.
%
% This function uses:
%
%   - P      : pressure waveform
%   - mgrid  : mSOUND set_grid object
%   - medium : mSOUND medium structure
%   - xdcr   : transducer parameter structure, made with msound_xdcr.m
%   - exci   : transmit parameter structure, made with msound_excite.m

% 2019-10-26 - Keita Yokoyama (UNC/NCSU)
%              initial version

    axial=medium.c0 * mgrid.t(size(exci.t,1)+1:end)';
    chanData=P( size(exci.t,1)+knnsearch(axial,0)+1 : end,   : );
    
    clearvars -except mgrid medium xdcr P axial chanData
    for line=1:size(chanData,2)

        % get positions (presumes elements have same spatial frequency as sim)
        pitch=mgrid.dx;
        x_orig=(-xdcr.dim(1)/2 : pitch : xdcr.dim(1)/2)+... % lateral vector
               pitch*floor(line-size(chanData,2)/2);        % shift to line pos

        % below is assumption-free
        axi2D=repmat(axial,  [1, size(chanData,2)]);
        lat2D=repmat(x_orig, [size(chanData,1), 1]);

        % dynamic-focus receive
        delay=sqrt( axi2D.^2 + lat2D.^2 )./medium.c0;
        delay_idx=round( (delay-min(delay,[],2))./mgrid.dt );

        % shift channel data by above indexc, for all positions
        idx2D=repmat( (1:size(chanData,1))', [1, size(chanData,2)]);
        idx2D=idx2D-delay_idx;

        % correct for over-shifted values
        idx2D( idx2D>size(chanData,1) )=size(chanData,1);
        idx2D( idx2D<1 )=1;

        RF(:,line)=sum(chanData(idx2D), 2);
    end
end