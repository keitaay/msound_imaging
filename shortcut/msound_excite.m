function exci=msound_excite(mgrid, medium, xdcr, varargin)
% exci = MSOUND_EXCITE( mgrid, medium, xdcr )
% exci = MSOUND_EXCITE( mgrid, medium, xdcr, focus, Fnum )
% 
% Create a structure that encodes information about a transmit pressure
% waveform in mSOUND forward simulations.
%
% REQUIRED INPUT:
%          mgrid = mSOUND set_grid object
%         medium = mSOUND medium structure
%           xdcr = transducer settings, created using "msound_xdcr.m"
%
% OPTIONAL INPUT:
%          focus = focal depth
%                  [m] (scalar)
%                  DEFAULT: focal depth of transducer
%                  CAUTION: inf = plane-wave
%           Fnum = F-number
%                  (1x2 vector - [lat,ele])
%                  DEFAULT: full transducer used for aperture
%                  CAUTION: overridden by "focus" for plane wave transmits

% 2019-10-26 - Keita Yokoyama (UNC/NCSU)
%              initial version
% 2019-10-29 - Keita Yokoyama (UNC/NCSU)
%              generalized function for 1/2/3-D environments;
%              added more helpful help texts
    
    % define number of dimensions in simulation
    nD=msound_nDim(mgrid);
    
    % define focus and F-number
    switch length(varargin)
        case 0 % inputted: n/a
            focus=xdcr.focus;
            exci.Fnum=[focus/xdcr.dim(1), focus/xdcr.dim(2)];
            
        case 1 % inputted: focus
            focus=varargin{1};
            if isinf(focus)
                exci.Fnum=[0, 0];
            else
                exci.Fnum=[focus/xdcr.dim(1), focus/xdcr.dim(2)];
            end
            
        case 2 % inputted: focus, F-number
            focus=varargin{1};
            if isinf(focus)
                warning('Plane-wave transmit detected! Ignoring manual F-number input.');
                exci.Fnum=[0, 0];
            else
                exci.Fnum=[varargin{2}(1), varargin{2}(2)];
            end
    end
    
    % using F-number, derive size of active aperture
    if isinf(focus), aperSize=xdcr.dim;
    else,            aperSize=exci.Fnum./focus;
    end
    
    % create time vector
    exci.t=(-4/xdcr.fc : mgrid.dt : 4/xdcr.fc)';

    % define impulse response function based on time delay "dly"
    impresp=@(dly) sin( 2*pi * xdcr.fc * (exci.t+dly) ) .* ...
                   exp(-(exci.t+dly).^2 * (xdcr.fc)^2 / 2);
    
    % define incident pressure
    exci.P0=1e6;
    
    switch nD
        case 1
            % focusing, time delays etc. are trivial in a 1-D case
            exci.Pi=repmat(exci.P0, [exci.t,1]);
            
        case 2
            % calculate element- and depth-wise time delays
            if isinf(focus)
                delay=zeros(1,mgrid.num_x);
            else
                delay=sqrt((mgrid.x).^2 + (focus).^2)./medium.c0;
            end
            delay=repmat(delay-min(delay(:)), [length(exci.t) 1]);
            
            % convert impulse time vector to spatial matrix
            exci.t=repmat(exci.t, [1 mgrid.num_x]);
            
            % convert excitation to pressure waveform
            exci.Pi=exci.P0 .* impresp(delay);
            exci.Pi(:, abs(mgrid.x) > aperSize(1)/2)=0;
            
        case 3
    end
end
    