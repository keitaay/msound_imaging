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
    else,            aperSize=focus./exci.Fnum;
    end
    
    % create time vector
    exci.t=(-4/xdcr.fc : mgrid.dt : 4/xdcr.fc)';

    % define impulse response function. the impulse response function is:
    %   h = sin( ka t ) * exp( kb t^2 )
    %
    % formerly an array function:
    % impresp=@(t) sin(2*pi*xdcr.fc*t) .* exp(-t.^2 * xdcr.fc^2 / 2);
    ka=2*pi*xdcr.fc;   kb=-xdcr.fc^2/2;
    
    % define incident pressure
    exci.P0=1e6;
    
    switch nD
        case 1
            % calculate impulse response
            t_imp=exci.t;
            h=sin(ka.*t_imp).*exp(kb.*t_imp.^2);
            
            % focusing, time delays etc. are trivial in a 1-D case
            exci.Pi=exci.P0 .* h;
            
            % apply dynamic range limit to excitation
            exci.Pi( abs(exci.Pi) < max(exci.Pi(:))*10.^(-xdcr.dynrange/20) )=0;
            
        case 2
            % calculate time delays in each element for focusing
            if isinf(focus)
                delay=zeros(1,mgrid.num_x);
            else
                delay=sqrt((mgrid.x).^2 + (focus).^2)./medium.c0;
            end
            delay=delay-min(delay(:));
            
            % extend time delay vector for duration of impulse
            num_t=knnsearch(mgrid.t', max(delay));
            delay_ext=( exci.t(1)-mgrid.dt*floor(num_t/2) : mgrid.dt : exci.t(end)+mgrid.dt*floor(num_t/2) )';
            
            % use matrix of time delay values to calculate
            % waveform for transmitted pulse sequence
            t_imp=repmat(delay_ext, [1 length(delay)]) + repmat(delay, [length(delay_ext) 1]);
            
            % calculate impulse response at all of above timepoints
            % and element positions
            h=sin(ka.*t_imp).*exp(kb.*t_imp.^2);
            
            % convert excitation to pressure waveform
            exci.Pi=exci.P0 .* h;
            exci.Pi(:, abs(mgrid.x) > aperSize(1)/2)=0;
            
            % mask out pressure waveforms that are not coming from an
            % actual, physically-real element position
            exci.Pi(:, xdcr.chanmap==0)=0;
            
            % apply dynamic range limit to excitation
            exci.Pi( abs(exci.Pi) < max(exci.Pi(:))*10.^(-xdcr.dynrange/20) )=0;
            
        case 3
            % calculate time delays in each element for focusing
            if isinf(focus)
                delay=zeros(1,mgrid.num_x,mgrid.num_y);
            else
                x=repmat(mgrid.x', 1,mgrid.num_y);
                y=repmat(mgrid.y,  mgrid.num_x,1);
                delay=sqrt(x.^2 + y.^2 + focus.^2)./medium.c0;
            end
            delay=delay-min(delay(:));
            
            % extend time delay vector for duration of impulse
            num_t=knnsearch(mgrid.t', max(delay(:)));
            delay_ext=( exci.t(1)  - mgrid.dt*floor(num_t/2) : mgrid.dt :...
                        exci.t(end)+ mgrid.dt*floor(num_t/2))';
            
            % change "delay" into 3-D matrix
            % (this makes time delays easier to find)
            delay=permute(delay, [3,1,2]);
            
            % use matrix of time delay values to calculate
            % waveform for transmitted pulse sequence
            t_imp=repmat( delay_ext, [1, size(delay,2), size(delay,3)]) +...
                  repmat( delay,     [length(delay_ext), 1, 1]);
            
            % calculate impulse response at all of above timepoints
            % and element positions
            h=sin(ka.*t_imp).*exp(kb.*t_imp.^2);
            
            % convert excitation to pressure waveform
            exci.Pi=exci.P0 .* h;
            exci.Pi(:, abs(mgrid.x) > aperSize(1)/2)=0;
            exci.Pi(:, abs(mgrid.y) > aperSize(2)/2)=0;
            
            % mask out pressure waveforms that are not coming from an
            % actual, physically-real element position
            maskNow=false(size(xdcr.chanmap));    % initialize
            maskNow(xdcr.chanmap==0)=true;        % select non-element
            maskNow=repmat( permute(maskNow,[3,1,2]),...
                           [size(exci.Pi,1),1,1]);% stretch to 3D, over time
            exci.Pi(maskNow)=0;                   % zero-pressure
            
            % apply dynamic range limit to excitation
            exci.Pi( abs(exci.Pi) < max(exci.Pi(:))*10.^(-xdcr.dynrange/20) )=0;
    end
end
    