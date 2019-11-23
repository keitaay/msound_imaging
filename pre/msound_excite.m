function exci=msound_excite(mgrid, medium, xdcr, lineloc, focus, varargin)
% exci = MSOUND_EXCITE( mgrid, medium, xdcr, lineloc, focus )
% exci = MSOUND_EXCITE( mgrid, medium, xdcr, lineloc, focus, Fnum )
% 
% Create a structure that encodes information about a transmit pressure
% waveform in mSOUND forward simulations.
%
% REQUIRED INPUT:
%          mgrid = mSOUND set_grid object
%
%         medium = mSOUND medium structure
%
%           xdcr = transducer settings, created using "msound_xdcr.m"
%
%        lineloc = positions/orientations for origin of lines
%                  [m/deg] (Nx4 vector - [lat, elev., theta, phi])
%
%          focus = focal depth
%                  [m] (scalar)
%                  DEFAULT: focal depth of transducer
%                  CAUTION: inf = plane-wave
% OPTIONAL INPUT:
%           Fnum = F-number
%                  (1x2 vector - [lat,ele])
%                  DEFAULT: full transducer used for aperture
%                  CAUTION: overridden by "focus" for plane wave transmits

% 2019-10-26 - Keita Yokoyama (UNC/NCSU)
%              initial version
% 2019-10-29 - Keita Yokoyama (UNC/NCSU)
%              (1) generalized function for 1/2/3-D environments
%              (2) added more helpful help texts
% 2019-11-14 - Keita Yokoyama (UNC/NCSU)
%              modified output scheme to allow for multiple transmits in
%              one imaging sequence (e.g. multiline conventional imaging)
% 2019-11-17 - Keita Yokoyama (UNC/NCSU)
%              (1) added required input for line positions
%              (2) added steering capabilities for plane/phased-array
%              (3) separated time-delay index finder into separate func.
    
% define number of dimensions in simulation
    nD=msound_nDim(mgrid);
    
% define line positions
    linepos=cell( size(lineloc,1),1 );
    for line=1:length(linepos), linepos{line}=lineloc( line,: ); end
    
% initialize flag for warning
    flagShown=false;
    
% define excitation object
    exci=struct('linepos', cell(length(linepos),1),...
                'focus',[], 'Fnum',[], 't',[], 'P0',[], 'Pchan',[], 'Pi',[],...
                'apodTx',[],'apodRx',[]);
    
    if nD==1, linepos={[0,0,0,0]}; end
    for line=1:length(linepos)
        
% define position of line + focal depth
        posNow=linepos{line};  exci(line).linepos=posNow;
        exci(line).focus=focus;
        if isinf(focus), wavetype='plane'; else, wavetype='focused'; end
        
% define F-number
        switch length(varargin)
            case 0 % inputted: n/a
                if strcmp(wavetype,'plane')
                    warning(['Plane wave detected, but no F-numbers were defined! ',...
                             'The entire array will be used.']);
                    exci(line).Fnum=[1, 1];
                else
                    exci(line).Fnum=[exci(line).focus/xdcr.dim(1),...
                                     exci(line).focus/xdcr.dim(2)];
                end

            case 1 % inputted: F-number
                if strcmp(wavetype,'plane') && ~flagShown
                    warning(['Plane wave detected! Interpreting F-numbers as the inverse of',...
                             'the fraction of the transducer used for active elements.']);
                    flagShown=true;
                end
                exci(line).Fnum=[varargin{1}(1), varargin{1}(2)];
        end

% using F-number, derive size of active aperture
        if strcmp(wavetype,'plane')
            aperSize=[ xdcr.dim(1)/exci(line).Fnum(1),  xdcr.dim(2)/exci(line).Fnum(2) ];
        else
            aperSize=exci(line).focus./exci(line).Fnum .* cosd([-lineloc(3) -lineloc(4)]);
        end

% define parameters for impulse response function (Gaussian-weighted sinusoid)
        %ka=2*pi*xdcr.fc;   kb=-xdcr.fc^2/2;
        fc=xdcr.fc;  bw=0.5;  dr=xdcr.dynrange;
        tLim=gauspuls('cutoff', fc , bw,[], -dr); 
        exci(line).t=(-tLim:mgrid.dt:tLim)';
        %exci(line).t=(-4/xdcr.fc : mgrid.dt : 4/xdcr.fc)'; % create time vector
        
% define incident pressure
        exci(line).P0=1e6;

        if nD==1
            t_imp=exci(line).t;                  % define impulse response's time vector
            %h=sin(ka.*t_imp).*exp(kb.*t_imp.^2); % calculate base transmit waveform
            h=gauspuls( t_imp, fc, bw );
            exci(line).Pchan=exci(line).P0 .* h; % set transmit pres. for all channels
            exci(line).Pi=exci(line).Pchan;      % set transmit across init. cond. mask

        else
            [delay, chLat, chEle]=getTimeDelays( xdcr, posNow, exci(line).focus, medium.c0, wavetype );
            chLat=permute(chLat, [2,3,1]);       % channels' lateral position wrt beam origin
            chEle=permute(chEle, [2,3,1]);       % channels' elevational position wrt beam origin
            
            % create time vector for duration of impulse
            num_t=knnsearch(mgrid.t', max(delay(:)));
            tAdj1=-mgrid.dt*floor(num_t/2);
            tAdjN= mgrid.dt*floor(num_t/2);
            tLim1=min([exci(line).t(1),  -max(abs(delay(:)))])+tAdj1;
            tLimN=max([exci(line).t(end), max(abs(delay(:)))])+tAdjN;
            delay_ext=( tLim1 : mgrid.dt : tLimN)';
            
            % use matrix of time delay values to calculate waveform for transmitted pulse
            t_imp=repmat( delay_ext, [1, size(delay,2), size(delay,3)]) +...
                  repmat( delay,     [length(delay_ext), 1, 1]);

% calculate impulse response at all of above timepoints/element positions
            %h=sin(ka.*t_imp).*exp(kb.*t_imp.^2);
            h=gauspuls( t_imp, fc, bw );

% scale excitation by peak pressure out
            exci(line).Pchan=exci(line).P0 .* h;
            
% mask out channels that fall outside of the desired aperture
            aperHL=( aperSize(1)/2 )*( cosd(posNow(3)) + sind(posNow(3))*tand(posNow(3)) );
            aperEL=( aperSize(2)/2 )*( cosd(posNow(4)) + sind(posNow(4))*tand(posNow(4)) );
            active=( abs(chLat) <= aperHL )&...
                   ( abs(chEle) <= aperEL ); % change to middle of 'active'
            exci(line).Pchan( :, ~active )=0;
            
% define apodization (default: 2D Hamming)
            apodTx=zeros(size(active));
            apodL=hamming( max(sum(active,1)) );
            apodE=hamming( max(sum(active,2)) )';
            apodTx(active)=apodL*apodE;
            apodRx=apodTx;  % make receive apod. same as transmit
            exci(line).apodTx=apodTx;   exci(line).apodRx=apodRx;
            
            % extend apodization over time, and apply to transmit
            apodExt=permute( repmat(apodTx,[ 1,1,length(delay_ext) ]), [3,1,2]);
            exci(line).Pchan=exci(line).Pchan.*apodExt;
            
% apply dynamic range limit to transmit (i.e. suppress very tiny pres. amplitudes)
            exci(line).Pi( abs(exci(line).Pchan) <...
                           max(exci(line).Pchan(:))*10.^(-xdcr.dynrange/20) )=0;

            % prepare storage matrix for transmit pressure waveforms
            if nD==2, exci(line).Pi=permute( zeros(size(xdcr.mask,1),...
                                                   size(t_imp,1)), [2,1]);
            else,     exci(line).Pi=permute( zeros(size(xdcr.mask,1),...
                                                   size(xdcr.mask,2),...
                                                   size(t_imp,1)), [3,1,2]);
            end

% apply channel-based pressure waveform to simulation grid
            for chE=1:size(xdcr.chanID,2)
            for chL=1:size(xdcr.chanID,1)
                chanNow=xdcr.chanID(chL,chE);

                % get mask for current channel in transmit waveform
                maskNow=xdcr.chanmap==chanNow;
                maskNow=permute( repmat(maskNow,[1,1,size(t_imp,1)]), [3,1,2]);

                % get number of lat/ele. gridspaces in this channel
                nStepL=mean(nonzeros( sum(maskNow(1,:,:),2 )));
                nStepE=mean(nonzeros( sum(maskNow(1,:,:),3 )));

                % copy transmit waveform for all subele. in current channel
                outNow=squeeze(exci(line).Pchan(:,chL,chE));
                outNow=repmat(outNow, [1, nStepL, nStepE]);

                % store pressure waveform to transmit matrix
                exci(line).Pi(maskNow)=outNow;
            end
            end
        end
    end
end
    