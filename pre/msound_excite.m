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
% 2019-11-14 - Keita Yokoyama (UNC/NCSU)
%              modified output scheme to allow for multiple transmits in
%              one imaging sequence (e.g. multiline conventional imaging)
    
% define number of dimensions in simulation
    nD=msound_nDim(mgrid);
    
% define line positions
    lineL=-0.01 : 1e-3 : 0.01;   lineE=0;
    linepos=cell( length(lineL)*length(lineE),1 );
    nLines=0; % counter for allocating line positions as 1D array of cells
    for lat=1:length(lineL)
    for ele=1:length(lineE)
        nLines=nLines+1;  linepos{nLines}=[lineL(lat), lineE(ele)];
    end
    end
    
% define excitation object
    exci=struct('linepos', cell(length(linepos),1),...
                'Fnum',[], 't',[], 'P0',[], 'Pchan',[], 'Pi',[]);
    
% define details of transmit/excitation for each line
    if nD==1, linepos={[0,0]}; end
    for line=1:length(linepos)
        posNow=linepos{line};  exci(line).linepos=posNow;
        
% define focus and F-number
        switch length(varargin)
            case 0 % inputted: n/a
                exci(line).focus=xdcr.focus;
                exci(line).Fnum=[exci(line).focus/xdcr.dim(1),...
                                 exci(line).focus/xdcr.dim(2)];

            case 1 % inputted: focus
                exci(line).focus=varargin{1};
                if isinf( exci(line).focus )
                    exci(line).Fnum=[0, 0];
                else
                    exci(line).Fnum=[exci(line).focus/xdcr.dim(1),...
                                     exci(line).focus/xdcr.dim(2)];
                end

            case 2 % inputted: focus, F-number
                exci(line).focus=varargin{1};
                if isinf( exci(line).focus )
                    warning('Plane-wave transmit detected! Ignoring F-number input.');
                    exci(line).Fnum=[0, 0];
                else
                    exci(line).Fnum=[varargin{2}(1), varargin{2}(2)];
                end
        end

% using F-number, derive size of active aperture
        if isinf(exci(line).focus)
            aperSize=xdcr.dim;
        else
            aperSize=exci(line).focus./exci(line).Fnum;
        end

% create time vector
        exci(line).t=(-4/xdcr.fc : mgrid.dt : 4/xdcr.fc)';

% define coefficients for impulse response function, taking the form:
% h = sin( ka t ) * exp( kb t^2 )
        ka=2*pi*xdcr.fc;   kb=-xdcr.fc^2/2;

% define incident pressure
        exci(line).P0=1e6;

        if nD==1
            t_imp=exci(line).t;                  % define impulse response's time vector
            h=sin(ka.*t_imp).*exp(kb.*t_imp.^2); % calculate base transmit waveform
            exci(line).Pchan=exci(line).P0 .* h; % set transmit pres. for all channels
            exci(line).Pi=exci(line).Pchan;      % set transmit across init. cond. mask

        else
            % calculate time delays in each element for focusing
            if isinf( exci(line).focus )
                delay=zeros(1, size(xdcr.chanID,1),size(xdcr.chanID,2));
            else
                latPos=nan(size(xdcr.chanID,1), 1);  % prepare to find lat. pos. of elem.
                elePos=nan(size(xdcr.chanID,2), 1);  % prepare to find ele. pos. of elem.
                for lat=1:length(latPos),  latPos(lat)=xdcr.chanLoc{lat,1}(1);  end
                for ele=1:length(elePos),  elePos(ele)=xdcr.chanLoc{1,ele}(2);  end
                latPos=latPos-posNow(1);             % shift focus by lateral line pos.
                elePos=elePos-posNow(2);             % shift focus by elevat. line pos.
                x=repmat(latPos,  1,length(elePos)); % 2-D matrix of lat position
                y=repmat(elePos', length(latPos),1); % 2-D amtrix of elev. posit.
                delay=sqrt(x.^2 + y.^2 + exci(line).focus.^2)./medium.c0;
            end
            delay=delay-min(delay(:));

            % create time vector for duration of impulse
            num_t=knnsearch(mgrid.t', max(delay(:)));
            delay_ext=( exci(line).t(1)  - mgrid.dt*floor(num_t/2) : mgrid.dt :...
                        exci(line).t(end)+ mgrid.dt*floor(num_t/2))';

            % change element-wise delays into 3-D matrix for easier finding
            delay=permute(delay, [3,1,2]);

            % use matrix of time delay values to calculate waveform for transmitted pulse
            t_imp=repmat( delay_ext, [1, size(delay,2), size(delay,3)]) +...
                  repmat( delay,     [length(delay_ext), 1, 1]);

            % calculate impulse response at all of above timepoints/element positions
            h=sin(ka.*t_imp).*exp(kb.*t_imp.^2);

            % scale excitation by peak pressure out
            exci(line).Pchan=exci(line).P0 .* h;

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

            % zero out elements that fall outside of current F-number
            if nD==2
                exci(line).Pi(:, abs(mgrid.x) > aperSize(1)/2)=0;
            else
                exci(line).Pi(:, abs(mgrid.x) > aperSize(1)/2, :)=0;
                exci(line).Pi(:, :, abs(mgrid.y) > aperSize(2)/2)=0;
            end
        end

% apply dynamic range limit to excitation
        exci(line).Pi( abs(exci(line).Pi) < max(exci(line).Pi(:))*10.^(-xdcr.dynrange/20) )=0;
    end
end
    