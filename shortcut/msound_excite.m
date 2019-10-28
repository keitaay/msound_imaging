function exci=msound_excite(mgrid, medium, xdcr, focus)
% exci = MSOUND_EXCITE( mgrid, medium, xdcr, focus )
% 
% Using an mSOUND set_grid object "mgrid", a "medium" structure, a "xdcr"
% structure as generated from msound_xdcr.m, and a scalar focal depth
% "focus", generate another structure that encodes the excitation pulse to
% be emitted by an mSOUND forward simulation.

% 2019-10-26 - Keita Yokoyama (UNC/NCSU)
%              initial version

    exci.t=(-4/xdcr.fc : mgrid.dt : 4/xdcr.fc)'; %time vector

    % define impulse response function based on time delay "dly"
    impresp=@(dly) sin( 2*pi * xdcr.fc * (exci.t+dly) ) .* ...
                   exp(-(exci.t+dly).^2 * (xdcr.fc)^2 / 2); 

    % calculate element- and depth-wise time delays
    if isinf(focus)
        delay=zeros(1,mgrid.num_x);
    else
        delay=sqrt((mgrid.x).^2 + (xdcr.focus).^2)./medium.c0;
    end
    delay=repmat(delay-min(delay(:)), [length(exci.t) 1]);

    % convert impulse time vector to spatial matrix
    exci.t=repmat(exci.t, [1 mgrid.num_x]);

    % convert excitation to pressure waveform
    exci.P0=1e6;
    exci.Pi=exci.P0 .* impresp(delay);
    exci.Pi(:, abs(mgrid.x) > xdcr.dim(1)/2)=0; % need to rewrite for 3-D