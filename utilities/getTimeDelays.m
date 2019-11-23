function varargout=getTimeDelays( xdcr, linePos, axial, c, wavetype )
% delay=GETTIMEDELAYS( xdcr, linePos, axial, c, wavetype )
%
% [delay, chLat, chEle]=GETTIMEDELAYS( xdcr, linePos, axial, c, wavetype )
%
% Apply a standard (focused or plane) ultrasound beamformer to channel
% data, and get a matrix of relevant indices.
%
% REQUIRED INPUT:
%           xdcr = transducer settings, created using "msound_xdcr.m"
%
%        linePos = 1x4 matrix of the desired wavefront's position and
%                  orientation at its origin
%                  [m/deg] (1x4 vector - [lat,ele,theta,phi])
%
%          axial = vector of focal depths
%                  [m]
%
%              c = speed of sound
%                  [m/s]
%
%       wavetype = string to denote wavefront type; this function supports:
%                    - 'focused' = focused
%                    - 'plane'   = plane wave

% 2019-11-19 - Keita Yokoyama (UNC/NCSU)
%              initial version

% PLANE-WAVE ONLY: nullify axial values to prevent focusing
    if strcmp(wavetype,'plane')
        axial=zeros( length(axial), 1 );
    end

% convert line-focus positions/tilts to grid (Cartesian) system
    Eo=axial.*sind(linePos(3)).*cosd(-90+linePos(4));
    Lo=axial.*sind(linePos(3)).*sind(-90+linePos(4));
    A =axial.*cosd(linePos(3));
    
% shift lateral/elevational position of focus by location of beam origin
    Lo=repmat( Lo+linePos(1), [1,size(xdcr.chanID)] );
    Eo=repmat( Eo+linePos(2), [1,size(xdcr.chanID)] );

% initialize matrix with position of elements, by Cartesian dimensions
    Lx=nan([length(axial), size(xdcr.chanID)]);
    Ex=nan([length(axial), size(xdcr.chanID)]);
    for lat=1:size(Lx,2),  Lx(:,lat,:)=xdcr.chanLoc{lat,1}(1);  end
    for ele=1:size(Ex,3),  Ex(:,:,ele)=xdcr.chanLoc{1,ele}(2);  end
    
% shift channel positions based on line positions
    L=Lx+Lo;
    E=Ex+Eo;
    
% calculate radial distance of sound propagated from elements to focus
    R=sqrt( A.^2 + L.^2 + E.^2 );
    
% calculate time delays in each element for focusing
    delay=R./c;
    
% PLANE-WAVE ONLY: recover negative time-delays
    if strcmp(wavetype,'plane')
        makeNeg=( sign(L)==-1 | sign(E)==-1 );
        delay(~makeNeg)=-delay(~makeNeg);
    end
    
% adjust time delays to relative value (minimum per depth = 0)
    delay=delay-min(delay,[],[2 3]);
    
% format outputs
    for n=1:( max(nargout,1) )
        switch n
            case 1, varargout{1}=delay;
            case 2, varargout{2}=Lx;
            case 3, varargout{3}=Ex;
        end
    end
end