function plotChanTraces( chanData, mgrid, medium, xdcr, exci, lineID, scaler )
% PLOTCHANTRACES(chanData,mgrid,medium,xdcr,exci,lineID,scaler)
% 
% Plot the raw, pre-beamformed voltage traces received in each element of
% an mSOUND-simulated transducer, as well as maps of corresponding
% channels (including apodization weights and transmit pres. waveforms).
%
% REQUIRED INPUT:
%       chanData = mSOUND result, re-organized as channel data for the
%                  transducer "xdcr" using "msound_get_chan.m"
%          mgrid = mSOUND set_grid object
%         medium = mSOUND medium structure
%           xdcr = transducer settings, created using "msound_xdcr.m"
%           exci = excitation settings, created using "msound_excite.m"
%         lineID = index for line (4th dim of "chanData") to plot
%         scaler = magnify amplitudes of channel data by this scalar value

% 2019-11-15 - Keita Yokoyama (UNC/NCSU)
%              initial version

eleSlice=1; % find better way to slice "chanData" elevationally in the future

% initialize label vector for lateral position, with respect to
% element closest to position of inputted line
    latChan=zeros(xdcr.Nelem(1), 1);
    eleChan=zeros(xdcr.Nelem(2), 1);
    for chan=1:xdcr.Nelem(1)
        latChan( chan )=xdcr.chanLoc{chan,eleSlice}(1);
    end
    for chan=1:xdcr.Nelem(2)
        eleChan( chan )=xdcr.chanLoc{1,chan}(2);
    end
    latChan=1000*latChan;  eleChan=1000*eleChan;
    
% initialize label vectors
    axial=1000*medium.c0.*mgrid.t';  t=1e6*mgrid.t;

% isolate relevant channel waveforms into 2-D matrix
    chan2D=squeeze( chanData(:,:,eleSlice,lineID) );
    
% extract transmit pressure/time vector from excitation structure
    Pout=exci(lineID).Pchan;    tout=1e6*exci(lineID).t;
    
% extract transmit/receive apodization from excitation structure
    apodTx=exci(lineID).apodTx; apodRx=exci(lineID).apodRx;
    
% get relevant spatial vectors + channel/apodization maps
    switch msound_nDim(mgrid)
        case 1, latGrid=0;             eleGrid=0;
        case 2, latGrid=1000*mgrid.x'; eleGrid=0;
        case 3, latGrid=1000*mgrid.x'; eleGrid=1000*mgrid.y';
    end
    chanmap=xdcr.chanmap;
    
% configure figure
    figure('Position',[0 0 1200 800])
    
% display diagnostic/supplemental images
    subplot(4,2,1)       % channel map
        imagesc(latGrid,eleGrid,chanmap');
        colormap(gca,'copper');    caxis([0 1])
        xlabel('Lateral (mm)');    ylabel('Elevational (mm)')
        title('Channel Locations')
    
    subplot(4,2,3)       % transmit waveform
        wavemap=[zeros(1,41),      linspace(1,0,20),        ones(1,20);...
                 linspace(0.5,0,20),  zeros(1,41),  linspace(0,0.5,20);...
                 ones(1,20),       linspace(1,0,20),      zeros(1,41)]';
        imagesc(latChan,tout,Pout);
        colormap(gca,wavemap)
        xlabel('Lateral (mm)');    ylabel('Time (us)')
        title('Transmitted Pressure Waves')
        ax=gca; ax.YDir='normal'; ax.YTickLabel=fliplr(ax.YTickLabel);
    
    subplot(4,2,5)       % transmit apodization
        imagesc(latChan,eleChan,apodTx');
        colormap(gca,'bone');     caxis([0 1])
        xlabel('Lateral (mm)');   ylabel('Elevational (mm)')
        title('Transmit Apodization')
    
    subplot(4,2,7)       % receive apodization
        imagesc(latChan,eleChan,apodRx');
        colormap(gca,'hot');      caxis([0 1])
        xlabel('Lateral (mm)');   ylabel('Elevational (mm)')
        title('Receive Apodization')
    
% display voltage trace
    subplot(4,2,[2 4 6 8]); hold on
    for chan=1:xdcr.Nelem(1)
        
        % get voltage trace of channel (ampl. is peak-normalized)
        chanTrace=chan2D(:,chan)./max(abs(chan2D(:)));
        
        % scale amplitude of plotted voltage traces by this amount
        chanTrace=scaler*chanTrace;
        
        % plot voltage trace
        plot( chan+chanTrace, t, 'k-');
    end
    
% set axis labels for raw data
    ax1=gca;
        ax1.XAxisLocation='top';
        ax1.XLabel.String='Channel ID';
        ax1.YLabel.String='Time of Flight (us)';
        ax1.XLim=[0.2 size(chan2D,2)+0.8];%[1-scaler size(chan2D,2)+scaler];
        ax1.YLim=[t(1) t(end)];
        ax1.YDir='reverse';
    
% set axis labels for spatial reference
    ax2=axes('Position',ax1.Position, 'Color','none',...
             'XAxisLocation','bottom',...
             'YAxisLocation','right');
        ax2.XLabel.String='Position with respect to Center of Line (mm)';
        ax2.YLabel.String='Distance of Echo Propagation (mm)';
        ax2.XLim=[latChan(1) latChan(end)];
        ax2.YLim=[axial(1) axial(end)];
        ax2.YDir='reverse';

% set contextual image parameters
    title({[  'Line #',num2str(lineID),...
              ', Elevational Row #',num2str(eleSlice),...
              ' - Amplitude x',num2str(scaler)],...
           '',''}); % extra rows are to ensure space for 2nd x-label
end