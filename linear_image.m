close all; clear variables; clc
addpath(genpath('C:\Users\kayokoya\Documents\mSOUND'));

%% grid
% wave/transducer definition
fc=4E6;   focus=0.01;   fs=fc*10;

% worldgrid (wave propagates in Y for 1/2D; Z for 3D)
simSize=[0.02  0.03  0.01  3e-5];
simStep=[1e-4  1e-4  1e-4  1/fs];

%% define useful variables
mgrid=set_grid(simStep(4),simSize(4),... % time
               simStep(1),simSize(1),... % lateral
               simStep(2),simSize(2));   % axial

medium=msound_medium(mgrid, fc, 'point');
xdcr=msound_xdcr(mgrid, fc, focus);
exci=msound_excite(mgrid, medium, xdcr, inf);
%testMedium(mgrid,medium,xdcr)

clearvars -except mgrid medium exci xdcr depth

%% simulation
reflOrder=1;

P=Forward2D(mgrid, medium, exci.Pi, xdcr.mask, reflOrder,...
            'NRL', 'animation', 'correction');

RF=msound_beamform(P,mgrid,medium,xdcr,exci);

%% result
bimg=20*log10(abs(hilbert(RF)));  x=1000*mgrid.x;  y=1e6*mgrid.t;
close all; figure

subplot(1,2,1); imagesc(x, y, RF);   title('RF'); caxis([-1e4 1e4])
    axis image; colormap(gca,'gray'); xlabel('Lateral (mm)');
    
subplot(1,2,2); imagesc(x, y, bimg); title('B-mode')
    axis image; colormap(gca,'gray'); xlabel('Lateral (mm)');