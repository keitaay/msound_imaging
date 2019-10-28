close all; clear variables; clc; addpath(genpath(pwd));

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

%% beamforming
RF=msound_beamform(P,mgrid,medium,xdcr,exci);

%% Imaging
msound_images(RF, mgrid);