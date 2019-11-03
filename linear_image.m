close all; clear variables; clc; addpath(genpath(pwd));

%% grid
% wave/transducer definition
fc=4E6;   fs=fc*10;

% worldgrid
size_axi=0.02;  step_axi=5e-5; % direction of wave transmission
size_lat=0.02;  step_lat=1e-4;
size_ele=0.01;  step_ele=1e-4;
size_t = 4e-5;  step_t = 1/fs;

%% define useful variables
%mgrid=set_grid(step_t, size_t, step_axi, size_axi);
mgrid=set_grid(step_t, size_t, step_lat, size_lat, step_axi, size_axi);
%mgrid=set_grid(step_t, size_t, step_lat, size_lat, step_ele, size_ele, step_axi, size_axi);

medium=msound_medium(mgrid, fc, 'point');
%xdcr=msound_xdcr(mgrid, fc);
xdcr=msound_xdcr(mgrid, fc, [], [], [], [mgrid.dx 0], [0 0]);
exci=msound_excite(mgrid, medium, xdcr);

clearvars -except mgrid medium exci xdcr depth

%% simulation
reflOrder=1;

P=Forward2D(mgrid, medium, exci.Pi, xdcr.mask, reflOrder,...
            'NRL', 'correction', 'animation');

P=msound_reform(P, mgrid, xdcr); close all

%% beamforming
chanData=msound_get_chan(P, mgrid, xdcr);

[RF,axial]=msound_beamform(chanData,mgrid,medium,xdcr,exci);

%% Imaging
testMedium(mgrid,medium,xdcr)
msound_images(RF, axial, mgrid)