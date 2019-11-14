close all; clear variables; clc; addpath(genpath(pwd));

%% grid
% wave/transducer definition
fc=6E6;   fs=fc*10;

% worldgrid
size_axi=0.03;   step_axi=1e-4; % direction of wave transmission
size_lat=0.02;   step_lat=1e-4;
size_ele=0.005;  step_ele=1e-4;
size_t = 4e-5;   step_t = 1/fs;

%% define useful variables
%mgrid=set_grid(step_t, size_t, step_axi, size_axi);
mgrid=set_grid(step_t, size_t, step_lat, size_lat, step_axi, size_axi);
%mgrid=set_grid(step_t, size_t, step_lat, size_lat, step_ele, size_ele, step_axi, size_axi);

medium=msound_medium(mgrid, fc, 'point');
xdcr=msound_xdcr(mgrid, fc);
%xdcr=msound_xdcr(mgrid, fc, [], [], [], [mgrid.dx 0], [0 0]);
exci=msound_excite(mgrid, medium, xdcr);

clearvars -except mgrid medium exci xdcr depth

%% simulation/beamforming
order=1;

Pcell=cell(length(exci),1);

for line=1:length(exci)
    
    P=Forward2D(mgrid,medium,exci(line).Pi,xdcr.mask,order,'NRL','correction');
    
    close all;clc;disp(['Completed line ',num2str(line),' of ',num2str(length(exci))]);
    
    Pcell{line}=msound_reform(P,mgrid,xdcr);
end

%% beamforming
chanData=msound_get_chan(Pcell, mgrid, xdcr);

[RF,axial]=msound_beamform(chanData,mgrid,medium,xdcr,exci);

%% Imaging
%testMedium(mgrid,medium,xdcr)
msound_images(RF, axial, mgrid)