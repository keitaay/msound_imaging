close all; clear variables; clc; addpath(genpath(pwd));

%% grid
% wave/transducer definition
fc=4E6;   fs=fc*10;

% worldgrid
size_axi=0.02;   step_axi=1e-4; % direction of wave transmission
size_lat=0.016;  step_lat=1e-4;
size_ele=0.01;   step_ele=1e-4;
size_t = 5e-5;   step_t = 1/fs;

% line positions
lineL=( -3 : 0.25 : 3 ).*1e-3;
lineE=zeros(1,length(lineL));
angTh=zeros(1,length(lineL));
angPh=zeros(1,length(lineL));
lineloc=[lineL', lineE', angTh', angPh'];

%% define useful variables
%mgrid=set_grid(step_t, size_t, step_axi, size_axi);
mgrid=set_grid(step_t, size_t, step_lat, size_lat, step_axi, size_axi);
%mgrid=set_grid(step_t, size_t, step_lat, size_lat, step_ele, size_ele, step_axi, size_axi);

medium=msound_medium(mgrid, fc, 'point');
xdcr=msound_xdcr(mgrid, fc);
%xdcr=msound_xdcr(mgrid, fc, [], [], [], [mgrid.dx 0], [0 0]);

%exci=msound_excite(mgrid, medium, xdcr, lineloc, inf, [2 2]);  % plane
exci=msound_excite(mgrid, medium, xdcr, lineloc, size_axi/2, [2 2]); % focused

clearvars -except mgrid medium exci xdcr lineloc
%testMedium(mgrid,medium,xdcr)

%% simulation
order=1;

Pcell=cell(length(exci),1);
parpool(4);
parfor line=1:length(exci)
    P=Forward2D(mgrid, medium, exci(line).Pi, xdcr.mask, order, 'NRL', 'correction');%#ok<PFBNS>
    close all; clc; disp(['Completed line ',num2str(line),' of ',num2str(length(exci))]);%#ok<PFBNS>
    Pcell{line}=msound_reform(P,mgrid,xdcr);
end
delete(gcp('nocreate'));
chanData=msound_get_chan(Pcell, mgrid, xdcr);

save('simSample.mat', 'mgrid','medium','xdcr','exci','order','Pcell','chanData');
plotChanTraces( chanData, mgrid, medium, xdcr, exci, ceil(size(lineloc,1)/2), 100 )

%% beamforming
[RF,axial]=msound_beamform(chanData,mgrid,medium,xdcr,exci);

%% Imaging
msound_images(RF, axial, mgrid)