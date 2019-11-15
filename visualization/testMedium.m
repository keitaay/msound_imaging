function testMedium(mgrid,medium,xdcr)
% TESTMEDIUM( mgrid, medium, xdcr )
% 
% Using an mSOUND set_grid object "mgrid" and structures "medium" and
% "xdcr" (see respective functions), generate figures that map out
% parameters that characterize the simulation field.
%
% Note that mSOUND uses five parameters (speed of sound, density,
% nonlinearity coefficient, attenuation coefficient, attenuation exponent);
% the sixth figure normalizes each of the five parameters, sums them
% together, and creates a composite image that depicts changes in ANY of
% the five values for each grid point.

% 2019-10-28 - Keita Yokoyama (UNC/NCSU)
%              initial version
% 2019-10-29 - Keita Yokoyama (UNC/NCSU)
%              added indicators for transducer position for
%              users with the Image Processing Toolbox installed
% 2019-11-14 - Keita Yokoyama (UNC/NCSU)
%              added imaging of absorption layer cross-section
    
% define number of dimensions in simulation
    nD=msound_nDim(mgrid);

% create universal plot information
    if nD==1
        error('This function does not support plots in 1D simulations.');
    elseif nD==2
        vec_lat=1000*mgrid.x;
        vec_axi=1000*(mgrid.y' - xdcr.plane_depth);
    elseif nD==3
        vec_lat=1000*mgrid.x;
        vec_axi=1000*(mgrid.z' - xdcr.plane_depth);
    end
    x_annot='Lateral (mm)';   y_annot='Axial (mm)';
   
    
% if image processing toolbox exists, prepare to draw image
    if license('test', 'image_toolbox')
        xdcLoc=@(ax,color) ...
            drawrectangle(ax,...
                          'Position',1000*[-xdcr.dim(1)/2, 0, xdcr.dim(1) xdcr.thickness],...
                          'Color',color,...
                          'LineWidth',5,...
                          'Deletable',false,'FaceAlpha',0.1,...
                          'InteractionsAllowed','none',...
                          'LabelVisible','off');
    else
        xdcLoc=@(ax,color) false;
    end
    
    figure('Position',[0 0 1500 600]);
% speed of sound (acoustic scatterers modeled here)
    subplot(2,5,1)
        imagesc(vec_lat, vec_axi, medium.c);        
        axis image; xlabel(x_annot); ylabel(y_annot); cb=colorbar;
        title('Speed of Sound');   colormap(gca,'parula'); ylabel(cb,'m/s')
        xdcLoc(gca,'red');
        
% density
    subplot(2,5,2)
        imagesc(vec_lat, vec_axi, medium.rho)
        axis image; xlabel(x_annot); ylabel(y_annot); cb=colorbar;
        title('Density');          colormap(gca,'cool');   ylabel(cb,'kg/m^{3}')
        xdcLoc(gca,'black');
        
% nonlinearity coefficient
    subplot(2,5,3)
        imagesc(vec_lat, vec_axi, medium.beta)
        axis image; xlabel(x_annot); ylabel(y_annot); cb=colorbar;
        title('Nonlinearity');     colormap(gca,'autumn'); ylabel(cb,'nonlinearity')
        xdcLoc(gca,'cyan');
        
% attenuation coefficient
    subplot(2,5,6)
        imagesc(vec_lat, vec_axi, medium.ca)
        axis image; xlabel(x_annot); ylabel(y_annot); cb=colorbar;
        title('Attenuation');      colormap(gca,'bone');   ylabel(cb,'dB/cmMHz')
        xdcLoc(gca,'yellow');
        
% power law exponent (2 = thermoviscous; no dispsersion)
    subplot(2,5,7)
        imagesc(vec_lat, vec_axi, medium.cb)
        axis image; xlabel(x_annot); ylabel(y_annot); cb=colorbar;
        title('Power Law Exponent');colormap(gca,'pink');  ylabel(cb,'Power Law Exp.')
        xdcLoc(gca,'magenta');
    
% composite of above five parameters
    subplot(2,5,8)
        Z=getZ(medium);
        imagesc(vec_lat, vec_axi, Z)
        axis image; xlabel(x_annot); ylabel(y_annot); cb=colorbar;
        title('Multiparametric');  colormap(gca,'parula');  ylabel(cb,'a.u.')
        xdcLoc(gca,'red');
        
        
    subplot(2,5,[4 5 9 10])
    % based on letter to editor in JASA 2012 vol 131, p999 (Yun Jing, NCSU)
    % "On the use of an absorption layer for the angular spectrum approach"
        if isfield(medium,'NRL_gamma')
            gamma = medium.NRL_gamma./...
                    (cosh(medium.NRL_alpha.*mgrid.abx_vec)).^2;
        else
            gamma=zeros(mgrid.num_x,1);
        end
        plot(vec_lat, gamma, 'linewidth',2); axis square
        xlabel(x_annot); ylabel('Absorption Magnitude (a.u.)');
        title('Absorption Strength');
        
end
function Z=getZ(med)
    % normalize each of the 5 features...
    nmlz_c=(med.c  - mean(med.c(:)))   ./ std(med.c(:));
    nmlz_r=(med.rho -mean(med.rho(:))) ./ std(med.rho(:));
    nmlz_b=(med.beta-mean(med.beta(:)))./ std(med.beta(:));
    nmlz_a=(med.ca - mean(med.ca(:)))  ./ std(med.ca(:));
    nmlz_p=(med.cb - mean(med.cb(:)))  ./ std(med.cb(:));

    % ...impute invalid values as needed...
    if isnan(sum(nmlz_c)),   nmlz_c(:,:)=0;  end
    if isnan(sum(nmlz_r)),   nmlz_r(:,:)=0;  end
    if isnan(sum(nmlz_b)),   nmlz_b(:,:)=0;  end
    if isnan(sum(nmlz_a)),   nmlz_a(:,:)=0;  end
    if isnan(sum(nmlz_p)),   nmlz_p(:,:)=0;  end

    % ...concatenate them together...
    Z(:,:,1)=nmlz_c;  Z(:,:,2)=nmlz_r;  Z(:,:,2)=nmlz_b;
    Z(:,:,2)=nmlz_a;  Z(:,:,2)=nmlz_p;
    
    % ...and sum
    Z=sum(abs(Z),3);
end