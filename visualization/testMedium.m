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

    % create universal plot information
    vec_lat=1000*mgrid.x;
    vec_axi=1000*(mgrid.y' - mgrid.y(1) - xdcr.plane_depth);
    x_annot='Lateral (mm)';   y_annot='Axial (mm)';
    
    figure;
% speed of sound (acoustic scatterers modeled here)
    subplot(2,3,1)
        imagesc(vec_lat, vec_axi, medium.c);        
        axis image; xlabel(x_annot); ylabel(y_annot); cb=colorbar;
        title('Speed of Sound');   colormap(gca,'parula'); ylabel(cb,'m/s')
        
% density
    subplot(2,3,2)
        imagesc(vec_lat, vec_axi, medium.rho)
        axis image; xlabel(x_annot); ylabel(y_annot); cb=colorbar;
        title('Density');          colormap(gca,'cool');   ylabel(cb,'kg/m^{3}')
        
% nonlinearity coefficient
    subplot(2,3,3)
        imagesc(vec_lat, vec_axi, medium.beta)
        axis image; xlabel(x_annot); ylabel(y_annot); cb=colorbar;
        title('Nonlinearity');     colormap(gca,'autumn'); ylabel(cb,'nonlinearity')
        
% attenuation coefficient
    subplot(2,3,4)
        imagesc(vec_lat, vec_axi, medium.ca)
        axis image; xlabel(x_annot); ylabel(y_annot); cb=colorbar;
        title('Attenuation');      colormap(gca,'bone');   ylabel(cb,'dB/cmMHz')
        
% power law exponent (2 = thermoviscous; no dispsersion)
    subplot(2,3,5)
        imagesc(vec_lat, vec_axi, medium.cb)
        axis image; xlabel(x_annot); ylabel(y_annot); cb=colorbar;
        title('Power Law Exponent');colormap(gca,'pink');  ylabel(cb,'Power Law Exp.')
    
    subplot(2,3,6)
        % normalize each of the 5 features...
        nmlz_c=(medium.c  - mean(medium.c(:)))   ./ std(medium.c(:));
        nmlz_r=(medium.rho -mean(medium.rho(:))) ./ std(medium.rho(:));
        nmlz_b=(medium.beta-mean(medium.beta(:)))./ std(medium.beta(:));
        nmlz_a=(medium.ca - mean(medium.ca(:)))  ./ std(medium.ca(:));
        nmlz_p=(medium.cb - mean(medium.cb(:)))  ./ std(medium.cb(:));
        
        % ...impute invalid values as needed...
        if isnan(sum(nmlz_c)),   nmlz_c(:,:)=0;  end
        if isnan(sum(nmlz_r)),   nmlz_r(:,:)=0;  end
        if isnan(sum(nmlz_b)),   nmlz_b(:,:)=0;  end
        if isnan(sum(nmlz_a)),   nmlz_a(:,:)=0;  end
        if isnan(sum(nmlz_p)),   nmlz_p(:,:)=0;  end
        
        % ...concatenate them together...
        Z(:,:,1)=nmlz_c;  Z(:,:,2)=nmlz_r;  Z(:,:,2)=nmlz_b;
        Z(:,:,2)=nmlz_a;  Z(:,:,2)=nmlz_p;
        
        % ...and plot
        imagesc(vec_lat, vec_axi, sum(abs(Z),3))
        axis image; xlabel(x_annot); ylabel(y_annot)
        title('Multiparametric');  colormap(gca,'parula')
end