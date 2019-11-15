function medium=msound_medium(mgrid, fc, phtType)
% medium = MSOUND_MEDIUM( mgrid, fc, phtType )
%
% Create a structure that encodes information about the material
% where an mSOUND-simulated sound wave will propagate through.
%
% This function is only intended as a way to quickly create a "medium"
% structure. If you are aiming to develop a simulation environment with
% more complicated features, please create your structure manually
% or only use this function as a first step (to easily generate scatterer
% fields).
%
% REQUIRED INPUT:
%          mgrid = mSOUND set_grid object
%             fc = center frequency of transducer
%                  [Hz] (scalar)
%        phtType = type of digital phantom to be generated; choose from:
%                  - point (single point target in center of grid space)
%                  - block (small block of scatterers in grid space)
%                  - field (entire simulation area contains scatterers)

% 2019-10-26 - Keita Yokoyama (UNC/NCSU)
%              initial version
% 2019-10-29 - Keita Yokoyama (UNC/NCSU)
%              generalized function for 1/2/3-D environments;
%              added more helpful help texts

%% Initialization
    % define number of dimensions in simulation
    nD=msound_nDim(mgrid);
    
    % environmental settings (parameter + std deviation)
    c0=1540;    d_c=20;
    rho0=1000;
    beta0=0;
    ca0=0;
    cb0=2;
    
    % non-reflecting layer settings
    NRLgamma=1;
    NRLalpha=0.2;
    
    % define template for matrices of material properties
    if nD==1, map_base=zeros(mgrid.num_x+1, 1); end
    if nD==2, map_base=zeros(mgrid.num_x, mgrid.num_y+1); end
    if nD==3, map_base=zeros(mgrid.num_x, mgrid.num_y,   mgrid.num_z+1); end
    
    % create map of all parameters
    map_c    =  c0   + map_base;   % speed of sound
    map_rho  = rho0  + map_base;   % density
    map_beta = beta0 + map_base;   % nonlinearity coefficient
    map_ca   =  ca0  + map_base;   % attenuation coefficient
    map_cb   =  cb0  + map_base;   % power law exponent
    
%% Scatterer Placement
    % get resolution cell size + number of scatterers in simulation
    n_inCell=11;        % usually, 11 is enough for speckle SNR = 1.91
    
    % use scatterer density per res cell to find overall scatterer count
    res_axi=c0/fc;      % axial res needed regardless of sim dimensionality
    N_axi=round(mgrid.y_length/res_axi);
    switch nD
        case 1
            Nsctr=n_inCell*N_axi;
        case 2
            res_lat=c0/fc;   N_lat=round(mgrid.x_length/res_lat);
            Nsctr=n_inCell*N_axi*N_lat;
        case 3
            res_lat=c0/fc;   N_lat=round(mgrid.x_length/res_lat);
            res_ele=c0/fc;   N_ele=round(mgrid.y_length/res_ele);
            Nsctr=n_inCell*N_axi*N_lat*N_ele;
    end
    
    % define scatterer positions
    seedPos=0;    rng(seedPos,'twister');
    switch nD
        case 1
            sctrPos_1D=randi( uint64(mgrid.num_x),[Nsctr,1] );
        case 2
            sctrPos_1D=randi( uint64(mgrid.num_x*mgrid.num_y),[Nsctr,1] );
        case 3
            sctrPos_1D=randi( uint64(mgrid.num_x*mgrid.num_y*mgrid.num_z),[Nsctr,1] );
    end
    
    % define scatterer amplitudes
    seedAmp=0;    rng(seedAmp,'twister');
    sctrAmp_1D=randn( Nsctr, 1);
    
    % define scatterers as speed-of-sound inhomogeneities
    switch phtType
        case {'field','block'}
            % define scatterer field as column vector...
            fields_c=reshape( map_c, [],1);
            fields_c(sctrPos_1D)=d_c*sctrAmp_1D + fields_c(sctrPos_1D);
            
            % ...then revert into matrix form, as appropriate
            switch nD
                case 1
                    fields_c=reshape(fields_c,...
                              mgrid.num_x+1, 1);
                case 2
                    fields_c=reshape(fields_c,...
                              mgrid.num_x, mgrid.num_y+1);
                case 3
                    fields_c=reshape(fields_c,...
                              mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
            end

            if strcmp(phtType,'block')
                % create binary mask for localized area with scatterers
                height=0.006;    width=0.005;    depth=0.005;
                sctr_present=true(size(map_base));
                switch nD
                    case 1
                        % define region to create scatterers
                        sctr_axi=[1 : knnsearch(mgrid.x',-height/2),...
                                      knnsearch(mgrid.x', height/2) : mgrid.num_x];
                        
                        % impose binary mask onto scatterer field
                        sctr_present(sctr_axi)=false;
                        fields_c(~sctr_present)=c0;
                        
                    case 2
                        % define region to create scatterers
                        sctr_axi=[1 : knnsearch(mgrid.y',-height/2),...
                                      knnsearch(mgrid.y', height/2) : mgrid.num_y];
                        sctr_lat=[1 : knnsearch(mgrid.x',-width/2),...
                                      knnsearch(mgrid.x', width/2)  : mgrid.num_x];
                        
                        % impose binary mask onto scatterer field
                        sctr_present(:,sctr_axi)=false;
                        sctr_present(sctr_lat,:)=false;
                        fields_c(~sctr_present)=c0;
                        
                    case 3
                        % define region to create scatterers
                        sctr_axi=[1 : knnsearch(mgrid.z',-height/2),...
                                      knnsearch(mgrid.z', height/2) : mgrid.num_z];
                        sctr_lat=[1 : knnsearch(mgrid.x',-width/2),...
                                      knnsearch(mgrid.x', width/2)  : mgrid.num_x];
                        sctr_ele=[1 : knnsearch(mgrid.y',-depth/2),...
                                      knnsearch(mgrid.y', depth/2)  : mgrid.num_y];
                        
                        % impose binary mask onto scatterer field
                        sctr_present(:,:,sctr_axi)=false;
                        sctr_present(sctr_lat,:,:)=false;
                        sctr_present(:,sctr_ele,:)=false;
                        fields_c(~sctr_present)=c0;
                end
            end
            
        case 'point'
            fields_c=map_c;
            c_point=1600;   % speed of sound at point target
            switch nD
                case 1
                    fields_c(knnsearch(mgrid.x',0) )=c_point;
                case 2
                    fields_c(knnsearch(mgrid.x',0),...
                             knnsearch(mgrid.y',0) )=c_point;
                case 3
                    fields_c(knnsearch(mgrid.x',0),...
                             knnsearch(mgrid.y',0),...
                             knnsearch(mgrid.z',0) )=c_point;
            end
    end
    
    % consolidate variables into "medium" structure
    medium.c=fields_c;          medium.c0=min(medium.c(:));
    medium.rho=map_rho;         medium.beta=map_beta;
    medium.ca=map_ca;           medium.cb=map_cb;
    medium.NRL_gamma=NRLgamma;  medium.NRL_alpha=NRLalpha;
end