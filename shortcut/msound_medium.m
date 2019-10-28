function medium=msound_medium(mgrid, fc, phtType)
% medium = MSOUND_MEDIUM( mgrid, fc, phtType )
% 
% Using an mSOUND set_grid object "mgrid", and transducer center frequency
% "fc", create a structure that encodes mSOUND-formatted information about
% the material where a simulated sound wave will propagate through.
%
% Then phantom type should be denoted with a string "phtType" out of one 
% of the following:
%   - point  (single point is generated in the center of the sim area)
%   - block  (a block within the simulation area contains scatterers)
%   - field  (entire simulation area contains scatterers)

% 2019-10-26 - Keita Yokoyama (UNC/NCSU)
%              initial version

    % define template
    map_base=zeros(mgrid.num_x, mgrid.num_y+1);
    
    % environmental settings (parameter + std deviation)
    c0=1540;    d_c=20;
    rho0=1000;  %d_rho=
    beta0=0;    %d_beta=
    ca0=0;      %d_ca=
    cb0=2;      %d_cb=
    
    % get resolution cell size + number of scatterers in simulation
    n_inCell=11;
    Fnum=[1.5 1.5];
    res_axi=c0/fc;           N_axi=round(mgrid.y_length/res_axi);
    res_lat=c0/fc*Fnum(1);   N_lat=round(mgrid.x_length/res_lat);
    %res_ele=c0/fc*Fnum(2);   N_ele=round(mgrid.y_length/res_ele);
    Nsctr=n_inCell*N_axi*N_lat;
    
    % create map of all parameters
    map_c    =  c0   + map_base;   % speed of sound
    map_rho  = rho0  + map_base;   % density
    map_beta = beta0 + map_base;   % nonlinearity coefficient
    map_ca   =  ca0  + map_base;   % attenuation coefficient
    map_cb   =  cb0  + map_base;   % power law exponent
    
    % define scatterer positions
    seedPos=0;    rng(seedPos,'twister');
    sctrPos_1D=randi( uint64(mgrid.num_x*mgrid.num_y),[Nsctr,1] );
    
    % define scatterer amplitudes
    seedAmp=0;    rng(seedAmp,'twister');
    sctrAmp_1D=randn( Nsctr, 1);
    
    switch phtType
        case {'field','block'}
            % define scatterers as sub-resolution changes in speed of sound
            fields_c=reshape( map_c, [],1);
            fields_c(sctrPos_1D)=d_c*sctrAmp_1D + fields_c(sctrPos_1D);
            fields_c=reshape(fields_c, mgrid.num_x, mgrid.num_y+1);

            if strcmp(phtType,'block')
                sctr_axi=[1:knnsearch(mgrid.y',-0.003),...
                          knnsearch(mgrid.y',0.003):mgrid.num_y];
                sctr_lat=[1:knnsearch(mgrid.x',-0.003),...
                          knnsearch(mgrid.x',0.003):mgrid.num_x];
                
                % create binary mask
                X=true(size(map_base));
                X(:,sctr_axi)=false;   X(sctr_lat,:)=false;

                % impose binary mask onto scatterer field
                fields_c(~X)=c0;
            end
            
        case 'point'
            fields_c=map_c;
            fields_c( knnsearch(mgrid.x',0) ,knnsearch(mgrid.y',0) )=1600;
    end
    
    % consolidate variables into "medium" structure
    medium.c=fields_c;     medium.c0=min(medium.c(:));
    medium.rho=map_rho;    medium.beta=map_beta;
    medium.ca=map_ca;      medium.cb=map_cb;
    medium.NRL_gamma=0.2;  medium.NRL_alpha=0.02;
end