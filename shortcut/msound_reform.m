function P=msound_reform(P, mgrid, xdcr)
% P = MSOUND_REFORM( P, mgrid, xdcr )
% 
% Re-structure the output from an mSOUND forward simulation into dimensions
% that are more intuitive.
%
% Regardless of the simulation type used, the output "P" will ALWAYS have
% dimensions of:
%                  time x axial x lateral x elevational
%
% REQUIRED INPUT:
%              P = output from an mSOUND forward simulation
%          mgrid = mSOUND set_grid object
%           xdcr = transducer settings, created using "msound_xdcr.m"

% 2019-10-26 - Keita Yokoyama (UNC/NCSU)
%              initial version

    % define number of dimensions in simulation
    nD=msound_nDim(mgrid);
    
    % reshape pressure readings into [time x spatial ... x axial]
    switch nD
        case 1
            P=reshape(P, [ mgrid.num_t,...
                           round(xdcr.thickness/mgrid.dx) ]);
            P=permute(P, [1,2,3,4]);
            
        case 2
            P=reshape(P, [ mgrid.num_t,...
                           round(xdcr.dim(1)/mgrid.dx)+1,...
                           round(xdcr.thickness/mgrid.dy) ]);
            P=permute(P, [1,3,2,4]);
            
        case 3
            P=reshape(P, [ mgrid.num_t,...
                           round(xdcr.dim(1)/mgrid.dx)+1,...
                           round(xdcr.dim(2)/mgrid.dy)+1,...
                           round(xdcr.thickness/mgrid.dz) ]);
            P=permute(P, [1,4,2,3]);
    end
end