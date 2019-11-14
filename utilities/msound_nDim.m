function nD=msound_nDim(mgrid)
% define number of dimensions in simulation
% nD = MSOUND_NDIM( mgrid )
% 
% Using an mSOUND set_grid object "mgrid", identify the number of spatial
% dimensions to be simulated.
%
% This function, itself, is very simple, but it becomes very helpful to
% reduce numbers of lines of code written for more complicated functions.

% 2019-10-28 - Keita Yokoyama (UNC/NCSU)
%              initial version

    if     mgrid.num_y==0,  nD=1;
    elseif mgrid.num_z==0,  nD=2;
    else,                   nD=3;
    end
end