function m=modvec(x,y,z,lattice)
%===================================================================================
% function m=modvec(x,y,z,lattice)
%  ResLib v.3.1
%===================================================================================
%
%  Calculates the modulus of a vector, defined by its fractional cell coordinates or
%  Miller indexes.
%
% A. Zheludev, 1999-2001
% Oak Ridge National Laboratory
%====================================================================================

m=sqrt( scalar(x,y,z,x,y,z,lattice) );

