function [ Tx,Ty,Tz ] = TractionVectorCartesianComponents3d(Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,CosAx,CosAy,CosAz)
% TractionVectorCartesianComponents3d: Calculates the Cartesian components
%                   of the traction vector for 3D data from the tensors and
%                   direction cosines.
%                   Equation 6.40 and 6.41 in Pollard and Fletcher, 2005. 
%   
% usage #1:
% [ Tx,Ty,Tz ] = TractionVectorCartesianComponents3d(Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,CosAx,CosAy,CosAz)
%
% Arguments: (input)
% Pxx,Pyy,Pzz
% Pxy,Pxz,Pyz       - The stress tensors (col vects).
%
% CosAx,CosAy,CosAz - Direction cosines of the plane at the tensor
%                     locations (planes normal vector). 
%
% Arguments: (output)
% Tx,Ty,Tz          - Cartesian components of the 3D traction vector
%
% Example usage:
% 
% %Cartesian components on a dipping plane (45 deg) subject to a
% %extensional stress.
% Pxx=1; Pyy=0; Pzz=0;
% Pxy=0; Pxz=0; Pyz=0;
% CosAx=deg2rad(45);
% CosAy=0;
% CosAz=deg2rad(45);
% [ Tx,Ty,Tz ] = TractionVectorCartesianComponents3d(Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,CosAx,CosAy,CosAz)
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen
	
%Tx(n) [sxx sxy sxz] * [nx]
%Ty(n) [sxy syy syz] * [ny]
%Tz(n) [sxz syz szz] * [nz]
Tx = (bsxfun(@times,Pxx,CosAx))+(bsxfun(@times,Pxy,CosAy))+(bsxfun(@times,Pxz,CosAz)); 
Ty = (bsxfun(@times,Pxy,CosAx))+(bsxfun(@times,Pyy,CosAy))+(bsxfun(@times,Pyz,CosAz));
Tz = (bsxfun(@times,Pxz,CosAx))+(bsxfun(@times,Pyz,CosAy))+(bsxfun(@times,Pzz,CosAz));

end

