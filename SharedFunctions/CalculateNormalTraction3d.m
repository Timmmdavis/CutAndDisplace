function [ NormalTraction ] = CalculateNormalTraction3d( Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,CosAx,CosAy,CosAz )
% CalculateNormalTraction3d: Calculates normal traction on a given 3D plane
%                   with Cartesian stress tensor components defined.
%                   Based on Equation 6.49 in Pollard, D.D. and Fletcher,
%                   R.C., 2005. Fundamentals of structural geology.
%                   Cambridge University Press.
%               
% usage #1:
% [ NormalTraction ] = CalculateNormalTraction3d( Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,CosAx,CosAy,CosAz )
%
% Arguments: (input)
% Pxx,Pyy,Pzz 
% Pxy,Pxz,Pyz 		- The stress tensor components on this plane.
%                    (column vectors)
%
% CosAx,CosAy,CosAz - The seperated direction cosines of the surface 
%                    (column vectors)
%
% Arguments: (output)
% NormalTraction    - Column vector of normal traction magnitudes on this
%                     plane. Notation in scripts is typically: Tn.
%
% Example usage 1:
%
% % Get normal traction for a plane dipping 45 degrees facing east:
% CosAx=cosd(45);
% CosAy=0;
% CosAz=cosd(45);
% Sxx=1; Syy=0; Szz=0; Sxy=0; Sxz=0; Syz=0;
% [ NormalTraction ] = CalculateNormalTraction3d( Sxx,Syy,Szz,Sxy,Sxz,Syz,CosAx,CosAy,CosAz )
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


% Normal traction on the planes. 
% Equation 6.49 Pollard and Fletcher Fundamentals
AA=(bsxfun(@times,Pxx,CosAx.^2));
BB=(bsxfun(@times,Pyy,CosAy.^2));
CC=(bsxfun(@times,Pzz,CosAz.^2));
DD=(bsxfun(@times,2*Pxy,CosAx.*CosAy));
EE=(bsxfun(@times,2*Pyz,CosAy.*CosAz));
FF=(bsxfun(@times,2*Pxz,CosAz.*CosAx));
NormalTraction=AA+BB+CC+DD+EE+FF;

end

