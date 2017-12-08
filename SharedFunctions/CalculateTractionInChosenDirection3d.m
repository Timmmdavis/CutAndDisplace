function [ TractionInDirection ] = CalculateTractionInChosenDirection3d( Tx,Ty,Tz,CosAx,CosAy,CosAz,ChosenDirectionCos )
% CalculateTractionInChosenDirection3d: Calculates traction on a plane in 
%                   any direction the user chooses. The user just needs to
%                   supply the direction cosines of the direction they want
%                   this in and the Cartesian traction components on this
%                   plane. If needed these can be calculated from a full
%                   stress tensor on the plane using function:
%                   "TractionVectorCartesianComponents3d.m".
%                   Based on Equation 6.52 in Pollard, D.D. and Fletcher,
%                   R.C., 2005. Fundamentals of structural geology.
%                   Cambridge University Press.
%               
% usage #1:
% [ TractionInDirection ] = CalculateTractionInChosenDirection3d( Tx,Ty,Tz,CosAx,CosAy,CosAz,ChosenDirectionCos )
%
% Arguments: (input)
% Tx,Ty,Tz          - The Cartesian traction components on the plane.
%
% CosAx,CosAy,CosAz - The seperated direction cosines of the surface 
%                    (column vectors)
%
% ChosenDirectionCos- 3*n array of the direction the user desired the
%                    traction in. This is direction cosines
%                    [CosAx,CosAy,CosAz] of a vector that faces in any
%                    direction.
%
% Arguments: (output)
% TractionInDirection- Column vector of traction magnitudes in the chosen
%                      direction.
%
% Example usage 1:
%
% % Calc dip slip traction for a plane dipping 45 degrees facing east:
% CosAx=cosd(45);
% CosAy=0;
% CosAz=cosd(45);
% Sxx=1; Syy=0; Szz=0; Sxy=0; Sxz=0; Syz=0;
% [ Tx,Ty,Tz ] = TractionVectorCartesianComponents3d(  Sxx,Syy,Szz,Sxy,Sxz,Syz,CosAx,CosAy,CosAz );
% [ StrikeSlipCosine,DipSlipCosine ] = CalculateDSandSSDirs( [CosAx,CosAy,CosAz] ) ;  
% [ DipSlipTraction ] = CalculateTractionInChosenDirection3d( Tx,Ty,Tz,CosAx,CosAy,CosAz,DipSlipCosine )
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Split up the equation into seperate parts (each line in the book):
A=(bsxfun(@times,(bsxfun(@times,Tx,(1-(CosAx.^2))))-(bsxfun(@times,Ty,CosAx.*CosAy))  -(bsxfun(@times,Tz,CosAx.*CosAz))  ,ChosenDirectionCos(:,1)));
B=(bsxfun(@times,(bsxfun(@times,Tx,-CosAx.*CosAy)) +(bsxfun(@times,Ty,(1-(CosAy.^2))))-(bsxfun(@times,Tz,CosAy.*CosAz))  ,ChosenDirectionCos(:,2)));
C=(bsxfun(@times,(bsxfun(@times,Tx,-CosAz.*CosAx)) -(bsxfun(@times,Ty,CosAz.*CosAy))  +(bsxfun(@times,Tz,(1-(CosAz.^2)))),ChosenDirectionCos(:,3)));
%Sum the parts of the equation
TractionInDirection=A+B+C; 

end

