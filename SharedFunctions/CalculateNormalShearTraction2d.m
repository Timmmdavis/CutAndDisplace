function [ Tn,Ts ] = CalculateNormalShearTraction2d( Tx,Ty,CosAx,CosAy)
% CalculateNormalShearTraction2d: Calculates the normal and shear 
%                   traction on a 2D plane.
%                   The user just needs to supply the direction cosines of
%                   the plane and the Cartesian traction components at this
%                   plane. Use function
%                   "TractionVectorCartesianComponents2d.m" to find the
%                   Cartesian traction components if needed. 
%                   Based on Eq. 6.54. in Pollard, D.D. and Fletcher, R.C.,
%                   2005. Fundamentals of structural geology. Cambridge
%                   University Press.
%                   Positive shear traction is counter clockwise from the
%                   normal vector and normal traction is positive when
%                   facing in the direction of the normal vector (tension).
%   
% usage #1:
% [ Tn,Ts ] = CalculateNormalShearTraction2d( Tx,Ty,CosAx,CosAy)
%
% Arguments: (input)
% Tx,Ty             - Traction vector cartesian components
%
% CosAx,CosAy		- Direction cosines of plane. 
%
% Arguments: (output)
% Tn,Ts             - 3 column vectors of the traction magnitude in the 
%				      normal (nn), dip (ds) and strike (ss) direction.
%
% Example usage 1:
%
% % Calculating tractions for a 2D plane dipping 45 degrees under tension
% % in Sxx. The planes normal points upwards and faces in the positive X-axis. 
% CosAx=cosd(45);
% CosAy=cosd(45);
% Sxx=1; Syy=0; Sxy=0; 
% [ Tx,Ty ] = TractionVectorCartesianComponents2d( Sxx,Syy,Sxy,CosAx,CosAy );
% [ Tn,Ts ] = CalculateNormalShearTraction2d( Tx,Ty,CosAx,CosAy)
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Converts the traction XY to normal and shear traction components.
%Eq. 6.54. in Pollard, D.D. and Fletcher
Tn = (bsxfun(@times,Tx,CosAx))+(bsxfun(@times,Ty,CosAy));
Ts = (bsxfun(@times,-Tx,CosAy))+(bsxfun(@times,Ty,CosAx));

end

