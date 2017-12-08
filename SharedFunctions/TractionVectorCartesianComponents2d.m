function [ Tx,Ty ] = TractionVectorCartesianComponents2d( Pxx,Pyy,Pxy,CosAx,CosAy )
% TractionVectorCartesianComponents2d: Calculates the Cartesian components
%                   of the traction vector for 2D data from the tensors and
%                   direction cosines.
%                   Equation 6.40 and 6.41 in Pollard and Fletcher, 2005. 
%   
% usage #1:
% [ Tx,Ty ] = TractionVectorCartesianComponents2d( Pxx,Pyy,Pxy,CosAx,CosAy )
%
% Arguments: (input)
% Pxx,Pyy,Pxy       - The stress tensors (col vects).
%
% CosAx,CosAy       - Direction cosines of the plane at the tensor
%                     locations (planes normal vector). 
%
% Arguments: (output)
% Tx,Ty             - Cartesian components of the 2D traction vector
%
% Example usage:
%
% %Cartesian components on a dipping plane (45 deg) subject to a
% %extensional stress.
% Pxx=1; Pyy=0; Pxy=0;
% CosAx=deg2rad(45);
% CosAy=deg2rad(45);
% [ Tx,Ty ] = TractionVectorCartesianComponents2d( Pxx,Pyy,Pxy,CosAx,CosAy )
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Tx(n) = [sxx sxy] * [nx]
%Ty(n) = [sxy syy] * [ny]
Tx = (bsxfun(@times,Pxx,CosAx))+(bsxfun(@times,Pxy,CosAy));
Ty = (bsxfun(@times,Pxy,CosAx))+(bsxfun(@times,Pyy,CosAy));

end

