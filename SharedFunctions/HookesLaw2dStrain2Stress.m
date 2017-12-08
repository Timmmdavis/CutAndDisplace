function [ Sxx,Syy,Sxy ] = HookesLaw2dStrain2Stress( Exx,Eyy,Exy,E,nu,mu )
% HookesLaw2dStrain2Stress: Using Plane strain 2d Hooke's law
%                   Converting the strain tensors to stress tensors using
%                   the Young's modulus, shear modulus and Poisson's ratio.
%                   Works with column vectors or single values.
%
%                   Equation from Pollard and Fletcher, 2005, Eq. 8.34 
%
%               
% usage #1:
% [ Sxx,Syy,Sxy ] = HookesLaw2dStrain2Stress( Exx,Eyy,Exy,E,nu,mu )
%
% Arguments: (input)
% Exx,Eyy,Exy       - Strain tensor components. (Can be column vectors). 
%
% E                 - Young's modulus
%
% nu                - The Poisson's ratio.
%
% mu                - The shear modulus
%
% Arguments: (output)
% Sxx,Syy,Sxy       - Stress tensor components. 
%
% Example usage:
%
% mu=5; %[Gpa]
% nu=0.2;
% E = mu*(2*(1+nu)) ; %Young's Mod
% Exx=0.2;
% Eyy=0;
% Exy=0;
% [ Sxx,Syy,Sxy ] = HookesLaw2dStrain2Stress( Exx,Eyy,Exy,E,nu,mu )
%
% Additional elastic constant conversions if needed:
% nu =lambda/(2*(mu+lamda);         Equation 8.28 Pollard       
% mu = E/(2*(1+nu));                Equation 8.26 Pollard
% E = mu*(2*(1+nu)) ;               Equation 8.26 Pollard (rearranged)
% lambda= E*nu/((1+nu)*(1-2*nu));   Equation 8.27 Pollard
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Lame's  constant,  Eq. 8.27 Pollard
lambda= E*nu/((1+nu)*(1-2*nu));   
%Hooke's law Plane Strain Eq. 8.34 
Sxx=(2*mu+lambda)*Exx+lambda*Eyy;
Syy=(2*mu+lambda)*Eyy+lambda*Exx;	
clear Exx Eyy
%We assume here Exy is not engineering strain.
Sxy= 2*mu*Exy; 
	

end

