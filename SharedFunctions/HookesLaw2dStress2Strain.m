function [ Exx,Eyy,Exy ] = HookesLaw2dStress2Strain( Sxx,Syy,Sxy,E,nu )
% HookesLaw2dStress2Strain: Using Plane strain 2d Hooke's law
%                   Converting the stress tensors to strain tensors using
%                   Lame's constant and the shear modulus.
%                   Works with column vectors or single values.
%
%                   Equation from Pollard and Fletcher, 2005, Eq. 8.48
%
%               
% usage #1:
% [ Exx,Eyy,Exy ] = HookesLaw2dStress2Strain( Sxx,Syy,Sxy,E,nu )
%
% Arguments: (input)
% Sxx,Syy,Sxy       - Stress tensor components. (Can be column vectors). 
%
% E                 - Young's modulus
%
% nu                - The Poisson's ratio.
%
% Arguments: (output)
% Exx,Eyy,Exy       - Strain tensor components. 
%
% Example usage:
%
% mu=5; %[Gpa]
% nu=0.2;
% E = mu*(2*(1+nu)) ; %Young's Mod
% Sxx=0.2;
% Syy=0;
% Sxy=0;
% [ Exx,Eyy,Exy ] = HookesLaw2dStress2Strain( Sxx,Syy,Sxy,E,nu )
%
% Additional elastic constant conversions if needed:
% nu =lambda/(2*(mu+lamda);         Equation 8.28 Pollard       
% mu = E/(2*(1+nu));                Equation 8.26 Pollard
% E = mu*(2*(1+nu)) ;               Equation 8.26 Pollard (rearranged)
% lambda= E*nu/((1+nu)*(1-2*nu));   Equation 8.27 Pollard
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Old function took in lambda and mu as the input args, for now I have this
%check so when using old scripts you are alerted. I will remove later.
if inputname(4)=='lambda'
    error('Input arguments for this function have changed')
end
	  
%Equation 8.48, Pollard and fletcher. 
Exx = 1./E.*(Sxx.*(1-nu.^2)-(Syy.*nu.*(1+nu)));
Eyy = 1./E.*(Syy.*(1-nu.^2)-(Sxx.*nu.*(1+nu)));
clear Sxx Syy
Exy = ((1+nu)./E).*Sxy;


end

