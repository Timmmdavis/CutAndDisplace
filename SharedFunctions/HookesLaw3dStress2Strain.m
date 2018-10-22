function [ Exx,Eyy,Ezz,Exy,Exz,Eyz ] = HookesLaw3dStress2Strain( Sxx,Syy,Szz,Sxy,Sxz,Syz,lambda,mu )
% HookesLaw3dStress2Strain: Using 3D Hooke's law
%                   Converting the stress tensors to strain tensors using
%                   Lame's constant and the shear modulus.
%                   Works with column vectors or single values.
%
%                   Eq.7.150 from Pollard and Fletcher, 2005 (See Errata!)
%                   %https://pangea.stanford.edu/projects/structural_geology/errata/errata.pdf
%
%               
% usage #1:
% [ Exx,Eyy,Ezz,Exy,Exz,Eyz ] = HookesLaw3dStress2Strain( Sxx,Syy,Szz,Sxy,Sxz,Syz,lambda,mu )
%
% Arguments: (input)
% Sxx,Syy,Szz       
% Sxy,Sxz,Syz       - Stress tensor components. (Can be column vectors). 
%
% lambda            - Lame's constant
%
% mu                - The shear modulus
%
% Arguments: (output)
% Exx,Eyy,Ezz       
% Exy,Exz,Eyz       - Strain tensor components. 
%
% Example usage:
%
% mu=5; %[Gpa]
% nu=0.2;
% E = mu*(2*(1+nu)) ;
% lambda= E*nu/((1+nu)*(1-2*nu));
% Sxx=0.2; Syy=0; Szz=0;
% Sxy=0;   Sxz=0; Syz=0
% [ Exx,Eyy,Ezz,Exy,Exz,Eyz ] = HookesLaw3dStress2Strain( Sxx,Syy,Szz,Sxy,Sxz,Syz,lambda,mu )
%
% % Additional elastic constant conversions if needed:
% % nu =lambda/(2*(mu+lamda);         Equation 8.28 Pollard & Fletcher       
% % mu = E/(2*(1+nu));                Equation 8.26 Pollard & Fletcher
% % E = mu*(2*(1+nu)) ;               Equation 8.26 Pollard & Fletcher (rearranged)
% % lambda= E*nu/((1+nu)*(1-2*nu));   Equation 8.27 Pollard & Fletcher
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

[~,E,~,nu,~]=ElasticConstantsCheck( lambda,mu );

Exx = 1/E.*(Sxx-(nu.*(Syy+Szz)));
Eyy = 1/E.*(Syy-(nu.*(Sxx+Szz)));
Ezz = 1/E.*(Szz-(nu.*(Sxx+Syy)));
clear Sxx Syy Szz
Exy = ((1+nu)/E).*Sxy; clear Sxy
Exz = ((1+nu)/E).*Sxz; clear Sxz
Eyz = ((1+nu)/E).*Syz;

end

