function [ Sxx,Syy,Szz,Sxy,Sxz,Syz ] = HookesLaw3dStrain2Stress( Exx,Eyy,Ezz,Exy,Exz,Eyz,lambda,mu )
% HookesLaw3dStrain2Stress: Using 3D Hooke's law
%                   Converting the strain tensors to stress tensors using
%                   the Young's modulus, shear modulus and Poisson's ratio.
%                   Works with column vectors or single values.
%
%                   Eq. 7.131 and 7.132 of Pollard and Fletcher, 2005, 
%
%               
% usage #1:
% [ Sxx,Syy,Szz,Sxy,Sxz,Syz ] = HookesLaw3dStrain2Stress( Exx,Eyy,Ezz,Exy,Exz,Eyz,lambda,mu )
%
% Arguments: (input)
% Exx,Eyy,Ezz       
% Exy,Exz,Eyz       - Strain tensor components. (Can be column vectors). 
%
% lambda            - Lame's constant
%
% mu                - The shear modulus
%
% Arguments: (output)
% Sxx,Syy,Szz       
% Sxy,Sxz,Syz        - Stress tensor components. 
%
% Example usage:
%
% mu=5; %[Gpa]
% nu=0.2;
% E = mu*(2*(1+nu)) ;
% lambda= E*nu/((1+nu)*(1-2*nu));
% Exx=0.2; Eyy=0; Ezz=0;
% Exy=0;   Exz=0; Eyz=0
% [ Sxx,Syy,Szz,Sxy,Sxz,Syz ] = HookesLaw3dStrain2Stress( Exx,Eyy,Ezz,Exy,Exz,Eyz,lambda,mu )
%
% Additional elastic constant conversions if needed:
% nu =lambda/(2*(mu+lamda);         Equation 8.28 Pollard       
% mu = E/(2*(1+nu));                Equation 8.26 Pollard
% E = mu*(2*(1+nu)) ;               Equation 8.26 Pollard (rearranged)
% lambda= E*nu/((1+nu)*(1-2*nu));   Equation 8.27 Pollard
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Equation 7.131 and 7.132 in Pollard and Fletcher 2005 Book. 
Sxx = 2*mu.*Exx+(lambda.*(Exx+Eyy+Ezz)); 
Syy = 2*mu.*Eyy+(lambda.*(Exx+Eyy+Ezz));
Szz = 2*mu.*Ezz+(lambda.*(Exx+Eyy+Ezz));
clear Exx Eyy Ezz
Sxy = 2*mu.*Exy;     clear Exy %We assume here Exy is not engineering strain. 
Sxz = 2*mu.*Exz;     clear Exz
Syz = 2*mu.*Eyz;     

end

