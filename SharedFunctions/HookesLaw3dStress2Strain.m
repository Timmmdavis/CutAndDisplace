function [ Exx,Eyy,Ezz,Exy,Exz,Eyz ] = HookesLaw3dStress2Strain( Sxx,Syy,Szz,Sxy,Sxz,Syz,lambda,mu )
%Hooke'sLaw3dStrain2Stress using Hooke's law
%Converting the strain tensors to stress tensors using the shear
%modulus and lambda
%Works with column vectors single values

%   Copyright 2017, Tim Davis, The University of Aberdeen
%Equation 7.150 Pollard and fletcher. (see errata of book -
%https://pangea.stanford.edu/projects/structural_geology/errata/errata.pdf)
nu =lambda/(2*(mu+lambda));            %Poisson's ratio, Equation 8.28 Pollard
E=(lambda*(1+nu)*(1-2*nu))/nu; 		 %Young's Mod,	Equation 8.27 Pollard. 
Exx = 1/E.*(Sxx-(nu.*(Syy+Szz)));
Eyy = 1/E.*(Syy-(nu.*(Sxx+Szz)));
Ezz = 1/E.*(Szz-(nu.*(Sxx+Syy)));
Exy = ((1+nu)/E).*Sxy;
Exz = ((1+nu)/E).*Sxz;
Eyz = ((1+nu)/E).*Syz;

end

