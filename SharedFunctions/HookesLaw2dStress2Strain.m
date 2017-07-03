function [ Exx,Eyy,Exy ] = HookesLaw2dStress2Strain( Sxx,Syy,Sxy,lambda,mu )
%Hooke'sLaw2dStress2Strain using Plane strain 2d Hooke's law
%Converting the strain tensors to stress tensors using the shear
%modulus and lambda
%Works with column vectors single values

%   Copyright 2017, Tim Davis, The University of Aberdeen
%Equation 8.48 Pollard and fletcher. 
nu =lambda./(2.*(mu+lambda));            %Poisson's ratio, Equation 8.28 Pollard
E=(lambda.*(1+nu).*(1-2*nu))./nu; 		  %Young's Mod,      Equation 8.27 Pollard. 
Exx = 1./E.*(Sxx.*(1-nu.^2)-(Syy.*nu.*(1+nu)));
Eyy = 1./E.*(Syy.*(1-nu.^2)-(Sxx.*nu.*(1+nu)));
Exy = ((1+nu)./E).*Sxy;


end

