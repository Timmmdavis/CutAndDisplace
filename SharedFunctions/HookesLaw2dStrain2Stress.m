function [ Sxx,Syy,Sxy ] = HookesLaw2dStrain2Stress( Exx,Eyy,Exy,E,nu,mu )
%Hooke'sLaw2dStrain2Stress using Plane strain 2d Hooke's law
%Converting the strain tensors to stress tensors using the Young's
%modulus, shear modulus and Poisson's ratio
%Works with column vectors or single values

%   Copyright 2017, Tim Davis, The University of Aberdeen

	%Hooke's law Plane Strain Eq 8.34 Pollard
	lambda= E*nu/((1+nu)*(1-2*nu));   %Lamé's  constant,  Equation 8.27 Pollard
	Sxx=(2*mu+lambda)*Exx+lambda*Eyy;
	Syy=(2*mu+lambda)*Eyy+lambda*Exx;	
	Sxy= 2*mu*Exy; %We assume here Exy is not engineering strain.
	

end

