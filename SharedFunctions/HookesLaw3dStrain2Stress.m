function [ Sxx,Syy,Szz,Sxy,Sxz,Syz ] = HookesLaw3dStrain2Stress( Exx,Eyy,Ezz,Exy,Exz,Eyz,lambda,mu )
%Hooke'sLaw3dStrain2Stress using Hooke's law
%Converting the strain tensors to stress tensors using the shear
%modulus and lambda
%Works with column vectors single values

%   Copyright 2017, Tim Davis, The University of Aberdeen
%Equation 7.131 and 7.132 in Pollard and Fletcher 2005 Book. 
Sxx = 2*mu.*Exx+(lambda.*(Exx+Eyy+Ezz)); 
Syy = 2*mu.*Eyy+(lambda.*(Exx+Eyy+Ezz));
Szz = 2*mu.*Ezz+(lambda.*(Exx+Eyy+Ezz));
clear Exx Eyy Ezz
Sxy = 2*mu.*Exy;     clear Exy %We assume here Exy is not engineering strain. 
Sxz = 2*mu.*Exz;     clear Exz
Syz = 2*mu.*Eyz;     clear Eyz

end

