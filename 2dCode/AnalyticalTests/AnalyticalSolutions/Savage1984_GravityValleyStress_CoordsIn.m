function [Sxx,Syy,Sxy,x,y] = Savage1984_GravityValleyStress_CoordsIn(tectxx,rg,nu,a,b,x,y)
% Savage1984_GravityValleyStress_CoordsIn: Returns Cartesian stresses
%              at grid points under a incised valley surface. Stress are
%              due to unloading due to erosion and take into account
%              gravitational loading. The valley surface is free of the
%              traction that was at this surface before it was incised at
%              the end of the solution. 
%
%              The solution is a MATLAB Version of the Fortran code from:
%              Savage, W.Z., Powers, P.S. and Swolfs, H.S., 1984. In situ
%              geomechanics of crystalline and sedimentary rocks; Part V,
%              RVT, a Fortran program for an exact elastic solution for
%              tectonics and gravity stresses in isolated symmetric ridges
%              and valleys (No. 84-827). US Geological Survey,.
%
% usage #1: 
% [Sxx,Syy,Sxy,x,y] = Savage1984_GravityValleyStress_CoordsIn(tectxx,rg,nu,a,b,x,y)
%
% Arguments: (input)
%     tectxx  - Far field remote stress Sxx either pulling or pushing the
%              ridge. Tension positive convention.
%
%       mu    - Shear Modulus.
%
%       nu    - Poisson's Ratio.
%
%       rg    - Weight, the density of material times the acceleration due to
%              gravity.
%
%       a,b   - Parameters describing the valleys slope, see the paper for 
%              a diagram.
%
%       x,y   - Observation points, actual x and y locations
%
%
% Arguments: (output)
%
%       x,y    -  same as inputs
%
%  Sxx,Syy,Sxy - 2D stress tensor components returned on a grid. This is
%               the total stress (including remote and gravitation
%               loading).
%
%
% Example usage:
%
%  nu = 0.25;
%  mu = 500;
%  tectx=1;
%  rg=1*9.81;
%  a=2;
%  b=-1;
%  u= linspace(0,4,50);
%  v = linspace(-0.3265,-4,46);
%  [u,v] = meshgrid(u,v);
%  
%  [Sxx,Syy,Sxy,X,Y] =...
%  Savage1984_GravityValleyStress_CoordsIn(tectx,rg,nu,a,b,x,y);
% 
%  DrawScatterPlots2d( X,Y,[],Sxx,Syy,Sxy )
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen
%  Modified from Steve Martel's fracture mechanics homework

z=complex(x,y);

%rad is the small radius from point w = -ia.
rad = 1.0e-04;
so = rg*b;
po = complex(0.0, 1.0);
ai = po*a;

% Phi at w = -ia for gravity solution from equation 4.
phia = -so*((4.*a+b)*(1-nu)+b)/(8.*(1-nu)*(2.*a+b));

% phi at is phi at w = -ia for tectonic solution from equation 12.
phiat = -(tectxx*b)/(4.*(2.*a+b));

%d2phia is the 2nd derivative of phi at w = -ia for tectonic sol
d2phia = -(tectxx*b*(4.*a+b)*(b-12.*a));
d2phia = d2phia/(2.*a*(2.*a+b)*(4.*a+b)^3);

w=zeros(size(z));
w=complex(w);
for i=1:numel(z)
	if x(i)<=0
		w(i)=(ai+z(i)-sqrt(z(i)^2 -2 .*ai .*z(i)-4 .*a.*b+ai.^2))./2;
	else
		w(i)=(ai+z(i)+sqrt(z(i)^2 -2 .*ai .*z(i)-4 .*a.*b+ai.^2))./2;
	end
end

u=real(w);
v=imag(w);

r = sqrt(u.^2+(v+a).^2);

X = x;
Y = y;

dz = ((w-ai).^2-a.*b)./((w-ai).^2);
%aw is a(w) given by equation 3 in text.
awl = po*(4.*a+b)./(8.*(w-ai));
aw2 = a.*b.*(w-3.*ai)./(8.*(1-nu).*(w-ai).^3);
aw = -awl-aw2;

%phi is phi(w) for the gravity solution given by equation 6.
phi = -aw.*so./dz+a.*b.*phia./(dz.*(w-ai).^2);

%From equation 11.
awt = -(a.*b.*tectxx)./(2.*(w-ai).^2);

%From equation 10.
phit = -awt./dz+a.*b.*phiat./(dz.*(w-ai).^2);

%From equation 8.
sumt = 4.*real(phit)+tectxx;
sum = sumt+4.*real(phi)+rg.*Y./(1-nu);

%First derivative of phi(w) in equation 5.
dlawl = po.*(4.*a+b)./(8.*(w-ai).^2);
dlaw2 = 2.*a.*b./(8.*(1-nu).*(w-ai).^3);
dlaw3 = 6.*po.*a.^2.*b./(8.*(1-nu).*(w-ai).^4);
dlaw = dlawl+dlaw2-dlaw3;
dlphi = -so.*dlaw./dz-(2.*a.*b.*(phi+phia))./(dz.*((w-ai).^3));

%Second derivative of phi(w) to be used in equation 5
%when w = -ia

d2phil = 2.*phi./((w-ai).^2);
d2phi2 = 4.*dlphi./(w-ai);
d2phi3 = so.*po.*a.^2.*b./(2.*(1-nu).*((w-ai).^5));
d2phi = -(d2phil+d2phi2+d2phi3)./dz;

%From equation 7.
bwl = po.*(4.*a+b)./(8.*(w-ai));
bw2 = (1-2.*nu).*a.*b.*(w-3.*ai)./(8.*(1-nu).*(w-ai).^3);
bw = -so.*(bwl+bw2);

psil = w.*dlphi+bw+phi;
dlphil = -(a.*b.*tectxx.*(4.*a+b).*(w-ai));
dlphi2 = 2.*(2.*a+b).*(((w-ai).^2-a.*b).^2);
dlphit = dlphil./dlphi2;
bwt = -awt;

%Part of equation 9.
psilt = w.*dlphit+bwt+phit;

%Test on closeness of w to -ia. If w is near -ia, the Taylor
% expansion about -ia is used
if (r<rad)
	psi2 = .5.*a.*b.*d2phi;
	psi2t = .5.*a.*b.*d2phia;
else
	psi2 = a.*b.*dlphi./(w+ai)-a.*b.*(phi-phia)./((w+ai).^2);
	psi2t = a.*b.*dlphit./(w+ai)-a.*b.*(phit-phiat)./((w+ai).^2);
end

psi = -(psil+psi2)./dz;
psit = -(psilt+psi2t)./dz;
zbar = conj(z);

%From equation 5.
strl = 2.*(zbar.*dlphi./dz+psi);
%Added
strlt = 2.*(zbar.*dlphit./dz+psit);

%From equation 9.
%The equations in this block form the sums and differences
%of equations 2, 5, 8, and 9, to obtain stresses.

dift = real (strlt)-tectxx;
dif = dift+real (strl)+rg.*Y.*(1-2.*nu)./(1-nu);
Sxy = imag(strl)+imag(strlt);
Sxy = Sxy./2;
Sxx = (sum-dif)./2;
Syy = (sum+dif)./2;
