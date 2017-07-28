
function [sigx,sigy,sigxy,x,y] = SavageSlope1984FortranCodeMATLABFunc(ts,rg,pr,a,b,u,v)

%MATLAB Version of the Fortran code:

%   Copyright 2017, Tim Davis, The University of Aberdeen
%Savage, William Z., Philip S. Powers, and Henri S. Swolfs. In situ
%geomechanics of crystalline and sedimentary rocks; Part V, RVT, a Fortran
%program for an exact elastic solution for tectonics and gravity stresses
%in isolated symmetric ridges and valleys. No. 84-827. US Geological
%Survey,, 1984.

%Parameters to change below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tectonic Stress
% ts=0; 
% 
% %Weight and elastic constants
% rg= 1; %Weight (Dens*AccelGrav)
% pr=0.25;
% 
% %A and B describing slope, See paper for diagram
% a=2;
% b=-1;
% 
% %Setting up grid points U and V, these are mapped to a different location
% %later
% u = linspace(-0,4,50);
% v = linspace(0,-4,50); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%rad is the small radius from point w = -ia.
rad = 1.0e-04;
so = rg*b;
po = complex(0.0, 1.0);
ai = po*a;

% Phi at w = -ia for gravity solution from equation 4.
phia = -so*((4.*a+b)*(1-pr)+b)/(8.*(1-pr)*(2.*a+b));

% phi at is phi at w = -ia for tectonic solution from equation 12.
phiat = -(ts*b)/(4.*(2.*a+b));

%d2phia is the 2nd derivative of phi at w = -ia for tectonic sol
d2phia = -(ts*b*(4.*a+b)*(b-12.*a));
d2phia = d2phia/(2.*a*(2.*a+b)*(4.*a+b)^3);

w = complex(u, v);
r = sqrt(u.^2+(v+a).^2);

%From equation 1.
z = w+(a.*b)./(w-ai);

x = real(z);
y = imag(z);
dz = ((w-ai).^2-a.*b)./((w-ai).^2);
%aw is a(w) given by equation 3 in text.
awl = po*(4.*a+b)./(8.*(w-ai));
aw2 = a.*b.*(w-3.*ai)./(8.*(1-pr).*(w-ai).^3);
aw = -awl-aw2;

%phi is phi(w) for the gravity solution given by equation 6.
phi = -aw.*so./dz+a.*b.*phia./(dz.*(w-ai).^2);

%From equation 11.
awt = -(a.*b.*ts)./(2.*(w-ai).^2);

%From equation 10.
phit = -awt./dz+a.*b.*phiat./(dz.*(w-ai).^2);

%From equation 8.
sumt = 4.*real(phit)+ts;
sum = sumt+4.*real(phi)+rg.*y./(1-pr);

%First derivative of phi(w) in equation 5.
dlawl = po.*(4.*a+b)./(8.*(w-ai).^2);
dlaw2 = 2.*a.*b./(8.*(1-pr).*(w-ai).^3);
dlaw3 = 6.*po.*a.^2.*b./(8.*(1-pr).*(w-ai).^4);
dlaw = dlawl+dlaw2-dlaw3;
dlphi = -so.*dlaw./dz-(2.*a.*b.*(phi+phia))./(dz.*((w-ai).^3));

%Second derivative of phi(w) to be used in equation 5
%when w = -ia

d2phil = 2.*phi./((w-ai).^2);
d2phi2 = 4.*dlphi./(w-ai);
d2phi3 = so.*po.*a.^2.*b./(2.*(1-pr).*((w-ai).^5));
d2phi = -(d2phil+d2phi2+d2phi3)./dz;

%From equation 7.
bwl = po.*(4.*a+b)./(8.*(w-ai));
bw2 = (1-2.*pr).*a.*b.*(w-3.*ai)./(8.*(1-pr).*(w-ai).^3);
bw = -so.*(bwl+bw2);

psil = w.*dlphi+bw+phi;
dlphil = -(a.*b.*ts.*(4.*a+b).*(w-ai));
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

dift = real (strlt)-ts;
dif = dift+real (strl)+rg.*y.*(1-2.*pr)./(1-pr);
sigxy = imag(strl)+imag(strlt);
sigxy = sigxy./2;
sigx = (sum-dif)./2;
sigy = (sum+dif)./2;

% %Plotting figures
% %Sxx
% figure,[C,h] =contourf(x,y,sigx);
% hold on
% contour(x,y,sigx,'k');
% xlabel('x'); ylabel('y'); axis('equal'); title('Sxx');h.LevelStep=0.2;set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1);
% hold off
% %Syy
% figure,[C,h] =contourf(x,y,sigy);
% hold on
% contour(x,y,sigy,'k');
% xlabel('x'); ylabel('y'); axis('equal'); title('Syy');h.LevelStep=0.5;set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1);
% hold off
% %Sxy
% figure,[C,h] =contourf(x,y,sigxy);
% hold on
% contour(x,y,sigxy,'k');
% xlabel('x'); ylabel('y'); axis('equal'); title('Sxy');h.LevelStep=0.02;set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1);
% hold off
