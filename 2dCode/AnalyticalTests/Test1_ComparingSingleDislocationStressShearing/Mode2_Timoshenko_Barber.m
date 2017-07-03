function [sxx,syy,sxy,ux,uy]=Mode2_Timoshenko_Barber(k,U,minx,maxx,spacing,a,B,pr)

%   Copyright 2017, Tim Davis, The University of Aberdeen
%Modified from Steve martels fracture homeworks

% Returns Cartesian displacements and stresses at grid points
% for a displacement discontinuity
% with a unit Burger's vector (B=1) and unit half length (a=1).
% The displacement discontintuity extends along the x-axis
% from x = -a to x = +a.
% The solution is based on equations for a glide dislocation from Barber 2010
%
%k Kolosovs constant for plane strain
%U ShearMod
%a Unit half-length displacement discontinuity
%B Length of the Burger's vector.



% ux = displacement in the z-direction
% uy = displacement in the y-direction
% sxx = sigma xx
% sxy = sigma xy
% syy = sigma yy
% Parameters k, and PR are elastic parameters (see Barber, 2010)
% x,y = coordinates of observation gridpoints,
%
% IMPORTANT: In Barber's solutions, a dislocation (d) cut extends
% from the origin to the right along the x-axis (not the left). 
% % 
% % %%%%
% % % Example:
% % PR = 0.25;
% % k = (3-PR)/(1+PR); %Kolosov's constant for plane stress
% % U = 500; %ShearMod
% % spacing=0.1;
% % minx=-4; maxx=4;
% % [x,y] = meshgrid(minx:spacing:maxx); %large grid
% % a = 1;  % Unit half-length displacement discontinuity
% % B=0.0001; % Length of the Burger's vector.

%Define larger grid that covers where both dislocations will be, this is
%used to calculate stress once. 
[x2,y2] = meshgrid(minx-a:spacing:maxx+a);
r = sqrt(x2.^2 + y2.^2);
sint = y2./r; 
cost = x2./r;  
t = atan2(y2,x2);
logr = log(r);
axr = cost; ayr = sint; axt = -ayr; ayt = axr;

% Calculate polar rt STRESS components due to unit displacement discontinuity
% Positive end Negative end
% Barber page 201 equations 13.19-13.20, these are for a climb dislocation, not glide. 
srr = -((2*U)*B*sint)./(pi*(k+1)*r);
stt = srr; 
srt = ((2*U)*B*cost)./(pi*(k+1)*r);
str = srt; 


% % % Convert stresses from polar coordinates to cart coordinates
% [x2,y2] = pol2cart(t,r) %Use if we didnt have x and y already
dim=size(srr);
srr=srr(:);stt=stt(:);srt=srt(:); %all col vecs
t=t(:); 
sxx=zeros(size(srr)); sxy=zeros(size(srr)); 
syy=zeros(size(srr)); syx=zeros(size(srr)); 
%http://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Polar_Coords/Polar_Coords.htm
%Section 2.7
%Loop to convert
for i=1:numel(srr)
    Cyli=[srr(i),srt(i);
          srt(i),stt(i)];
    RV1= [cos(t(i)),-sin(t(i));
          sin(t(i)), cos(t(i))];
    RV2= [cos(t(i)), sin(t(i));
         -sin(t(i)), cos(t(i))]; 
    Cart=RV1*Cyli*RV2;
    sxx(i)=Cart(1,1);
    syy(i)=Cart(2,2);
    sxy(i)=Cart(1,2);
    syx(i)=Cart(2,1);
end   
sxx=reshape(sxx,dim); sxy=reshape(sxy,dim);
syy=reshape(syy,dim); syx=reshape(syx,dim);

%Now separating the stresses from a Barber dislocation into two arrays that define the positive
%and negative ends of the dislocation in real coordinates 
sxxp = ClipToPosNeg(sxx,a,spacing,1);
sxxn = ClipToPosNeg(sxx,a,spacing,0);
sxyp = ClipToPosNeg(sxy,a,spacing,1);
sxyn = ClipToPosNeg(sxy,a,spacing,0);
syxp = ClipToPosNeg(syx,a,spacing,1);
syxn = ClipToPosNeg(syx,a,spacing,0);
syyp = ClipToPosNeg(syy,a,spacing,1);
syyn = ClipToPosNeg(syy,a,spacing,0);
% Normal discontintuity, Both ends (note minus sign)
% Superpose dislocations extending right from x = -a and x = +a.
sxx = sxxp-sxxn; sxy = sxyp-sxyn; syx= syxp-syxn; syy = syyp-syyn; 

%Displacement equations, see Pollard and Fletcher 2005 Eq 8.36-8.37. 
%Equations from Fig 08_10
b=-B;
mu = U; lm = (2*mu*pr)/(1-2*pr); %pr=PR % Elastic moduli
c1 = (0.5*mu)/(lm+2*mu); c2 = (lm+mu)/(lm+2*mu);
Ux1 = -(b/(2*pi))*atan2(y2,x2);
Ux2 = -(b/(2*pi))*c2*(x2.*y2)./(x2.^2+y2.^2);
ux=Ux1+Ux2; bbb=abs(ux)==inf; ux(bbb)=0 ;
Uy1 = -(b/(2*pi))*(-c1).*log(x2.^2+y2.^2);
Uy2 = -(b/(2*pi))*c2*(y2.^2)./(x2.^2+y2.^2);
uy=Uy1+Uy2; bbb=abs(uy)==inf; uy(bbb)=0 ;

uxp = ClipToPosNeg(ux,a,spacing,1);
uxn = ClipToPosNeg(ux,a,spacing,0);
uyp = ClipToPosNeg(uy,a,spacing,1);
uyn = ClipToPosNeg(uy,a,spacing,0);

ux=uxp-uxn;
uy=uyp-uyn;

function Var = ClipToPosNeg(Var,a,spacing,PosNeg)

clip=a/spacing;
%length=size(Var);


if PosNeg==1 %Positive
    Var=Var(1+clip:(size(Var,1))-clip,1:(size(Var,2))-(2*clip));
else
    Var=Var(1+clip:(size(Var,1))-clip,(2*clip+1):(size(Var,2)));
end



