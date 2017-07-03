%   Copyright 2017, Tim Davis, The University of Aberdeen

%eq from Analytical volcano deformation sourcemodels
%Chapter 8 Michael Lisowski
close all;clear



%Depth of chamber
D=50; %depth
%Obs points
X=linspace(0,100,100);
Y=zeros(1,numel(X));

R=sqrt(D.^2+X.^2); %Distance from centre of chamber to point
Radius=1; %Radius of the magma chaber 
P=-1; %Pressure (Negative, a compressive force)
nu=0.25; %Poisson's ratio
G=1; %Shear mod

lsi=Radius^3*-P*((1-nu)/G);%leftside of eq 8.15

%Mctigue eq from the book, doesn't work correctly
%lsi=lsi*(1+((Radius/D)^3))*(((1+nu)/(2*(-7+5*nu)))+((15*(D^2)*(-2+nu))./(4*(R.^2)*(-7+5*nu))));

u=(X./(R.^3)).*lsi;
v=(Y./(R.^3)).*lsi;
w=(D./(R.^3)).*lsi;

plot(X,u)
hold on
plot(X,w)

%not truly analyical but its a nice approximation, check paul segalls book
%for the derivationa and why its an approx
