function [X,Y,Sxx,Syy,Sxy,Ux,Uy]=Barber1992_GlideDislocation(k,mu,X,Y,a,b,nu)
% Barber1992_GlideDislocation: Returns Cartesian displacements and stresses
%              at grid points for a displacement discontinuity with a unit
%              Burger's vector (B=1) and unit half length (a=1). The
%              displacement discontintuity extends along the x-axis from 
%              x = -a to x = +a. 
%              Note a positive Burgers vector here drives a right lateral
%              shearing. 
%
%              The solution is based on equations for a glide
%              dislocation from Barber: Elasticity (1992) and Pollard and
%              Fletcher (2005). 
%
%
%
% usage #1: [X,Y,Sxx,Syy,Sxy,Ux,Uy] =
% Barber1992_GlideDislocation(k,mu,X,Y,a,b,nu);
%
% Arguments: (input)
%       k     - Kolosov's constant for plane strain.
%
%       mu    - Shear Modulus.
%
%  minx,maxx  - The bounds of the grid that stresses and displacements are
%              found on. (Y and X).
%
%  spacing    - The spacing between the points on this grid.
%
%       a     - The half length of the dislocation.
%
%       b     - The Burgers vector (separation of dislocation walls).
%
%       nu    - Poisson's Ratio.
%
%
% Arguments: (output)
%
%       X,Y    - X and Y locations of the observation points.
%
%  Sxx,Syy,Sxy - 2D stress tensor components returned on a grid.
%
%       Ux,Uy  - Displacement of stress tensor components returned on a
%               grid.
%
%
% Example usage:
%
%  nu = 0.25;
%  k = 3-4*nu; 
%  mu = 500;
%  spacing=0.1;
%  minx=-4; maxx=4;
%  [X,Y] = meshgrid(minx:spacing:maxx);
%  a = 1;  
%  b=0.0001; 
%  
%  [X,Y,Sxx,Syy,Sxy,Ux,Uy] =...
%   Barber1992_GlideDislocation(k,mu,X,Y,a,b,nu);
% 
%  quiver(X(:),Y(:),Ux(:),Uy(:))
%  DrawContourFPlots2d( X,Y,[], Sxx,Syy,Sxy )
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen
%  Modified from Steve Martel's fracture mechanics homework



%Define larger grid that covers where both dislocations will be, this is
%used to calculate stress once. 
[Sxxp,Syyp,Sxyp,Uxp,Uyp]=CompD(mu,b,k,nu,X-a,Y);
[Sxxn,Syyn,Sxyn,Uxn,Uyn]=CompD(mu,b,k,nu,X+a,Y);

Ux=Uxp-Uxn;
Uy=Uyp-Uyn;

Sxx=Sxxp-Sxxn;
Syy=Syyp-Syyn;
Sxy=Sxyp-Sxyn;

function [Sxx,Syy,Sxy,Ux,Uy]=CompD(mu,b,k,nu,X,Y)

%Set up vars
r = sqrt(X.^2 + Y.^2);
sint = Y./r; 
cost = X./r;  
theta = atan2(Y,X);
    
    
% Calculate polar rt stress components due to unit displacement discontinuity
% Positive end Negative end
% Barber page 201 equations 13.24-13.25
Srr = ((2*mu)*b*sint)./(pi*(k+1)*r);
Stt = Srr; 
Srt = -((2*mu)*b*cost)./(pi*(k+1)*r);

%Converting the components into Cartesian tensor
dimx = size(Srr,1);
dimy = size(Srr,2);
%Calling external function
[Sxx,Syy,Sxy]=StressTensorTransformation2d(Srr(:),Stt(:),Srt(:),cos(theta(:)),cos((pi/2)-theta(:)));
[Sxx,Syy,Sxy]=ReshapeData2d(dimx,dimy,Sxx,Syy,Sxy);

%Displacement equations, see Pollard and Fletcher 2005 Eq 8.36-8.37. 
%Equations from Fig 08_10
lambda = (2*mu*nu)/(1-2*nu); 
%Constants
c1 = (0.5*mu)/(lambda+2*mu); 
c2 = (lambda+mu)/(lambda+2*mu);
%Equations
%Ux
Ux1 = -(b/(2*pi))*atan2(Y,X);
Ux2 = -(b/(2*pi))*c2*(X.*Y)./(X.^2+Y.^2);
Ux=Ux1+Ux2; bbb=abs(Ux)==inf; Ux(bbb)=0 ;
%Uy
Uy1 = -(b/(2*pi))*(-c1).*log(X.^2+Y.^2);
Uy2 = -(b/(2*pi))*c2*(Y.^2)./(X.^2+Y.^2);
Uy=Uy1+Uy2; bbb=abs(Uy)==inf; Uy(bbb)=0 ;


end

end

