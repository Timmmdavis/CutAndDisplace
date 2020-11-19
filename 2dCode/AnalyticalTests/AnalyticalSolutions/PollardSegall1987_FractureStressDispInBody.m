function [Sxx,Syy,Sxy,Ux,Uy]=PollardSegall1987_FractureStressDispInBody(mu,nu,Pyy,Pxy,X,Y)
% PollardSegall1987_FractureStressDispInBody: Returns Cartesian
%              displacements and stresses at grid points surronding a line
%              crack loaded by a far field stress, the crack has a unit
%              half length (a=1) and extends along the x-axis from x = -a
%              to x = +a.
%
%               See: Pollard, D.D. and Segall, P., 1987. Theoretical
%               displacements and stresses near fractures in rock: with
%               applications to faults, joints, veins, dikes, and solution
%               surfaces. Fracture mechanics of rock, 277(349), pp.277-349.
%
% usage #1:
% [SxxChange,SyyChange,SxyChange,Ux,Uy]...
%  =PollardSegall1987_FractureStressDispInBody...
%   (mu,nu,Pyy,Pxy,X,Y)
%
% Arguments: (input)
%       nu    - Poisson's Ratio.
%
%       mu    - Shear Modulus.
%
%       Pyy   - The normal stress at infinity. (traction Tn)
%
%       Pxy   - The shear stress at infinity.  (traction -Ts) 
%
%  minx,maxx  - The bounds of the grid that stresses and displacements are
%              found on. (Y and X).
%
%  spacing    - The spacing between the points on this grid.
%
%
% Arguments: (output)
%
%  Sxx,Syy,Sxy - 2D stress tensor components returned on a grid. These are
%               the stress change due to the dislocation displacement and
%               the driving remote stress components would need to be added
%               back in to get the final solution.
%
%       Ux,Uy  - Displacement of stress tensor components returned on a
%               grid.
%
%
% Example usage:
%
%  mu = 500;
%  nu = 0.25;
%  spacing=0.1;
%  minx=-4; maxx=4;
%  [X,Y] = meshgrid(minx:spacing:maxx);
%  Syy= 1;
%  Sxy= 0;
%  
%  [SxxChange,SyyChange,SxyChange,Ux,Uy]...
%  =PollardSegall1987_FractureStressDispInBody...
%   (mu,nu,Syy,Sxy,X,Y);
% 
%  DrawScatterPlots2d( X,Y,[],SxxChange,SyyChange,SxyChange,Ux,Uy )
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

x=X(:);
y=Y(:);

% Unit half-length displacement discontinuity. This is going to have its 
% centre at 0,0. 
a = 1;  
%Defining stress. Extension positive "Engineering convention"
Pyy=Pyy;
Pxy=Pxy;
Syz=0; %Out of Plane stresses all done at the bottom of code
Sxz=0; %Out of Plane stresses all done at the bottom of code

% Now calculating the locations and angles of each point relative to the
% two ends (positive and negative) and the centre of the fracture. 
% See figure 8.3

% Centre of the crack
r = sqrt(x.^2 + y.^2);      %   Array giving each points radial distance from the centre. 
theta = atan2(y,-x);        %   Angle measured from the x axis. 
theta=abs(theta-pi);        %   making this start from x axis and finish at 2piR.


% Positive end of the crack
xp = x-a;
rp = sqrt(xp.^2 + y.^2);
tp = atan2(y,-xp);          %   theta here is not measured away from positve X
tpp=abs(tp-pi);             %   making this start from x axis and finish at 2piR as described in paper. 

xn = x+a;
rn = sqrt(xn.^2 + y.^2);
tn = atan2(y,-xn);          %   theta here is not measured away from positve X
tnn=abs(tn-pi);             %   making this start from x axis and finish at 2piR. 


% Creating the well used variables that are often used
R=sqrt(rp.*rn);     %as described in pollard 8.31a
THETA=(tpp+tnn)/2;  %as described in pollard 8.31a
R1=R.^-1;R3=R.^-3;  %Radius


% Calculate displacements eq 8.33a and 8.33b Pollard and Segall 1984
% These are seperated into seperate components based on related driving stress so the equation as a whole is easier to read 
vsyy=Pyy.*(2*(1-nu)*(R.*sin(THETA)-r.*sin(theta))-r.*sin(theta).*(r.*R1.*cos(theta-THETA)-1));
usyy=Pyy.*((1-(2*nu))*(R.*cos(THETA)-r.*cos(theta))-r.*sin(theta).*(r.*R1.*sin(theta-THETA)));
usxy=Pxy.*(2*(1-nu)*(R.*sin(THETA)-r.*sin(theta))+r.*sin(theta).*(r.*R1.*cos(theta-THETA)-1));
vsxy=(Pxy.*((1-(2*nu))*(R.*cos(THETA)-r.*cos(theta))+r.*sin(theta).*(r.*R1.*sin(theta-THETA)))).*-1; %Why do I need to do this?

Uy=(vsyy+vsxy)/2*mu;%The 2/G corresponds to the left hand side of the eq 8.33 that i have pushed over
Ux=(usyy+usxy)/2*mu;

% Calculate the stresses from eq 8.44 onwards from Pollard and Segall 1984
% Note this contains no elastic constants related to rigidity. The
% resultant stress related to a fracture and is only dependant on input stress magntitude
% on the fracture surface (this stress must already alude to the strength of the elastic solid) 
SyyChange_syy= Pyy.*(r.*R1.*cos(theta-THETA)-1+a^2.*r.*R3.*sin(theta).*sin(3*THETA));   %8.44a
SyyChange_sxy= Pxy.*(a^2.*r.*R3.*sin(theta).*cos(3*THETA));
Syy= SyyChange_syy + SyyChange_sxy;%Can add remote here if want to find total stress, see eq

SxyChange_sxy= Pxy.*(r.*R1.*cos(theta-THETA)-1-a^2.*r.*R3.*sin(theta).*sin(3*THETA));   %8.44b
SxyChange_syy= Pyy.*((a^2).*r.*R3.*sin(theta).*cos(3*THETA));                             
Sxy=SxyChange_syy + SxyChange_sxy; %Can add remote here if want to find total stress, see eq

SxxChange_syy= Pyy.*(r.*R1.*cos(theta-THETA)-1-a^2.*r.*R3.*sin(theta).*sin(3*THETA));   %8.44c 
SxxChange_sxy= Pxy.*(2.*r.*R1.*sin(theta-THETA)-a^2.*r.*R3.*sin(theta).*cos(3*THETA));
Sxx= SxxChange_syy + SxxChange_sxy; %Can add remote here if want to find total stress, see eq


if Syz>0 %Doing if user has chosen to run with out of plane stress
    % Calculate the stresses from eq 8.44 onwards from Pollard and Segall 1984
    % Note this contains no elastic constants related to rigidity. The
    % resultant stress related to a fracture and is only dependant on input stress magntitude
    % on the fracture surface (this stress must already alude to the strength of the elastic solid) 
    SyzChange_syz= Syz.*(r.*R1.*cos(theta-THETA)-1);                                        %8.44d
    Syz= SyzChange_syz;   %Can add remote here if want to find total stress, see eq


end

if Syz>0 %Doing if user has chosen to run with out of plane stress
    % Calculate the stresses from eq 8.44 onwards from Pollard and Segall 1984
    % Note this contains no elastic constants related to rigidity. The
    % resultant stress related to a fracture and is only dependant on input stress magntitude
    % on the fracture surface (this stress must already alude to the strength of the elastic solid) 
    SxzChange_syz= Syz.*(r.*R1.*sin(theta-THETA));                                        %8.44d
    Sxz= SxzChange_syz;   %Can add remote here if want to find total stress, see eq

end