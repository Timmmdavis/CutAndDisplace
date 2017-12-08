function [rOvb,SrrNorm,SttNorm]=CrouchStar1983_InhomogeneousAnnulus(mu1,nu1,mu2,a,b,r)
% CrouchStar1983_InhomogeneousAnnulus: Returns the normalised radial and
%               tangential stresses at points at a radial distance
%               within and away from an annulus with 2 material constants
%               and an internal pressure.
%
%               See: Crouch, S.L. and Starfield, A.M., 1983. Boundary
%               element methods in solid mechanics: With applications in
%               rock mechanics and geological engineering: Winchester,
%               Massachusetts.
%
% usage #1:
% [rOvb,SrrNorm,SttNorm]=CrouchStar1983_InhomogeneousAnnulus(mu1,nu1,mu2,a,b,r)
%
% Arguments: (input)
%       nu1    - Poisson's Ratio for inner elastic.
%
%       mu1    - Shear Modulus for inner elastic.
%
%       mu2    - Shear Modulus for outer elastic.
%
%       a      - The central holes radius
%
%       b      - The distance from the centre to the Annulus's outer edge
%               (elastic interface)
%
%       r      - The radial distance of an observation point to the
%               problems centre
%
%
% Arguments: (output)
%       rOvb  - The radial distance of the observation point from the
%              problems centre, normalised with respect to distance b.
%
%    SrrNorm  - The radial stress at the observation points, values 
%              normalised to the driving pressure.
%
%    SttNorm  - The tangential stress at the observation points, values 
%              normalised to the driving pressure
%
%
% Example usage:
%
%  mu1 = 500; 
%  nu1 = 0.25;
%  mu2 = 1000; 
%  a=1;  
%  b=2;       
%  r=linspace(0,b+2,50);
%  r(r<a)=[]; %Removing points in the hole
%  
%  [rOvb,SrrNorm,SttNorm]=CrouchStar1983_InhomogeneousAnnulus(mu1,nu1,mu2,a,b,r);
% 
%  plot(r,SrrNorm);hold on;plot(r,SttNorm);
%  xlabel('radius'); ylabel('Stress'); legend('Srr','Stt');
%  
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Parameters
P=1;    %Internal hole pressure    

%R between a and b (1st elastic)
Rad_a_r_bFlg=r<b & r>a; 
%R above b (2nd elastic)
Rad_b_rFlg=r>b;
%Grabbing these
Rad_a_r_b=r(Rad_a_r_bFlg);
Rad_b_r=r(Rad_b_rFlg);


%Equations for PPrime (assigned pr in script)
one=2*(1-nu1)*P*(a^2)/(b^2);
two=2*(1-nu1)+(mu1/mu2-1)*(1-(a^2)/(b^2));
pr= one/two; %C&S 7.5.24

%Equations for calculating Srr (Eq.1 Crouch and Star 7.5.23)
one=(1/(1-(a^2)/(b^2)));
two=(((P*(a^2))/(b^2)-pr)-(P-pr)*(a^2)./(Rad_a_r_b.^2));
Srr_a_r_b=one*two;
Srr_b_r=-pr*(b^2)./(Rad_b_r.^2);


%Equations for calculating Stt (Eq.1 Crouch and Star 7.5.23)
one=(1/(1-(a^2)/(b^2)));
two=(((P*(a^2))/(b^2)-pr)+(P-pr)*(a^2)./(Rad_a_r_b.^2));
Stt_a_r_b=one*two;
Stt_b_r=+pr*(b^2)./(Rad_b_r.^2);


%Normalizing
SrrNorm_a_r_b=Srr_a_r_b/P;
SrrNorm_b_r=Srr_b_r/P;

SttNorm_a_r_b=Stt_a_r_b/P;
SttNorm_b_r=Stt_b_r/P;

rOvb_a_r_b=Rad_a_r_b/b;
rOvb_b_r=Rad_b_r/b;


%Appending
SrrNorm=[SrrNorm_a_r_b,SrrNorm_b_r];
SttNorm=[SttNorm_a_r_b,SttNorm_b_r];
rOvb=[rOvb_a_r_b,rOvb_b_r];


%Plotting

% plot(rOvb,SrrNorm);
% hold on
% plot(rOvb,SttNorm);
% title('Srr (blue) Stt (orange)'), xlabel('rad/b')
% ylabel('SrrSttNorm');