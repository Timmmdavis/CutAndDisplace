function [Stress] = LDstressFS(x,y,xe,ye,a,Beta,Ds,Dn,nu,E)
% LDstressFS: LineDislocationInducedStressFullSpace Computes the influence 
%               (stress) of a single planar line crack shearing and/or
%               opening on the surronding points. (full space).
%
% usage #1:
%[Stress] = LDstressFS(x,y,xe,ye,a,Beta,Ds,Dn,nu,E)
%
% Arguments: (input)
%      x & y  - The observation points locations in real Cartesian coords. 
%
%     xe & ye - The element midpoint location in real coords. 
%
%       a     - The elements half length
%
%       Beta  - The angle of the element away from the x-axis (radians).
%               When the normal in the -y axis down this is 0. In degrees
%               when the normal points east this is 90, west -90 and north
%               180.
%
%     Dn,Ds   - The defined displacement of each element.(normal and shear)
%               Dn+ is opening, Ds+ is left lateral shearing. 
%
%       nu    - The Poisson's ratio
%
%       E     - The Young's modulus
%
% Arguments: (output)
% Stress - Is the stress caused by the movement of the dislocation at the 
%               observataion points. [Sxx,Syy,Sxy].
%
% Example usage:
%
%  [x,y]=meshgrid(-2:0.1:2,-2:0.1:2);
%  dimx = length(x(:,1)); 
%  dimy = length(x(1,:));   
%  Ds=0; Dn=1;
%  [Stress] = LDstressFS(x(:),y(:),0,0,1,0,Ds,Dn,0.25,1);
%  [Sxx,Syy,Sxy ] = ExtractCols( Stress );
%  [Sxx,Syy,Sxy]=ReshapeData2d( dimx,dimy,Sxx,Syy,Sxy);
%  DrawContourFPlots2d(x,y,[],Sxy,Syy,Sxx);
%
%  Author: Tim Davis/Steve Martel
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen
%  Modified form of Steve Martels BEM scripts from:
%  http://www.soest.hawaii.e.du/martel/Martel.BEM_dir/
%  Also includes some notation from Ritz & Pollard 2012.


% The shear modulus, sm, is related to the prescribed elastic constants.
sm = E/(2*(1+nu));
% Define material constant used in calculating influence coefficients.
con = 1/(4*pi*(1-nu)); 
cons = 2*sm;
H=a;
Dxb = Ds; Dyb = -Dn; 
sb = sin(Beta); cb = cos(Beta);
s2b = sin(2*Beta); c2b = cos(2*Beta);
% Define array of local coordinates for the observation grid relative to
%   the midpoint and orientation of the ith element.
% Refer to (Figure 5.6, C&S, p. 91) and eqs. 4.5.1 of C&S, p. 57. 
XB = (x-xe)*cb + (y-ye)*sb;
YB = -(x-xe)*sb + (y-ye)*cb;

% Calculate derivatives of the function f(x,y), eq. 5.2.5 of C&S, p. 81. 
%   which are used to calculate the displacement and stress components. 
% It is understood that X and Y refer to XB and YB.
% First abbreviate repeated terms in the derivatives of f(x,y):
Y2 = YB.^2;
XMH = XB-H; XPH = XB+H; 
XMH2 = XMH.^2; XPH2 = XPH.^2; 
R1S = XMH2 + Y2; R1S2 = R1S.^2; 
R2S = XPH2 + Y2; R2S2 = R2S.^2;
% % The following derivatives are eqs. 4.5.5a thru d of C&S, p. 58.
% FF4 = con*(YB./R1S - YB./R2S); 
% FF5 = con*(XMH./R1S - XPH./R2S);
% 
% % The following derivatives are eqs. 5.5.3a and b of C&S, p. 91.
% FF6 = con*((XMH2 - Y2)./R1S2 - (XPH2 - Y2)./R2S2);
% FF7 = 2*con*YB.*(XMH./R1S2 - XPH./R2S2);
% 
% % Calculate the stress components using eqs. 5.5.5 of C&S, p. 92.
% Sxx =  cons*Dxb*(2*(cb*cb)*FF4 + s2b*FF5 + YB.*(c2b*FF6-s2b*FF7))...
%      +cons*Dyb*(-FF5 + YB.*(s2b*FF6 + c2b*FF7));
% Syy =  cons*Dxb*(2*(sb*sb)*FF4 - s2b*FF5 - YB.*(c2b*FF6-s2b*FF7))...
%      +cons*Dyb*(-FF5 - YB.*(s2b*FF6 + c2b*FF7));
% Sxy =  cons*Dxb*(s2b*FF4 - c2b*FF5 + YB.*(s2b*FF6+c2b*FF7))...
%      +cons*Dyb*(-YB.*(c2b*FF6 - s2b*FF7)); 

% The following derivatives are eqs. 4.5.5a thru d of C&S, p. 58.
FF4 = (YB./R1S - YB./R2S); 
FF5 = (XMH./R1S - XPH./R2S);

% The following derivatives are eqs. 5.5.3a and b of C&S, p. 91.
FF6 = ((XMH2 - Y2)./R1S2 - (XPH2 - Y2)./R2S2);
FF7 = 2*YB.*(XMH./R1S2 - XPH./R2S2);

% Calculate the stress components using eqs. 5.5.5 of C&S, p. 92.
Sxx =  Dxb*(2*(cb*cb)*FF4 + s2b*FF5 + YB.*(c2b*FF6-s2b*FF7))...
     +Dyb*(-FF5 + YB.*(s2b*FF6 + c2b*FF7));
Syy =  Dxb*(2*(sb*sb)*FF4 - s2b*FF5 - YB.*(c2b*FF6-s2b*FF7))...
     +Dyb*(-FF5 - YB.*(s2b*FF6 + c2b*FF7));
Sxy =  Dxb*(s2b*FF4 - c2b*FF5 + YB.*(s2b*FF6+c2b*FF7))...
     +Dyb*(-YB.*(c2b*FF6 - s2b*FF7)); 

Stress=[Sxx(:),Syy(:),Sxy(:)]*cons*con;

