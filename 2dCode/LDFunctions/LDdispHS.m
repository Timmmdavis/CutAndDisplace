function [Disp] = LDdispHS(x,y,xe,ye,a,Beta,Ds,Dn,nu)
% LDdispHS: LineDislocationInducedDisplacementHalfSpace Computes the 
%               influence (displacement) of a single planar line crack
%               shearing and/or opening on the surronding points. Half
%               space at 0.
%
% usage #1:
%[Disp] = LDdispHS(x,y,xe,ye,a,Beta,Ds,Dn,nu)
%
% Arguments: (input)
%      x & y  - The observation points locations in real Cartesian coords.
%              (y must be below 0)
%
%     xe & ye - The element midpoint location in real coords. 
%              (must have a ye below 0). 
%
%       a     - The elements half length.
%
%       Beta  - The angle of the element away from the x-axis (radians).
%               When the normal in the -y axis down this is 0. In degrees
%               when the normal points east this is 90, west -90 and north
%               180.
%
%     Dn,Ds   - The defined displacement of each element.(normal and shear).
%               Dn+ is opening, Ds+ is left lateral shearing. 
%
%       nu    - The Poisson's ratio.
%
%       E     - The Young's modulus.
%
% Arguments: (output)
% Disp        - Is the displacement caused by the movement of the
%             dislocation at the observataion points. [Ux,Uy].
%
% Example usage:
% 
%  [x,y]=meshgrid(-2:0.05:2,-4:0.05:0);
%  dimx = length(x(:,1)); 
%  dimy = length(x(1,:)); 
%  Ds=0; Dn=1;
%  [Disp] = LDdispHS(x(:),y(:),0,-2,1,0,Ds,Dn,0.25);
%  [ Ux,Uy ] = ExtractCols( Disp );
%  [Ux,Uy]=ReshapeData2d( dimx,dimy,Ux,Uy );
%  figure,quiver(x(:),y(:),Ux(:),Uy(:));
%
%  Author: Tim Davis/Steve Martel
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen
%  Modified form of Steve Martels BEM scripts from:
%  http://www.soest.hawaii.e.du/martel/Martel.BEM_dir/
%  Also includes some notation from Ritz & Pollard 2012.


if any(y>0)
    error('Half-space solution: Z coordinates must be negative!')
end

% Define material constant used in calculating influence coefficients.
con = 1/(4*pi*(1-nu)); 
Dxb = Ds; Dyb = -Dn;
sb = sin(Beta); cb = cos(Beta);
s2b = sin(2*Beta); c2b = cos(2*Beta);
s3b = sin(3*Beta); c3b = cos(3*Beta);


% Define array of local coordinates for the observation grid relative to
%   the midpoint and orientation of the ith element.
% Refer to (Figure 5.6, C&S, p. 91) and eqs. 4.5.1 of C&S, p. 57. 
XB = (x-xe)*cb + (y-ye)*sb;
YB = -(x-xe)*sb + (y-ye)*cb;

length=size(x);

% Coordinates of the image dislocation
XBi = (x-xe)*cb - (y+ye)*sb;		%equation 7.4.6 C&S
YBi = (x-xe)*sb + (y+ye)*cb ;

% Fix roundoff errors in Ybi and Yb from trig function problems
bad = abs(YBi)<1e-10; %Steve Martels Fix
bad2 = abs(YB)<1e-10;
YBi(bad) = 0;
YB(bad2) = 0;
clear bad bad2


% Calculate derivatives of the function f(x,y), eq. 5.2.5 of C&S, p. 81. 
%   which are used to calculate the displacement and stress components. 
% It is understood that X and Y refer to XB and YB.
% First abbreviate repeated terms in the derivatives of f(x,y):
Y2 = YB.^2;
XMa = XB-a; XPa = XB+a; 
XMa2 = XMa.^2; XPa2 = XPa.^2; 
R1S = XMa2 + Y2; 
R2S = XPa2 + Y2; 

% Same thing for the image dislocation
Y2i = YBi.^2;
XMai = XBi-a; XPai = XBi+a; 
XMa2i = XMai.^2; XPa2i = XPai.^2; 
R1Si = XMa2i + Y2i; R1S2i = R1Si.^2; 
R2Si = XPa2i + Y2i; R2S2i = R2Si.^2;

% The following derivatives are eqs. 4.5.5a thru d of C&S, p. 58.
FF2 = con*(log(sqrt(R1S)) - log(sqrt(R2S)));

% Flag all observation points for which Yb is not 0    
i1 = find(YB);
% Flag any observation points on the element
i2 = find(YB == 0 & abs(XB) < a);
i3 = find( YB==0 & abs(XB) > a );
% Steve Martels Solution to elements lying on same plane
% FB3 = 0 for pts colinear with element, CON*pi for pts. on element 
% FB3 = difference of arc tangents for all other pts.
FF3=zeros(size(x));
FF3(i1) = atan2(YB(i1),XMa(i1)) - atan2(YB(i1),XPa(i1)); 
FF3(i2) = pi.*ones(size(i2));
FF3(i3) = zeros(size(i3));
FF3 = -con.*(FF3)';	
FF3 = reshape(FF3,length);

FF4 = con*(YB./R1S - YB./R2S); 
FF5 = con*(XMa./R1S - XPa./R2S);

% Calculate intermediate functions Fnifor the image dislocation
FF2i = con*(log(sqrt(R1Si)) - log(sqrt(R2Si)) ); %Equations 4.5.5 C&S

% Flag all observation points for which Yb is not 0    
i1i = find(YBi);
% Flag any observation points on the element
i2i = find(YBi == 0 & abs(XBi) < a);
i3i = find( YBi==0 & abs(XBi) > a );
% Steve Martels Solution to elements lying on same plane
% FB3 = 0 for pts colinear with element, CON*pi for pts. on element 
% FB3 = difference of arc tangents for all other pts.
FF3i=zeros(size(x));
FF3i(i1i) = atan2(YBi(i1i),XMai(i1i)) - atan2(YBi(i1i),XPai(i1i)); 
FF3i(i2i) = pi.*ones(size(i2i));
FF3i(i3i) = zeros(size(i3i));
FF3i = -con.*(FF3i)';	
FF3i = reshape(FF3i,length);

FF4i = con*(YBi./R1Si - YBi./R2Si); 
FF5i = con*(XMai./R1Si - XPai./R2Si);


% The halfspace examples of eqs. 5.5.3a and b of C&S, p. 91.
% See Appendix A of: Martel, S.J. and Langley, J.S., 2006. Propagation of
% normal faults to the surface in basalt, Koae fault system, Hawaii.
% Journal of Structural Geology, 28(12), pp.2123-2143.
FF6i = con*((XMa2i - Y2i)./R1S2i - (XPa2i - Y2i)./R2S2i);
FF7i = 2*con*YBi.*(XMai./R1S2i - XPai./R2S2i);


% Define material constants used in calculating displacements.
pr1 = 1-2*nu; pr2 = 2.*(1-nu); pr3=3-4.*nu;%pr3 = 1-pr;
% Calculate the displacement components using eqs. 5.5.4 of C&S, p. 91.
Ux = Dxb*(-pr1*sb*FF2 + pr2*cb*FF3 + YB.*(sb*FF4 - cb*FF5))...
    +Dyb*(-pr1*cb*FF2 - pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5));
Uy = Dxb*(+pr1*cb*FF2 + pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5))...
    +Dyb*(-pr1*sb*FF2 + pr2*cb*FF3 - YB.*(sb*FF4 - cb*FF5));


%  See equations 7.4.8 and 7.4.9 in Crouch and Starfield
%  Calculate IMAGE AND SUPPLEMENTAL DISPLACEMENT components due to unit SHEAR 
%  displacement discontinuity
Uxi_s =Dxb* (pr1.*sb.*FF2i - pr2.*cb.*FF3i +....
    (pr3.*(y.*s2b - YB.*sb) + 2.*y.*s2b).*FF4i +...
    (pr3.*(y.*c2b - YB.*cb) - y.*(1-2.*c2b)).*FF5i +...
    2.*y.*(y.*s3b - YB.*s2b).*FF6i -....
    2.*y.*(y.*c3b - YB.*c2b).*FF7i);

Uyi_s =Dxb* (-pr1.*cb.*FF2i - pr2.*sb.*FF3i -....
    (pr3.*(y.*c2b - YB.*cb) + y.*(1-2.*c2b)).*FF4i +...  
    (pr3.*(y.*s2b - YB.*sb) - 2.*y.*s2b).*FF5i +...		 
    2.*y.*(y.*c3b - YB.*c2b).*FF6i +....
    2.*y.*(y.*s3b - YB.*s2b).*FF7i);     

%  Calculate IMAGE AND SUPPLEMENTAL DISPLACEMENT components due to unit NORMAL 
%  displacement discontinuity
Uxi_n =Dyb* ( pr1.*cb.*FF2i + pr2.*sb.*FF3i -....
    (pr3.*(y.*c2b - YB.*cb) - y).*FF4i +...
    pr3.*(y.*s2b - YB.*sb).*FF5i -...
    2.*y.*(y.*c3b - YB.*c2b).*FF6i -....
    2.*y.*(y.*s3b - YB.*s2b).*FF7i);    

Uyi_n =Dyb* (pr1.*sb.*FF2i - pr2.*cb.*FF3i -....
    pr3.*(y.*s2b - YB.*sb).*FF4i -...
    (pr3.*(y.*c2b - YB.*cb) + y).*FF5i +...
    2.*y.*(y.*s3b - YB.*s2b).*FF6i -....
    2.*y.*(y.*c3b - YB.*c2b).*FF7i);


Uxi=Uxi_s+Uxi_n;
Uyi=Uyi_s+Uyi_n;

Ux=Ux+Uxi;
Uy=Uy+Uyi;

Disp=[Ux(:),Uy(:)];
