function [Disp] = LDdispFS(x,y,xe,ye,a,Beta,Ds,Dn,nu)
% LDdispFS: LineDislocationInducedDisplacementFullSpace, Computes the 
%               influence (displacement) of a single planar line crack
%               shearing and/or opening on the surronding points. (full
%               space).
%
% usage #1:
%[Disp] = LDdispFS(x,y,xe,ye,a,Beta,Ds,Dn,nu)
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
% Disp        - Is the displacement caused by the movement of the
%             dislocation at the observataion points. [Ux,Uy].
%
% Example usage:
%
%  [x,y]=meshgrid(-2:0.1:2,-2:0.1:2);
%  dimx = length(x(:,1)); 
%  dimy = length(x(1,:));   
%  Ds=0; Dn=1;
%  [Disp] = LDdispFS(x(:),y(:),0,0,1,0,Ds,Dn,0.25);
%  [ Ux,Uy ] = ExtractCols( Disp );
%  [Ux,Uy]=ReshapeData2d( dimx,dimy,Ux,Uy );
%  figure,quiver(x(:),y(:),Ux(:),Uy(:));
% 
%  Author: Tim Davis/Steve Martel
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen
%  Modified form of Steve Martels BEM scripts from:
%  http://www.soest.hawaii.e.du/martel/Martel.BEM_dir/
%  Also includes some notation from Ritz & Pollard 2012.


%Lengths for reshaping
length=size(x);
lengthrow=length(1,1);
lengthcol=length(1,2);
%Define material constant used in calculating influence coefficients.
con = 1/(4*pi*(1-nu)); 
H=a;
Dxb = Ds; Dyb = -Dn; 
sb = sin(Beta); cb = cos(Beta);
%Define array of local coordinates for the observation grid relative to
%the midpoint and orientation of the ith element.
%Refer to (Figure 5.6, C&S, p. 91) and eqs. 4.5.1 of C&S, p. 57. 
XB = (x-xe)*cb + (y-ye)*sb;
YB = -(x-xe)*sb + (y-ye)*cb;

%Flag all observation points for which Yb is not 0    
i1 = find(YB);
%Flag any observation points on the element
i2 = find(YB == 0 & abs(XB) < H);
i3 = find(YB == 0 & abs(XB) > H ); 

% %Calculate derivatives of the function f(x,y), eq. 5.2.5 of C&S, p. 81. 
% %which are used to calculate the displacement and stress components. 
% %It is understood that X and Y refer to XB and YB.
% %First abbreviate repeated terms in the derivatives of f(x,y):
% Y2 = YB.^2;
% XMH = XB-H; XPH = XB+H; 
% XMH2 = XMH.^2; XPH2 = XPH.^2; 
% R1S = XMH2 + Y2; 
% R2S = XPH2 + Y2; 
% %The following derivatives are eqs. 4.5.5a thru d of C&S, p. 58.
% FF2 = con*(log(sqrt(R1S)) - log(sqrt(R2S)));
% 
% %Steve Martels Solution to elements lying on same plane
% %FB3 = 0 for pts colinear with element, CON*pi for pts. on element 
% %FB3 = difference of arc tangents for all other pts.
% FF3(i1) = atan2(YB(i1),XMH(i1)) - atan2(YB(i1),XPH(i1)); 
% FF3(i2) = pi.*ones(size(i2));
% FF3(i3) = zeros(size(i3));
% FF3 = -con.*(FF3)';	
% FF3 = reshape(FF3,lengthrow,lengthcol);
% FF4 = con*(YB./R1S - YB./R2S); 
% FF5 = con*(XMH./R1S - XPH./R2S);
% 
% 
% %Define material constants used in calculating displacements.
% pr1 = 1-2*nu; pr2 = 2-2*nu;
% %Calculate the displacement components using eqs. 5.5.4 of C&S, p. 91.
% Ux = Dxb*(-pr1*sb*FF2 + pr2*cb*FF3 + YB.*(sb*FF4 - cb*FF5))...
%     +Dyb*(-pr1*cb*FF2 - pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5));
% Uy = Dxb*(+pr1*cb*FF2 + pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5))...
%     +Dyb*(-pr1*sb*FF2 + pr2*cb*FF3 - YB.*(sb*FF4 - cb*FF5));

%Calculate derivatives of the function f(x,y), eq. 5.2.5 of C&S, p. 81. 
%which are used to calculate the displacement and stress components. 
%It is understood that X and Y refer to XB and YB.
%First abbreviate repeated terms in the derivatives of f(x,y):
Y2 = YB.^2;
XMH = XB-H; XPH = XB+H; 
XMH2 = XMH.^2; XPH2 = XPH.^2; 
R1S = XMH2 + Y2; 
R2S = XPH2 + Y2; 
%The following derivatives are eqs. 4.5.5a thru d of C&S, p. 58.
FF2 = (log(sqrt(R1S)) - log(sqrt(R2S)));

%Steve Martels Solution to elements lying on same plane
%FB3 = 0 for pts colinear with element, CON*pi for pts. on element 
%FB3 = difference of arc tangents for all other pts.
FF3(i1) = atan2(YB(i1),XMH(i1)) - atan2(YB(i1),XPH(i1)); 
FF3(i2) = pi.*ones(size(i2));
FF3(i3) = zeros(size(i3));
FF3 = -(FF3)';	
FF3 = reshape(FF3,lengthrow,lengthcol);
FF4 = (YB./R1S - YB./R2S); 
FF5 = (XMH./R1S - XPH./R2S);


%Define material constants used in calculating displacements.
pr1 = 1-2*nu; pr2 = 2*(1-nu);
%Calculate the displacement components using eqs. 5.5.4 of C&S, p. 91.
Ux = Dxb*(-pr1*sb*FF2 + pr2*cb*FF3 + YB.*(sb*FF4 - cb*FF5))...
    +Dyb*(-pr1*cb*FF2 - pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5));
Uy = Dxb*(+pr1*cb*FF2 + pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5))...
    +Dyb*(-pr1*sb*FF2 + pr2*cb*FF3 - YB.*(sb*FF4 - cb*FF5));



Disp=[Ux(:),Uy(:)]*con;

