function [StressDisp] = ...
coeffhs_func(x,y,xe,ye,a,Beta,Sd,Nd,pr,E)
%   Slightly Modified form of Steve Martels BEM scripts from:
%   http://www.soest.hawaii.e.du/martel/Martel.BEM_dir/
%   Also includes some notation from Ritz & Pollard 2012: 
%   Copyright 2017, Tim Davis, The University of Aberdeen


if any(y>0)
    error('Half-space solution: Z coordinates must be negative!')
end
% The shear modulus, sm, is related to the prescribed elastic constants.
    sm = E/(2*(1+pr));
% Define material constant used in calculating influence coefficients.
    con = 1/(4*pi*(1-pr)); 

    H=a;
    Dxb = Sd; Dyb = Nd;
    sb = sin(Beta); cb = cos(Beta);
    s2b = sin(2*Beta); c2b = cos(2*Beta);
    s3b = sin(3*Beta); c3b = cos(3*Beta);
    s4b = sin(4*Beta); c4b = cos(4*Beta);
    
% Define array of local coordinates for the observation grid relative to
%   the midpoint and orientation of the ith element.
% Refer to (Figure 5.6, C&S, p. 91) and eqs. 4.5.1 of C&S, p. 57. 
    XB = (x-xe)*cb + (y-ye)*sb;
    YB = -(x-xe)*sb + (y-ye)*cb;
    
    length=size(x);
    lengthrow=length(1,1);
    lengthcol=length(1,2);

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
    XMH = XB-H; XPH = XB+H; 
    XMH2 = XMH.^2; XPH2 = XPH.^2; 
    R1S = XMH2 + Y2; R1S2 = R1S.^2; 
    R2S = XPH2 + Y2; R2S2 = R2S.^2;
    
    % Same thing for the image dislocation
    Y2i = YBi.^2;
    XMHi = XBi-H; XPHi = XBi+H; 
    XMH2i = XMHi.^2; XPH2i = XPHi.^2; 
    R1Si = XMH2i + Y2i; R1S2i = R1Si.^2; 
    R2Si = XPH2i + Y2i; R2S2i = R2Si.^2;
    
% The following derivatives are eqs. 4.5.5a thru d of C&S, p. 58.
    FF2 = con*(log(sqrt(R1S)) - log(sqrt(R2S)));
    
    % Flag all observation points for which Yb is not 0    
    i1 = find(YB);
    % Flag any observation points on the element
    i2 = find(YB == 0 & abs(XB) < H);
    i3 = find( YB==0 & abs(XB) > H );
    % Steve Martels Solution to elements lying on same plane
    % FB3 = 0 for pts colinear with element, CON*pi for pts. on element 
    % FB3 = difference of arc tangents for all other pts.
    FF3(i1) = atan2(YB(i1),XMH(i1)) - atan2(YB(i1),XPH(i1)); 
    FF3(i2) = pi.*ones(size(i2));
    FF3(i3) = zeros(size(i3));
    FF3 = -con.*(FF3)';	
    FF3 = reshape(FF3,lengthrow,lengthcol);
    
    FF4 = con*(YB./R1S - YB./R2S); 
    FF5 = con*(XMH./R1S - XPH./R2S);
% The following derivatives are eqs. 5.5.3a and b of C&S, p. 91.
    FF6 = con*((XMH2 - Y2)./R1S2 - (XPH2 - Y2)./R2S2);
    FF7 = 2*con*YB.*(XMH./R1S2 - XPH./R2S2);
    
% Calculate intermediate functions Fnifor the image dislocation
    FF2i = con*(log(sqrt(R1Si)) - log(sqrt(R2Si)) ); %Equations 4.5.5 C&S

    % Flag all observation points for which Yb is not 0    
    i1i = find(YBi);
    % Flag any observation points on the element
    i2i = find(YBi == 0 & abs(XBi) < H);
    i3i = find( YBi==0 & abs(XBi) > H );
    % Steve Martels Solution to elements lying on same plane
    % FB3 = 0 for pts colinear with element, CON*pi for pts. on element 
    % FB3 = difference of arc tangents for all other pts.
    FF3i(i1i) = atan2(YBi(i1i),XMHi(i1i)) - atan2(YBi(i1i),XPHi(i1i)); 
    FF3i(i2i) = pi.*ones(size(i2i));
    FF3i(i3i) = zeros(size(i3i));
    FF3i = -con.*(FF3i)';	
    FF3i = reshape(FF3i,lengthrow,lengthcol);
    
    FF4i = con*(YBi./R1Si - YBi./R2Si); 
    FF5i = con*(XMHi./R1Si - XPHi./R2Si);
% The halfspace examples of eqs. 5.5.3a and b of C&S, p. 91.
    FF6i = con*((XMH2i - Y2i)./R1S2i - (XPH2i - Y2i)./R2S2i);
    FF7i = 2*con*YBi.*(XMHi./R1S2i - XPHi./R2S2i);
    FF8i = 2*con*YBi.*((Y2i-3*XMH2i)./R1Si.^3-(Y2i - 3.*XPH2i)./R2Si.^3);
    FF9i = 2*con*((XMHi.^3-3.*Y2i)./R1Si.^3-(XPHi.^3-3.*Y2i)./R2Si.^3);

    
% Calculate the stress components using eqs. 5.5.5 of C&S, p. 92.
    Sxx =  2*sm*Dxb*(2*cb*cb*FF4 + s2b*FF5 + YB.*(c2b*FF6-s2b*FF7))...
         +2*sm*Dyb*(-FF5 + YB.*(s2b*FF6 + c2b*FF7));
    Syy =  2*sm*Dxb*(2*sb*sb*FF4 - s2b*FF5 - YB.*(c2b*FF6-s2b*FF7))...
         +2*sm*Dyb*(-FF5 - YB.*(s2b*FF6 + c2b*FF7));
    Sxy =  2*sm*Dxb*(s2b*FF4 - c2b*FF5 + YB.*(s2b*FF6+c2b*FF7))...
         +2*sm*Dyb*(-YB.*(c2b*FF6 - s2b*FF7)); 
     
%  Calculate IMAGE AND SUPPLEMENTAL STRESS components due to unit SHEAR 
%  displacement discontinuity
    Sxxi_s = 2*sm*Dxb*(FF4i - 3.*(c2b.*FF4i - s2b.*FF5i) +....
	(2.*y.*(cb - 3.*c3b) + 3.*YB.*c2b).*FF6i +.....
	(2.*y.*(sb - 3.*s3b) + 3.*YB.*s2b).*FF7i -.....
	2.*y.*(y.*c4b - YB.*c3b).*FF8i -....
	2.*y.*(y.*s4b - YB.*s3b).*FF9i);

    Syyi_s = 2*sm*Dxb*(FF4i - (c2b.*FF4i - s2b.*FF5i) -....
	(4.*y.*sb.*s2b - YB.*c2b).*FF6i +.....
 	(4.*y.*sb.*c2b + YB.*s2b).*FF7i +....
 	2.*y.*(y.*c4b - YB.*c3b).*FF8i +....
 	2.*y.*(y.*s4b - YB.*s3b).*FF9i);

    Sxyi_s = 2*sm*Dxb*(s2b.*FF4i + c2b.*FF5i +.... 
	(2.*y.*s2b.*(1+4.*c2b) - YB.*s2b).*FF6i +.....
	(2.*y.*cb.*(3-4.*c2b) + YB.*c2b).*FF7i +..... 
	2.*y.*(y.*s4b - YB.*s3b).*FF8i -.....
	2.*y.*(y.*c4b - YB.*c3b).*FF9i);


%  Calculate IMAGE AND SUPPLEMENTAL STRESS components due to unit NORMAL
%  displacement discontinuity
    Sxxi_n = 2*sm*Dyb*(FF5i + (2.*y.*(sb - 2.*s3b) +.... 
	3.*YB.*s2b).*FF6i - (2.*y.*(cb - 2.*c3b) +....
	3.*YB.*c2b).*FF7i - 2.*y.*(y.*s4b - YB.*s3b).*FF8i +....
	2.*y.*(y.*c4b - YB.*c3b).*FF9i);

    Syyi_n = 2*sm*Dyb*(FF5i - (2.*y.*sb - YB.*s2b).*FF6i +.... 
	(2.*y.*cb - YB.*c2b).*FF7i +....
	2.*y.*(y.*s4b - YB.*s3b).*FF8i -.... 
	2.*y.*(y.*c4b - YB.*c3b).*FF9i);

    Sxyi_n = 2*sm*Dyb*((4.*y.*sb.*s2b + YB.*c2b).*FF6i -.... 
	(4.*y.*sb.*c2b - YB.*s2b).*FF7i -.....
	2.*y.*(y.*c4b - YB.*c3b).*FF8i -.... 
	2.*y.*(y.*s4b - YB.*s3b).*FF9i);   


    Sxxi=Sxxi_s+Sxxi_n;
    Syyi=Syyi_s+Syyi_n;
    Sxyi=Sxyi_s+Sxyi_n;
    
    Sxx=Sxx+Sxxi;
    Syy=Syy+Syyi;
    Sxy=Sxy+Sxyi;
     
% Define material constants used in calculating displacements.
    pr1 = 1-2*pr; pr2 = 2.*(1-pr); pr3=3-4.*pr;%pr3 = 1-pr;
% Calculate the displacement components using eqs. 5.5.4 of C&S, p. 91.
    UX = Dxb*(-pr1*sb*FF2 + pr2*cb*FF3 + YB.*(sb*FF4 - cb*FF5))...
        +Dyb*(-pr1*cb*FF2 - pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5));
    UY = Dxb*(+pr1*cb*FF2 + pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5))...
        +Dyb*(-pr1*sb*FF2 + pr2*cb*FF3 - YB.*(sb*FF4 - cb*FF5));

    
    %  See equations 7.4.8 and 7.4.9 in Crouch and Starfield
    %  Calculate IMAGE AND SUPPLEMENTAL DISPLACEMENT components due to unit SHEAR 
    %  displacement discontinuity
    UXi_s =Dxb* (pr1.*sb.*FF2i - pr2.*cb.*FF3i +....
            (pr3.*(y.*s2b - YB.*sb) + 2.*y.*s2b).*FF4i +...
            (pr3.*(y.*c2b - YB.*cb) - y.*(1-2.*c2b)).*FF5i +...
            2.*y.*(y.*s3b - YB.*s2b).*FF6i -....
            2.*y.*(y.*c3b - YB.*c2b).*FF7i);
        
    UYi_s =Dxb* (-pr1.*cb.*FF2i - pr2.*sb.*FF3i -....
            (pr3.*(y.*c2b - YB.*cb) + y.*(1-2.*c2b)).*FF4i +...  
            (pr3.*(y.*s2b - YB.*sb) - 2.*y.*s2b).*FF5i +...		 
            2.*y.*(y.*c3b - YB.*c2b).*FF6i +....
            2.*y.*(y.*s3b - YB.*s2b).*FF7i);     
        
    %  Calculate IMAGE AND SUPPLEMENTAL DISPLACEMENT components due to unit NORMAL 
    %  displacement discontinuity
    UXi_n =Dyb* ( pr1.*cb.*FF2i + pr2.*sb.*FF3i -....
            (pr3.*(y.*c2b - YB.*cb) - y).*FF4i +...
            pr3.*(y.*s2b - YB.*sb).*FF5i -...
            2.*y.*(y.*c3b - YB.*c2b).*FF6i -....
            2.*y.*(y.*s3b - YB.*s2b).*FF7i);    
        
    UYi_n =Dyb* (pr1.*sb.*FF2i - pr2.*cb.*FF3i -....
            pr3.*(y.*s2b - YB.*sb).*FF4i -...
            (pr3.*(y.*c2b - YB.*cb) + y).*FF5i +...
            2.*y.*(y.*s3b - YB.*s2b).*FF6i -....
            2.*y.*(y.*c3b - YB.*c2b).*FF7i);
    
    
    UXi=UXi_s+UXi_n;
    UYi=UYi_s+UYi_n;
    
    UX=UX+UXi;
    UY=UY+UYi;
        
 StressDisp=[Sxx(:),Syy(:),Sxy(:),UX(:),UY(:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
