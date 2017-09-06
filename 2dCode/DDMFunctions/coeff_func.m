function [StressDisp] = ...
coeff_func(x,y,xe,ye,a,Beta,Sd,Nd,pr,E)
%   Slightly Modified form of Steve Martels BEM scripts from:
%   http://www.soest.hawaii.e.du/martel/Martel.BEM_dir/
%   Also includes some notation from Ritz & Pollard 2012: 
%   Copyright 2017, Tim Davis, The University of Aberdeen


% Lengths for reshaping
    length=size(x);
    lengthrow=length(1,1);
    lengthcol=length(1,2);
% The shear modulus, sm, is related to the prescribed elastic constants.
    sm = E/(2*(1+pr));
% Define material constant used in calculating influence coefficients.
    con = 1/(4*pi*(1-pr)); 
    cons = 2*sm;%E/(1+pr);


    H=a;
    Dxb = Sd; Dyb = -Nd; 
    sb = sin(Beta); cb = cos(Beta);
    s2b = sin(2*Beta); c2b = cos(2*Beta);
% Define array of local coordinates for the observation grid relative to
%   the midpoint and orientation of the ith element.
% Refer to (Figure 5.6, C&S, p. 91) and eqs. 4.5.1 of C&S, p. 57. 
    XB = (x-xe)*cb + (y-ye)*sb;
    YB = -(x-xe)*sb + (y-ye)*cb;

    % Flag all observation points for which Yb is not 0    
    i1 = find(YB);
    % Flag any observation points on the element
    i2 = find(YB == 0 & abs(XB) < H);
    i3 = find(YB == 0 & abs(XB) > H ); 
      
% Calculate derivatives of the function f(x,y), eq. 5.2.5 of C&S, p. 81. 
%   which are used to calculate the displacement and stress components. 
% It is understood that X and Y refer to XB and YB.
% First abbreviate repeated terms in the derivatives of f(x,y):
    Y2 = YB.^2;
    XMH = XB-H; XPH = XB+H; 
    XMH2 = XMH.^2; XPH2 = XPH.^2; 
    R1S = XMH2 + Y2; R1S2 = R1S.^2; 
    R2S = XPH2 + Y2; R2S2 = R2S.^2;
% The following derivatives are eqs. 4.5.5a thru d of C&S, p. 58.
    FF2 = con*(log(sqrt(R1S)) - log(sqrt(R2S)));
    
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
    

% Calculate the stress components using eqs. 5.5.5 of C&S, p. 92.
    SXX =  cons*Dxb*(2*(cb*cb)*FF4 + s2b*FF5 + YB.*(c2b*FF6-s2b*FF7))...
         +cons*Dyb*(-FF5 + YB.*(s2b*FF6 + c2b*FF7));
    SYY =  cons*Dxb*(2*(sb*sb)*FF4 - s2b*FF5 - YB.*(c2b*FF6-s2b*FF7))...
         +cons*Dyb*(-FF5 - YB.*(s2b*FF6 + c2b*FF7));
    SXY =  cons*Dxb*(s2b*FF4 - c2b*FF5 + YB.*(s2b*FF6+c2b*FF7))...
         +cons*Dyb*(-YB.*(c2b*FF6 - s2b*FF7)); 



     
% Define material constants used in calculating displacements.
    pr1 = 1-2*pr; pr2 = 2-2*pr;
% Calculate the displacement components using eqs. 5.5.4 of C&S, p. 91.
    UX = Dxb*(-pr1*sb*FF2 + pr2*cb*FF3 + YB.*(sb*FF4 - cb*FF5))...
        +Dyb*(-pr1*cb*FF2 - pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5));
    UY = Dxb*(+pr1*cb*FF2 + pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5))...
        +Dyb*(-pr1*sb*FF2 + pr2*cb*FF3 - YB.*(sb*FF4 - cb*FF5));
    

 %Stress=[SXX,SYY,SXY];
 StressDisp=[SXX(:),SYY(:),SXY(:),UX(:),UY(:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
