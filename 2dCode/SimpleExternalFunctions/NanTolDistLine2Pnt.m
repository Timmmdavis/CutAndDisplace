function [X,Y]=NanTolDistLine2Pnt( X,Y,P1,P2,LineNormalVector,Distance )
% NanTolDistLine2Pnt: Nans points within a set distance of a line. Works for
%                   multiple line segments and input points.
%
%  .out
%    -----------------
% |       .in          |
%|    ~~~~~~~~~~~       |
% |                    |
%    -----------------
%               
% usage #1:
% [X,Y]=NanTolDistLine2Pnt( X,Y,P1,P2,MidPoint,LineNormalVector,Distance )
%
%
% Arguments: (input)
%    X,Y           - The query points that will be set to NAN if within
%                    the distance.
%
%  MidPoint        - The x and y locations of each elements midpoint [x,y]. 
%
%  P1,P2           - The start (P1) and end (P2) points of each separate 
%                   element in [x,y] coordiantes. First col is x. 
%
% LineNormalVector - The direction cosines, CosAx (Nx) and CosAy in a list
%                   for each element. 
%
% Distance         - Set distance around each line segment that is removed.
%                   
%
%
%
% Arguments: (output)
%    X,Y           - The query points that will have been set to NAN if within
%                    the distance.
%
%
% Example usage:
%
% N/A.
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


for i = 1:numel(P1(:,1))
    
    %Get extra X and Y bits:
    ExX=(LineNormalVector(i,1)*Distance);
    ExY=(LineNormalVector(i,2)*Distance);

    %Init bounds:
    P1Plus=P1;  P1Neg=P1;
    P2Plus=P2;  P2Neg=P2;

    %Move Points to create rectangular bounds:
    %X    
    P1Plus(i,1)=P1(i,1)+ExX;
    P2Plus(i,1)=P2(i,1)+ExX;
    P1Neg(i,1) =P1(i,1)-ExX;    
    P2Neg(i,1) =P2(i,1)-ExX;
    %Y
    P1Plus(i,2)=P1(i,2)+ExY;
    P2Plus(i,2)=P2(i,2)+ExY;
    P1Neg(i,2) =P1(i,2)-ExY;    
    P2Neg(i,2) =P2(i,2)-ExY;

    %Accumulate as rectanglei
    XRect=[P1Plus(i,1),P1Neg(i,1),P2Neg(i,1),P2Plus(i,1)];
    YRect=[P1Plus(i,2),P1Neg(i,2),P2Neg(i,2),P2Plus(i,2)];

    %To draw the current retangle:
    %hold on
    %fill(XRect,YRect,'r','FaceAlpha',0.5)
    
    %Now check if point is within the bound: 
    in = inpolygon(X,Y,XRect,YRect); 
    
    %Now check if point is within radius of first tip
    Xr=X-P1(i,1);
    Yr=Y-P1(i,2);
    [~,r1] = cart2pol(Xr,Yr); 

    %Now check if point is within radius of second tip
    Xr=X-P2(i,1);
    Yr=Y-P2(i,2);
    [~,r2] = cart2pol(Xr,Yr); 
    
    %Nan Points inside
    X(in)=nan;
    Y(in)=nan;
    %Nan Points within distance of tip
    X(r1<Distance)=nan;
    Y(r1<Distance)=nan;
    X(r2<Distance)=nan;
    Y(r2<Distance)=nan;    
    
end

if sum(isnan(X(:)))==numel(X) 
    disp('Function DistLine2Pnt has set all points to NaN, the distance tolerance used it too high.')
elseif sum(isnan(X(:)))>0
    disp('Function DistLine2Pnt has set some points to NaN')    
end

end

