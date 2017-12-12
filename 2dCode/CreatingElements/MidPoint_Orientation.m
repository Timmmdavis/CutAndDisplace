function [ MidPoint,HalfLength,P1,P2,LineNormalVector ] = MidPoint_Orientation( XBEG,XEND,YBEG,YEND )
% MidPoint_Orientation: From the start and end points of each element in X
%                   and Y this function calculates the midpoint locations,
%                   halflength and the direction cosines of each element.
%               
% usage #1:
% [ MidPoint,HalfLength,P1,P2,LineNormalVector]...
%  = MidPoint_Orientation( XBEG,XEND,YBEG,YEND )
%
%
% Arguments: (input)
% XBEG,XEND        - The start and end points in x of each
%                   separate element
%
% YBEG,YEND        - The start and end points in y of each
%                   separate element
%
% Arguments: (output)
%  MidPoint        - The x and y locations of each elements midpoint [x,y]. 
%
%  HalfLength      - The half length of each element
%
%  P1,P2           - The start (P1) and end (P2) points of each separate 
%                   element in [x,y] coordiantes. First col is x. 
%
% LineNormalVector - The direction cosines, CosAx (Nx) and CosAy in a list
%                   for each element. 
%
%
% Example usage:
%
%[ xe,ye,HalfLength,Points,LineNormalVector ] = MidPoint_Orientation( XBEG,XEND,YBEG,YEND );
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Creating midpoints and creating element orientations from start and
%end points of each element
xe = (XBEG + XEND)/2;				% element midpoints
ye = (YBEG + YEND)/2;				% element midpoints
MidPoint=[xe,ye];
%Collating
P1=[XBEG,YBEG];  
P2=[XEND,YEND]; 

xd = XEND-XBEG;
yd = YEND-YBEG;
HalfLength = sqrt(xd.*xd +yd.*yd)/2;			% half-length of each element
%In the lines below, B = beta = orientation of element relative to global x-axis
Beta = atan2(yd,xd);                        % This is the Beta Angle used in the Crouch and Starfield functions
LineNormalVector=[sin(Beta),-cos(Beta)];  


  
end

