function [MidPoint,HalfLength,P1,P2,LineNormalVector]...
    = RemovingFixedEls2d(Fdisp,P1,P2)
% RemovingFixedEls2d: Removes the 'fixed' elements as the
%               displacement on these is not calculated. Simply the
%               calculation assumed thier midpoints did not displace.
%
% usage #1:
% [MidPoint,HalfLength,Points,LineNormalVector]...
%     = RemovingFixedEls2d(Fdisp,Points)
%
% Arguments: (input)
%     Fdisp     - Flag where ones represent elements that have midpoints
%                 that do not displace. 
%
%  P1,P2           - The start (P1) and end (P2) points of each separate 
%                   element in [x,y] coordiantes. First col is x. 
%
%
% Arguments: (output)
%    MidPoint      - The x and y locations of each elements midpoint. 
%
%  HalfLength      - The half length of each element
%
%  P1,P2           - The start (P1) and end (P2) points of each separate 
%                   element in [x,y] coordiantes. First col is x. 
%
% LineNormalVector - The direction cosines, CosAx (Nx) and CosAy in a list
%                   for each element. 
%
% Example usage:
%
% [MidPoint,HalfLength,P1,P2,LineNormalVector]...
%     = RemovingFixedEls2d(Fdisp,P1,P2)
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Extract the beginning and end points of each element.
XBEG=P1(:,1);
XEND=P2(:,1);
YBEG=P1(:,2);
YEND=P2(:,2);

%Make sure Fdisp is a logical
FD=logical(Fdisp);

%Remove the points at elements that are fixed
XBEG=XBEG(~FD,:);
XEND=XEND(~FD,:);
YBEG=YBEG(~FD,:);
YEND=YEND(~FD,:);

%Recompute the element orientations for the elements that do displace. 
[ MidPoint,HalfLength,P1,P2,LineNormalVector ] = MidPoint_Orientation( XBEG,XEND,YBEG,YEND );

end

