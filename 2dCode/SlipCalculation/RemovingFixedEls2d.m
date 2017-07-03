function [XBEG,XEND,YBEG,YEND,NUM,x,y,xe,ye,HalfLength,Beta,CosB,Points,NormAng]...
    = RemovingFixedEls2d(Fdisp,XBEG,XEND,YBEG,YEND,NUM)
%% Removing any fixed elements
%Removes using Fdisp logical flag. The element can be removed 
% as the solution as we have used these for thier purpose
% We know they have no displacement so they will not contribute to anything from now on.
% Therefore its visually better if these are not around. 

%   Copyright 2017, Tim Davis, The University of Aberdeen
FD=logical(Fdisp);
XBEG=XBEG(~FD,:);
XEND=XEND(~FD,:);
YBEG=YBEG(~FD,:);
YEND=YEND(~FD,:);
Sm = sum(Fdisp);
NUM=NUM-Sm;
[ x,y,xe,ye,HalfLength,Beta,CosB,Points,NormAng ] = MidPoint_Orientation( XBEG,XEND,YBEG,YEND,NUM );

end

