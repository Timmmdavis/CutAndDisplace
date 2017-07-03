function [Triangles,FaceNormalVector,MidPoint,P1,P2,P3]...
    = RemovingFixedEls3d(Triangles,FaceNormalVector,MidPoint,P1,P2,P3,Fdisp)
%% Removing any fixed elements
%Removes using Fdisp logical flag. The element can be removed 
% as the solution as we have used these for thier purpose
% We know they have no displacement so they will not contribute to anything from now on.
% Therefore its visually better if these are not around. 

%   Copyright 2017, Tim Davis, The University of Aberdeen
FD=logical(Fdisp);
Triangles=Triangles(~FD,:);
FaceNormalVector=FaceNormalVector(~FD,:);
MidPoint=MidPoint(~FD,:);
P1=P1(~FD,:);
P2=P2(~FD,:);
P3=P3(~FD,:);

end

