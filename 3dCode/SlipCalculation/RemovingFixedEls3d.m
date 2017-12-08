function [Triangles,FaceNormalVector,MidPoint,P1,P2,P3]...
    = RemovingFixedEls3d(Triangles,FaceNormalVector,MidPoint,P1,P2,P3,Fdisp)
% RemovingFixedEls3d: Removes the 'fixed' elements as the
%               displacement on these is not calculated. Simply the
%               calculation assumed the midpoints did not displace.
%
% usage #1:
% [Triangles,FaceNormalVector,MidPoint,P1,P2,P3]...
%     = RemovingFixedEls3d(Triangles,FaceNormalVector,MidPoint,P1,P2,P3,Fdisp)
%
% Arguments: (input)
%     Fdisp     - Flag where ones represent elements that have midpoints
%                 that do not displace. 
%
% Triangles     -  Triangles is a list where each row contains 3 index
%                    locations in "Points" which contains the XYZ location
%                    of each corner of the triangle.
%
% FaceNormalVector  - The normal vectors n*3, (CosAx,CosAy,CosAz) of each
%                    triangle.
%
% MidPoint    - The midpoints n*3, (XYZ) of each triangle.
%
% P1,P2,P3    - n*3 Column vectors where each 'P' represents the
%               different corner points of one of the triangles (XYZ).
%
% Arguments: (output)
%       As Inputs but elements that were 'fixed' have been removed. 
%
% Example usage:
%
% [Triangles,FaceNormalVector,MidPoint,P1,P2,P3]...
%     = RemovingFixedEls3d(Triangles,FaceNormalVector,MidPoint,P1,P2,P3,Fdisp)
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Get logical flag
FD=logical(Fdisp);

%Remove bits. 
Triangles=Triangles(~FD,:);
FaceNormalVector=FaceNormalVector(~FD,:);
MidPoint=MidPoint(~FD,:);
P1=P1(~FD,:);
P2=P2(~FD,:);
P3=P3(~FD,:);

end

