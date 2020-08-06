function [X,Y,Z]=NanTolDistTri2Pnt( X,Y,Z,P1,P2,P3,MidPoint,FaceNormalVector,Distance )
% DistTri2Pnt: Nans points within a set distance of a triangle. Works for
%                   multiple triangles and input points.
%
%               
% usage #1:
% [X,Y,Z]=NanTolDistTri2Pnt( X,Y,Z,P1,P2,P3,MidPoint,FaceNormalVector,Distance )
%
%
% Arguments: (input)
%    X,Y,Z         - The query points that will be set to NAN if within
%                    the distance.
%
%  MidPoint        - The x, y and Z locations of each triangles midpoint
%                   [x,y,z]. 
%
% P1,P2,P3          - n*3 Column vector where each 'P' represents the
%                    different corner points of one of the triangles (XYZ).
%                    Not as efficient in terms of storage but easier to
%                    understand.
%
% FaceNormalVector  - The normal vectors n*3, (CosAx,CosAy,CosAz) of each
%                    triangle.
%
% Distance         - Set distance around each line segment that is removed.
%                   
%
%
%
% Arguments: (output)
%    X,Y,Z           - The query points that will have been set to NAN if within
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
    
    %Vector from midpoint to Tris P1
    Mid2P1=normr([(P1(i,1)-MidPoint(i,1)),...
                  (P1(i,2)-MidPoint(i,2)),...
                  (P1(i,3)-MidPoint(i,3))]);
    %Amount we move P1 out:          
    NewP1Offset=[Mid2P1(1).*Distance,Mid2P1(2).*Distance,Mid2P1(3).*Distance];
    
    P1(i,:)=P1(i,:)+NewP1Offset;
    
    %Vector from midpoint to Tris P2
    Mid2P2=normr([(P2(i,1)-MidPoint(i,1)),...
                  (P2(i,2)-MidPoint(i,2)),...
                  (P2(i,3)-MidPoint(i,3))]);
    %Amount we move P2 out:          
    NewP2Offset=[Mid2P2(1).*Distance,Mid2P2(2).*Distance,Mid2P2(3).*Distance];
    
    P2(i,:)=P2(i,:)+NewP2Offset;    
    
    PointsLst=[P1(i,:),P2(i,:),P3(i,:)];
    Seperation=[FaceNormalVector(i,1).*Distance,FaceNormalVector(i,2).*Distance,FaceNormalVector(i,3).*Distance];

        %Vector from midpoint to Tris P3
    Mid2P3=normr([(P3(i,1)-MidPoint(i,1)),...
                  (P3(i,2)-MidPoint(i,2)),...
                  (P3(i,3)-MidPoint(i,3))]);
    %Amount we move P3 out:          
    NewP3Offset=[Mid2P3(1).*Distance,Mid2P3(2).*Distance,Mid2P3(3).*Distance];
    
    P3(i,:)=P3(i,:)+NewP3Offset;
    
    %Start axis lims
    %The triangles after being seperated by chonsen amount
    [PointsNew2P]=DisplaceTris2(PointsLst,Seperation);
    [PointsNew2N]=DisplaceTris2(PointsLst,-Seperation);

    %Create an alpha shape around the points. 
    shp = alphaShape([PointsNew2P(:,1);PointsNew2N(:,1)],...
        [PointsNew2P(:,2);PointsNew2N(:,2)],...
        [PointsNew2P(:,3);PointsNew2N(:,3)],1E9);
     %Drawing if wanted: 
     %figure;
     %plot(shp);  hold on
     %scatter3(X(:),Y(:),Z(:),'.k')
    
    X=X(:);
    Y=Y(:);
    Z=Z(:);
    in=zeros(numel(X),1);
    for k=1:numel(X)
        if isnan(X(k))
        else
        %Now check if point is within the bound: 
        if inShape(shp,X(k),Y(k),Z(k)) 
            in(k)=1;
        end
        end
    end
    %Now check if point is within radius of 1st corner
    Xr1=X-P1(i,1);
    Yr1=Y-P1(i,2);
    Zr1=Z-P1(i,3);    
    [~,~,r1] = cart2sph(Xr1,Yr1,Zr1); 
    %Now check if point is within radius of 2nd corner
    Xr2=X-P2(i,1);
    Yr2=Y-P2(i,2);
    Zr2=Z-P2(i,3);    
    [~,~,r2] = cart2sph(Xr2,Yr2,Zr2);     
    %Now check if point is within radius of 3rd corner
    Xr3=X-P3(i,1);
    Yr3=Y-P3(i,2);
    Zr3=Z-P3(i,3);    
    [~,~,r3] = cart2sph(Xr3,Yr3,Zr3);     

    %Nan Points inside shp
    X(logical(in))=nan;
    Y(logical(in))=nan;
    Z(logical(in))=nan;
    %Nan Points within distance of edge points
    X(r1<Distance)=nan;
    Y(r1<Distance)=nan;
    Z(r1<Distance)=nan; 
    X(r2<Distance)=nan;
    Y(r2<Distance)=nan;
    Z(r2<Distance)=nan;
    X(r3<Distance)=nan;
    Y(r3<Distance)=nan;
    Z(r3<Distance)=nan;          
    
end

if sum(isnan(X(:)))==numel(X) 
    disp('Function DistTri2Pnt has set all points to NaN, the distance tolerance used it too high.')
elseif sum(isnan(X(:)))>0
    disp('Function DistTri2Pnt has set some points to NaN')    
end


end


%Internal func
function [PointsNewLoc]=DisplaceTris2(Points,Movement)
    %Takes inputs and moves each triangles points by the amount specified in the 3*3 vector 'Movement'
    N=3;
    PointsNewMvX=bsxfun(@plus,Points(:,1:N:end),(Movement(:,1))); 
    PointsNewMvY=bsxfun(@plus,Points(:,2:N:end),(Movement(:,2))); 
    PointsNewMvZ=bsxfun(@plus,Points(:,3:N:end),(Movement(:,3))); 
    PointsNewLoc=[reshape(PointsNewMvX.',1,[])',reshape(PointsNewMvY.',1,[])',reshape(PointsNewMvZ.',1,[])'];
end
