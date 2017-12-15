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
    
    
    %Putting triangle flat    
    V2=[0,0,1 ]; %Pointing up
    V1=FaceNormalVector(i,:); %Pointing East
    Xtri=[P1(i,1),P2(i,1),P3(i,1)];
    Ytri=[P1(i,2),P2(i,2),P3(i,2)];
    Ztri=[P1(i,3),P2(i,3),P3(i,3)];
    %Now flatten tri:
    [Xtri,Ytri,~] = RotateObject3dAllignVectors(V1,V2,Xtri,Ytri,Ztri,MidPoint(i,1),MidPoint(i,2),MidPoint(i,3));

    %Edge 1
    [X1e1,X2e1,Y1e1,Y2e1]=MoveLine(Xtri(1),Xtri(2),Ytri(1),Ytri(2),Distance);
    %Edge 2 
    [X1e2,X2e2,Y1e2,Y2e2]=MoveLine(Xtri(2),Xtri(3),Ytri(2),Ytri(3),Distance);
    %Edge 3
    [X1e3,X2e3,Y1e3,Y2e3]=MoveLine(Xtri(3),Xtri(1),Ytri(3),Ytri(1),Distance);

    Xtri=[X1e1,X2e1,X1e2,X2e2,X1e3,X2e3];
    Ytri=[Y1e1,Y2e1,Y1e2,Y2e2,Y1e3,Y2e3];
    Ztri=[1,   1,   1,   1,   1,   1];

    %Now rotate back to real coords 
    V1=[0,0,1]; %Pointing up
    V2=FaceNormalVector(i,:); %Pointing East
    %Moving each by distance in relation to tri norm vector
    [Xup,Yup,Zup] = RotateObject3dAllignVectors(V1,V2,Xtri,Ytri,Ztri*Distance,0,0,0);
    [Xdwn,Ydwn,Zdwn] = RotateObject3dAllignVectors(V1,V2,Xtri,Ytri,Ztri*-Distance,0,0,0);
     

    %Create an alpha shape around the points. 
    shp = alphaShape([Xup;Xdwn],[Yup;Ydwn],[Zup;Zdwn],1E9);
    %Drawing if wanted: 
    %plot(shp)    
    
    %Now check if point is within the bound: 
    in = inShape(shp,Xtri,Ytri,Ztri) ; 

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
    X(in)=nan;
    Y(in)=nan;
    Z(in)=nan;
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

function [X1,X2,Y1,Y2]=MoveLine(X1,X2,Y1,Y2,Dist)
    
    %Call Func
    [ ~,~,P1F,P2F,LineNormalVectorF ]...
    = MidPoint_Orientation( X1,X2,Y1,Y2 );
    %Get extra X and Y bits:
    ExX=(LineNormalVectorF(:,1)*Dist);
    ExY=(LineNormalVectorF(:,2)*Dist);

    %Move Points along normal
    %X    
    P1F(:,1)=P1F(:,1)+ExX;
    P2F(:,1)=P2F(:,1)+ExX;
    %Y
    P1F(:,2)=P1F(:,2)+ExY;
    P2F(:,2)=P2F(:,2)+ExY;
    
    %Outputs
    X1=P1F(:,1);
    X2=P2F(:,1);
    Y1=P1F(:,2);
    Y2=P2F(:,2);     

end

end


