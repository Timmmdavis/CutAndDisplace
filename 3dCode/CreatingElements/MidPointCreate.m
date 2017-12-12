function [MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles,draw)
% MidPointCreate: Creates Midpoints for the fault triangles and also the 
%               face normal vector. (not using MATLABs triangulation
%               functions). 
%               
% usage #1:
% [MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles)
%
% usage #2:  No drawing of the surface
% [MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles,0)
%
% Arguments: (input)
% Points            - Columns 2 3 and 4 are the XYZ locations of one the
%                    corner points of a triangle. Column 1 is the index.
%                    Not needed unless you want to draw a surface too. 
%
% Triangles         -  Triangles is a list where each row contains 3 index
%                    locations in "Points" which contains the XYZ location
%                    of each corner of the triangle.
%                    Not needed unless you want to draw a surface too. 
%
% Draw              -  Draw is a flag where 1 means the surface will be
%                     drawn as a figure (showing its normals).
%
%
% Arguments: (output)
% MidPoint          - The midpoints n*3, (XYZ) of each triangle. 
%
% FaceNormalVector  - The normal vectors n*3, (CosAx,CosAy,CosAz) of each
%                    triangle.
%
%
% Example usage:
%
% [x,y] = meshgrid(-2:.2:2);                                
% z = x .* exp(-x.^2 - y.^2);
% Triangles = delaunay(x(:),y(:));
% Points=[[1:numel(x)]',x(:),y(:),z(:)];
% [MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

if nargin==2
    draw=1;
end    

%Calculate mesh surface, midpoints and normals on the triangles. 
Points2=Points(:,2:4);

%Could use commented lines below, these relies on matlabs triangulation and
%additional toolboxes. 
% % TR = triangulation(Triangles,Points2);
% % MidPointOld = incenter(TR);
% % FaceNormalVectorOld = faceNormal(TR);

%Prepping for loop
MidPoint=zeros(size(Triangles));
FaceNormalVector=zeros(size(Triangles));
for i = 1:numel(Triangles(:,1))
    
    %Grabbing the current triangle points for the loop
    CurrentT1=(Triangles(i,1));
    CurrentT2=(Triangles(i,2));
    CurrentT3=(Triangles(i,3));
    %Getting XYZ list for each vertex on the triangle
    Pa=Points2(CurrentT1,:);
    Pb=Points2(CurrentT2,:);
    Pc=Points2(CurrentT3,:);
    %Distance formula between two points.   
    c=sqrt((Pa(1)-Pb(1))^2+(Pa(2)-Pb(2))^2+(Pa(3)-Pb(3))^2);
    a=sqrt((Pb(1)-Pc(1))^2+(Pb(2)-Pc(2))^2+(Pb(3)-Pc(3))^2);
    b=sqrt((Pa(1)-Pc(1))^2+(Pa(2)-Pc(2))^2+(Pa(3)-Pc(3))^2);
    if a==0 || b==0 || c==0
        error('Slither triangles or triangles with no area exist on your surface')
    end
    %Calculating midpoint using
    %http://mathworld.wolfram.com/Incenter.html
    %Gives the same result as MATLABS incenter calc
    MidPointx=((Pa(1)*a)+(Pb(1)*b)+(Pc(1)*c))/(a+b+c);
    MidPointy=((Pa(2)*a)+(Pb(2)*b)+(Pc(2)*c))/(a+b+c);
    MidPointz=((Pa(3)*a)+(Pb(3)*b)+(Pc(3)*c))/(a+b+c);
    MidPoint(i,:)=[MidPointx,MidPointy,MidPointz];

    %Now calculating the normal orientation
    %New vectors (see calculating normals online
    U=Pb-Pa;
    V=Pc-Pa;
    %Cross product of the vectors
    Nx = (U(2)*V(3)) - (U(3)*V(2));
    Ny = (U(3)*V(1)) - (U(1)*V(3));
    Nz = (U(1)*V(2)) - (U(2)*V(1));
    %Vector Magnitude
    aMag=sqrt((Nx * Nx) + (Ny * Ny) + (Nz * Nz));
    %Norm values
    Ax=Nx/aMag;
    Ay=Ny/aMag;
    Az=Nz/aMag;
    FaceNormalVector(i,:)=[Ax,Ay,Az];
    
end

if draw==1
    % Drawing the calculated trianglulation
    figure ('name','Step 2.1'); trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',(.8),'facecolor', 'cyan');
    axis equal;
    hold on;
    quiver3(MidPoint(:,1),MidPoint(:,2),MidPoint(:,3), ...
         FaceNormalVector(:,1),FaceNormalVector(:,2),FaceNormalVector(:,3),0.5, 'color','r');
    xlabel('x'); ylabel('y');
    title('Loaded surface showing normals') 
    hold off;
    WhiteFigure
end

end
