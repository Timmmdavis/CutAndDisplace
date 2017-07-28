function [MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles,draw)
%Creates MidPoints for the fault triangles and the face normal vector that is used to calculate traction. 
%This flips all normals that point down on the imported surface upwards.
%This fixes TDE calculation issues.
%MidPoint is the midpoint of each triangle
%FaceNormalVector is the normal vector or 'direction cosines' ax ay az
%TR is the MATLAB triangulation class
%Triangles is a list where each row in the list of three rows in 'Points2' that
%is the XYZ of one of the points of the triangle. 

%   Copyright 2017, Tim Davis, The University of Aberdeen
if nargin==2
    draw=1;
end    


%Calculate mesh surface, midpoints and normals on the triangles. 
Points2=Points(:,2:4);

%Loop to create triangle midpoint
%Could use commented MATLAB func that relies on triangulation below
% % TR = triangulation(Triangles,Points2);
% % MidPointOld = incenter(TR);
% % FaceNormalVectorOld = faceNormal(TR);
%Prepping for loop
MidPoint=zeros(size(Triangles));
FaceNormalVector=zeros(size(Triangles));
for i = 1:numel(Triangles(:,1))
    
    %Calculating midpoint (the incenter)
    
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


%Adding warning for surfaces with flat triangles, when calculating shear or
%normal traction on these both are calculated where really only one should
%be. 
flag=FaceNormalVector(:,1)==0 & FaceNormalVector(:,2)==0; %only Az exists (ie normal points up)
if any(flag)==1
    disp('WARNING!!')
    disp('Your surface has flat triangles, this does not work correctly for friction... maybe other parts of the code.');
    disp('Test properly or remesh surface');
end    

if draw==1
% Drawing the calculated trianglulation
figure ('name','Step 2.1'); trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',(.8),'facecolor', 'cyan');
%trisurf(TR,'FaceColor', 'cyan', 'faceAlpha',0.8);
axis equal;
hold on;
quiver3(MidPoint(:,1),MidPoint(:,2),MidPoint(:,3), ...
     FaceNormalVector(:,1),FaceNormalVector(:,2),FaceNormalVector(:,3),0.5, 'color','r');
xlabel('x'); ylabel('y');
title('Loaded surface showing normals') 
hold off;
end

end
