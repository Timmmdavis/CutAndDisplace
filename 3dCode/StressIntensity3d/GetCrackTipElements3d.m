function [FeP1P2S,FeP1P3S,FeP2P3S]=GetCrackTipElements3d...
    (MidPoint,P1,P2,P3,FaceNormalVector)
% GetCrackTipElements3d: Creates some structures the size of P1P2P3 that 
%                   correspond to connections between these points. This
%                   contains flags saying if the two points are an edge.
%                   The length of this edge, the vector from the Triangles
%                   Midpoint to the edge and the vector along this edge.
%                   This function was originally created to get variables
%                   needed to approximate stress intensity factors at edge
%                   triangles. 
%               
% usage:
% [FeP1P2S,FeP1P3S,FeP2P3S]=GetCrackTipElements3d...
% (MidPoint,P1,P2,P3);
%
% Arguments: (input)
% MidPoint          - The [X,Y,Z] of each triangles MidPoint, same length
%                    as Triangles so index's correspond.
%
% P1,P2,P3          - The corner point of each triangle in 'Triangles'.
%                    Arranged so the row index's correspond exactly to
%                    triangles and MidPoint. 
%
% FaceNormalVector  - The direction cosines of the triangle normal 
%				     [CosAx,CosAy,CosAz]
%
% Arguments: (output)
% FePaPbS           - A structure array containing information between
%                    Point-a and Point-b. The arrays inside will be NaNs
%                    unless this triangle connector is an edge element.
%                    Rows correspond to index in Triangles/P1,P2,P3
%                    Inside these structures are:
%
%                       FeLe    - FreeEdgeLength (Length between Pa and Pb)
%                       FeMd    - FreeEdgeMidPoint (MidPoint XYZ of the 
%                                   free edge connection).     
%                       FeEv    - FreeEdgeEdgeVector (Dir cosines 
%                                   CosAx,CosAy,CosAz that point along the edge).
%                       FeM2Ev  - FreeEdgeMidPointToEdgeVector (Dir cosines
%                                   CosAx,CosAy,CosAz that point from the 
%                                   MidPoint of the traingle to the MidPoint
%                                   of the free edge connection).
%                       FreeFlg - FreeEdgeFlag (A flag that says if the
%                                   connection is a free edge).
%
%                       FeM2ELe - The length of the Tris 'midpoint' to the
%                                   free edge. 
%                       IntAng -  The angle in degrees of the triangles
%                                   corner facing the free edge.
%                             
% Example usage:
%
% [x,y] = meshgrid(-2:.2:2);                                
% z = x .* exp(-x.^2 - y.^2);
% Triangles = delaunay(x(:),y(:));
% Points=[[1:numel(x)]',x(:),y(:),z(:)];
% [MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
% [P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 
% hold on
% [FeP1P2S,FeP1P3S,FeP2P3S]=GetCrackTipElements3d...
% (MidPoint,P1,P2,P3,FaceNormalVector);
% %As an example drawing the vector going from the midpoint of the edge 
% % triangles to the edges MidPoint:
% FreeEdMdX=[FeP1P2S.FeMd(FeP1P2S.FreeFlg,1);FeP1P3S.FeMd(FeP1P3S.FreeFlg,1);FeP2P3S.FeMd(FeP2P3S.FreeFlg,1)];
% FreeEdMdY=[FeP1P2S.FeMd(FeP1P2S.FreeFlg,2);FeP1P3S.FeMd(FeP1P3S.FreeFlg,2);FeP2P3S.FeMd(FeP2P3S.FreeFlg,2)];
% FreeEdMdZ=[FeP1P2S.FeMd(FeP1P2S.FreeFlg,3);FeP1P3S.FeMd(FeP1P3S.FreeFlg,3);FeP2P3S.FeMd(FeP2P3S.FreeFlg,3)];
% FreeEdM3EV=[FeP1P2S.FeM2Ev(FeP1P2S.FreeFlg,:);FeP1P3S.FeM2Ev(FeP1P3S.FreeFlg,:);FeP2P3S.FeM2Ev(FeP2P3S.FreeFlg,:)];
% quiver3(FreeEdMdX,FreeEdMdY,FreeEdMdZ,FreeEdM3EV(:,1),FreeEdM3EV(:,2),FreeEdM3EV(:,3),'b')
% FreeEd=sum([FeP1P2S.FreeFlg,FeP1P3S.FreeFlg,FeP2P3S.FreeFlg],2)>0; %Logical of free edges
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University


%Get edge triangles
[P1P2FreeFlg,P2P3FreeFlg,P1P3FreeFlg]=EdgeCons(P1,P2,P3,MidPoint);

%Do for P1 P2: (Function at base of file)
[Flg,FeLe,FeMd,FeEv,FeM2Ev,FeM2ELe,IntAng]=GetValues(P1P2FreeFlg,P1,P2,P3,MidPoint,FaceNormalVector);
%Put in structure
FeP1P2S.FeLe=FeLe;
FeP1P2S.FeMd=FeMd;
FeP1P2S.FeEv=FeEv;
FeP1P2S.FeM2Ev=FeM2Ev;
FeP1P2S.FreeFlg=Flg;
FeP1P2S.FeM2ELe=FeM2ELe;
FeP1P2S.IntAng=IntAng;

%Do for P1 P3: (Function at base of file)
[Flg,FeLe,FeMd,FeEv,FeM2Ev,FeM2ELe,IntAng]=GetValues(P1P3FreeFlg,P3,P1,P2,MidPoint,FaceNormalVector);
%Put in structure
FeP1P3S.FeLe=FeLe;
FeP1P3S.FeMd=FeMd;
FeP1P3S.FeEv=FeEv;
FeP1P3S.FeM2Ev=FeM2Ev;
FeP1P3S.FreeFlg=Flg;
FeP1P3S.FeM2ELe=FeM2ELe;
FeP1P3S.IntAng=IntAng;

%Do for P2 P3: (Function at base of file)
[Flg,FeLe,FeMd,FeEv,FeM2Ev,FeM2ELe,IntAng]=GetValues(P2P3FreeFlg,P2,P3,P1,MidPoint,FaceNormalVector);
%Put in structure
FeP2P3S.FeLe=FeLe;
FeP2P3S.FeMd=FeMd;
FeP2P3S.FeEv=FeEv;
FeP2P3S.FeM2Ev=FeM2Ev;
FeP2P3S.FreeFlg=Flg;
FeP2P3S.FeM2ELe=FeM2ELe;
FeP2P3S.IntAng=IntAng;

function [I,FeLe,FeMd,FeEv,FeM2Ev,FeM2ELe,IntAng]=GetValues(Flg,Pa,Pb,Pc,MidPoint,FaceNormalVector)
%Fills array with values if the connection is a free edge. 

%Initialise some arrays:

%Lengths of Free edges
FeLe=nan(numel(MidPoint)/3,1);
%Length from InnerPoint to edge Midpoint
FeM2ELe=nan(numel(MidPoint)/3,1);
%MidPoints of Free edges
FeMd=nan(numel(MidPoint)/3,3);
%Vector parallel to Free edges (edge vector Ev)
FeEv=nan(numel(MidPoint)/3,3);
%Vector pointing from midpoint to midpoint of Free edge (mid to edge vector
%M2Ev)
FeM2Ev=nan(numel(MidPoint)/3,3);
%Internal angles
IntAng=nan(numel(MidPoint)/3,1);

%Create index thats says if the connection is a free edge. 
I=logical(Flg); 
%Length of edge
FeLe(I)=sqrt(((Pa(I,1)-Pb(I,1)).^2)+((Pa(I,2)-Pb(I,2)).^2)+((Pa(I,3)-Pb(I,3)).^2));
%Midpoint of edge
FeMd(I,:)=  [((Pa(I,1)+Pb(I,1))/2), ((Pa(I,2)+Pb(I,2))/2), ((Pa(I,3)+Pb(I,3))/2)];
%Length of mid to edge dist
FeM2ELe(I)=sqrt(((FeMd(I,1)-MidPoint(I,1)).^2)+((FeMd(I,2)-MidPoint(I,2)).^2)+((FeMd(I,3)-MidPoint(I,3)).^2));
%FePc2ELe(I)=sqrt(((FeMd(I,1)-Pc(I,1)).^2)+((FeMd(I,2)-Pc(I,2)).^2)+((FeMd(I,3)-Pc(I,3)).^2));
%Vector from midpoint to edge midpoint. 
FeM2Ev(I,:)=normr([FeMd(I,1)-MidPoint(I,1),FeMd(I,2)-MidPoint(I,2),FeMd(I,3)-MidPoint(I,3)]);


%Vector pointing along edge
FeEv(I,:)=normr([(Pa(I,1)-Pb(I,1)),  (Pa(I,2)-Pb(I,2)),     (Pa(I,3)-Pb(I,3))]);

%Internal angle of the triangle (angle between edges that are not the free
%edge in question). 
v=normr([(Pa(I,1)-Pc(I,1)),  (Pa(I,2)-Pc(I,2)),     (Pa(I,3)-Pc(I,3))]);
w=normr([(Pb(I,1)-Pc(I,1)),  (Pb(I,2)-Pc(I,2)),     (Pb(I,3)-Pc(I,3))]);

Indx=find(I);
for i = 1:numel(v(:,1))
    IntAng(Indx(i))=rad2deg(acos(dot(v(i,:),w(i,:))));
end
%|a| is the magnitude (length) of vector a
%|b| is the magnitude (length) of vector b
%? is the angle between a and b
%a · b = |a| × |b| × cos(?) 
%so acos the dot of the normalised gives theta

% figure;
% quiver3(0,0,0,v(1,1),v(1,2),v(1,3))
% hold on; quiver3(0,0,0,w(1,1),w(1,2),w(1,3))
% axis('equal')

%%
%Fix to make sure that the the Edge vector is counter clockwise from the
%mid2edge vector when looking in the normal direction:
for i=1:numel(Indx)
    %First rotate to flat:
    V1=[0,0,1]; %Pointing up
    %Get Nx Ny Nz for the vectors and put in here. 
    X=[FeM2Ev(Indx(i),1),FeEv(Indx(i),1)];
    Y=[FeM2Ev(Indx(i),2),FeEv(Indx(i),2)];
    Z=[FeM2Ev(Indx(i),3),FeEv(Indx(i),3)];
    %Rotate so vectors are flat:
    [X,Y,~] = RotateObject3dAllignVectors(FaceNormalVector(Indx(i),:),V1,X,Y,Z,0,0,0);
    %Now get the two vectors as 2D coords (2nd we rotate by 90 counter
    %clockwise)
    VM2Ev=[X(1),Y(1)]; 
    VEv=[Y(2),-X(2)]; 
    %Calculate the dot product
    AllignFlag=dot(VM2Ev,VEv);
    %See if these allign or not:
    if AllignFlag>0
        FeEv(Indx(i),:)=-FeEv(Indx(i),:);
    end
end

%Then the direction of this vector: 
%Vector from midpoint to edge midpoint (cross Pro)
FeM2Ev2=cross(FaceNormalVector(I,:),FeEv(I,:),2);
%Check if edge vector and slip vector match in sign
Vect=(dot(FeM2Ev(I,:)',FeM2Ev2'))<=0;
%Flip if not
FeM2Ev2(Vect==1)=-FeM2Ev2(Vect==1); 
FeM2Ev(I,:)=FeM2Ev2;




end

end