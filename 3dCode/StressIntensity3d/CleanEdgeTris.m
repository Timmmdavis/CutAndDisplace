function [P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector]=CleanEdgeTris(MidPoint,P1,P2,P3,TR,FaceNormalVector)

%% Get edge triangles
disp('Grabbing edge tris') 
[P1P2FreeFlg,P2P3FreeFlg,P1P3FreeFlg]=EdgeCons(P1,P2,P3,MidPoint);

%% PART 1: Get edge triangles that share an inner point
disp('Now find those that share an inner point') 
FreeTris=(P1P2FreeFlg+P2P3FreeFlg+P1P3FreeFlg)>0;
[P1P2TET,P2P3TET,P1P3TET]=GetConnectedTrianglesOnEdge(P1(FreeTris,:),P2(FreeTris,:),P3(FreeTris,:),MidPoint(FreeTris,:));
%TET for 'triangle edge touching'
%indexs to the points connecting the two tris

% % Collate
FreeTrisIndx=find(FreeTris);
hold on

%Find the tris with this connected edge (Will be in order of connected ones): 
BadTris=[FreeTrisIndx(P1P2TET);FreeTrisIndx(P2P3TET);FreeTrisIndx(P1P3TET)];
%Find the points on edges (that are not part of the connected edge)
P3Locs=FreeTrisIndx(P1P2TET); scatter3(P3(P3Locs,1),P3(P3Locs,2),P3(P3Locs,3),'filled','blue')
P1Locs=FreeTrisIndx(P2P3TET); scatter3(P1(P1Locs,1),P1(P1Locs,2),P1(P1Locs,3),'filled','blue')
P2Locs=FreeTrisIndx(P1P3TET); scatter3(P2(P2Locs,1),P2(P2Locs,2),P2(P2Locs,3),'filled','blue')
%Finding the inner shared point
P3InBad=ismember(BadTris,find(P1P2FreeFlg)); %Each indx is where P1P2 is a free edge on the bad tri (i.e. P3 is the good one)
P2InBad=ismember(BadTris,find(P1P3FreeFlg));
P1InBad=ismember(BadTris,find(P2P3FreeFlg));
P3Locs2=BadTris(P3InBad); scatter3(P3(P3Locs2,1),P3(P3Locs2,2),P3(P3Locs2,3),'filled','red')
P2Locs2=BadTris(P2InBad); scatter3(P2(P2Locs2,1),P2(P2Locs2,2),P2(P2Locs2,3),'filled','red')
P1Locs2=BadTris(P1InBad); scatter3(P1(P1Locs2,1),P1(P1Locs2,2),P1(P1Locs2,3),'filled','red')

%Assuming just two connections
[~,SortedTriangles,~] = ConnectedTrianglesFinder(TR,MidPoint);
SortedTriangles=SortedTriangles(BadTris,1:4);
Logic=ismember(SortedTriangles(:,:),BadTris); %Logical, if 1&2 or 3&4 are flagged its two we are looking for
Grab=(Logic(:,1)+Logic(:,2))==2;
Connections=[SortedTriangles(Grab,1:2);SortedTriangles(~Grab,3:4)];
Connections= unique(sort(Connections,2), 'rows');

   
%And the edge points:
PointC=[];PointB=[];PointA=[];
for ii=1:numel(Connections(:,1))
    
    A=P3Locs(find(ismember(P3Locs,Connections(ii,1))));
    B=P2Locs(find(ismember(P2Locs,Connections(ii,1))));
    C=P1Locs(find(ismember(P1Locs,Connections(ii,1))));
    if ~isempty(A)
        PointA=[PointA;P3(A,:)];
    end
    if ~isempty(B)
        PointA=[PointA;P2(B,:)];    
    end
    if ~isempty(C)        
        PointA=[PointA;P1(C,:)];   
    end
    
    D=P3Locs(find(ismember(P3Locs,Connections(ii,2))));
    E=P2Locs(find(ismember(P2Locs,Connections(ii,2))));    
    F=P1Locs(find(ismember(P1Locs,Connections(ii,2))));
    if ~isempty(D)
        PointB=[PointB;P3(D,:)];
    end
    if ~isempty(E)
        PointB=[PointB;P2(E,:)];    
    end
    if ~isempty(F)
        PointB=[PointB;P1(F,:)];   
    end
    
    G=P3Locs2(find(ismember(P3Locs2,Connections(ii,1))));
    H=P2Locs2(find(ismember(P2Locs2,Connections(ii,1))));    
    I=P1Locs2(find(ismember(P1Locs2,Connections(ii,1))));
    if ~isempty(G)
        PointC=[PointC;P3(G,:)];
    end
    if ~isempty(H)
        PointC=[PointC;P2(H,:)];    
    end
    if ~isempty(I)
        PointC=[PointC;P1(I,:)];   
    end    
    
end

%% Make sure we are facing the correct way:
disp('Fix normal direction') 
for JJ=1:numel(Connections(:,1))
    AvgVect=(FaceNormalVector(Connections(JJ,1),:)+FaceNormalVector(Connections(JJ,2),:))/2;
    NewTriNormalVector = CalculateTriangleNormal( PointA(JJ,:),PointB(JJ,:),PointC(JJ,:) );
    C = dot(AvgVect,NewTriNormalVector);
    if C<0
        PcTmp=PointA(JJ,:);
        PointA(JJ,:)=PointC(JJ,:);
        PointC(JJ,:)=PcTmp;
    end
end

%Remove the bad tris
P1(BadTris,:)=[];
P2(BadTris,:)=[]; 
P3(BadTris,:)=[];
disp('Removing rows here, would be a desireable place to keep track of indexing')

%Add new tris
P1=[P1;PointA];
P2=[P2;PointB];
P3=[P3;PointC];

%%
disp('Collate') 

Points=[zeros(size(P1));zeros(size(P1));zeros(size(P1))];
Points(1:3:end)=P1;
Points(2:3:end)=P2;
Points(3:3:end)=P3;
Triangles=1:1:numel(Points)/3;
Triangles=reshape(Triangles,3,[])';

figure;
trisurf(Triangles,Points(:,1),Points(:,2),Points(:,3));
Points=[(1:1:numel(Points(:,1)))',Points];

%% PART 2: Now Rotate so always eq lat tris on edge
%[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles,0);


%% Get edge triangles (Of cleaned tri)
[P1P2FreeFlg,P2P3FreeFlg,P1P3FreeFlg]=EdgeCons(P1,P2,P3,MidPoint);

hold on
scatter3(MidPoint(P1P2FreeFlg,1),MidPoint(P1P2FreeFlg,2),MidPoint(P1P2FreeFlg,3),'filled','red')
scatter3(MidPoint(P2P3FreeFlg,1),MidPoint(P2P3FreeFlg,2),MidPoint(P2P3FreeFlg,3),'filled','green')
scatter3(MidPoint(P1P3FreeFlg,1),MidPoint(P1P3FreeFlg,2),MidPoint(P1P3FreeFlg,3),'filled','blue')
title('CleanedDupEdges - FreeEdgeIndx - Going into Equi eqs')

%% Get new point if all tris are now isos
%Do for P1 P2: (Function at base of file)
[P1,P2,P3]=MakeEqEdgeTris(P1P2FreeFlg,P1,P2,P3,MidPoint,FaceNormalVector);

%Do for P1 P3: (Function at base of file)
[P1,P3,P2]=MakeEqEdgeTris(P1P3FreeFlg,P1,P3,P2,MidPoint,FaceNormalVector);

%Do for P2 P3: (Function at base of file)
[P2,P3,P1]=MakeEqEdgeTris(P2P3FreeFlg,P2,P3,P1,MidPoint,FaceNormalVector);


Points=[zeros(size(P1));zeros(size(P1));zeros(size(P1))];
Points(1:3:end)=P1;
Points(2:3:end)=P2;
Points(3:3:end)=P3;
Triangles=1:1:numel(Points)/3;
Triangles=reshape(Triangles,3,[])';

figure;
trisurf(Triangles,Points(:,1),Points(:,2),Points(:,3));
axis('equal')

Points=[(1:1:numel(Points(:,1)))',Points];
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles,0);

function [Pa,Pb,Pc]=MakeEqEdgeTris(Flg,Pa,Pb,Pc,MidPoint,FaceNormalVector)
%Fills array with values if the connection is a free edge. 

%Initialise some arrays:

%Lengths of Free edges
FeLe=nan(numel(MidPoint)/3,1);
%Length from innerpoint to edge Midpoint
FeIn2ELe=nan(numel(MidPoint)/3,1);
%MidPoints of Free edges
FeMd=nan(numel(MidPoint)/3,3);
%Vector parallel to Free edges (edge vector Ev)
FeEv=nan(numel(MidPoint)/3,3);
%Vector pointing from midpoint to midpoint of Free edge (mid to edge vector
%M2Ev)
FeM2Ev=nan(numel(MidPoint)/3,3);
%Internal angles
IntAng=nan(numel(MidPoint)/3,1);
%Vector from inner point to edge midpoint
FeIn2Ev=nan(numel(MidPoint)/3,3);

%Create index thats says if the connection is a free edge. 
Flag=logical(Flg); 
%Midpoint of edge
FeMd(Flag,:)=  [((Pa(Flag,1)+Pb(Flag,1))/2), ((Pa(Flag,2)+Pb(Flag,2))/2), ((Pa(Flag,3)+Pb(Flag,3))/2)];
%Vector from midpoint to edge midpoint. 
FeM2Ev(Flag,:)=normr([FeMd(Flag,1)-MidPoint(Flag,1),FeMd(Flag,2)-MidPoint(Flag,2),FeMd(Flag,3)-MidPoint(Flag,3)]);
%Vector from inner point to edge midpoint
FeIn2Ev(Flag,:)=normr([FeMd(Flag,1)-Pc(Flag,1),FeMd(Flag,2)-Pc(Flag,2),FeMd(Flag,3)-Pc(Flag,3)]);

%Vector pointing along edge
FeEv(Flag,:)=normr([(Pa(Flag,1)-Pb(Flag,1)),  (Pa(Flag,2)-Pb(Flag,2)),     (Pa(Flag,3)-Pb(Flag,3))]);

%Internal angle of the triangle (angle between edges that are not the free
%edge in question). 
v=normr([(Pa(Flag,1)-Pc(Flag,1)),  (Pa(Flag,2)-Pc(Flag,2)),     (Pa(Flag,3)-Pc(Flag,3))]);
w=normr([(Pb(Flag,1)-Pc(Flag,1)),  (Pb(Flag,2)-Pc(Flag,2)),     (Pb(Flag,3)-Pc(Flag,3))]);

Indx=find(Flag);
for i = 1:numel(v(:,1))
    IntAng(Indx(i))=rad2deg(acos(dot(v(i,:),w(i,:))));
end
%|a| is the magnitude (length) of vector a
%|b| is the magnitude (length) of vector b
%? is the angle between a and b
%a · b = |a| × |b| × cos(?) 
%so acos the dot of the normalised gives theta


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

%% 
%If triangles are not completly equilateral then the mid2Ed vec calculated
%above is not perpendicular to the edge (meaning poor calculation of K2).
%This fixes this for non eq tris.


%First we recompute the mid2edvec length (perp):
%Angle between vectors: 
Ang=pi/2-(acos(dot(FeIn2Ev(Flag,:)',FeEv(Flag,:)')));
Upsidedown=(FaceNormalVector(Flag,3)<0)';
Ang(Upsidedown)=-Ang(Upsidedown);
% %Length of R
% FePc2ELe(Flag,:);
% %Centre of Rotation
% FeMd(Flag,:);
%Default axis
Vect=FaceNormalVector(Flag,:);
%Place Point to be rotated in correct pos
PcCoords=Pc(Flag,:)-FeMd(Flag,:);
% EXAMPLE:
%     Rotate point (1;2;3) around vector (4;5;6) by an angle of pi/2
%     P = [1;2;3];  % create the point
%     V = [4;5;6];  % create vector around which rotation is performed
%     Qrot = qGetRotQuaternion( pi/2, V );
%     P2 = qRotatePoint( P, Qrotate ); ];
%
Pc2=zeros(numel(Ang),3);
for i=1:numel(Ang)
    Qrot = qGetRotQuaternion( Ang(i), Vect(i,:)' );
    Pc2(i,:) = qRotatePoint( PcCoords(i,:)', Qrot ); 
end

%Move back to orig coords:
Pc2=Pc2+FeMd(Flag,:);


% %% Only if you want eq tris
% %Length of edge
% FeLe(Flag,:)=sqrt(((Pa(Flag,1)-Pb(Flag,1)).^2)+((Pa(Flag,2)-Pb(Flag,2)).^2)+((Pa(Flag,3)-Pb(Flag,3)).^2));
% %Length of innerpoint to edge innerpoint
% FeIn2ELe(Flag,:)=sqrt(((FeMd(Flag,1)-Pc(Flag,1)).^2)+((FeMd(Flag,2)-Pc(Flag,2)).^2)+((FeMd(Flag,3)-Pc(Flag,3)).^2));
% Shdbe=(FeLe*sqrt(3))/2;
% ToMove=Shdbe-FeIn2ELe;
% FeIn2Ev(Flag,:)=normr([FeMd(Flag,1)-Pc2(:,1),FeMd(Flag,2)-Pc2(:,2),FeMd(Flag,3)-Pc2(:,3)]);
% Pc2=Pc2+(-ToMove(Flag,:).*FeIn2Ev(Flag,:));
% %%


hold on
scatter3(Pc2(:,1),Pc2(:,2),Pc2(:,3),'g','filled')

%And clean up connected tris
for i=1:numel(Indx)
    InPa=ismember(Pa(:,:),Pc(Indx(i),:),'rows');
    InPb=ismember(Pb(:,:),Pc(Indx(i),:),'rows');
    InPc=ismember(Pc(:,:),Pc(Indx(i),:),'rows');
    InPaIndx=find(InPa);
    for j=1:numel(InPaIndx)
        Pa(InPaIndx(j),:)=Pc2(i,:);
    end
    
    InPbIndx=find(InPb);
    for j=1:numel(InPbIndx)      
        Pb(InPbIndx(j),:)=Pc2(i,:);
    end    
    
    InPcIndx=find(InPc);
    for j=1:numel(InPcIndx)    
        Pc(InPcIndx(j),:)=Pc2(i,:);
    end    
end

Pc(Flag,:)=Pc2;



end

end