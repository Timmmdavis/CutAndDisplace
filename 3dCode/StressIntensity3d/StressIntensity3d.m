function [FeP1P2S,FeP1P3S,FeP2P3S] = StressIntensity3d...
    (Dn,Dss,Dds,mu,nu,FaceNormalVector,FeP1P2S,FeP1P3S,FeP2P3S)
%StressIntensity3d: Returns the value of stress intensity at elements at
%               the edge of a mesh. Outputs are structures as in inputs but
%               these contain the stress intensity approximation for the
%               connection (if its a free edge). 
%
%               Do not use in a half-space. Cracks open differently in this case
%				and the function would need to be changed.
%
% usage:
% [FeP1P2S,FeP1P3S,FeP2P3S] = StressIntensity3d...
% (Dn,Dss,Dds,E,nu,FaceNormalVector,FeP1P2S,FeP1P3S,FeP2P3S)
%
% Arguments: (input)
%  Dss,Dds,Dn       - Vectors that describe how much the elements displace in the
%                     normal (Dn) and strike slip (Dss) and dipslip (Dds)
%                     directions on the elements.
%
%       mu          - Shear modulus.
%
%       nu          - The Poisson's ratio.
%
% FaceNormalVector  - The direction cosines, CosAx (Nx), CosAy and CosAz in 
%                     a list for each element. 
%
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
%                                   MidPoint of the triangle to the MidPoint
%                                   of the free edge connection).
%                       FreeFlg - FreeEdgeFlag (A flag that says if the
%                                   connection is a free edge).
%                       FeM2ELe - The length of the Tris 'midpoint' to the
%                                   free edge. 
%                       IntAng -  The angle in degrees of the triangles
%                                   corner facing the free edge.
%
% Arguments: (output)
%                    
% FePaPbS            - As in inputs but with additional arrays in the
%                      structure that are the stress intensities if the
%                      connection is a free edge.
%
% Example usage:
%
% [x,y] = meshgrid(-2:.2:2);                                
% z = x .* exp(-x.^2 - y.^2);
% Triangles = delaunay(x(:),y(:));
% Points=[[1:numel(x)]',x(:),y(:),z(:)];
% [MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
% [P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 
% [FeP1P2S,FeP1P3S,FeP2P3S]=GetCrackTipElements3d...
% (Triangles,Points,MidPoint,P1,P2,P3);
% %Just open crack:
% mu=5;           	
% nu = 0.25;     
% Dn=ones(size(P1(:,1))); 
% Dss=zeros(size(Dn));
% Dds=zeros(size(Dn));
% [FeP1P2S,FeP1P3S,FeP2P3S] = StressIntensity3d...
%    (Dn,Dss,Dds,mu,nu,FaceNormalVector,FeP1P2S,FeP1P3S,FeP2P3S);
% K1=sum([Dds,FeP1P2S.K1,FeP1P3S.K1,FeP2P3S.K1]','omitnan')';
% PlotSlipDistribution3d(Triangles,Points,[],K1);
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University

%Call internal function (base of file)
[K1_P1P2,K2_P1P2,K3_P1P2]=K1K2K3TriDislocation(FeP1P2S.FreeFlg,mu,nu,Dn,FeP1P2S.FeM2Ev,FeP1P2S.FeEv,FeP1P2S.FeM2ELe,FeP1P2S.IntAng,FaceNormalVector,Dss,Dds);
%Call internal function (base of file)
[K1_P1P3,K2_P1P3,K3_P1P3]=K1K2K3TriDislocation(FeP1P3S.FreeFlg,mu,nu,Dn,FeP1P3S.FeM2Ev,FeP1P3S.FeEv,FeP1P3S.FeM2ELe,FeP1P3S.IntAng,FaceNormalVector,Dss,Dds);
%Call internal function (base of file)
[K1_P2P3,K2_P2P3,K3_P2P3]=K1K2K3TriDislocation(FeP2P3S.FreeFlg,mu,nu,Dn,FeP2P3S.FeM2Ev,FeP2P3S.FeEv,FeP2P3S.FeM2ELe,FeP2P3S.IntAng,FaceNormalVector,Dss,Dds);

%Put results into the strucs:
%P1P2
FeP1P2S.K1=K1_P1P2;
FeP1P2S.K2=K2_P1P2;
FeP1P2S.K3=K3_P1P2;
%P1P3
FeP1P3S.K1=K1_P1P3;
FeP1P3S.K2=K2_P1P3;
FeP1P3S.K3=K3_P1P3;
%P2P3
FeP2P3S.K1=K1_P2P3;
FeP2P3S.K2=K2_P2P3;
FeP2P3S.K3=K3_P2P3;



function [K1,K2,K3]=K1K2K3TriDislocation(LocFlg,mu,nu,Dn,FeM2Ev,FeEv,FeM2ELe,IntAng,FaceNormalVector,Dss,Dds)
%Calculates K1 K2 and K3 on connections that are free edges using simple
%formulas based on the displacement discontinuity of the triangle the
%connection borders. 

%Calculate direction cosines for the displacements
[ StrikeSlipCosine,DipSlipCosine ] = CalculateDSandSSDirs( FaceNormalVector(LocFlg,:) );
%Cart vector components of the strike slip displacement
DssCart=bsxfun(@times,StrikeSlipCosine,Dss(LocFlg));
%Cart vector components of the dip slip displacement 
DdsCart=bsxfun(@times,DipSlipCosine,Dds(LocFlg));
%Total displacement vector lying on triangle plane
DPlaneCart=DssCart+DdsCart; 
    
%Here we rotate all the different cosines so they are 2D in the plane of
%the triangle (using the normal vector). 

%Init some vars:
Num=sum(LocFlg);
V1=zeros(Num,3); %Vector pointing up
angle=zeros(Num,1);
DssFlt=zeros(Num,2);
DdsFlt=DssFlt;
FeM2EvFlt=DssFlt;
FeEvFlt=DssFlt;
V1(:,3)=1;
%Find the index for the loop
Indx=find(LocFlg);

%Loop to get the directions of the vectors in the coordiantes of the
%triangle plane. 
for i=1:Num
    
    Fnv=FaceNormalVector(Indx(i),:);
    %Rotate the slip directions to flat
    [DssFlt(i,1),DssFlt(i,2),~] = RotateObject3dAllignVectors(Fnv,V1(i,:),StrikeSlipCosine(i,1),StrikeSlipCosine(i,2),StrikeSlipCosine(i,3),0,0,0);
    [DdsFlt(i,1),DdsFlt(i,2),~] = RotateObject3dAllignVectors(Fnv,V1(i,:),DipSlipCosine(i,1),DipSlipCosine(i,2),DipSlipCosine(i,3),0,0,0);
    %Rotate the stress intensity directions to flat    
    [FeM2EvFlt(i,1),FeM2EvFlt(i,2),~] = RotateObject3dAllignVectors(Fnv,V1(i,:),FeM2Ev(Indx(i),1),FeM2Ev(Indx(i),2),FeM2Ev(Indx(i),3),0,0,0);
    [FeEvFlt(i,1),FeEvFlt(i,2),~] = RotateObject3dAllignVectors(Fnv,V1(i,:),FeEv(Indx(i),1),FeEv(Indx(i),2),FeEv(Indx(i),3),0,0,0);    
    
    %Angle between the strike slip vector and mid to edge vector in the
    %plane of the triangle
    angle(i) =abs(acos(DssFlt(i,1))-acos(FeM2EvFlt(i,1)));
    
end

% DefineDirection
Dir1=(pi/2)+angle; 	%direction (az) of component 1 (degrees)
% Calculating direction cosines
CosAx=sin(Dir1); %Angle away from X axis for Dir 1
CosAy=cos(Dir1); %Angle away from Y axis for Dir 1
%Calculating new vectors (using 'traction' function)
[ DAlongEd,DMid2Ed ] = CalculateNormalShearTraction2d( Dds(Indx),Dss(Indx),CosAx,CosAy);
%Ignoring sign for time being
DAlongEd=abs(DAlongEd);
DMid2Ed=abs(DMid2Ed);

%Init some vars
K1=nan(numel(Dn),1);
K2=K1;
K3=K2;




%% Correct sign of shear components
% This means these match the drawings in Fig 9.30 of Pollard and Fletcher
% assuming our end element normal corresponds to the y-axis in this figure.

%Assuming out vectors for the edge triangles point in the correct
%directions (MidPoint2Edge for (KII) and counter-clockwise from this vector
%when looking along the normal direction (KIII)).
%We check if the slip vector also points in this direction
%or not. We adjust sign accordingly. 
%Check if mid2edge vector and slip vector match in sign
Vect=(dot(FeM2Ev(LocFlg,:)',DPlaneCart'))<=0;
%Flip if not
DMid2Ed(Vect==1)=-DMid2Ed(Vect==1);
%Check if edge vector and slip vector match in sign
Vect=(dot(FeEv(LocFlg,:)',DPlaneCart'))<=0;
%Flip if not
DAlongEd(Vect==1)=-DAlongEd(Vect==1); 

%% Now calculating stress intensities

%Constants:
h=FeM2ELe(LocFlg); %Currently whole tri length FePc2ELe(LocFlg)/2
disp('Danger, Changed above for testing')

%Approximate stress intensities. 
K1(LocFlg)=(mu*sqrt(pi)./(2*sqrt(h).*(1-nu))).*Dn(LocFlg);
K2(LocFlg)=(mu*sqrt(pi)./(2*sqrt(h).*(1-nu))).*DMid2Ed;
K3(LocFlg)=(mu*sqrt(pi)./(2*sqrt(h).*(1-nu))).*DAlongEd.*(1-nu);       

%% And correcting errors in formulas

%Correction factors based on errors in formula above when we reduce 'h'
%towards zero and compare to results of the analytical solution and the
%overestimation from the constant movement of elements in 3D. 
C_K1=1.834;
C_K2=1.834;
C_K3=1.834;
K1=K1./C_K1;
K2=K2./C_K2;
K3=K3./C_K3;


end

end