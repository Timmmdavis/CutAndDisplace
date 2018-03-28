% This test checks the estimated stress intensity for two 3D problems,
% firstly a curved crack subject to biaxial stress (2D plane strain
% problem) and then a inclined penny shaped crack that is subject to
% uniaxial tension. This checks that both problems are within a tolerance
% (assuming the elastic parmeters have not changed).
%
% Problem 1 in: Shiah, Y.C. and Lin, Y.J., 2004. An infinitely large plate
% weakened by a circular-arc crack subjected to partially distributed
% loads. Journal of engineering mathematics, 49(1), pp.1-18.
%
% Problem 2 in: Nejati, M., Paluszny, A. and Zimmerman, R.W., 2015. A
% disk-shaped domain integral method for the computation of stress
% intensity factors using tetrahedral meshes. International Journal of
% Solids and Structures, 69, pp.230-251.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you do not need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

%Conventions:
%
% This follows the same unit convention as Poly3d so ''assumes
% dimensionally consistent units for physical quantities with
% dimensions of length or stress. Thus if you use kilometers for any
% quantity with dimensions of length (e.g. a coordinate position), you
% must use kilometers for all quantities with dimensions of length
% (e.g. displacement discontinuity).'' Andrew Lyell Thomas, Masters
% thesis, Stanford University.
%
% Engineering stress convention is used in this script for both tensors
% and principal stresses. S1 is the LEAST compressive stress. 
% See figure 6.13 & 6.14 in Pollard and Fletcher 2005. 

%Clearing old figures and arrays. 
clear;close all

%Loading a colourmap to be used for figures.
cmap = colormap_cpt('Ccool-warm');
cmap2 = colormap_cpt('Ccool-warm2');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: Import the fault surface.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option C = User defined surface
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r=3; %Circle rad
HalfLength=1; %Crack half length
Density=10;

%Yscl
Yscl=4;
X = linspace(-HalfLength,HalfLength,Density); 
numX=numel(X);
Y = linspace(-HalfLength*Yscl,HalfLength*Yscl,Density*Yscl);
numY=numel(Y);
[X,Y] = meshgrid(X,Y); 
spacing=HalfLength/(Density/2);
%Adding tris on eastern edge
X=[X,(ones(1,numY)*1+spacing)']; 
Y=[Y,((linspace(-HalfLength*Yscl,HalfLength*Yscl,Density*Yscl))+spacing/2)']; 
%Adding tris on western edge
X=[X,-(ones(1,numY)*1+spacing)']; 
Y=[Y,((linspace(-HalfLength*Yscl,HalfLength*Yscl,Density*Yscl))+spacing/2)']; 

Y(Y>Yscl)=Yscl;

[ Triangles,Points ] = MeshSurfaceXYPnts( X,Y );
%Curving crack
Points(:,4)=-(sqrt(r.^2-(Points(:,2)).^2));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Appending and flagging the fixed data points.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Checks and Creates Fdisp
chk=exist('Fdisp','var');
if chk==0 %exists in workspace
    Fdisp=[];
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define full/halfspace and elastic constants.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Defining halfspace. 1 is on, 0 is off. The code will run much faster with
% this off as the calculations are simpler in a full space.
halfspace = 0; 

% Defining the height of the freesurface relative to the value of 0 of
% your imported surfaces. 
freesurface_height = 0;

% Defining elastic constants
mu=5;           	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
[ ~,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu );

% Getting surface ready for script. 
Points=[Points(:,1),Points(:,2),Points(:,3),(Points(:,4)-freesurface_height)];
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, remote (with
%friction), remote with no opening components, tractions on the elements
%(internal pressure) etc.
%
%If Cartesian stress tensor components are not known but principal stresses
%are use the function 'TensorsFromPrincipal' to calculate the tensors in a
%Cartesian reference frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[Sxx,Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,Mu,Sf,strain,SecondSurface ] = CreateBlankVars;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. 
    %Choose to define in stress or strain.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	%%%%%%%%%%%%%%
  %StrainInput       
	%%%%%%%%%%%%%%

%Put to 1 to define the stresses defined in 'stress input' as strain values    
strain=0;                   
    
	%%%%%%%%%%%%%% 
  %StressInput
	%%%%%%%%%%%%%%
    

Sxx = 0;                    
Szz = 1;
Syy = (Sxx+Szz)*nu; 
Sxy = 0;        			
Sxz = 0;
Syz = 0; 

Option='B'; 




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(Option,'A')    
    [Dss,Dds,Dn,Sxx,Syy,Szz,Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,...
     Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
     FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points);
end


    % Removing any fixed elements
if any(Fdisp)==1
    [Triangles,FaceNormalVector,MidPoint,P1,P2,P3]...
    = RemovingFixedEls3d(Triangles,FaceNormalVector,MidPoint,P1,P2,P3,Fdisp);
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing figures of slip distribution on the crack.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,Dss,Dds,Dn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Animate Fault movement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To do this look at file Animate3d.m

%As its a 2D problem we remove values out of range:
bad=abs(MidPoint(:,2))>1;
Dss(bad)=nan;
Dds(bad)=nan;
Dn(bad)=nan;


%% The analytical test: 

[FeP1P2S,FeP1P3S,FeP2P3S]=GetCrackTipElements3d(MidPoint,P1,P2,P3,FaceNormalVector);
[FeP1P2S,FeP1P3S,FeP2P3S]=StressIntensity3d(Dn,Dss,Dds,mu,nu,FaceNormalVector,FeP1P2S,FeP1P3S,FeP2P3S);

%Location of free edge midpoints.  
FreeEdMdY=[FeP1P2S.FeMd(FeP1P2S.FreeFlg,2);FeP1P3S.FeMd(FeP1P3S.FreeFlg,2);FeP2P3S.FeMd(FeP2P3S.FreeFlg,2)];
%Accumulate from structure into big vectors:
K1=[FeP1P2S.K1(FeP1P2S.FreeFlg);FeP1P3S.K1(FeP1P3S.FreeFlg);FeP2P3S.K1(FeP2P3S.FreeFlg)];
K2=[FeP1P2S.K2(FeP1P2S.FreeFlg);FeP1P3S.K2(FeP1P3S.FreeFlg);FeP2P3S.K2(FeP2P3S.FreeFlg)];

%Analytical formula for stress intensity at curved crack tips
CrkHeight=max(Points(:,4))-min(Points(:,4));
[K1an,K2an] = ShiahLin2004_StrIntCurvedCrack(Szz(1),CrkHeight,HalfLength);


%Data is the perimeter
figure;
scatter(FreeEdMdY,ones(size(K1))*K1an)
hold on
scatter(FreeEdMdY,ones(size(K1))*K2an)
scatter(FreeEdMdY,K1,15,'k','filled')
scatter(FreeEdMdY,K2,15,'k','filled')
title('Stress intensity for curved crack subjected to biaxial tension (in xz-plane)')
legend('Analytical solution K1','Analytical solution K2','Numerical result K1','Numerical result K2')
xlabel('Length along crack in 3d (y), num values past 1 filtered')
ylabel('K1/K2 at crack tips')
WhiteFigure;
ylim([0 max(K1an)+max(K1an)/4])
hold off

Residuala=ones(size(K1))*K1an-K1;
Residualb=ones(size(K2))*K2an-K2;
Residual=[Residuala;Residualb];
Residual(isnan(Residual))=[]; %remove nans
if any(abs(Residual)>0.3)
    error('Numerical value of K1 deviates far from analytical solution for curved line crack')
else
    disp('Good match to analytical solution')
end

disp('Starting Part 2')
		
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: Import the fault surface.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  

%  string='CircleMesh_1a_500Faces.ts';
%  [ Points,Triangles ] = GoCadAsciiReader( string );
% Radius=1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option C = User defined surface
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Density=10; (defined above)
InnerDist=(1/(Density));
OuterEdgeSmpl=((2*pi)/InnerDist);

% OuterEdgeSmpl=round(rad2deg((2*pi)/(Density/2)));
OuterEdgeSmpl=round(OuterEdgeSmpl/2)*2; %Making sure its divible by two. 
Radius=1; 
X = linspace(-Radius,Radius,Density);
Y = linspace(-Radius,Radius,Density); 
[X,Y] = meshgrid(X,Y); 
%Inner circle of points
%[Th,R]=cart2pol(X,Y); X(R>0.8)=[];Y(R>0.8)=[]; %Circle
%EdgeSkew
Rx=1;%0.8;
Ry=1;%1.05;
X=X(:);Y=Y(:);
[ External ] = InsideEllipse(X, Y, Rx, Ry, 0.01);
X(find(External))=[];Y(find(External))=[];
%Outer Edge isosceles length h, compared to one for a eq tri. 
h=1; 

[ xe,ye ] = CreateCircleXY( OuterEdgeSmpl,1,0 );
xe=xe.*Rx;
ye=ye.*Ry;
X=[X;xe]; Y=[Y;ye];
%Outer circle of points
[ xe,ye ] = CreateCircleXY( OuterEdgeSmpl/2,1+((2*h)*InnerDist),0 ); 
%Rotating so we get equilateral tris. 
Theta=deg2rad(360/OuterEdgeSmpl)/4;
[Xrot,Yrot] = RotateObject2d(X,Y,Theta);
X=[X;xe]; Y=[Y;ye];
%Scaling so radius is one
X=X/(1+((2*h)*InnerDist));
Y=Y/(1+((2*h)*InnerDist));

[ Triangles,Points ] = MeshSurfaceXYPnts( X,Y );

%%
 string='CircleMesh_1a_500Faces.ts';
 [ Points,Triangles ] = GoCadAsciiReader( string );
%%

%Pennys angle away from Z. 
Beta=25; 
%Rotate this (XZ)
[Points(:,2),Points(:,4)] = RotateObject2d(Points(:,2),Points(:,4),deg2rad(90-Beta));



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Appending and flagging the fixed data points.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Checks and Creates Fdisp
chk=exist('Fdisp','var');
if chk==0 %exists in workspace
    Fdisp=[];
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define full/halfspace and elastic constants.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Defining halfspace. 1 is on, 0 is off. The code will run much faster with
% this off as the calculations are simpler in a full space.
halfspace = 0; 

% Defining the height of the freesurface relative to the value of 0 of
% your imported surfaces. 
freesurface_height = 0;

% Defining elastic constants
mu=5;           	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
[ ~,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu );

% Getting surface ready for script. 
Points=[Points(:,1),Points(:,2),Points(:,3),(Points(:,4)-freesurface_height)];
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, remote (with
%friction), remote with no opening components, tractions on the elements
%(internal pressure) etc.
%
%If Cartesian stress tensor components are not known but principal stresses
%are use the function 'TensorsFromPrincipal' to calculate the tensors in a
%Cartesian reference frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[Sxx,Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,Mu,Sf,strain,SecondSurface ] = CreateBlankVars;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. 
    %Choose to define in stress or strain.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	%%%%%%%%%%%%%%
  %StrainInput       
	%%%%%%%%%%%%%%

%Put to 1 to define the stresses defined in 'stress input' as strain values    
strain=0;                   
    
	%%%%%%%%%%%%%% 
  %StressInput
	%%%%%%%%%%%%%%
    

Sxx = 0;                    
Syy = 0; 
Szz = 1;
Sxy = 0;        			
Sxz = 0;
Syz = 0; 

Option='B'; 




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(Option,'A')    
    [Dss,Dds,Dn,Sxx,Syy,Szz,Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,...
     Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
     FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points);
end


    % Removing any fixed elements
if any(Fdisp)==1
    [Triangles,FaceNormalVector,MidPoint,P1,P2,P3]...
    = RemovingFixedEls3d(Triangles,FaceNormalVector,MidPoint,P1,P2,P3,Fdisp);
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing figures of slip distribution on the crack.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,Dss,Dds,Dn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Animate Fault movement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To do this look at file Animate3d.m


%% The analytical test: 


[FeP1P2S,FeP1P3S,FeP2P3S]=GetCrackTipElements3d(MidPoint,P1,P2,P3,FaceNormalVector);
[FeP1P2S,FeP1P3S,FeP2P3S]=StressIntensity3d(Dn,Dss,Dds,mu,nu,FaceNormalVector,FeP1P2S,FeP1P3S,FeP2P3S);

%Calculate theta (location around the crack
FreeEdMdX=[FeP1P2S.FeMd(FeP1P2S.FreeFlg,1);FeP1P3S.FeMd(FeP1P3S.FreeFlg,1);FeP2P3S.FeMd(FeP2P3S.FreeFlg,1)];
FreeEdMdY=[FeP1P2S.FeMd(FeP1P2S.FreeFlg,2);FeP1P3S.FeMd(FeP1P3S.FreeFlg,2);FeP2P3S.FeMd(FeP2P3S.FreeFlg,2)];
FreeEdMdZ=[FeP1P2S.FeMd(FeP1P2S.FreeFlg,3);FeP1P3S.FeMd(FeP1P3S.FreeFlg,3);FeP2P3S.FeMd(FeP2P3S.FreeFlg,3)];
% Flatten back to calc theta
[FreeEdMdX,FreeEdMdZ] = RotateObject2d(FreeEdMdX,FreeEdMdZ,deg2rad(-(90-Beta)));
%scatter3(FreeEdMdX,FreeEdMdY,FreeEdMdZ); axis('equal'); 
Theta=rad2deg(cart2pol(FreeEdMdX,FreeEdMdY));
%Accumulate from structure into big vectors:
K1=[FeP1P2S.K1(FeP1P2S.FreeFlg);FeP1P3S.K1(FeP1P3S.FreeFlg);FeP2P3S.K1(FeP2P3S.FreeFlg)];
K2=[FeP1P2S.K2(FeP1P2S.FreeFlg);FeP1P3S.K2(FeP1P3S.FreeFlg);FeP2P3S.K2(FeP2P3S.FreeFlg)];
K3=[FeP1P2S.K3(FeP1P2S.FreeFlg);FeP1P3S.K3(FeP1P3S.FreeFlg);FeP2P3S.K3(FeP2P3S.FreeFlg)];

%%
FreeEdIndx=[find(FeP1P2S.FreeFlg);find(FeP1P3S.FreeFlg);find(FeP2P3S.FreeFlg)];
FreeEdMdX=[FeP1P2S.FeMd(FeP1P2S.FreeFlg,1);FeP1P3S.FeMd(FeP1P3S.FreeFlg,1);FeP2P3S.FeMd(FeP2P3S.FreeFlg,1)];
FreeEdMdY=[FeP1P2S.FeMd(FeP1P2S.FreeFlg,2);FeP1P3S.FeMd(FeP1P3S.FreeFlg,2);FeP2P3S.FeMd(FeP2P3S.FreeFlg,2)];
FreeEdMdZ=[FeP1P2S.FeMd(FeP1P2S.FreeFlg,3);FeP1P3S.FeMd(FeP1P3S.FreeFlg,3);FeP2P3S.FeMd(FeP2P3S.FreeFlg,3)];
FreeEdEV=[FeP1P2S.FeEv(FeP1P2S.FreeFlg,:);FeP1P3S.FeEv(FeP1P3S.FreeFlg,:);FeP2P3S.FeEv(FeP2P3S.FreeFlg,:)];
FreeEdM2EV=[FeP1P2S.FeM2Ev(FeP1P2S.FreeFlg,:);FeP1P3S.FeM2Ev(FeP1P3S.FreeFlg,:);FeP2P3S.FeM2Ev(FeP2P3S.FreeFlg,:)];
quiver3(FreeEdMdX,FreeEdMdY,FreeEdMdZ,FreeEdEV(:,1),FreeEdEV(:,2),FreeEdEV(:,3),'r'); hold on
quiver3(FreeEdMdX,FreeEdMdY,FreeEdMdZ,FreeEdM2EV(:,1),FreeEdM2EV(:,2),FreeEdM2EV(:,3),'b'); hold on
quiver3(FreeEdMdX,FreeEdMdY,FreeEdMdZ,FaceNormalVector(FreeEdIndx,1),FaceNormalVector(FreeEdIndx,2),FaceNormalVector(FreeEdIndx,3),'g'); hold on
%%


[K1an,K2an,K3an] = NejatiEtAl2015_StrIntInclinedPennyTension(Szz,Beta,Theta,Radius,nu);
K1an=(K1an);
K2an=(K2an);
K3an=(K3an);

%Data is the perimeter
figure;
scatter(Theta,ones(size(K1)).*K1an,[],[0.4,0.4,0.4])
hold on
scatter(Theta,ones(size(K1)).*K2an,'r')
scatter(Theta,ones(size(K1)).*K3an,'b')
scatter(Theta,K1,15,'k','filled')
scatter(Theta,K2,15,'r','filled')
scatter(Theta,K3,15,'b','filled')
title('Stress intensity for inclined penny crack subjected to tension')
legend('Analytical solution K1','Analytical solution K2','Analytical solution K3'...
    ,'Numerical result K1','Numerical result K2','Numerical result K3')
xlabel('Theta (angle around crack)')
ylabel('Stress intensity at crack tips (ignoring sign)')
WhiteFigure;
ylim([0 max(K1an(1))+max(K1an(1))/4])
hold off

ylim([-1 1])

Residuala=(ones(size(K1)).*K1an)-K1;
Residualb=K2an-K2;
Residualc=K3an-K3;
Residual=[Residuala;Residualb;Residualc];
Residual(isnan(Residual))=[]; %remove nans
if any(abs(Residual)>0.15)
    error('Numerical value of K deviates far from analytical solution for penny-shaped crack')
else
    disp('Good match to analytical solution')
end

%Max Residual
MxResidual=[MxResidual;[max(Residuala),max(Residualb),max(Residualc)]];
