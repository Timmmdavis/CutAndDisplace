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
    %Option A = Loading Gocad ascii data.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  string='SphereUniformDistributionON2_500Faces.ts';
%  [ Points,Triangles ] = GoCadAsciiReader( string );
%  string='SphereFix4.ts';
%  [ PointsF,TrianglesF ] = GoCadAsciiReader( string );  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Loading Stl data.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   string='Sphere5120.stl';    
%   [ Points,Triangles ] = STLReader( string );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option C = User defined surface
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Density=10;
a=1; 
X = linspace(-a,a,Density);
Y = linspace(-a,a,Density); 
[X,Y] = meshgrid(X,Y); 
[ Triangles,Points ] = MeshSurfaceXYPnts( X,Y );
%Rotate this
[Points(:,2),Points(:,4)] = RotateObject2d(Points(:,2),Points(:,4),deg2rad(200));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Appending and flagging the fixed data points.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Appending the triangles and points to the places we are going to fix. 
%[Points,Triangles,Fdisp] = DataAppender3d( Points,PointsF,Triangles,TrianglesF);


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
mu=1;           	%Shear Mod, mu or G. Relates shear stress to shear strain. 
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
    %Option A = Constant Slip across the fault surface.
	%Runs a constant slip on every element. Output 'stresses' are blank
	%arrays of 0's as this option uses displacement boundary conditions.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
% sz=zeros(numel(P1(:,1)),1);
% Dss =1;      Dss    = sz+Dss;   %Positive = LeftLatMovement 
% Dds =0;      Dds    = sz+Dds;   %Positive = Reverse movement
% Dn  =0;      Dn     = sz+Dn;    %Positive = Opening movement
% clear sz 
% 
% Sxx = 0; 				
% Syy = 0; 
% Szz = 0; 
% Sxy = 0; 					
% Sxz = 0; 					
% Syz = 0;
% Option='A'; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. 
    %Choose to define in stress or strain.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% 	%%%%%%%%%%%%%%
%   %StrainInput       
% 	%%%%%%%%%%%%%%
% 
% %Put to 1 to define the stresses defined in 'stress input' as strain values    
% strain=0;                   
%     
% 	%%%%%%%%%%%%%% 
%   %StressInput
% 	%%%%%%%%%%%%%%
%     
% 
% Sxx = 0;                    
% Syy = 0; 
% Szz = 0.01;
% Sxy = 0;        			
% Sxz = 0;
% Syz = 0; 
% 
% Option='B'; 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option C = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. This option includes frictional contact
    %properties on the fault surface, elements cannot interpenetrate and
    %slip is reduced by the frictional parameters. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
Sxx = 0;       			
Syy = 0; 
Szz = -1;
Sxy = 0;        		
Sxz = 0;
Syz = 0; 
% %Finding distance of points along X from 0,0. 
X=(MidPoint(:,1));
ne=size(MidPoint(:,1));
Mu  = ones((ne))*0.1;    %Coefficient of friction
Sf  = ones(size(Mu))*0;  %Sliding friction 
PlotSlipDistribution3d(Triangles,Points,cmap2,Sf)
Option='C'; 


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option D = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. No opening components, this option solves
    %the boundary value problem through shearing of elements alone. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% 	%%%%%%%%%%%%%%
%     %StrainInput        
% 	%%%%%%%%%%%%%%
%     
% % %Put to 1 to define the stresses defined in 'stress input' as strain values    
% % strain=0;    
% 
% 
% 	%%%%%%%%%%%%%%
%     %StressInput
% 	%%%%%%%%%%%%%%
%      
% Sxx = 0;         			
% Syy = 0; 
% Szz = 1; 
% Sxy = 0;                    
% Sxz = 0;
% Syz = 0;
% 
% Option='D';

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Option E = Defining boundary conditions as traction at the element
%     %centres. 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
% Tss = 0;     %positive - right lateral
% Tds = 0.1;     %positive (normals point to sky) - drives slip down dip 
% Tn = 0;     %positive is pressure that loads the boundary in tension
% 
% Option='E';


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



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XYZ observation points to calculate stress strain
%and displacement on. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Regular/randomly spaced grid of points bounding the faults.
    %Define how far it extends away and its density.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%Define the number of cells on the square grid
cells=10;   
%How many extra cells you want to add away from the faults, 
%reduce to 0 if having half space issues 
padding=0;  
[maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ,sz]=...
    MinMaxDataExtents3d(Points(:,2:4),cells,padding);

%Uniform grid spacing 3D
[X,Y,Z] = meshgrid(mingriX:sz:maxgriX,mingriY:sz:maxgriY,mingriZ:sz:maxgriZ);    
[dimx,dimy,dimz] = size(X);

%Uniform grid spacing 2D
% [X,Y] = meshgrid(mingriX:sz:maxgriX,mingriY:sz:maxgriY);    
% [dimx,dimy] = size(X);  
% Z=zeros(size(X));

%%Random points within the calculated bounds
% NumPnts=1000; %number of points, do not have cells , 10000
% lengthx=maxgriX-mingriX;
% lengthy=maxgriY-mingriY;
% lengthz=maxgriZ-mingriZ;
% xmv=(maxgriX+mingriX)/2; 
% ymv=(maxgriY+mingriY)/2; 
% zmv=(maxgriZ+mingriZ)/2; 
% x=rand(1,NumPnts)*(lengthx*2);
% y=rand(1,NumPnts)*(lengthy*2);
% z=rand(1,NumPnts)*(lengthz*2);
% X=x-lengthx+xmv;
% Y=y-lengthy+ymv;
% Z=z-lengthz+zmv;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = User defined XYZ grid.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Creating grid with user defined sampling
% Density=10;
% X = linspace(679000,681000,Density);
% Y = linspace(1504000,1506000,Density); 
% Z = linspace(-2000,0,Density); 
% [X,Y,Z] = meshgrid(X,Y,Z); 
% [dimx,dimy,dimz] = size(X);

% %Creating 2D grid with user defined sampling
%Density=10;
%X = linspace(675000,685000,Density);
%Y = linspace(1498000,1508000,Density); 
%[X,Y] = meshgrid(X,Y);
%[dimx,dimy] = size(X);
%Z=zeros(size(X));

% %Observation points at midpoints
% X=MidPoint(:,1);
% Y=MidPoint(:,2);
% Z=MidPoint(:,3);


% %Creating a flat line
%  X = linspace(1,15,100); 
%  Y = zeros(size(X)); 
%  Z = zeros(size(X));

% %Random Data
% % Random data at on n*n square grid centred at xmv,ymv
% Density=5000;
% n=4;
% xmv=0; ymv=0;zmv=0;
% x=rand(1,Density)*n;
% y=rand(1,Density)*n;
% z=rand(1,Density)*n;
% X=x-n+xmv;
% Y=y-n+ymv;
% Z=z-n+zmv;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option C = Import XYZ from a .text file that is on the loaded paths. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% S = load('Data.dat');	
% X = S(:,1);                   
% Y = S(:,2);                         
% Z = S(:,3);  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option D = Load a secondary Fault Surface
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Secondary Surface
% SecondSurface=1; %Flag set
% string='SphereUniformDistributionON2_500Faces.ts';
% [ PointsObs,TrianglesObs ] = GoCadAsciiReader( string );
% PointsObs=[PointsObs(:,1),PointsObs(:,2),PointsObs(:,3),(PointsObs(:,4)-freesurface_height)];
% [MidPointObs,FaceNormalVectorObs] = MidPointCreate(PointsObs,TrianglesObs);
% X=MidPointObs(:,1);
% Y=MidPointObs(:,2);
% Z=MidPointObs(:,3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing the obs Point data just created
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
Dist=0.2;
[X,Y,Z]=NanTolDistTri2Pnt( X,Y,Z,P1,P2,P3,MidPoint,FaceNormalVector,Dist );
%Drawing figure of surface and the observation points
figure;
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',(.2),'FaceColor', [0.5 0 0.9 ]);
hold on
scatter3(X(:),Y(:),Z(:),5,'k')  %Showing surface and obs points
xlabel('x'); ylabel('z'); axis('equal'); title('Surface and Obs Points');
WhiteFigure;



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses/Disps at defined observation points in XYZ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
[StressTTotal,StrainTTotal,StressTChg,StrainTChg,StressTReg,StrainTReg]=...
CalculateStressOnSurroundingPoints3d(Dss,Dds,Dn,mu,lambda,X,Y,Z,Sxx,...
Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Calculate Displacements.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ux,Uy,Uz] = CalculateDisplacementOnSurroundingPoints3d...
(Dss,Dds,Dn,nu,X,Y,Z,P1,P2,P3,halfspace);

 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 6: Finite strain computation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    

% [exx,eyy,ezz,exy,exz,eyz,ExxInfErrorPerc,EyyInfErrorPerc,...
% EzzInfErrorPerc,ExyInfErrorPerc,ExzInfErrorPerc,EyzInfErrorPerc]...
% =FiniteStrainLagrangian3d(Ux,Uy,Uz,X,Y,Z);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: Visualisation and analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% figure,quiver3(X(:),Y(:),Z(:),Ux(:),Uy(:),Uz(:));
% xlabel('x'); ylabel('y'); axis('equal');
% title('Vectors showing displacement of points') ;

[Sxx,Syy,Szz,Sxy,Sxz,Syz ] = ExtractCols( StressTTotal );
[Exx,Eyy,Ezz,Exy,Exz,Eyz ] = ExtractCols( StrainTTotal );


[S1,S2,S3,S1dir,S2dir,S3dir] = EigCalc3d(Sxx(:),Syy(:),Szz(:),Sxy(:),Sxz(:),Syz(:));

%If the stresses don't draw well use this function
%FilterValue=2.5;
%[S1,S2,Sxx,Syy,Sxy] = NanOrRemoveBadPoints( FilterValue,0,S1,S2,Sxx,Syy,Sxy );

%Function that draws stress ellipsoids
DrawStressEllipsoidsPrincipal(StressTTotal,StrainTTotal,X,Y,Z,'Triangles',Triangles,'Points',Points);

%Function that draws principal stress directions
DrawS1S2S3Directions(StressTChg,X,Y,Z,'Triangles',Triangles,'Points',Points);

%Check if data is on a uniform or non uniform grid
[ Uniform ] = UniformGridCheck3d( X,Y,Z );

%If data is random just just scatter. 
if Uniform==0
    DrawScatterPlots3d( X,Y,Z,cmap, S1,S2,S3 )
    %Drawing the stress on the 2nd surface if it exists
if SecondSurface==1
    Mu=ones(size(Sxx))*0.6; %Coeff Friction
    Cohesion=zeros(size(Sxx)); %Coeff Friction
    [CSS] = CalculateMaximumCoulombStress(MidPointObs,FaceNormalVectorObs,...
        Sxx,Syy,Szz,Sxy,Sxz,Syz,Mu,Cohesion,PointsObs,TrianglesObs,cmap2 );
end

%If data is uniform plot as gridded data.
else
    %Reshaping stresses to grid dimensions
    [X,Y,Z,S1,S2,S3,Sxx,Syy,Szz,Sxy,Sxz,Syz ]=ReshapeData3d...
        ( dimx,dimy,dimz,X,Y,Z,S1,S2,S3,Sxx,Syy,Szz,Sxy,Sxz,Syz );
    %Drawing isocontours of the stress
    IsoContoursPrincipalStressPercentiles...
        ( S1,S2,S3,X,Y,Z,'Triangles',Triangles,'Points',Points,'Alpha',0.8);
end


