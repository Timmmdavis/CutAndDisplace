%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you do not need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

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

   
% %      %MediumSampleSetup
  string='BurgmannFreeBndE1.ts'; %Loading inner annlus, free bnd elastic 1
 [ PointsFBE1,TrianglesFBE1 ] = GoCadAsciiReader( string );		
   string='BurgmannInterfaceE1E2.ts'; %Loading interface elastic 1
 [ PointsIFE1E2,TrianglesIFE1E2 ] = GoCadAsciiReader( string );		
   string='BurgmannInterfaceE2E1.ts'; %Loading interface elastic 2
 [ PointsIFE2E1,TrianglesIFE2E1 ] = GoCadAsciiReader( string );
    string='BurgmannFreeBndE2.ts'; %Loading interface elastic 2
 [ PointsFBE2,TrianglesFBE2 ] = GoCadAsciiReader( string );


%Adding interface E1E2, Points towards E2 
[Points,Triangles,BoundaryFlag] = DataAppender3d( PointsFBE1,PointsIFE1E2,TrianglesFBE1,TrianglesIFE1E2);

Cng=BoundaryFlag==1;
BoundaryFlag(Cng)=2; %1 represents fixed elements that do not exist in this problem

%Adding interface E2E1, Points towards E1
[Points,Triangles,BoundaryFlag] = DataAppender3d( Points,PointsIFE2E1,Triangles,TrianglesIFE2E1,BoundaryFlag,3 );

%Adding FreeBnd E2E1, Points towards E1
[Points,Triangles,BoundaryFlag] = DataAppender3d( Points,PointsFBE2,Triangles,TrianglesFBE2,BoundaryFlag,4 );


%Quick Low Down on BoundaryFlag
%0=free boundary E1
%1=fixed bits of the free boundary of E1
%2=E1-E2 interface, E1 elastic properties, normals point towards E2
%3=E2-E1 interface, E2 elastic properties, normals point towards E1
%4=If existed would be free boundary E2, 
%5=fixed bits of the free boundary of E2

 

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
mu=100;           	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
[ ~,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu );

% Getting surface ready for script. 
Points=[Points(:,1),Points(:,2),Points(:,3),(Points(:,4)-freesurface_height)];
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 

title('All Surfaces Loaded');

%Locked els
Fdisp=(BoundaryFlag == 1 | BoundaryFlag == 5);



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

[Sxx,Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,Mu,Sf,strain ] = CreateBlankVars;

	%%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%

%Elastic 2 properties
nu2=0.25;
mu2=1000; %0.5
lambda2=mu2*((2*nu2)/(1-2*nu2));

%Appending these to the back of E1
nu=[nu;nu2];
mu=[mu;mu2];
lambda=[lambda;lambda2];
  
Sxx = 0;       			
Syy = 0; 
Szz = 0;
Sxy = 1;        		
Sxz = 0;
Syz = 0; 
Option='F';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Dss,Dds,Dn,Sxx,Syy,Szz,Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,...
 Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
 FaceNormalVector,BoundaryFlag,strain,Mu,Sf,Option,Triangles,Points);

%Removing fixed disps from BoundaryFlag
BoundaryFlag(Fdisp)=[]; %removing

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


IFFlag=(BoundaryFlag == 2 | BoundaryFlag == 3 );

ShearDispNoint=Dss(~IFFlag);

PointsX=MidPoint(:,1);
PointsZ=MidPoint(:,3);
PointsXNoint=PointsX(~IFFlag);
PointsZNoint=PointsZ(~IFFlag);

MiddlePartNoInt=(PointsZNoint<100 & PointsZNoint>-100);

PointsXNointMid=PointsXNoint(MiddlePartNoInt);
ShearDispNointMid=ShearDispNoint(MiddlePartNoInt);


[Rsorted, SortIndex] = sort(PointsXNointMid);
ShearDispNointMidSort = ShearDispNointMid(SortIndex);

figure;plot(Rsorted,abs(ShearDispNointMidSort))
xlabel('x'); ylabel('y');  title('SlipDistribution'),


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XYZ observation points to calculate stress strain
%and displacement on. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = User defined XY grid.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Uniform grid spacing 2D
Density=50;
X = linspace(-100,100,Density);
Y = linspace(-650,-350,Density); 
[X,Y] = meshgrid(X,Y); 
Z=zeros(size(X));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing and fixing Obs Point data just created
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dist=5;
[X,Y,Z]=NanTolDistTri2Pnt( X,Y,Z,P1,P2,P3,MidPoint,FaceNormalVector,Dist );


E1Flag=X<0;
XE1=X(E1Flag);
YE1=Y(E1Flag);
ZE1=Z(E1Flag);

XE2=X(~E1Flag);
YE2=Y(~E1Flag);
ZE2=Z(~E1Flag);

%Octave can't handle transparent objects
figure;trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',(.2),'FaceColor', [0.5 0 0.9 ]);
hold on
title('fractures and obs points'), xlabel('x'), ylabel('y');axis equal
hold on
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
    scatter(XE1(:),XE1(:),'b');
    scatter(XE2(:),YE2(:),'r');
elseif isOctave==0
    scatter(XE1(:),YE1(:),'b','filled','MarkerFaceAlpha',1/8);
    scatter(XE2(:),YE2(:),'r','filled','MarkerFaceAlpha',1/8);
end
title('Observation Points, Blue E1 Red E2 and displacement of E1 els'), xlabel('x'), ylabel('y')
hold off         
        

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses/Disps at defined observation points in XYZ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    
%E1 Elastic 1
E1Bits =(BoundaryFlag == 0 | BoundaryFlag == 2);
DnE1=Dn(E1Bits );
DssE1=Dss(E1Bits );
DdsE1=Dds(E1Bits );
Part=1; %Which elastic we want to extract
[MidPointE1,P1E1,P2E1,P3E1,muE1,lambdaE1,FaceNormalVectorE1,nuE1,NUME1,FdispE1,FB_E1,IF_E1]...
= ExtractElasticParts3d( Part,BoundaryFlag,MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,Fdisp );
    
%Calculating the stress on the points that make up the annulus (elastic 1).
%0 input stress as the original analytical calc is internal pressure not
%boundary stress. its the same Tn anyway for this geometry
 [StressTTotalE1,StrainTTotalE1,StressTChgE1,StrainTChgE1,StressTRegE1,StrainTRegE1]...
     =CalculateStressOnSurroundingPoints3d(DssE1,DdsE1,DnE1,muE1,lambdaE1,XE1,YE1,ZE1,0,0,0,0,0,0,P1E1,P2E1,P3E1,halfspace,nuE1);
 
[UxE1,UyE1,UzE1] = CalculateDisplacementOnSurroundingPoints3d( DssE1,DdsE1,DnE1, nuE1,XE1,YE1,ZE1, P1E1,P2E1,P3E1,halfspace);
 
%Extract stresses
 SxxE1=StressTTotalE1(:,1)';SyyE1=StressTTotalE1(:,2)';SxyE1=StressTTotalE1(:,4)';
 sz=numel(Sxx);
 
 
 %E2 Elastic 2
E2Bits =(BoundaryFlag == 3 | BoundaryFlag == 4);
DnE2=Dn(E2Bits );
DssE2=Dss(E2Bits );
DdsE2=Dds(E2Bits );
Part=2; %Which elastic we want to extract
[MidPointE2,P1E2,P2E2,P3E2,muE2,lambdaE2,FaceNormalVectorE2,nuE2,NUME2,FdispE2,FB_E2,IF_E2]...
= ExtractElasticParts3d( Part,BoundaryFlag,MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,Fdisp );
    
%Calculating the stress on the points that make up the annulus (elastic 1).
%0 input stress as the original analytical calc is internal pressure not
%boundary stress. its the same Tn anyway for this geometry
 [StressTTotalE2,StrainTTotalE2,StressTChgE2,StrainTChgE2,StressTRegE2,StrainTRegE2]...
     =CalculateStressOnSurroundingPoints3d(DssE2,DdsE2,DnE2,muE2,lambdaE2,XE2,YE2,ZE2,0,0,0,0,0,0,P1E2,P2E2,P3E2,halfspace,nuE2);

 [UxE2,UyE2,UzE2] = CalculateDisplacementOnSurroundingPoints3d( DssE2,DdsE2,DnE2, nuE2,XE2,YE2,ZE2, P1E2,P2E2,P3E2,halfspace);
 
 %Extract stresses
SxxE2=StressTTotalE2(:,1)';SyyE2=StressTTotalE2(:,2)';SxyE2=StressTTotalE2(:,4)';


X=[XE1;XE2];
Y=[YE1;YE2];
Ux=[UxE1;UxE2];
Uy=[UyE1;UyE2];
%Draw data with quiver, doesn't matter if its gridded or not
figure,quiver(X(:),Y(:),Ux(:),Uy(:));
xlabel('x'); ylabel('y'); axis('equal'); title('Disp'); 


Sxx=[StressTChgE1(:,1);StressTChgE2(:,1)];
Syy=[StressTChgE1(:,2);StressTChgE2(:,2)];
Sxy=[StressTChgE1(:,3);StressTChgE2(:,3)];

%Calclating 2d EigenValues
[S1,S2,S1dir,S2dir]=EigCalc2d(Sxx,Syy,Sxy);

DrawScatterPlots2d( X,Y,cmap2,Sxy,Syy,Sxx,S1,S2 );
