% Test: That midpoints of surfaces are traction free at the end of each loop. 
% 
% Proof: This tests the implementation of the equation system is correct
% for both friction and non friction


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
 
Density=15;
a=1; 
X = linspace(-a,a,Density);
Y = linspace(-a,a,Density); 
[X,Y] = meshgrid(X,Y); % define Cartesian grid
[ Triangles,Points ] = MeshSurfaceXYPnts( X,Y );


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
mu=5000;           	%Shear Mod, mu or G. Relates shear stress to shear strain. 
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


	%%%%%%%%%%%%%%
    %StrainInput        If halfspace = 1, Z stress components are removed
	%%%%%%%%%%%%%%
     
strain=0;                   %Put to 1 to define the stresses defined in 'stress input' as strain values
    
	%%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%

Mu=0;
Sf=0;

Option='B'; %Satisfying Sxz
Sxx = 0;Syy = 0;Szz = 0;Sxy = 0;Sxz = 1;Syz = 0; 
TestBoundaryCondIsSatisfied(Points,Triangles,cmap2,MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option)
disp('Sxz is satisfied')

Option='C'; %Satisfying Sxz Fric solver
TestBoundaryCondIsSatisfied(Points,Triangles,cmap2,MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option)
disp('Sxz is satisfied Friction')

Option='B'; %Satisfying Szz
Sxx = 0;Syy = 0;Szz = 1;Sxy = 0;Sxz = 0;Syz = 0; 
TestBoundaryCondIsSatisfied(Points,Triangles,cmap2,MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option)
disp('Szz is satisfied')

Option='C'; %Satisfying Szz Fric solver
TestBoundaryCondIsSatisfied(Points,Triangles,cmap2,MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option)
disp('Szz is satisfied Friction')

%Rotating so normals face in xx
Rotatation=-90; 
%Rotating around XZ
[Points(:,2),Points(:,4)] = RotateObject2d(Points(:,2),Points(:,4),deg2rad(Rotatation) );
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 
Option='B'; %Satisfying Sxx
Sxx = 1;Syy = 0;Szz = 0;Sxy = 0;Sxz = 0;Syz = 0; 
TestBoundaryCondIsSatisfied(Points,Triangles,cmap2,MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option)
disp('Sxx is satisfied')

Option='C'; %Satisfying Sxx Fric solver
TestBoundaryCondIsSatisfied(Points,Triangles,cmap2,MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option)
disp('Sxx is satisfied Friction')

Option='B'; %Satisfying Sxy
Sxx = 0;Syy = 0;Szz = 0;Sxy = 1;Sxz = 0;Syz = 0; 
TestBoundaryCondIsSatisfied(Points,Triangles,cmap2,MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option)
disp('Sxy is satisfied')

Option='C'; %Satisfying Sxy Fric solver
TestBoundaryCondIsSatisfied(Points,Triangles,cmap2,MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option)
disp('Sxy is satisfied Friction')

Option='B'; %Satisfying Syz
Sxx = 0;Syy = 0;Szz = 0;Sxy = 0;Sxz = 0;Syz = 1; 
TestBoundaryCondIsSatisfied(Points,Triangles,cmap2,MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option)
disp('Syz is satisfied')

Option='C'; %Satisfying Syz Fric solver
TestBoundaryCondIsSatisfied(Points,Triangles,cmap2,MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option)
disp('Syz is satisfied Friction')

%Rotating so normals face in yy
Rotatation=-90; 
%Rotating around XZ
[Points(:,2),Points(:,3)] = RotateObject2d(Points(:,2),Points(:,3),deg2rad(Rotatation) );
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 
Option='B'; %Satisfying Syy
Sxx = 0;Syy = 1;Szz = 0;Sxy = 0;Sxz = 0;Syz = 0; 
TestBoundaryCondIsSatisfied(Points,Triangles,cmap2,MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option)
disp('Syy is satisfied')

Option='C'; %Satisfying Syy Fric solver
TestBoundaryCondIsSatisfied(Points,Triangles,cmap2,MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option)
disp('Syy is satisfied Friction')

disp('Surfaces are traction free for any orientation/tensor')

function TestBoundaryCondIsSatisfied(Points,Triangles,cmap2,MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option)


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
    

% %Observation points at midpoints
 X=MidPoint(:,1);
 Y=MidPoint(:,2);
 Z=MidPoint(:,3);

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

[Sxx,Syy,Szz,Sxy,Sxz,Syz ] = ExtractCols( StressTTotal );
 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: Visualisation and analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
[ Tn,Tds,Tss ] = CalculateNormalAndShearTractions3d( FaceNormalVector,Sxx,Syy,Szz,Sxy,Sxz,Syz );

if any(abs([Tn;Tds;Tss])>1e-8) %~ Good enough approximation of traction free
	error('The solution is not correctly satisfying constraints that the surface is traction free')
end	
SumOfResultantTnTdsTss=Tn+Tds+Tss;
PlotSlipDistribution3d(Triangles,Points,cmap2,SumOfResultantTnTdsTss)

end