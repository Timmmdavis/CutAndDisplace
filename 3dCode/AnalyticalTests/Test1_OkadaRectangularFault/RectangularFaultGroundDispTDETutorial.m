%For Tutorial 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you do not need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

clear;close all

% %===== add file paths ==========================
pathstring = pwd;                                   %Get the address of the current working directory
if ispc; parts = strsplit(mfilename('fullpath'), '\');       %Getting the address of the script working directory and splitting into cell array
else; parts = strsplit(mfilename('fullpath'), '/');  end
[~,n] = find(~cellfun(@isempty,strfind(parts,'BEM_DDM_MATLAB'))); %finding the scripts root directory and its location (n)
if ispc; addpath(genpath(strjoin(parts(1,1:n),'\'))); %adding all folders to the path that are in this dir 
else; addpath(genpath(strjoin(parts(1,1:n),'/'))); end;%mac/linux
cd(pathstring)                                       %jumping back to the current dir

	%Script to run dislocations on a fault surface

	% This follows the same unit convention as Poly3d so ''assumes dimensionally consistent units for physical quantities with
	% dimensions of length or stress. Thus if you use kilometers for any quantity with
	% dimensions of length (e.g. a coordinate position), you must use kilometers for all
	% quantities with dimensions of length (e.g. displacement discontinuity).''
	% Andrew Lyell Thomas, Masters thesis, Stanford University. 


cmap = colormap_cpt('Ccool-warm'); %Loads a colourmap to be used for figures. This is Kenneth Morelands diverging colourmap.
cmap2 = colormap_cpt('Ccool-warm2');

%Parameters:

%Slip:
StrikeSlipDisp  =0;     
DipSlipDisp     =-1;    
TensileSlipDisp =0;  
%Fault surface Dip=
Dip=45;
%Fault surface Azimuth=
Azimuth=45;
%Fault Depth
Depth=2.5;
%Fault width
a=2;
%Fault height
b=2;

		
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: Import the fault surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Creating our own fault surface
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Density=2;    
X = linspace(-a,a,Density);
Y = linspace(-b,b,Density); 
[X,Y] = meshgrid(X,Y); % define Cartesian grid
Z=ones(size(X));
%Triangulating
TRI = delaunay(X,Y);
%Making points and triangles.
Triangles=TRI; clear TRI

%Rotate around dip
[Y,Z] = RotateObject2d(Y,Z,deg2rad(Dip));
%Rotate around az
[X,Y] = RotateObject2d(X,Y,deg2rad(Azimuth));


%Appending points
Points=[(1:1:numel(X))',X(:),Y(:),Z(:)-Depth];


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define Full/Halfspace and Elastic constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %Defining halfspace. 1 is on, 0 is off. The code will run much faster with this off as the calculations are simpler in a full space. 
	%Use this if the algorithm is throwing errors for vertex's above 0m.
    
halfspace =0; 

	%Defining the height of the freesurface relative to the value of 0 of
	%your imported surfaces. Used when deformation is not at
	%sealevel. See the output MaxZ height if your half space clashes with
	%your points and adjust until this outputs as a negative value. 
freesurface_height = 0;

    %Defining elastic constants, All three are required to run this whole script.  See equations below if you need to calculate these 
	%from other constants. 
	%Note nu is needed in the displacement calculation.  
	
mu=1;%0.4;%4000 	%Shear Mod, mu or G. Relates shear stress to shear strain. 
lambda=1;%4000;  %Lame's constant  0 perfectly compressible like cork, Infinate for incompressible material like rubber.
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;

	%Equations you can use to calculate elastic parameters that you need above. 
	
% E = 10000;                        %Young's Modulus
% K = E/3*(1-2*nu);                 %Bulk Modulus.    Equation 8.25 Pollard
% mu = E/(2*(1+nu));                %Shear Modulus.   Equation 8.26 Pollard  which is Geo convention 
% lambda= E*nu/((1+nu)*(1-2*nu));   %Lamé's  constant,  Equation 8.27 Pollard
% nu =lambda/(2*(mu+lamda);            %Poisson's ratio, Equation 8.28 Pollard


	%Function to extract the three vertex's for each triangle and 
	%triangles defined by a combination of vertex Id's. This is used later by the algorithms 
    %The max Z value is output as text in the command window giving an idea
    %of minimum free surface height. 
Points=[Points(:,1),Points(:,2),Points(:,3),(Points(:,4)-freesurface_height)];
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 
	%Fictitious disp flag, if set to one then this forces the triangles with this to 0 displacement. 
Fdisp= zeros(size(Triangles(:,1))); 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, Remote, Remote
%no opening, tractions on the elements or user imported slip value for each face. 

%Stress tensor outputs from this part of the toolbox i.e. Sxx Syy etc are simply the user defined remote stress
%components (boundary conditions). These are used later in the stress calculations. If strain is input these are converted to stress 
%using Hooke's law. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Constant Slip across the fault surface
	%Runs a constant slip on every triangle. A mixture of two slips will give an inclined slip.
	%Output 'stresses' are blank arrays of 0's as this option is not driven by boundary conditions. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
a=size(P1);
b=a(1,1);
c=zeros(b,1); 
    
StrikeSlipDisp     = c+StrikeSlipDisp;     %Positive = LeftLatMovement 
DipSlipDisp		   = c+DipSlipDisp;        %Positive = Reverse movement
TensileSlipDisp	   = c+TensileSlipDisp;    %Positive = Opening movement
clear a b c             	%a b and c are used to create an array of slips the size of the number of triangles on the fault surface.

Sxx = 0; 					
Syy = 0; 
Szz = 0; 
Sxy = 0; 					
Sxz = 0; 					
Syz = 0;   



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XYZ observation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Simply define flat observation plane of points. XYZ with defined step size 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
[X,Y,Z] = meshgrid(-5:0.5:10,-5:0.5:10,0);
rowcount = length(X(:,1)); 
colcount = length(X(1,:));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses  on dispersed XYZ
%observation points. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%TotalStress,Stress and StressChange are Xx12 arrays of strain tensors and stress tensors. XX YY ZZ XY XZ YZ
	%Strain is the first 6 coloumns and stress the last 6. Stress is the regional stress 
	%StressChange is the stress on the points from the event.
	%TotalStress is the driving stress and stress change from the event added. 
	
 [TotalStrainStress,StressStrainChange,StressStrainRemoteSxx,Syy,Szz,Sxy,Sxz,Syz]=CalculateStressOnSurroundingPoints(StrikeSlipDisp,DipSlipDisp,TensileSlipDisp,mu,lambda,X,Y,Z,Sxx,Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu);
[exx,eyy,ezz,exy,exz,eyz,Sxx,Syy,Szz,Sxy,Sxz,Syz ] = ExtractCols( TotalStrainStress );


     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Calculate Displacements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  [Ux,Uy,Uz] = CalculateDisplacementOnSurroundingPoints( StrikeSlipDisp,DipSlipDisp,TensileSlipDisp, nu, X,Y,Z, P1, P2, P3,halfspace);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: VISUALISATION AND ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
titlesz=25;
fntsz=21;

%Tri dislocations
[E1,E2,E1dir,E2dir]=EigCalc2d(exx,eyy,exy);

%Plotting the displacement induced in the TDE code and its fault
figure;quiver3(X(:),Y(:),Z(:),Ux,Uy,Uz)
xlabel('x'); ylabel('y'); axis('equal'); title('Ground displacement and fault')
hold on
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceColor', 'cyan', 'faceAlpha',0.8);
hold off
ChangeFontSizes(fntsz,titlesz);


%Plotting Tri Dis displacements vs the Okada displacements
DrawDeformedGrid2d( X,Y,Ux,Uy,5,cmap2,(E1+E2));hax1=gca;no1=gcf;title('Ground deformation');

colorbar;colormap(cmap2);
xlabel('x'); ylabel('y');gca;% title('Triangles');
ChangeFontSizes(fntsz,titlesz);


%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Option D = Load a secondary Fault Surface
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     %Secondary Surface
%     SecondSurface=1; %Flag set
%     string='GoCadExportHector.ts';
%     [ PointsObs,TrianglesObs ] = GoCadAsciiReader( string,pathstring );
%     scl=3.0000e-04;
%     PointsObs=[PointsObs(:,1),(PointsObs(:,2).*scl)-mean(PointsObs(:,2).*scl)-3,(PointsObs(:,3).*scl)-mean(PointsObs(:,3).*scl),(PointsObs(:,4).*scl)];
%     [MidPointObs,FaceNormalVectorObs] = MidPointCreate(PointsObs,TrianglesObs);
%     hold on
%     trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',(.2),'FaceColor', [0.5 0 0.9 ]);
%     hold off
%     X=MidPointObs(:,1);
%     Y=MidPointObs(:,2);
%     Z=MidPointObs(:,3);
%     
%     %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %STEP 5: Calculate Stresses and/or displacements on dispersed XYZ
% %observation points. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Option A = Calculate Stresses and Strains
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	%TotalStrainStress,StrainStressChange and StrainStressRemote are Xx12 arrays of strain tensors and stress tensors. XX YY ZZ XY XZ YZ
% 	%Strain is the first 6 coloumns and stress the last 6. StressStrainRemoteSxx is the regional stress 
% 	%StressStrainChange is the stress on the points from the event.
% 	%TotalStressStrain is the driving stress and stress change from the event added. 
% 	
%  [TotalStrainStress,StrainStressChange,StrainStressRemote,Sxx,Syy,Szz,Sxy,Sxz,Syz]=CalculateStressOnSurroundingPoints(StrikeSlipDisp,...
%  DipSlipDisp,TensileSlipDisp,mu,lambda,X,Y,Z,Sxx,Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu);
% 
% if SecondSurface==1
% Mu=ones(size(Sxx))*0.6; %Coeff Friction
% Cohesion=zeros(size(Sxx)); %Coeff Friction
% [CSS] = CalculateCoulombStressOnPlane( X,Y,Z,FaceNormalVectorObs,Sxx,Syy,Szz,Sxy,Sxz,Syz,Mu,Cohesion,PointsObs,TrianglesObs,cmap2 );
% hold on
% trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',(.2),'FaceColor', [0.5 0 0.9 ]);
% hold off
% end