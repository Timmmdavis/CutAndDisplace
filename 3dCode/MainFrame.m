%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you do not need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   
clear;close all

%   Copyright 2017, Tim Davis, The University of Aberdeen

% %===== add file paths ==========================
pathstring = pwd;                                   %Get the address of the current working directory
if ispc; parts = strsplit(mfilename('fullpath'), '\');       %Getting the address of the script working directory and splitting into cell array
else; parts = strsplit(mfilename('fullpath'), '/'); end
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



	%Conventions: 
		%See figure 6.14 in David Pollards book for tensors. 
		%Geological convention is that a confining pressure is a positive stress. Geological convention is used throughout this script. 
		%For a normal stress: positive if it produces compression in the material.
		%For a shear stress: negative if, when acting on the + face identified by the first index, it points
		%in the + direction identified by the second index. Example: Exy is - if on the +x face it points
		%in the +y direction.
		%For example a positive value in shear stress Sxy will cause RightLateral movement.
		
		%Strain uses the same convention as engineering stress convention where extension is a positve value. When entering strain values remember this. 
		%For a normal stress: positive (negative) if it produces tension (compression) in the material.
		%For a shear stress: positive if, when acting on the + face identified by the first index, it points
		%in the + direction identified by the second index. Example: Exy is + if on the +x face it points
		%in the +y direction.
		%For example a positive value in shear strain Exy will cause Leftlateral movement.
		
		
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: Import the fault surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading formatted ascii data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    %  Loads the triangles that define the fault surface 
    %  Splits this into two files that define the vertex's and which vertex's
    %  make which triangle. Use Gocad ascii format and rejig like the
    %  example files. Note column 1 for the points must start at 1.
    %  Reexport Gocad ascii file to get this. 

% S = load('TwoTriangle.dat');          %Loads the dat file -	GocadAsciiexport formatted with no text
% Points = S([1:4],[1,2,3,4]);          %Points XYZ that are the fault. First numbers define which rows, second which columns.
% Triangles = S([5:6],[1,2,3]);         %which points make a triangle
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Loading Gocad ascii data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    string='GoCadExportHector.ts';
% % string='VolcanoTopo.ts';
% % string='SphereUniformDistributionON2_500Faces.ts';
 [ Points,Triangles ] = GoCadAsciiReader( string,pathstring );

   
%  string='SphereFix4.ts';
%  [ PointsF,TrianglesF ] = GoCadAsciiReader( string,pathstring );  


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option C = Loading Stl data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   string='Sphere5120.stl';    %WAS OTHER SURFACE BEFORE  2560
%   [ Points,Triangles ] = STLReader( string );

  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Appending and flagging the fixed data points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Appending the triangles and points to the places we are going to fix. 
%[Points,Triangles,Fdisp] = DataAppender3d( Points,PointsF,Triangles,TrianglesF);
	%Fictitious disp flag, if set to one then this forces those triangles to 0 displacement.
    %Checks and Creates Fdisp
chk=exist('Fdisp','var');
if chk==0 %exists in workspace
Fdisp=zeros(size(Points(:,1)));
end

% %% SCALING, REMOVE
Points(:,2)=Points(:,2).*1;
Points(:,3)=Points(:,3).*1;
Points(:,4)=Points(:,4).*1;
% %%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define Full/Halfspace and Elastic constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %Defining halfspace. 1 is on, 0 is off. The code will run much faster with this off as the calculations are simpler in a full space. 
	%Use this if the algorithm is throwing errors for vertex's above 0m.
    
halfspace = 0; 

	%Defining the height of the freesurface relative to the value of 0 of
	%your imported surfaces. Used when deformation is not at
	%sealevel. See the output MaxZ height if your half space clashes with
	%your points and adjust until this outputs as a negative value. 
freesurface_height =0;

    %Defining elastic constants, All three are required to run this whole script.  See equations below if you need to calculate these 
	%from other constants. 
	%Note nu is needed in the displacement calculation.  
	
mu=5000;%0.4;%4000 	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
lambda=mu*(2*nu/(1-2*nu)); %Lame's constant  0 perfectly compressible like cork, Infinite for incompressible material like rubber.
E = mu*(2*(1+nu)) ; %Young's Mod
	%Equations you can use to calculate elastic parameters that you need above. 
	
% E = 10000;                        %Young's Modulus
% K = E/3*(1-2*nu);                 %Bulk Modulus.    Equation 8.25 Pollard
% mu = E/(2*(1+nu));                %Shear Modulus.   Equation 8.26 Pollard  which is Geo convention 
% lambda= E*nu/((1+nu)*(1-2*nu));   %Lamé's  constant,  Equation 8.27 Pollard
% nu =lambda/(2*(mu+lambda));            %Poisson's ratio, Equation 8.28 Pollard


	%Function to extract the three vertex's for each triangle and 
	%triangles defined by a combination of vertex Id's. This is used later by the algorithms 
    %The max Z value is output as text in the command window giving an idea
    %of minimum free surface height. 
Points=[Points(:,1),Points(:,2),Points(:,3),(Points(:,4)-freesurface_height)];
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, Remote, Remote
%no opening, tractions on the elements or user imported slip value for each face. 

%Stress tensor outputs from this part of the toolbox i.e. Sxx Syy etc are simply the user defined remote stress
%components (boundary conditions). These are used later in the stress calculations. If strain is input these are converted to stress 
%using Hooke's law. 
%If Cartesian stress tensor components are not known but principal stresses
%are use the function 'TensorsFromPrincipal' to calculate the tensors in a
%Cartesian reference frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[Sxx,Syy,Szz,Sxy,Sxz,Syz,Tnn,Tss,Tds,Mu,Sf,strain,SecondSurface ] = CreateBlankVars;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Constant Slip across the fault surface
	%Runs a constant slip on every triangle. A mixture of two slips will give an inclined slip.
	%Output 'stresses' are blank arrays of 0's as this option is not driven by boundary conditions. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
% a=size(P1);
% b=a(1,1);
% c=zeros(b,1); 
%     
% StrikeSlipDisp  =0;      StrikeSlipDisp     = c+StrikeSlipDisp*-1;%Positive = RightLatMovement (any orientation)
% DipSlipDisp     =-1;      DipSlipDisp		= c+DipSlipDisp;        %Positive = reverse movement
% TensileSlipDisp =0;      TensileSlipDisp	= c+TensileSlipDisp;    %Positive = extensional movement
% clear a b c             	%a b and c are used to create an array of slips the size of the number of triangles on the fault surface.
% 
% Sxx = 0; 				
% Syy = 0; 
% Szz = 0; 
% Sxy = 0; 					
% Sxz = 0; 					
% Syz = 0;   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. Choose to define in stress or strain. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	%%%%%%%%%%%%%%
    %StrainInput        If halfspace = 1, Z stress components are removed
	%%%%%%%%%%%%%%
     
strain=0;                   %Put to 1 to define the stresses defined in 'stress input' as strain values
    
	%%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%
    

Sxx = 0;                    %Positive is tension
Syy = 0; 
Szz = 0;
Sxy = 0.5;        			%Positive = left lat when frac strikes NS and right lat when EW. See Pollard Fletcher Diagram - 6.13
Sxz = 0;
Syz = 0; 

Option='B'; 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option C = Friction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Sxx = 0;       			%Positive is tension 
% Syy = 0; 
% Szz = 0;
% Sxy = 0;        		%Positive = left lat when frac strikes NS and right lat when EW. See Pollard Fletcher Diagram - 6.13	
% Sxz = 0;
% Syz = 0; 
% % %Finding distance of points along X from 0,0. 
% X=(MidPoint(:,1));
% ne=size(MidPoint(:,1));
% Mu  = ones((ne))*0;    %Coefficient of friction
% Sf  = ones(size(Mu))*0;  
% figure;trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Sf);
% Option='C'; %slip from uniform remote stress with friction


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option D = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. This option has no opening components.  
	%Choose to define in stress or strain.  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%
    %StrainInput        If halfspace = 1, Z stress components are removed
	%%%%%%%%%%%%%%
    
% strain=0;                   %Put to 1 to define the stresses defined in 'stress input' as strain values
%     
% 
% 	%%%%%%%%%%%%%%
%     %StressInput
% 	%%%%%%%%%%%%%%
%      
% Sxx = -1;         			%Positive is tension
% Syy = 0; 
% Szz = 0; 
% Sxy = 0;                    %Positive = left lat when frac strikes NS and right lat when EW. See Pollard Fletcher Diagram - 6.13
% Sxz = 0;
% Syz = 0;
% 
% Option='D'; %slip no opening; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option D = Define tractions on triangles instead of slip. Inf code
    %works out the superposition of this pressure
	%Output 'stresses' are blank arrays of 0's as this options input stresses are not within the global CS    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
% Tss = 0;     %positive - right lateral
% Tds = 0;     %positive (normals point to sky) - drives slip down dip 
% Tnn = 1;     %positive is pressure that loads the boundary in tension
% 
% Option='E'; %Traction defined on elements (pressure)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[StrikeSlipDisp,DipSlipDisp,TensileSlipDisp,Sxx,Syy,Szz,...
Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tnn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removing any fixed triangles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removing any fixed triangles
if any(Fdisp)==1
[Triangles,FaceNormalVector,MidPoint,P1,P2,P3]...
    = RemovingFixedEls3d(Triangles,FaceNormalVector,MidPoint,P1,P2,P3,Fdisp);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Animate Fault movement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To do this just call 'Animate' in the cmd window. To change the
%parameters, save movie or the camera angle go into the 'Animate' function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting the calculated slips.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,StrikeSlipDisp,DipSlipDisp,TensileSlipDisp)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XYZ observation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Regular Grid that sits around the fault surfaces. Define how far it
    %extends away and its density
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%Calculate distances%%%%
cells=5;   %Define the number of observations points on the grid around the fault (not extended) %6
padding=1;  %How much extra bumph you want to add away from the grid. DISTANCE not cells.        %1
[maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ,sz]=MinMaxDataExtents3d(Points,cells,padding);

% %%%%3d grid%%%%
[X,Y,Z] = meshgrid(mingriX:sz:maxgriX,mingriY:sz:maxgriY,mingriZ:sz:maxgriZ);    
dimx = length(X(:,:,1)); %is this right?
dimy = length(X(:,1,:));
dimz = length(X(1,:,:));

% %%%%2d grid%%%%
% [X,Y] = meshgrid(mingriX:sz:maxgriX,mingriY:sz:maxgriY);    
% dimx = length(X(:,1));
% dimy = length(X(1,:));    
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = User defined XY grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%3d grid%%%%                 Regular and of your location and density choice
% Density=10;
% X = linspace(679000,681000,Density);
% Y = linspace(1504000,1506000,Density); 
% Z = linspace(-2000,0,Density); 
% [X,Y,Z] = meshgrid(X,Y,Z); 

%%%%2d grid%%%%                 Regular and of your location and density choice
%Density=10;
%X = linspace(675000,685000,Density);
%Y = linspace(1498000,1508000,Density); 
%[X,Y] = meshgrid(X,Y); % define Cartesian grid
%Z=zeros(size(X));

% %%%%Using midpoints%%%%
% X=MidPoint(:,1);
% Y=MidPoint(:,2);
% Z=MidPoint(:,3);


% %%%%Creating a flat line%%%% 
%  X = linspace(1,15,100); % add eps to avoid singularity at origin
%  Y = zeros(size(X)); 
%  Z = zeros(size(X));

%%%%Random Data%%%%
% % Random data at on n*n square grid centred at xmv,ymv
% Density=5000
% n=4;
% xmv=0; ymv=0;zmv=0;
% x=rand(1,Density)*n;
% y=rand(1,Density)*n;
% z=rand(1,Density)*n;
% X=x-n+xmv;
% Y=y-n+ymv;
% Z=z-n+zmv; %for 2d put me to 0;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option C = Import XYZ from a .text file within the folder you are
    %running this script from
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Loads the dat file, This contains 3 columns XYZ
	%XYZ are the points the stresses/straina and displements will be calculated on.

    
% S = load('CathalsIssueObs.dat');	
% X = S(:,1);                   
% Y = S(:,2);                         
% Z = S(:,3);  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option D = Load a secondary Fault Surface
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     %Secondary Surface
%     SecondSurface=1; %Flag set
%     string='RFaultEWStrike_1000Faces.ts';
%     [ PointsObs,TrianglesObs ] = GoCadAsciiReader( string,pathstring );
%     PointsObs=[PointsObs(:,1),PointsObs(:,2),PointsObs(:,3),(PointsObs(:,4)-freesurface_height)];
%     [MidPointObs,FaceNormalVectorObs] = MidPointCreate(PointsObs,TrianglesObs);
%     X=MidPointObs(:,1);
%     Y=MidPointObs(:,2);
%     Z=MidPointObs(:,3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing the obs Point data just created
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Drawing figure of surface and the observation points
figure;trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',(.2),'FaceColor', [0.5 0 0.9 ]);
hold on
scatter3(X(:),Y(:),Z(:),5,'k')  %Showing surface and obs points
xlabel('x'); ylabel('z'); axis('equal'); title('Surface and Obs Points');




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses and/or displacements on dispersed XYZ
%observation points. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%TotalStrainStress,StrainStressChange and StrainStressRemote are Xx12 arrays of strain tensors and stress tensors. XX YY ZZ XY XZ YZ
	%Strain is the first 6 coloumns and stress the last 6. StressStrainRemoteSxx is the regional stress 
	%StressStrainChange is the stress on the points from the event.
	%TotalStressStrain is the driving stress and stress change from the event added. 
	
 [TotalStrainStress,StrainStressChange,StrainStressRemote,Sxx,Syy,Szz,Sxy,Sxz,Syz]=CalculateStressOnSurroundingPoints(StrikeSlipDisp,...
 DipSlipDisp,TensileSlipDisp,mu,lambda,X,Y,Z,Sxx,Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu);

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Calculate Displacements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% [Ux,Uy,Uz] = CalculateDisplacementOnSurroundingPoints( StrikeSlipDisp,DipSlipDisp,TensileSlipDisp, nu, X,Y,Z, P1, P2, P3,halfspace);

 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 6: Finite strain computation, error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

    

% [exx,eyy,ezz,exy,exz,eyz,ExxInfErrorPerc,EyyInfErrorPerc,EzzInfErrorPerc,ExyInfErrorPerc,ExzInfErrorPerc,EyzInfErrorPerc]=FiniteStrainLagrangian3d(Ux,Uy,Uz,X,Y,Z);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: VISULISATION AND ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure,quiver3(X(:),Y(:),Z(:),Ux(:),Uy(:),Uz(:));
% xlabel('x'); ylabel('y'); axis('equal');
% title('Vectors showing displacement of points') ;

[exx,eyy,ezz,exy,exz,eyz,Sxx,Syy,Szz,Sxy,Sxz,Syz ] = ExtractCols( TotalStrainStress );


[S1,S2,S3,S1dir,S2dir,S3dir] = EigCalc3d(Sxx(:),Syy(:),Szz(:),Sxy(:),Sxz(:),Syz(:));

%If the stresses don't draw well use this function
%FilterValue=2.5;
%[S1,S2,Sxx,Syy,Sxy] = NanOrRemoveBadPoints( FilterValue,5,1,S1,S2,Sxx,Syy,Sxy );

%Function that draws stress ellipsoids
DrawStressEllipsoidsPrincipal(TotalStrainStress,X,Y,Z,Triangles,Points);

%Function that draws principal stress directions
DrawS1S2S3Directions(StrainStressChange,X,Y,Z,Triangles,Points);

%Check if data is on a uniform or non uniform grid
[ Uniform ] = UniformGridCheck3d( X,Y,Z );

%If data is random just just scatter. 
if Uniform==0
DrawScatterPlots3d( X,Y,Z,cmap, S1,S2,S3 )
%Drawing the stress on the 2nd surface if it exists
if SecondSurface==1
Mu=ones(size(Sxx))*0.6; %Coeff Friction
Cohesion=zeros(size(Sxx)); %Coeff Friction
[CSS] = CalculateCoulombStressOnPlane( X,Y,Z,FaceNormalVectorObs,Sxx,Syy,Szz,Sxy,Sxz,Syz,Mu,Cohesion,PointsObs,TrianglesObs,cmap2 );
end

%If data is uniform plot as gridded data.
else
%Reshaping stresses to grid dimensions
[X,Y,Z,S1,S2,S3,Sxx,Syy,Szz,Sxy,Sxz,Syz ]=ReshapeData3d( dimx,dimy,dimz,X,Y,Z,S1,S2,S3,Sxx,Syy,Szz,Sxy,Sxz,Syz );
%Drawing isocontours of the stress
[ hlink ] =IsoContoursPrincipalStressPercentiles( S1,S2,S3,X,Y,Z,Triangles,Points );
end


