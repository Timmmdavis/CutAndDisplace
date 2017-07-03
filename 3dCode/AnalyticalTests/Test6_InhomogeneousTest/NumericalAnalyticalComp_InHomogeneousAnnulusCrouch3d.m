%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you don’t need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

%   Copyright 2017, Tim Davis, The University of Aberdeen
clear;close all

%===== add file paths ==========================
pathstring = pwd;                                   %Get the address of the current working directory
parts = strsplit(mfilename('fullpath'), '\');       %Getting the address of the script working directory and splitting into cell array
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
    %The file must list the traingle vertex's then triangles. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  string='Tube_FixedInnnerTris.ts'; %Loading inner annlus, free bnd elastic 1
 [ PointsFBE1F,TrianglesFBE1F ] = GoCadAsciiReader( string,pathstring );
 FDsz=numel(PointsFBE1F(:,1));
    
        %VeryHighSampleSetup
%  string='Tube_1000mRad_4000FacesONVHS.ts'; %Loading inner annlus, free bnd elastic 1
% [ PointsFBE1,TrianglesFBE1 ] = GoCadAsciiReader( string,pathstring );		
%   string='Tube_2000mRad_4000FacesONVHS.ts'; %Loading interface elastic 1
% [ PointsIFE1E2,TrianglesIFE1E2 ] = GoCadAsciiReader( string,pathstring );		
%   string='Tube_2000mRad_4000FacesINVHS.ts'; %Loading interface elastic 2
% [ PointsIFE2E1,TrianglesIFE2E1 ] = GoCadAsciiReader( string,pathstring );		        
        
 string='Tube_1000mRad_1500FacesONVHS.ts'; %Loading inner annlus, free bnd elastic 1
[ PointsFBE1,TrianglesFBE1 ] = GoCadAsciiReader( string,pathstring );		
  string='Tube_2000mRad_1500FacesONVHS.ts'; %Loading interface elastic 1
[ PointsIFE1E2,TrianglesIFE1E2 ]= GoCadAsciiReader( string,pathstring );		
  string='Tube_2000mRad_1500FacesINVHS.ts'; %Loading interface elastic 2
[ PointsIFE2E1,TrianglesIFE2E1 ] = GoCadAsciiReader( string,pathstring );		


 %%%        HighSampleSetup 
%  string='Tube_1000mRad_1500FacesONHS.ts'; %Loading inner annlus, free bnd elastic 1
%  [ PointsFBE1,TrianglesFBE1 ] = GoCadAsciiReader( string,pathstring );		
%   string='Tube_2000mRad_1500FacesONHS.ts'; %Loading interface elastic 1
% [ PointsIFE1E2,TrianglesIFE1E2 ]  = GoCadAsciiReader( string,pathstring );		
%   string='Tube_2000mRad_1500FacesINHS.ts'; %Loading interface elastic 2
% [ PointsIFE2E1,TrianglesIFE2E1 ] = GoCadAsciiReader( string,pathstring );	

 
 
% %      %MediumSampleSetup
%   string='Tube_1000mRad_1500FacesONMS.ts'; %Loading inner annlus, free bnd elastic 1
%  [ PointsFBE1,TrianglesFBE1 ] = GoCadAsciiReader( string,pathstring );		
%    string='Tube_2000mRad_1500FacesONMS.ts'; %Loading interface elastic 1
%  [ PointsIFE1E2,TrianglesIFE1E2 ] = GoCadAsciiReader( string,pathstring );		
%    string='Tube_2000mRad_1500FacesINMS.ts'; %Loading interface elastic 2
%  [ PointsIFE2E1,TrianglesIFE2E1 ] = GoCadAsciiReader( string,pathstring );


% %    %LowSampleSetup
%  string='Tube_1000mRad_1500FacesONLS.ts'; %Loading inner annlus, free bnd elastic 1
% [ PointsFBE1,TrianglesFBE1 ] = GoCadAsciiReader( string,pathstring );		
%   string='Tube_2000mRad_1500FacesONLS.ts'; %Loading interface elastic 1
% [ PointsIFE1E2,TrianglesIFE1E2 ] = GoCadAsciiReader( string,pathstring );		
%   string='Tube_2000mRad_1500FacesINLS.ts'; %Loading interface elastic 2
% [ PointsIFE2E1,TrianglesIFE2E1 ] = GoCadAsciiReader( string,pathstring );	
 	
%Adding Fixed parts of E1 to E1 inner ann
[Points,Triangles,BoundaryFlag] = DataAppender3d( PointsFBE1,PointsFBE1F,TrianglesFBE1,TrianglesFBE1F );

%Adding interface E1E2, Points towards E2 
[Points,Triangles,BoundaryFlag] = DataAppender3d( Points,PointsIFE1E2,Triangles,TrianglesIFE1E2,BoundaryFlag,2 );

%Adding interface E2E1, Points towards E1
[Points,Triangles,BoundaryFlag] = DataAppender3d( Points,PointsIFE2E1,Triangles,TrianglesIFE2E1,BoundaryFlag,3 );

% %How to add Fixed parts of E2 to data
% [Points,Triangles,BoundaryFlag] = DataAppender3d( Points,PointsFBE1F,Triangles,TrianglesFBE1F,BoundaryFlag,5 );
% %BoundaryFlag(Cng)=5; %5 is the fixed els for elastic 2

%Quick Low Down on BoundaryFlag
%0=free boundary E1
%1=fixed bits of the free boundary of E1
%2=E1-E2 interface, E1 elastic properties, normals point towards E2
%3=E2-E1 interface, E2 elastic properties, normals point towards E1
%4=If existed would be free boundary E2, 
%5=fixed bits of the free boundary of E2
%Only structured for two elastics at the moment. Introducing E3 would mean there are
%potentially 6 interfaces. E1-E2 E1-E3, E2-E1 E2=E3, E3-E1 E3-E2


 

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
freesurface_height = 0;

    %Defining elastic constants, All three are required to run this whole script.  See equations below if you need to calculate these 
	%from other constants. 
	%Note nu is needed in the displacement calculation.  
	
mu=500;%0.4;%4000 	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
lambda=mu*((2*nu)/(1-2*nu));%4000;    %Lame's constant  0 perfectly compressible like cork, Infinate for incompressible material like rubber.

	%Equations you can use to calculate elastic parameters that you need above. 
	
% E = mu*(2*(1+nu)) ;               %Young's Modulus
% K = E/3*(1-2*nu);                 %Bulk Modulus.    Equation 8.25 Pollard
% mu = E/(2*(1+nu));                %Shear Modulus.   Equation 8.26 Pollard  which is Geo convention 
% lambda= E*nu/((1+nu)*(1-2*nu));   %Lamé's  constant,  Equation 8.27 Pollard
% nu =lambda/(2*(mu+lamda);            %Poisson's ratio, Equation 8.28 Pollard

	%Function to extract the three vertex's for each triangle and 
	%triangles defined by a combination of vertex Id's. This is used later by the algorithms 
    %The max Z value is output as text in the command window giving an idea
    %of minimum free surface height. 
Points=[Points(:,1),Points(:,2),Points(:,3),(Points(:,4)-freesurface_height)];
%Creating Data
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 

title('All Surfaces Loaded');

	%Fictitious disp flag, if set to one then this forces the triangles with this to 0 displacement. 
Fdisp=(BoundaryFlag == 1 | BoundaryFlag == 5);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, Remote, Remote
%no opening, tractions on the elements or user imported slip value for each face. 

%Stress tensor outputs from this part of the toolbox i.e. Sxx Syy etc are simply the user defined remote stress
%components (boundary conditions). These are used later in the stress calculations. If strain is input these are converted to stress 
%using Hooke's law. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[Sxx,Syy,Szz,Sxy,Sxz,Syz,Tnn,Tss,Tds,Mu,Sf,strain ] = CreateBlankVars;
	%%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%

%Elastic 2 properties
nu2=0.25;
mu2=100; %0.5
lambda2=mu2*((2*nu2)/(1-2*nu2));%4000;  %Lame's constant  0 perfectly compressible like cork, Infinate for incompressible material like rubber.

%Appending these to the back of E1
nu=[nu;nu2];
mu=[mu;mu2];
lambda=[lambda;lambda2];
  
    
Sxx = 0.002;       			%Positive is tension
Syy = 0.002; 
Szz = 0; %plane strain conditions
Sxy = 0;        			%Positive = left lat when frac strikes NS and right lat when EW. See Pollard Fletcher Diagram - 6.13
Sxz = 0;
Syz = 0; 
Option='F';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[StrikeSlipDisp,DipSlipDisp,TensileSlipDisp,Sxx,Syy,Szz,...
Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tnn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,BoundaryFlag,strain,Mu,Sf,Option);

%Removing fixed disps from BoundaryFlag
BoundaryFlag(Fdisp)=[]; %removing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removing any fixed triangles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removing any fixed triangles
if any(Fdisp)==1
[Triangles,FaceNormalVector,MidPoint,P1,P2,P3]...
    = RemovingFixedEls3d(Triangles,FaceNormalVector,MidPoint,P1,P2,P3,Fdisp);
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting the calculated slips.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,StrikeSlipDisp,DipSlipDisp,TensileSlipDisp)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XYZ observation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Simply define flat observation plane of points. XYZ with defined step size 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear size
[X,Y] = meshgrid(-3000:50:3000,-3000:50:3000);
% % % [X,Y] = meshgrid(1200:50:3500,-3581:50:-1500);
dimx = length(X(:,1));
dimy = length(X(1,:));
Y=Y-freesurface_height;


%Creating the radial coordinates for these obs points, these will be used
%to filter these. 
[TH,R] = cart2pol(X,Y);


%Removing points in hole
Bad=R<1000;
X=X(~Bad);Y=Y(~Bad);TH=TH(~Bad);R=R(~Bad);
Theta=radtodeg(TH);%figure;scatter(X(:),Y(:),20,Theta(:))

%Now only selecting a certain portion of points (Makes the figure neater)
%Good=Theta<135 & Theta>45; %90 deg rad
% Good=Theta<15 & Theta>-15;   %4 deg rad
% X=X(Good);Y=Y(Good);R=R(Good);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Removing values within in 100m of dislocations
disp('Any points closer than 100m from dislocations removed')
%Removing points close to freebnd1
Bad=R<1100;
X=X(~Bad);Y=Y(~Bad);Theta=Theta(~Bad);R=R(~Bad);
%Removing points close to interface where BEM solution breaks down
Bad=R<2100 & R>1900;
X=X(~Bad);Y=Y(~Bad);Theta=Theta(~Bad);R=R(~Bad);
%Removing points really far away so the graph is nicer
Bad=R>3000;
X=X(~Bad);Y=Y(~Bad);Theta=Theta(~Bad);R=R(~Bad);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Flagging points in 2nd elastic
OuterPart=R>2000;
XE2=X(OuterPart);YE2=Y(OuterPart);
XE1=X(~OuterPart);YE1=Y(~OuterPart);
ZE1=zeros(size(YE1));
ZE2=zeros(size(YE2));

%Drawing the observation points on one of the previous figures
%Octave can't handle transparent objects
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1;
scatter(XE1(:),YE1(:),'b');
scatter(XE2(:),YE2(:),'r');
elseif isOctave==0;
figure;scatter(XE1(:),YE1(:),'b','filled','MarkerFaceAlpha',1/8);hold on
scatter(XE2(:),YE2(:),'r','filled','MarkerFaceAlpha',1/8);hold off
end
title('Observation Points, Blue E1 Red E2'), xlabel('x'), ylabel('y')
        

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
	    
%E1 Elastic 1
E1Bits =(BoundaryFlag == 0 | BoundaryFlag == 2);
TensileSlipDispE1=TensileSlipDisp(E1Bits );
StrikeSlipDispE1=StrikeSlipDisp(E1Bits );
DipSlipDispE1=DipSlipDisp(E1Bits );
Part=1; %Which elastic we want to extract
[MidPointE1,P1E1,P2E1,P3E1,muE1,lambdaE1,FaceNormalVectorE1,nuE1,NUME1,FdispE1,FB_E1,IF_E1]...
= ExtractElasticParts3d( Part,BoundaryFlag,MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,Fdisp );
    
%Calculating the stress on the points that make up the annulus (elastic 1).
%0 input stress as the original analytical calc is internal pressure not
%boundary stress. its the same Tn anyway for this geometry
 [TotalStrainStressE1,StrainStressChangeE1,StrainStressRemoteE1,SyyE1,SzzE1,SxyE1,SxzE1,SyzE1]...
     =CalculateStressOnSurroundingPoints(StrikeSlipDispE1,DipSlipDispE1,TensileSlipDispE1,muE1,lambdaE1,XE1,YE1,ZE1,0,0,0,0,0,0,P1E1,P2E1,P3E1,halfspace,nuE1);

 %Extract stresses
 SxxE1=StrainStressChangeE1(:,7)';SyyE1=StrainStressChangeE1(:,8)';SxyE1=StrainStressChangeE1(:,10)';
 sz=numel(Sxx);
 
 
 %E2 Elastic 2
E2Bits =(BoundaryFlag == 3 | BoundaryFlag == 4);
TensileSlipDispE2=TensileSlipDisp(E2Bits );
StrikeSlipDispE2=StrikeSlipDisp(E2Bits );
DipSlipDispE2=DipSlipDisp(E2Bits );
Part=2; %Which elastic we want to extract
[MidPointE2,P1E2,P2E2,P3E2,muE2,lambdaE2,FaceNormalVectorE2,nuE2,NUME2,FdispE2,FB_E2,IF_E2]...
= ExtractElasticParts3d( Part,BoundaryFlag,MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,Fdisp );
    
%Calculating the stress on the points that make up the annulus (elastic 1).
%0 input stress as the original analytical calc is internal pressure not
%boundary stress. its the same Tn anyway for this geometry
 [TotalStrainStressE2,StrainStressChangeE2,StrainStressRemoteE2,SyyE2,SzzE2,SxyE2,SxzE2,SyzE2]...
     =CalculateStressOnSurroundingPoints(StrikeSlipDispE2,DipSlipDispE2,TensileSlipDispE2,muE2,lambdaE2,XE2,YE2,ZE2,0,0,0,0,0,0,P1E2,P2E2,P3E2,halfspace,nuE2);

 %Extract stresses
 SxxE2=StrainStressChangeE2(:,7)';SyyE2=StrainStressChangeE2(:,8)';SxyE2=StrainStressChangeE2(:,10)';
 
%Now appending the two sets of cartesian stress components
Sxx=[SxxE1,SxxE2];
Syy=[SyyE1,SyyE2];
Sxy=[SxyE1,SxyE2];

%Now appending the two sets X Y points, putting in the same orientation as
%the tensors (rows)

X=[XE1;XE2]';
Y=[YE1;YE2]';


%%%%%%%%%%%%%%%%%%
%Creating radial coordinates for each observation point and then converting
%tensors from cart to polar coordinates
[TH,R] = cart2pol(X,-Y);
[ SRR,STT,SRT ] = StressTensorTransformation2d(Sxx,Syy,Sxy,cos(TH),cos((pi/2)-TH));


Tn=0.002;

%%%%
%Normalizing
%%%%
SrrNorm_a_r_b=SRR/Tn(1,1);
SttNorm_a_r_b=STT/Tn(1,1);

rOvb_a_r_b=R/2000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating and Plotting analytical results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
[rOvb,SrrNorm,SttNorm]=AnnulusEq_Func(muE1,nuE1,muE2,nuE2);
figure;plot(rOvb,SrrNorm,'LineWidth',4,'color','blue');
hold on
plot(rOvb,SttNorm,'LineWidth',4,'color','m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Plotting numerical results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%Octave can't handle transparent objects
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1;
scatter(rOvb_a_r_b(:),SrrNorm_a_r_b(:),2,'k'); 
scatter(rOvb_a_r_b(:),SttNorm_a_r_b(:),2,'k');
elseif isOctave==0;
%scatter(rOvb_a_r_b(:),SrrNorm_a_r_b(:),12,'filled','MarkerFaceAlpha',1/8);
%scatter(rOvb_a_r_b(:),SttNorm_a_r_b(:),12,'filled','MarkerFaceAlpha',1/8);
scatter(rOvb_a_r_b(:),SrrNorm_a_r_b(:),24,'k','filled','MarkerFaceAlpha',1/8);
scatter(rOvb_a_r_b(:),SttNorm_a_r_b(:),24,'k','filled','MarkerFaceAlpha',1/8);
end
title('Analytical - lines, Srr (blue) Stt (pink), Numerical result as points'), xlabel('rad/b')
ylabel('Srr Stt Norm');
titlesz=25;
fntsz=21;
grid on
ChangeFontSizes(fntsz,titlesz);

%Calculating analytical solution stresses for points used in numerical calc
[TH,Rin] = cart2pol(XE1,YE1);[TH,ROut] = cart2pol(XE2,YE2);clear TH
%Calculating soultion for points same as points used in numerical calc
[rOvbPnts,SrrNormPnts,SttNormPnts]=AnnulusEq_Func_pnts(muE1,nuE1,muE2,nuE2,(Rin')/1e3,(ROut')/1e3);
%Calculating polar stress residuals.
ResidualSrr=SrrNormPnts-SrrNorm_a_r_b;
ResidualStt=SttNormPnts-SttNorm_a_r_b;
% plot(rOvbPnts,ResidualSrr,'--')
% plot(rOvbPnts,ResidualStt,'--')

%disp(max(ResidualSrr));
fprintf('MaxResidualSrr %i.\n',max(abs(ResidualSrr(:))))
%disp(max(ResidualStt));
fprintf('MaxResidualStt %i.\n',max(abs(ResidualStt(:))))
biggestres=max(ResidualStt)+max(ResidualSrr);

if  abs(biggestres) > 0.2
    error('Not a good match to analytical solution')
else
    disp('Everything looks good, tolerance checks max polar stress residuals and flags errors above 0.2')
end

