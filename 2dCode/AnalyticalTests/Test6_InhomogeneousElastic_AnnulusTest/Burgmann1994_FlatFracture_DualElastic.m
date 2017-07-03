% Test: Test shows how to run with multiple elastic bodies and that this is setup to work
% correctly. A hole with an internal pressure is sits in the centre of a
% ring with elastic property 1. This ring sits inside a matrix of elastic property 2. 
% A solution is calculated and the radial stresses at all points within an observation grid
% are calculated. 
% 
% Solution: Crouch, S.L. and Starfield, A.M., 1982. Boundary element methods in solid 
% mechanics: with applications in rock mechanics and geological engineering. Allen & Unwin.
%  
% Proof: This tests that this has been correctly implemented and is the
% correct methodology for multiple elastics. A slightly cleaner formulation
% (restructured) could be used but this is a good starting point. 

%   Copyright 2017, Tim Davis, The University of Aberdeen
% To Do: 
% Check if removing observation points a certain distance from the radial
% tips fixes spurious stress values. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you don’t need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

%===== add file paths ==========================
AdRs0=mfilename('fullpath'); %Directory containing this script
parts = strsplit(AdRs0, '\');
DirPart = parts(1,1:end-1); %Removing the file name
AdRs0 = strjoin(DirPart,'\');
cd(AdRs0) %Making sure we are in the place we are running from.

cd('..') % go up one level
cd('..') % go up one level
AdRs05 = pwd;
cd('..') % go up one level 
AdRs1 = pwd; % get the address of the new directory you're in
cd(AdRs0) % back to your original position
addpath(genpath(AdRs1))%Adds all code directorys and subfolders to the path

close all; clear 

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
%STEP 1: import the fault surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading XY ascii data or manually creating fractures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%First free boundary    
sz=400;     
Pointsxy=[linspace(-100,0,sz)',zeros(sz,1)];
mystruct.('line1')=(1:sz);
E1Fbnd=Pointsxy;

%Interface
PointsxyIFE1E2=[zeros(sz,1),linspace(-100,100,sz)'];

BoundaryFlag=zeros((sz-1),1);
%Adding the interface of E1-E2 in
[Pointsxy,mystruct,BoundaryFlag] = DataAppender2d( Pointsxy,PointsxyIFE1E2,mystruct,BoundaryFlag,2 );


%Free boundary elastic 1 is now 0's the interface is 2 in out var called
%BoundaryFlag

%Second interface
PointsxyIFE2E1=PointsxyIFE1E2;

%Adding the interface of E2-E1
[Pointsxy,mystruct,BoundaryFlag] = DataAppender2d( Pointsxy,PointsxyIFE2E1,mystruct,BoundaryFlag,3 );
names = fieldnames(mystruct);
numstructs=numel(names);
FlipNormalsFlag=zeros(numstructs,1);
FlipNormalsFlag(end,1)=1; %Flag that tells the line builder to flip the normals of the interface


%2nd elastic free boundary data
PointsxyFree2=[E1Fbnd(:,1)+100,E1Fbnd(:,2)];

[Pointsxy,mystruct,BoundaryFlag] = DataAppender2d( Pointsxy,PointsxyFree2,mystruct,BoundaryFlag,4 );
FlipNormalsFlag=[FlipNormalsFlag;0];% We don’t want to flip the last normals

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
    
halfspace = 0; 

	%Defining the height of the freesurface relative to the value of 0 of
	%your imported surfaces. Used when deformation is not at
	%sealevel. See the output MaxZ height if your half space clashes with
	%your points and adjust until this outputs as a negative value. 
    %Use this if the algorithm is throwing errors for vertex's above 0m.
    
freesurface_height = 0;


Pointsxy(:,2)=Pointsxy(:,2)-freesurface_height;

    %Defining elastic constants,  See equations below if you need to calculate these 
	%from other constants. 
	
mu=1; % 1%0.4;%4000 	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
E = mu*(2*(1+nu)) ;                     %Young's Modulus

	%Equations you can use to calculate elastic parameters that you need above. 
% K = E/3*(1-2*nu);                 %Bulk Modulus.    Equation 8.25 Pollard
% mu = E/(2*(1+nu));                %Shear Modulus.   Equation 8.26 Pollard  which is Geo convention 
% lambda= E*nu/((1+nu)*(1-2*nu));   %Lamé's  constant,  Equation 8.27 Pollard
% nu =lambda/(2*(mu+lamda);            %Poisson's ratio, Equation 8.28 Pollard

    %Reshaping the list of XY into usable variables
[ x,y,xe,ye,HalfLength,Beta,CosB,Points,NormAng,NUM,XBEG,XEND,YBEG,YEND]...
 = CreateElements2d( Pointsxy,mystruct,FlipNormalsFlag );

IFFlagE1=(BoundaryFlag == 2 );
XE1=x(IFFlagE1);
YE1=y(IFFlagE1);
figure;

subplot(1,2,1);scatter(XE1,YE1,15,(1:numel(XE1)));
title('Interface element numbering, must be the same'), xlabel('x'), ylabel('y');axis equal

IFFlagE2=(BoundaryFlag == 3 );
XE2=x(IFFlagE2);
YE2=y(IFFlagE2);
subplot(1,2,2);scatter(XE2,YE2,15,(1:numel(XE2)));
title('Interface element numbering, must be the same'), xlabel('x'), ylabel('y');axis equal


    %Plotting figure of the the loaded fracture/fractures
   
figure;
hold on 
line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','r')
nx=cos(NormAng);
ny=cos((pi/2)-NormAng);
title('fractures'), xlabel('x'), ylabel('y');axis equal
%Fixed elements shown as blue
Fdisp=(BoundaryFlag == 1 | BoundaryFlag == 5);
line([Points(Fdisp,1)';Points(Fdisp,2)'],[Points(Fdisp,3)';Points(Fdisp,4)'],'color','b')
%1st elastic normals shown as blue
E1Flg=(BoundaryFlag == 0 | BoundaryFlag == 1 | BoundaryFlag == 2);
quiver(xe(E1Flg),ye(E1Flg),nx(E1Flg),ny(E1Flg),'color','b')
%2nd elastic normals shown as green
E2Flg=(BoundaryFlag == 3 | BoundaryFlag == 4 | BoundaryFlag == 5);
quiver(xe(E2Flg),ye(E2Flg),nx(E2Flg),ny(E2Flg),'color','g')
hold off


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now passing these into the solver to find the unknown displacements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   
%
[SxxE1,SyyE1,SxyE1,Tn,Ts,Mu,Sf,strain ] = CreateBlankVars;

%Elastic 2 properties
nu2=0.25;
mu2=7; %0.5
E2 = mu2*(2*(1+nu2)) ;

%Appending these to the back of E1
nu=[nu;nu2];
mu=[mu;mu2];
E=[E;E2];

 

SxxE1 = 0;  SyyE1 = 0;  SxyE1 = 1; 
Option='F'; %Inhomogeneous Elastic


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ShearDisp,TensileDisp,SxxE1,SyyE1,SxyE1]=SlipCalculator2d(x,y,xe,ye,HalfLength,Beta,SxxE1,SyyE1,SxyE1,Tn,Ts,nu,E,halfspace,NUM,NormAng,strain,BoundaryFlag,Mu,Sf,Option);

%Removing fixed disps from BoundaryFlag
BoundaryFlag(Fdisp)=[]; %removing

%%
% Removing any fixed elements
if any(Fdisp)==1  
[XBEG,XEND,YBEG,YEND,NUM,x,y,xe,ye,HalfLength,Beta,CosB,Points,NormAng]...
    = RemovingFixedEls2d(Fdisp,XBEG,XEND,YBEG,YEND,NUM);
end

IFFlag=(BoundaryFlag == 2 | BoundaryFlag == 3 );

TensileDispNoint=TensileDisp(~IFFlag);
ShearDispNoint=ShearDisp(~IFFlag);
NumberOfFreeEls=(sum(~IFFlag));

% Plot displacement discontinuity magnitude for all elements
GraphSlipVsElNum( NumberOfFreeEls,TensileDispNoint,abs(ShearDispNoint) )


%Creating directions and magnitudes of slip for plotting
PlotOpeningVsShearOnEls( NormAng,TensileDisp,ShearDisp,Points,x,y,HalfLength )

%%
%Extracting elements for calculations
%E1 Elastic 1
E1Bits =(BoundaryFlag == 0 | BoundaryFlag == 2);
TensileDispE1=TensileDisp(E1Bits );
ShearDispE1=ShearDisp(E1Bits );
Part=1; %Which elastic we want to extract
[xE1,yE1,xeE1,yeE1,HalfLengthE1,BetaE1,nuE1,EE1,NormAngE1,NUME1,FdispE1,FB_E1,IF_E1] = ExtractElasticParts2d( Part,BoundaryFlag,x,y,xe,ye,HalfLength,Beta,nu,E,NormAng );
muE1 = EE1/(2*(1+nuE1));


%E2 Elastic 2
E2Bits =(BoundaryFlag == 3 | BoundaryFlag == 4);
TensileDispE2=TensileDisp(E2Bits );
ShearDispE2=ShearDisp(E2Bits );
Part=2; %Which elastic we want to extract
[xE2,yE2,xeE2,yeE2,HalfLengthE2,BetaE2,nuE2,EE2,NormAngE2,NUME2,FdispE2,FB_E2,IF_E2] = ExtractElasticParts2d( Part,BoundaryFlag,x,y,xe,ye,HalfLength,Beta,nu,E,NormAng );
muE2 = EE2/(2*(1+nuE2));


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating observation points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   
%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Regular Grid that sits around the fault surfaces. Define how far it
    %extends away and its density
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
cells=75;   %Define the number of observations points on the square grid
padding=25;  %How much extra bumph you want to add away from the fault sticks, reduce to 0 if having half space issues 
[maxgriX,mingriX,maxgriY,mingriY,sz]=MinMaxDataExtents(Points,cells,padding);
[X,Y] = meshgrid(mingriX:sz:maxgriX,mingriY+150:sz:maxgriY-150);
rowcount = size(X,1);
colcount = size(X,2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing and fixing Obs Point data just created
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Replacing XY of any obs points for Nan if these lie on top of the end points of a
    %dislocation, without this it can result in spurious values. 
[ X,Y  ] = NullPointsLyingOnElement( X,Y,XBEG,YBEG,XEND,YEND,NUM );

E1Flag=X<0;
XE1=X(E1Flag);
YE1=Y(E1Flag);

XE2=X(~E1Flag);
YE2=Y(~E1Flag);

%Octave can't handle transparent objects
figure;line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','r')
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating the stress on the observation points, elastic 1 then elastic 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

   
%Calculating the stress on the points that make up the annulus (elastic 1).
[StressChangeE1,TotalStressE1,StressRegE1,InducedDisplacementsE1]=StressDispOnSurroundingPoints(XE1,YE1,xeE1,yeE1,HalfLengthE1,SxxE1,SyyE1,SxyE1,nuE1,EE1,halfspace,ShearDispE1,TensileDispE1,BetaE1);

%Calculating the stress on the points that make up the external matrix (elastic 2).
[StressChangeE2,TotalStressE2,StressRegE2,InducedDisplacementsE2]=StressDispOnSurroundingPoints(XE2,YE2,xeE2,yeE2,HalfLengthE2,0,0,0,nuE2,EE2,halfspace,ShearDispE2,TensileDispE2,BetaE2);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: VISULISATION AND ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Drawing figures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Collating data
X=[XE1;XE2];
Y=[YE1;YE2];

Ux=[InducedDisplacementsE1(:,1);InducedDisplacementsE2(:,1)];
Uy=[InducedDisplacementsE1(:,2);InducedDisplacementsE2(:,2)];

Sxx=[StressChangeE1(:,1);StressChangeE2(:,1)];
Syy=[StressChangeE1(:,2);StressChangeE2(:,2)];
Sxy=[StressChangeE1(:,3);StressChangeE2(:,3)];


    %Elastic constants for each pnt
ObsE1=ones(numel(XE1),1);
ObsE2=ones(numel(XE2),1);

mu=[ObsE1*(mu(1));ObsE2*(mu(2))];
E=[ObsE1*(E(1));ObsE2*(E(2))];
nu=[ObsE1*(nu(1));ObsE2*(nu(2))];

lambda= E.*nu./((1+nu).*(1-2.*nu));   %Lamé's  constant,  Equation 8.27 Pollard

%Stress2strain
[ Exx,Eyy,Exy ] = Hooke'sLaw2dStress2Strain( Sxx,Syy,Sxy,lambda,mu );

%Calclating 2d EigenValues
[S1,S2,S1dir,S2dir]=EigCalc2d(Sxx,Syy,Sxy);

%Calclating 2d EigenValues
[E1,E2,E1dir,E2dir]=EigCalc2d(Exx,Eyy,Exy);
Dilatation2d=E1+E2;

%If the stresses don’t draw well use this function
FilterValue=2.5;
[S1,S2,Sxx,Syy,Sxy] = NanOrRemoveBadPoints( FilterValue,5,1,S1,S2,Sxx,Syy,Sxy );


%Draw data with quiver, doesn't matter if its gridded or not
figure,quiver(X(:),Y(:),Ux,Uy);
xlabel('x'); ylabel('y'); axis('equal'); title('Disp'); 

%Drawing S1 S2 directions
DrawS1S2Directions(X(:),Y(:),S1dir,S2dir,Points )


%Check if data is on a uniform or non uniform grid
[ Uniform ] = UniformGridCheck2d( X,Y );
%If data is not uniform plot as non grid data. 
if Uniform==0
    
DrawScatterPlots2d( X,Y,cmap2,Sxy,Syy,Sxx,S1,S2 );

%else reshapefrom col vectors to meshes then draw with nice grid functions
else

%Reshaping stresses to grid dimensions
[Sxy,Syy,Sxx,S1,S2]=ReshapeData2d( rowcount,colcount,Sxy,Syy,Sxx,S1,S2 );

DrawContourFPlots2d( X,Y,cmap2,Sxy,Syy,Sxx,S1,S2 );

Scl=1;
DrawDeformedGrid2d( X,Y,Ux,Uy,Scl,cmap2,(abs(Ux)+abs(Uy))); %Dilatation2d

end
