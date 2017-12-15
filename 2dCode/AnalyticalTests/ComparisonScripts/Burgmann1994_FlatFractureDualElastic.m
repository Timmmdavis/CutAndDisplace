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
%STEP 1: import the fault surface.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading XY ascii data or manually creating fractures.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
%First free boundary    
sz=100;     
Pointsxy=[linspace(-100,0,sz)',zeros(sz,1)];
mystruct.('line1')=(1:sz);
E1Fbnd=Pointsxy;

BoundaryFlag=zeros((sz-1),1);

%Interface
PointsxyIFE1E2=[zeros(sz,1),linspace(-250,250,sz)'];


%Adding the interface of E1-E2 in
[Pointsxy,mystruct,BoundaryFlag] = DataAppender2d( Pointsxy,PointsxyIFE1E2,mystruct,BoundaryFlag,2 );


%Free boundary elastic 1 is now 0's the interface is 2 in out var called
%'BoundaryFlag'

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
FlipNormalsFlag=[FlipNormalsFlag;0];% We do not want to flip the last normals

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
Pointsxy(:,2)=Pointsxy(:,2)-freesurface_height;

% Defining elastic constants
mu=1;           	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
[ ~,E,~,nu,mu ] = ElasticConstantsCheck( mu,nu );

%Reshaping the list of XY into usable variables
[MidPoint,HalfLength,P1,P2,LineNormalVector,Lines,Points,BoundaryFlag  ]...
 = CreateElements2d( Pointsxy,mystruct,BoundaryFlag, FlipNormalsFlag);

IFFlagE1=(BoundaryFlag == 2 );
XE1=MidPoint(IFFlagE1,1);
YE1=MidPoint(IFFlagE1,2);
figure;

subplot(1,2,1);scatter(XE1,YE1,15,(1:numel(XE1)));
title('Interface element numbering, must be the same'), xlabel('x'), ylabel('y');axis equal

IFFlagE2=(BoundaryFlag == 3 );
XE2=MidPoint(IFFlagE2,1);
YE2=MidPoint(IFFlagE2,2);
subplot(1,2,2);scatter(XE2,YE2,15,(1:numel(XE2)));
title('Interface element numbering, must be the same'), xlabel('x'), ylabel('y');axis equal


%Plotting figure of the the loaded fracture/fractures
figure;
hold on 
PlotFracture( P1,P2,'r' )
CosAx=LineNormalVector(:,1);
CosAy=LineNormalVector(:,2);
title('fractures'), xlabel('x'), ylabel('y');axis equal
%Fixed elements shown as blue
Fdisp=find(BoundaryFlag == 1 | BoundaryFlag == 5);
PlotFracture(P1(Fdisp,:),P2(Fdisp,:),'b' )
%1st elastic normals shown as blue
E1Flg=find(BoundaryFlag == 0 | BoundaryFlag == 1 | BoundaryFlag == 2);
quiver(MidPoint(E1Flg,1),MidPoint(E1Flg,2),CosAx(E1Flg),CosAy(E1Flg),'color','b')
%2nd elastic normals shown as green
E2Flg=find(BoundaryFlag == 3 | BoundaryFlag == 4 | BoundaryFlag == 5);
quiver(MidPoint(E2Flg,1),MidPoint(E2Flg,2),CosAx(E2Flg),CosAy(E2Flg),'color','g')
hold off

Fdisp=(BoundaryFlag == 1 | BoundaryFlag == 5);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now passing these into the solver to find the unknown displacements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   
%
[~,~,~,Tn,Ts,Mu,Sf,strain ] = CreateBlankVars;

%Elastic 2 properties
nu2=0.25;
mu2=7; %0.5
E2 = mu2*(2*(1+nu2)) ;

%Appending these to the back of E1
nu=[nu;nu2];
mu=[mu;mu2];
E=[E;E2];

SxxE1 = 0;  SyyE1 = 0;  SxyE1 = 1; 
Option='F'; 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ds,Dn,~,~,~]=SlipCalculator2d(...
    MidPoint,HalfLength,SxxE1,SyyE1,SxyE1,Tn,Ts,nu,E,halfspace,LineNormalVector,strain,BoundaryFlag,Mu,Sf,Option);

%Removing fixed disps from BoundaryFlag
BoundaryFlag(Fdisp)=[]; %removing

%%
% Removing any fixed elements
if any(Fdisp)==1  
    [MidPoint,HalfLength,P1,P2,LineNormalVector]...
    = RemovingFixedEls2d(Fdisp,P1,P2);
end


IFFlag=(BoundaryFlag == 2 | BoundaryFlag == 3 );

DnNoint=Dn(~IFFlag);
DsNoint=Ds(~IFFlag);
NumberOfFreeEls=(sum(~IFFlag));

% Plot displacement discontinuity magnitude for all elements
GraphSlipVsElNum( Dn,Ds )

%Creating directions and magnitudes of slip for plotting
PlotOpeningVsShearOnEls( LineNormalVector,Dn,Ds,P1,P2,MidPoint,HalfLength )


%%
%Extracting elements for calculations
%E1 Elastic 1
E1Bits =(BoundaryFlag == 0 | BoundaryFlag == 2);
DnE1=Dn(E1Bits );
DsE1=Ds(E1Bits );
Part=1; %Which elastic we want to extract
[MidPointE1,HalfLengthE1,nuE1,EE1,LineNormalVectorE1,NUME1,FdispE1,FB_E1,IF_E1]...
 = ExtractElasticParts2d( Part,BoundaryFlag,MidPoint,HalfLength,nu,E,LineNormalVector );
muE1 = EE1/(2*(1+nuE1));

%E2 Elastic 2
E2Bits =(BoundaryFlag == 3 | BoundaryFlag == 4);
DnE2=Dn(E2Bits );
DsE2=Ds(E2Bits );
Part=2; %Which elastic we want to extract
[MidPointE2,HalfLengthE2,nuE2,EE2,LineNormalVectorE2,NUME2,FdispE2,FB_E2,IF_E2]...
 = ExtractElasticParts2d( Part,BoundaryFlag,MidPoint,HalfLength,nu,E,LineNormalVector );
muE2 = EE2/(2*(1+nuE2));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating observation points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   
%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Regular/randomly spaced grid of points bounding the faults.
    %Define how far it extends away and its density.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%Define the number of cells on the square grid
cells=75;   
%How many extra cells you want to add away from the faults, 
%reduce to 0 if having half space issues 
padding=25;  
[maxgriX,mingriX,maxgriY,mingriY,sz]=MinMaxDataExtents2d(Points,cells,padding);
[X,Y] = meshgrid(mingriX:sz:maxgriX,mingriY:sz:maxgriY);
dimx = size(X,1);
dimy = size(X,2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing and fixing Obs Point data just created
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dist=4;
[X,Y]=NanTolDistLine2Pnt( X,Y,P1,P2,LineNormalVector,Dist );    
    
E1Flag=X<0;
XE1=X(E1Flag);
YE1=Y(E1Flag);

XE2=X(~E1Flag);
YE2=Y(~E1Flag);

%Octave can't handle transparent objects
figure;PlotFracture( P1,P2,'r' )
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
%Calculating the stress on the observation points, elastic 1 then elastic 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

   
%Calculating the stress on the points that make up the annulus (elastic 1).
[StressTTotalE1,StrainTTotalE1,StressTChgE1,StrainTChgE1,StressTRegE1,StrainTRegE1]...
    = CalculateStressesOnSurroundingPoints2d(XE1,YE1,MidPointE1,HalfLengthE1,0,0,0,nuE1,EE1,halfspace,DsE1,DnE1,LineNormalVectorE1);

[UxE1,UyE1] = CalculateDisplacementOnSurroundingPoints2d(XE1,YE1,MidPointE1,HalfLengthE1,nuE1,EE1,halfspace,DsE1,DnE1,LineNormalVectorE1);
	

%Calculating the stress on the points that make up the external matrix (elastic 2).
[StressTTotalE2,StrainTTotalE2,StressTChgE2,StrainTChgE2,StressTRegE2,StrainTRegE2]...
    = CalculateStressesOnSurroundingPoints2d(XE2,YE2,MidPointE2,HalfLengthE2,0,0,0,nuE2,EE2,halfspace,DsE2,DnE2,LineNormalVectorE2);

[UxE2,UyE2] = CalculateDisplacementOnSurroundingPoints2d(XE2,YE2,MidPointE2,HalfLengthE2,nuE2,EE2,halfspace,DsE2,DnE2,LineNormalVectorE2);
	

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: Visualisation and analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%Collating data
X=[XE1;XE2];
Y=[YE1;YE2];

Ux=[UxE1;UxE2];
Uy=[UyE1;UyE2];

Sxx=[StressTChgE1(:,1);StressTChgE2(:,1)];
Syy=[StressTChgE1(:,2);StressTChgE2(:,2)];
Sxy=[StressTChgE1(:,3);StressTChgE2(:,3)];

%Elastic constants for each pnt
ObsE1=ones(numel(XE1),1);
ObsE2=ones(numel(XE2),1);

mu=[ObsE1*(mu(1));ObsE2*(mu(2))];
E=[ObsE1*(E(1));ObsE2*(E(2))];
nu=[ObsE1*(nu(1));ObsE2*(nu(2))];

lambda= E.*nu./((1+nu).*(1-2.*nu));   

%Stress2strain
[ Exx,Eyy,Exy ] = HookesLaw2dStress2Strain( Sxx,Syy,Sxy,E,nu );

[X,Y,Sxx,Syy,Sxy,Ux,Uy,Exx,Eyy,Exy]=ReshapeData2d...
        ( dimx,dimy,X,Y,Sxx,Syy,Sxy,Ux,Uy,Exx,Eyy,Exy);

%Calclating 2d EigenValues
[S1,S2,S1dir,S2dir]=EigCalc2d(Sxx,Syy,Sxy);

%Calclating 2d EigenValues
[E1,E2,E1dir,E2dir]=EigCalc2d(Exx,Eyy,Exy);
Dilatation2d=E1+E2;

%If the stresses do not draw well use this function
FilterValue=2.5;
[S1,S2,Sxx,Syy,Sxy] = NanOrRemovePoints( FilterValue,S1,S2,Sxx,Syy,Sxy );

%Draw data with quiver, doesn't matter if its gridded or not
figure,quiver(X(:),Y(:),Ux(:),Uy(:));
xlabel('x'); ylabel('y'); axis('equal'); title('Disp'); 

%Drawing S1 S2 directions
DrawS1S2Directions(X(:),Y(:),S1dir,S2dir,'P1',P1,'P2',P2 )


%Check if data is on a uniform or non uniform grid
[ Uniform ] = UniformGridCheck2d( X,Y );

%If data is not uniform plot as non grid data. 
if Uniform==0
    DrawScatterPlots2d( X,Y,cmap2,Sxy,Syy,Sxx,S1,S2 );
else
    %Reshaping stresses to grid dimensions
    [X,Y,Sxx,Syy,Sxy,S1,S2,Exx,Eyy,Exy,E1,E2,Ux,Uy]=ReshapeData2d...
        ( dimx,dimy,X,Y,Sxx,Syy,Sxy,S1,S2,Exx,Eyy,Exy,E1,E2,Ux,Uy );
    DrawContourFPlots2d( X,Y,cmap2,Sxy,Syy,Sxx,S1,S2 );
    Dilatation2d=(E1+E2);
    DrawDeformedGrid2d( X,Y,Ux,Uy,cmap2,Dilatation2d);
end
