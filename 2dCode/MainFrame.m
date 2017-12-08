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
%STEP 1: import the fault surface.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading XY ascii data or manually creating fractures.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
% % % %Single flat user defined fracture
% x=linspace(-0.5,0.5,15);    
% y=zeros(1,numel(x)); 
% Pointsxy=[x;y]';
% mystruct.line1=(1:(length(Pointsxy(:,1))));
% Fnms=fieldnames(mystruct);
% nf=numel(Fnms);clear Fnms

% % %Circle
% rad=1;
% [ xe,ye ] = CreateCircleXY( 600,rad );
% sz=numel(xe);
% Pointsxy=[ye,xe];
% mystruct.line1=(1:sz);
% %Square of points inside circle
% rr=0.8;
% x=[-rr,0,rr,0,-rr];y=[0,-rr,0,rr,0];
% PointsxyF=[x;y]';
% [Pointsxy,mystruct,Fdisp]=DataAppender2d(Pointsxy,PointsxyF,mystruct);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Loading Shp
    %Loads Points that represent the file. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Pointsxy,mystruct]=m_shaperead('FracturesStoneHaven');
Fnms=fieldnames(mystruct);
nf=numel(Fnms);clear Fnms
Pointsxy=[Pointsxy(:,1),Pointsxy(:,3)]; 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Appending and flagging the fixed data points.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[Pointsxy,mystruct,nf,Fdisp]=DataAppender2d( Pointsxy,PointsxyF,mystruct);


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
Pointsxy(:,2)=Pointsxy(:,2)-freesurface_height;


% Defining elastic constants
mu=1;           	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
[ ~,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu );

% Reshaping the list of XY into usable variables
[MidPoint,HalfLength,P1,P2,LineNormalVector,Lines,Points,Fdisp  ]...
 = CreateElements2d( Pointsxy,mystruct,Fdisp );

% Drawing fracture
%Normal els as red.
PlotFracture( P1,P2,'r' )
%Fixed elements shown as blue
id=find(Fdisp); 
PlotFracture(P1(id,:),P2(id,:),'b' )
%Titles
title('fractures'), xlabel('x'), ylabel('y');axis equal
%Drawing Normals
quiver(MidPoint(:,1),MidPoint(:,2),LineNormalVector(:,1),LineNormalVector(:,2))



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, remote (with
%friction), remote with no opening components, tractions on the elements
%(internal pressure) etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[Sxx,Syy,Sxy,Tn,Ts,Mu,Sf,strain ] = CreateBlankVars;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Constant Slip across the fault surface.
	%Runs a constant slip on every element. Output 'stresses' are blank
	%arrays of 0's as this option uses displacement boundary conditions.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Option='A';     
% cc=zeros(NUM,1); 
% % Positive = left lateral displacement
% ShearDisp  = 0;       ShearDisp	  = cc+(ShearDisp); 
% % Positive = opening displacement
% TensileDisp =1;       TensileDisp = cc+(TensileDisp); 
% Sxx = 0;  Syy = 0;  Sxy = 0;   	
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. 
    %Choose to define in stress or strain, the stresses can be defined as
    %single values that are repeated for every element or defined as a
    %vector that is the stress at every element: i.e.  
    %Syy = linspace(0,1,NUM); %(linearly increasing for each element).
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
Sxy = 1;                   
Option='B'; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option C = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. This option includes frictional contact
    %properties on the fault surface, elements cannot interpenetrate and
    %slip is reduced by the frictional parameters. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	%%%%%%%%%%%%%%
    %StrainInput       
	%%%%%%%%%%%%%%

% %Put to 1 to define the stresses defined in 'stress input' as strain values    
% strain=0;    
    
	%%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%
	
%  Sxx = 0; 				
%  Syy = 0;
%  Sxy = 1;                  
% 
% % Frictional parameters, define as single value for all elements or vary
% % based on element number
%  Mu  = 0.6;       Mu=repmat(Mu,numel(xe),1);  %Coefficient of friction
%  Sf  = 0;         Sf=repmat(Sf,numel(xe),1);  %Frictional strength
%  
% Option='C'; 


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option D = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. No opening components, this option solves
    %the boundary value problem through shearing of elements alone. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%
    %StrainInput       
	%%%%%%%%%%%%%%

% %Put to 1 to define the stresses defined in 'stress input' as strain values    
% strain=0;    

	%%%%%%%%%%%%%%
    %StressInput
	%%%%%%%%%%%%%%
    
% Sxx = 0;         		
% Syy = 0; 
% Sxy = 1;                   
% Option='D';


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option E = Defining boundary conditions as traction at the element
    %centres. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%  cc=zeros(NUM,1); 
%     
%  Tn = 0;    Tn=cc+Tn;	    
%  Ts = 0.5;  Ts=cc+Ts;       
%  Option='E'; 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(Option,'A')
    [Ds,Dn,Sxx,Syy,Sxy]=SlipCalculator2d(MidPoint,HalfLength,Sxx,Syy,Sxy...
        ,Tn,Ts,nu,E,halfspace,LineNormalVector,strain,Fdisp,Mu,Sf,Option);
end

    % Removing any fixed elements
if any(Fdisp)==1  
    [MidPoint,HalfLength,P1,P2,LineNormalVector]...
    = RemovingFixedEls2d(Fdisp,P1,P2);
end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing figures of slip distribution on the crack.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Plot displacement discontinuity magnitude for all elements
GraphSlipVsElNum( Dn,Ds )

% Creating directions and magnitudes of slip for plotting
PlotOpeningVsShearOnEls( LineNormalVector,Dn,Ds,P1,P2,MidPoint,HalfLength )

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Animate Fault movement.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %To do this look at file Animate2d.m

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XY observation points to calculate stress strain
%and displacement on. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Regular/randomly spaced grid of points bounding the faults.
    %Define how far it extends away and its density.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define the number of cells on the square grid
cells=16;   
%How many extra cells you want to add away from the faults, 
%reduce to 0 if having half space issues 
padding=10;  
% Create a vector that is just all the points in X Y. 
PointsXY=[[P1(:,1);P2(:,1)],[P1(:,2);P2(:,2)]];
[maxgriX,mingriX,maxgriY,mingriY,sz]=MinMaxDataExtents2d(PointsXY,cells,padding);

%Uniform grid spacing
[X,Y] = meshgrid(mingriX:sz:maxgriX,mingriY:sz:maxgriY);
[dimx,dimy] = size(X);  

% %%Random points within the calculated bounds
% NumPnts=15000; %number of points, do not have cells 
% lengthx=maxgriX-mingriX;
% lengthy=maxgriY-mingriY;
% xmv=(maxgriX+mingriX)/2; ymv=(maxgriY+mingriY)/2;
% x=rand(1,NumPnts)*(lengthx*2);
% y=rand(1,NumPnts)*(lengthy*2);
% X=x-lengthx+xmv;
% Y=y-lengthy+ymv;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = User defined XY grid.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Creating grid with user defined sampling
% X = linspace(-10,15,25); 
% Y = linspace(-10,15,25); 
% [X,Y] = meshgrid(X,Y); 
% [dimx,dimy] = size(X);  

% %Creating a flat line
%  X = linspace(0,15,100); 
%  Y = zeros(size(X)); 
 
% % Random data at on n*n grid centred at xmv,ymv
% n=4;
% xmv=0; ymv=0;
% x=rand(1,10201)*n;
% y=rand(1,10201)*n;
% X=x-n+xmv;
% Y=y-n+ymv;

% %Observation points at midpoints
% X=x;
% Y=y;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option C = Import XY from a .text file that is on the loaded paths. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Loads the dat file, This contains 2 columns XY
% S = load('Points.dat');	
% X = S(:,1);                   
% Y = S(:,2);                         


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing and fixing Obs Point data just created.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Drawing these and the location compared to the Obs Points
figure;
PlotFracture( P1,P2,'r' )
title('Fractures and Observation points');
xlabel('x'), ylabel('y');axis equal
hold on
scatter(X(:),Y(:),'.')
hold off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses/Disps at defined observation points in XY.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
[StressTTotal,StrainTTotal,StressTChg,StrainTChg,StressTReg,StrainTReg]...
= CalculateStressesOnSurroundingPoints2d(X,Y,MidPoint,HalfLength...
,Sxx,Syy,Sxy,nu,E,halfspace,Ds,Dn,LineNormalVector);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Calculate Displacements.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ux,Uy] = CalculateDisplacementOnSurroundingPoints2d...
(X,Y,MidPoint,HalfLength,nu,E,halfspace,Ds,Dn,LineNormalVector);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 6: Finite strain computation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% [exx,eyy,exy,ExxInfErrorPerc,EyyInfErrorPerc,ExyInfErrorPerc]...
% =FiniteStrainLagrangian2d(Ux,Uy,X,Y);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: Visualisation and analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

  
%Remote stress
Pxx=Sxx(1,1);
Pyy=Syy(1,1);
Pxy=Sxy(1,1);

%Extracting stress
[ Sxx,Syy,Sxy ] = ExtractCols( StressTTotal );

%Stress2strain
[ Exx,Eyy,Exy ] = HookesLaw2dStress2Strain( Sxx,Syy,Sxy,E,nu );

%Calclating 2d EigenValues
[S1,S2,S1dir,S2dir]=EigCalc2d(Sxx,Syy,Sxy);

%Calclating 2d EigenValues
[E1,E2,E1dir,E2dir]=EigCalc2d(Exx,Eyy,Exy);
Dilatation2d=E1+E2;

%If the stresses don't draw well use this function
FilterValue=2.5;
Bad=abs(Sxx)>FilterValue;
[S1,S2,Sxx,Syy,Sxy] = NanOrRemovePoints( Bad,S1,S2,Sxx,Syy,Sxy );


%Draw data with quiver, doesn't matter if its gridded or not
figure,quiver(X(:),Y(:),Ux,Uy);WhiteFigure;
xlabel('x'); ylabel('y'); axis('equal'); title('Disp'); 

%Drawing S1 S2 directions.
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
    Scl=mu/3;
    Dilatation2d=(E1+E2);
    DrawDeformedGrid2d( X,Y,Ux,Uy,cmap2,Dilatation2d,'Scale',Scl);
end
