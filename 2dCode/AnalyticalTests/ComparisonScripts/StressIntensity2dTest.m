% This test checks the estimated stress intensity for two 2D, a
% line crack subject to tension then a curved crack subject to biaxial
% stress (2D plane strain problem)
%
% Problem 1 in: Shiah, Y.C. and Lin, Y.J., 2004. An infinitely large plate
% weakened by a circular-arc crack subjected to partially distributed
% loads. Journal of engineering mathematics, 49(1), pp.1-18.
%
% Problem 2 in: Pollard and Fletcher, 2005 (Fundamentals of structrual
% geology).


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
     
% %Single flat user defined fracture
a=1;
r=3;
x=linspace(-a,a,50);    
y=-(sqrt(r.^2-(x).^2));
CrkHeight=max(y)-min(y);


Pointsxy=[x;y]';
mystruct.line1=(1:(length(Pointsxy(:,1))));
Fnms=fieldnames(mystruct);
nf=numel(Fnms);clear Fnms


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
    %Option E = Defining boundary conditions as traction at the element
    %centres. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
 cc=zeros(size(MidPoint(:,1))); 
    
 Tn = 1;    Tn=cc+Tn;	    
 Ts = 0;    Ts=cc+Ts;       
 Option='E'; 


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


%% TEST.
[EndElsLoc] = GetCrackTipElements2d(P1,P2);
EndElsFlg=EndElsLoc>0;

hold on
%Draw end els
PlotFracture( P1(EndElsFlg,:),P2(EndElsFlg,:),'b' )
hold off

%Analytical formula for stress intensity at curved crack tips
[K1an,K2an] = ShiahLin2004_StrIntCurvedCrack(Tn(1),CrkHeight,a);

%%
[K1,K2] = StressIntensity2d(EndElsLoc,Dn,Ds,HalfLength,E,nu,LineNormalVector,MidPoint,P1,P2);
K2=(K2); %Ignoring sign

figure;
scatter([1,2],[K1an,K1an]);
hold on
scatter([1,2],[K2an,K2an]);
scatter([1,2],K1(EndElsFlg),15,'k','filled');
scatter([1,2],abs(K2(EndElsFlg)),15,'k','filled');
title('Stress intensity for curved line crack subjected to biaxial tension')
legend('Analytical solution K1','Analytical solution K2','Numerical result K1','Numerical result K2')
xlabel('~')
ylabel('K1/K2 at crack tips')
WhiteFigure;

hold off

Residual=[K1an;K1an;K2an;K2an]-[K1(EndElsFlg);abs(K2(EndElsFlg))];
if any(abs(Residual)>0.01)
    error('Numerical value of K1/K2 deviates far from analytical solution for curved line crack')
else
    disp('Good match to analytical solution')
end

disp('Starting Part 2')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: import the fault surface.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading XY ascii data or manually creating fractures.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
% %Single flat user defined fracture
a=1;
x=linspace(-a,a,50);  % x=fliplr(x); %Flips normals
y=zeros(1,numel(x)); 

Theta=deg2rad(20);
Beta=(0.5*pi)-Theta;
[x,y] = RotateObject2d(x,y,Theta);


Pointsxy=[x;y]';
mystruct.line1=(1:(length(Pointsxy(:,1))));
Fnms=fieldnames(mystruct);
nf=numel(Fnms);clear Fnms



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

figure;
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
Syy = 1;
Sxy = 0;
                
Option='B'; 


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

%% TEST.
[EndElsLoc] = GetCrackTipElements2d(P1,P2);
EndElsFlg=EndElsLoc>0;

hold on
%Draw end els
PlotFracture( P1(EndElsFlg,:),P2(EndElsFlg,:),'b' )
hold off

%%
[K1,K2] = StressIntensity2d(EndElsLoc,Dn,Ds,HalfLength,E,nu,LineNormalVector,MidPoint,P1,P2);
hold on
%Drawing values of stress intensity at the tips as points
scatter(MidPoint(EndElsFlg,1),MidPoint(EndElsFlg,2),15,K2(EndElsFlg,:),'filled');
hold off

%%
% %Stress intensity for a line crack loaded by tension:  
% K1an=Syy*sqrt(pi*a);

%Stress intensity for an inclined line crack loaded by tension
%Tada Stress analysis handbook of cracks, PP127, Section 5.2:  
Cons=sin(Beta)*sqrt(pi*a);
K1an=(Syy*sin(Beta))*Cons;
K2an=(Syy*cos(Beta))*Cons;
%Note this has no sign involved


figure;
scatter([1,2],[K1an,K1an]);
hold on
scatter([1,2],abs([K2an,K2an]));
scatter([1,2],K1(EndElsFlg),15,'k','filled');
scatter([1,2],abs(K2(EndElsFlg)),15,'k','filled');
title('Stress intensity for inclined line crack subjected to tension')
legend('Analytical solution K1','Analytical solution K2','Numerical result K1','Numerical result K2')
xlabel('~')
ylabel('K1/K2 at crack tips')
WhiteFigure;
%ylim([0 max(K1an)+max(K1an)/4])
hold off


Residual=[K1an;K1an;K2an;K2an]-[K1(EndElsFlg);abs(K2(EndElsFlg))];

if any(abs(Residual)>0.1)
    error('Numerical value of K1/K2 deviates far from analytical solution for curved line crack')
else
    disp('Good match to analytical solution')
end

