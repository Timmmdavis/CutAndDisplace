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


%%
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
     
% %Circle (inner ann)
rad=1000;
Smpling=200;
[ xe,ye ] = CreateCircleXY( Smpling,rad );
sz=numel(xe);
Pointsxy=[ye,xe];
mystruct.line1=(1:sz);

%Square of fixed displacement points inside elastic 1
%Fixed points inner circle (making a square)
rr=200;%size of the square
x=[-rr,0,rr,0,-rr];                        %List of x points
y=[0,-rr,0,rr,0];
PointsxyF=[x;y]';

%Fixed displacements for E1 Flagged.
[Pointsxy,mystruct,BoundaryFlag] = DataAppender2d( Pointsxy,PointsxyF,mystruct );
%Fdisp = 1 is the free boundary points of E1 %maxbndflg=1

% %Circle (Interface Elastic 1)
rad=2000;
[ xe,ye ] = CreateCircleXY( Smpling,rad );
PointsxyIFE1E2=[ye,xe];

%Adding the interface of E1-E2 in
[Pointsxy,mystruct,BoundaryFlag] = DataAppender2d( Pointsxy,PointsxyIFE1E2,mystruct,BoundaryFlag,2 );

%Free boundary elastic 1 is now 0's the interface is 2 in out var called
%BoundaryFlag
PointsxyIFE2E1=PointsxyIFE1E2;

%Adding the interface of E2-E1
[Pointsxy,mystruct,BoundaryFlag] = DataAppender2d( Pointsxy,PointsxyIFE2E1,mystruct,BoundaryFlag,3 );
%maxbndflg=3
names = fieldnames(mystruct);
numstructs=numel(names);
FlipNormalsFlag=zeros(numstructs,1);
FlipNormalsFlag(end,1)=1; %Flag that tells the line builder to flip the normals of the interface


%Fixing disps for 2nd elastic (in hole too)
[Pointsxy,mystruct,BoundaryFlag] = DataAppender2d( Pointsxy,PointsxyF,mystruct,BoundaryFlag,5 );
%5 is the fixed els in elastic 2
FlipNormalsFlag=[FlipNormalsFlag;0]; %not going to flip the normals on this one



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
mu=500;           	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
[ ~,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu );

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
[SxxE1,SyyE1,SxyE1,Tn,Ts,Mu,Sf,strain ] = CreateBlankVars;

%Elastic 2 properties
nu2=0.25;
mu2=100; %0.5
E2 = mu2*(2*(1+nu2)) ;

%Appending these to the back of E1
nu=[nu;nu2];
mu=[mu;mu2];
E=[E;E2];

SxxE1 = 0.002;  SyyE1 = 0.002;  SxyE1 = 0; 
Tn = 0.002;% 0.002; %This is same as Sxx and Syy with a constant value
Option='F'; %Inhomogeneous Elastic


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ds,Dn,Sxx,Syy,Sxy]=SlipCalculator2d(...
    MidPoint,HalfLength,SxxE1,SyyE1,SxyE1,Tn,Ts,nu,E,halfspace,LineNormalVector,strain,BoundaryFlag,Mu,Sf,Option);

%Removing fixed disps from BoundaryFlag
BoundaryFlag(Fdisp)=[]; %removing

%%
% Removing any fixed elements
if any(Fdisp)==1  
    [MidPoint,HalfLength,P1,P2,LineNormalVector]...
    = RemovingFixedEls2d(Fdisp,P1,P2);
end

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
Part=2; %Which elastic we want to extract[xeE2,yeE2,aE2,nuE2,EE2,LineNormalVectorE2,NUME2,FdispE2,FB_E2,IF_E2]...
[MidPointE2,HalfLengthE2,nuE2,EE2,LineNormalVectorE2,NUME2,FdispE2,FB_E2,IF_E2]...
 = ExtractElasticParts2d( Part,BoundaryFlag,MidPoint,HalfLength,nu,E,LineNormalVector );
muE2 = EE2/(2*(1+nuE2));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XY observation points to calculate stress strain
%and displacement on. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% cells=20;   %Define the number of observations points on the square grid
% padding=15;  %How much extra bumph you want to add away from the fault sticks, reduce to 0 if having half space issues 
% [maxgriX,mingriX,maxgriY,mingriY,sz]=MinMaxDataExtents2d(Points,cells,padding);
% [X,Y] = meshgrid(mingriX:sz:maxgriX,mingriY:sz:maxgriY);
% rowcount = size(X,1);
% colcount = size(X,2);

%Line of points
X=linspace(1000,3000,21)';
Y=zeros(size(X));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Removing values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Creating the radial coordinates for these obs points, these will be used
%to filter these. 
[theta,R] = cart2pol(X,Y);
Theta=radtodeg(theta);%figure;scatter(X(:),Y(:),20,Theta(:))

%Removing points in hole (Holewidth)
%Removing points near elements where BEM solution breaks down (tol width)
%Removing points really far away from the hole so the graph is nicer
Holewidth=1000; AnnWidth=2000; Toldis=10;
bad= (R<Holewidth+Toldis | R<AnnWidth+Toldis & R>AnnWidth-Toldis | R>Holewidth+AnnWidth)  ;
X(bad)=nan;Y(bad)=nan;Theta(bad)=nan;R(bad)=nan;
disp('Any points closer than 10m from dislocations removed') %toldis

%Now only selecting a certain portion of points (Makes the figure neater)
%Good=Theta<135 & Theta>45; %90 deg rad
% Good=Theta<2 & Theta>-2;   %4 deg rad
% X=X(Good);Y=Y(Good);R=R(Good);

%Flagging points in 2nd elastic
OuterPart=R>2000;
XE2=X(OuterPart);YE2=Y(OuterPart);
XE1=X(~OuterPart);YE1=Y(~OuterPart);


%Drawing the observation points on one of the previous figures
hold on

%Octave can't handle transparent objects
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
    scatter(X(:),Y(:),'b');
    scatter(XE2(:),YE2(:),'r');
elseif isOctave==0
    scatter(X(:),Y(:),'b','filled','MarkerFaceAlpha',1/8);
    scatter(XE2(:),YE2(:),'r','filled','MarkerFaceAlpha',1/8);
end
title('Observation Points, Blue E1 Red E2 '), xlabel('x'), ylabel('y')
         



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating the stress on the observation points, elastic 1 then elastic 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

%Calculating the stress on the points that make up the annulus (elastic 1).
%0 input stress as the original analytical calc is internal pressure not
%boundary stress. its the same Tn anyway for this geometry
[StressTTotalE1,StrainTTotalE1,StressTChgE1,StrainTChgE1,StressTRegE1,StrainTRegE1]...
    = CalculateStressesOnSurroundingPoints2d(XE1,YE1,MidPointE1,HalfLengthE1,0,0,0,nuE1,EE1,halfspace,DsE1,DnE1,LineNormalVectorE1);

[UxE1,UyE1] = CalculateDisplacementOnSurroundingPoints2d(XE1,YE1,MidPointE1,HalfLengthE1,nuE1,EE1,halfspace,DsE1,DnE1,LineNormalVectorE1);
	

%Calculating the stress on the points that make up the external matrix (elastic 2).
[StressTTotalE2,StrainTTotalE2,StressTChgE2,StrainTChgE2,StressTRegE2,StrainTRegE2]...
    = CalculateStressesOnSurroundingPoints2d(XE2,YE2,MidPointE2,HalfLengthE2,0,0,0,nuE2,EE2,halfspace,DsE2,DnE2,LineNormalVectorE2);

[UxE2,UyE2] = CalculateDisplacementOnSurroundingPoints2d(XE2,YE2,MidPointE2,HalfLengthE2,nuE2,EE2,halfspace,DsE2,DnE2,LineNormalVectorE2);
	
%Collating data
SxxE1=StressTChgE1(:,1)';
SyyE1=StressTChgE1(:,2)';
SxyE1=StressTChgE1(:,3)';
sz=numel(SxxE1);
SxxE2=StressTChgE2(:,1)';
SyyE2=StressTChgE2(:,2)';
SxyE2=StressTChgE2(:,3)';
%Now appending the two sets of cartesian stress components
Sxx=[SxxE1,SxxE2];
Syy=[SyyE1,SyyE2];
Sxy=[SxyE1,SxyE2];
X=[XE1;XE2];
Y=[YE1;YE2];
Ux=[UxE1;UxE2];
Uy=[UyE1;UyE2];


quiver(X,Y,Ux,Uy);
hold off
test=1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Running analytical solution, converting the numerical results 
%from cartesian to polar stress tensors and drawing comparative figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 



%Now appending the two sets X Y points, putting in the same orientation as
%the tensors (rows)
X=[XE1;XE2]';
Y=[YE1;YE2]';

%%%%%%%%%%%%%%%%%%
%Creating radial coordinates for each observation point and then converting
%tensors from cart to polar coordinates. Eq 6.96 Pollard and Fletcher 2005
[theta,R] = cart2pol(X,-Y);
CosAx=cos(theta);
CosAy=sin(theta);
[ SRR,STT,SRT ] = StressTensorTransformation2d(Sxx,Syy,Sxy,CosAx,CosAy);

% bad=abs(SRR)>10;SRR(bad)=0;%removing spurious values
% bad=abs(STT)>10;STT(bad)=0;
% bad=abs(SRT)>10;SRT(bad)=0;


%%%%
%Normalizing
%%%%
SrrNorm_a_r_b=SRR/Tn(1,1);
SttNorm_a_r_b=STT/Tn(1,1);
rOvb_a_r_b=R/2000;

%Removing nans
SrrNorm_a_r_b(isnan(SrrNorm_a_r_b))=[];
SttNorm_a_r_b(isnan(SttNorm_a_r_b))=[];
rOvb_a_r_b(isnan(rOvb_a_r_b))=[];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating and Plotting analytical results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
a=1;  
b=2;       
r= [linspace(a,b,50),linspace(b,3,50)];
[rOvb,SrrNorm,SttNorm]=CrouchStar1983_InhomogeneousAnnulus(muE1,nuE1,muE2,a,b,r);
figure;plot(rOvb,SrrNorm,'LineWidth',2,'color','blue');
hold on
plot(rOvb,SttNorm,'LineWidth',2,'color','g');


scatter(rOvb_a_r_b(:),SrrNorm_a_r_b(:),24,'k','filled');
scatter(rOvb_a_r_b(:),SttNorm_a_r_b(:),24,'k','filled');
title('Stress across elastic interface'), 
xlabel('Distance, interface at 1')
ylabel('Stress components')
grid on
legend('show')
legend('Stt','Srr','Numerical results')
plot([1 1],[-1 1.5],'k--') %interface
titlesz=25;
fntsz=21;
ChangeFontSizes(fntsz,titlesz);



%Calculating analytical solution stresses for points used in numerical calc
[~,Rin] = cart2pol(XE1,YE1);[~,ROut] = cart2pol(XE2,YE2);clear TH
%Calculating soultion for points same as points used in numerical calc     
r= [(Rin')/1e3,(ROut')/1e3];
[rOvbPnts,SrrNormPnts,SttNormPnts]=CrouchStar1983_InhomogeneousAnnulus(muE1,nuE1,muE2,a,b,r);
%Calculating polar stress residuals.
ResidualSrr=SrrNormPnts-SrrNorm_a_r_b;
ResidualStt=SttNormPnts-SttNorm_a_r_b;
PercentErrorStressrrInterpPnts=((100./SrrNormPnts').*SrrNorm_a_r_b)-100;
PercentErrorStressttInterpPnts=((100./SttNormPnts').*SttNorm_a_r_b)-100;
% plot(rOvbPnts,ResidualSrr,'--')
% plot(rOvbPnts,ResidualStt,'--')

%disp(max(ResidualSrr));
fprintf('MaxResidualSrr %i.\n',max(abs(ResidualSrr(:))))
%disp(max(ResidualStt));
fprintf('MaxResidualStt %i.\n',max(abs(ResidualStt(:))))
biggest_res=max(ResidualStt)+max(ResidualSrr);

if  abs(biggest_res) > 0.015
    error('Not a good match to analytical solution')
else
    disp('Everything looks good, tolerance checks max polar stress residuals and flags errors above 0.01')
end
