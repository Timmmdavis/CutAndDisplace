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
%STEP 0: Bits you do not need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

% %===== add file paths ==========================
pathstring = pwd;                                   %Get the address of the current working directory
if ispc; parts = strsplit(mfilename('fullpath'), '\');       %Getting the address of the script working directory and splitting into cell array
else; parts = strsplit(mfilename('fullpath'), '/');  end
[~,n] = find(~cellfun(@isempty,strfind(parts,'BEM_DDM_MATLAB'))); %finding the scripts root directory and its location (n)
if ispc; addpath(genpath(strjoin(parts(1,1:n),'\'))); %adding all folders to the path that are in this dir 
else; addpath(genpath(strjoin(parts(1,1:n),'/'))); end;%mac/linux
cd(pathstring)                                       %jumping back to the current dir

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
     
%disp('This setup will give rank deficient warning as interface is closed cont with no fixing')
    

ppp=300; %boundary sampling
% %Circle (Free boundary elastic 1)
ri = 1;
k = 0;
for theta = 0:pi/ppp:2*pi %600
    k = k+1;
    x(k) = ri*cos(theta);
    y(k) = ri*sin(theta);
end
x=x*1000;
y=y*1000;
sz=numel(x);
Pointsxy=[x',y'];
mystruct.line1=(1:sz);

%Square of fixed displacement points inside elastic 1
%Fixed points inner circle (making a square)
rr=0.2;%size of the square
x=[-rr,0,rr,0,-rr];                        %List of x points
y=[0,-rr,0,rr,0];
PointsxyF=[x;y]';
PointsxyF=PointsxyF.*4000;

%Fixed displacements for E1 Flagged.
[Pointsxy,mystruct,BoundaryFlag] = DataAppender2d( Pointsxy,PointsxyF,mystruct );
%Fdisp = 1 is the free boundary points of E1 %maxbndflg=1

clear theta x y    

% %Circle (Interface Elastic 1)
r = 1;
k = 0;
for theta = 0:pi/ppp:2*pi
    k = k+1;
    x(k) = r*cos(theta);
    y(k) = r*sin(theta);
end
x=x*2000;
y=y*2000;
PointsxyIFE1E2=[x;y]';

%Adding the interface of E1-E2 in
[Pointsxy,mystruct,BoundaryFlag] = DataAppender2d( Pointsxy,PointsxyIFE1E2,mystruct,BoundaryFlag,2 );
%maxbndflg=2

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
	
mu=500; % 1%0.4;%4000 	%Shear Mod, mu or G. Relates shear stress to shear strain. 
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
mu2=100; %0.5
E2 = mu2*(2*(1+nu2)) ;

%Appending these to the back of E1
nu=[nu;nu2];
mu=[mu;mu2];
E=[E;E2];

 

SxxE1 = 0.002;  SyyE1 = 0.002;  SxyE1 = 0; 
Tn =  0.002; %This is same as Sxx and Syy with a constant value
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


% % cells=20;   %Define the number of observations points on the square grid
% % padding=15;  %How much extra bumph you want to add away from the fault sticks, reduce to 0 if having half space issues 
% % [maxgriX,mingriX,maxgriY,mingriY,sz]=MinMaxDataExtents(Points,cells,padding);
% % [X,Y] = meshgrid(mingriX:sz:maxgriX,mingriY:sz:maxgriY);
% % rowcount = size(X,1);
% % colcount = size(X,2);

X=linspace(1000,3000,21)';
Y=zeros(size(X));


    %Replacing XY of any obs points for Nan if these lie on top of the end points of a
    %dislocation, without this it can result in spurious values. 
[ X,Y  ] = NullPointsLyingOnElement( X,Y,XBEG,YBEG,XEND,YEND,NUM );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Removing values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Creating the radial coordinates for these obs points, these will be used
%to filter these. 
[TH,R] = cart2pol(X,Y);
Theta=radtodeg(TH);%figure;scatter(X(:),Y(:),20,Theta(:))

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating the stress on the observation points, elastic 1 then elastic 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

%Calculating the stress on the points that make up the annulus (elastic 1).
%0 input stress as the original analytical calc is internal pressure not
%boundary stress. its the same Tn anyway for this geometry
[StressChangeE1,TotalStressE1,StressRegE1,InducedDisplacementsE1]=StressDispOnSurroundingPoints(XE1,YE1,xeE1,yeE1,HalfLengthE1,0,0,0,nuE1,EE1,halfspace,ShearDispE1,TensileDispE1,BetaE1);

SxxE1=StressChangeE1(:,1)';SyyE1=StressChangeE1(:,2)';SxyE1=StressChangeE1(:,3)';
sz=numel(SxxE1);


%Calculating the stress on the points that make up the external matrix (elastic 2).
[StressChangeE2,TotalStressE2,StressRegE2,InducedDisplacementsE2]=StressDispOnSurroundingPoints(XE2,YE2,xeE2,yeE2,HalfLengthE2,0,0,0,nuE2,EE2,halfspace,ShearDispE2,TensileDispE2,BetaE2);

SxxE2=StressChangeE2(:,1)';SyyE2=StressChangeE2(:,2)';SxyE2=StressChangeE2(:,3)';

% %%
% %Drawing def
% %Grabbing Displacements
% UxE1=InducedDisplacementsE1(:,1);
% UyE1=InducedDisplacementsE1(:,2);
% UxE2=InducedDisplacementsE2(:,1);
% UyE2=InducedDisplacementsE2(:,2);
% %Setting up drawing of displacement
% blank=zeros(rowcount,colcount);blank(bad)=nan;
% blank(~OuterPart)=blank(~OuterPart)+XE1;
% blank(OuterPart)=blank(OuterPart)+XE2;
% Xx=blank;
% blank(~OuterPart)=blank(~OuterPart)+YE1;
% blank(OuterPart)=blank(OuterPart)+YE2;
% Yy=blank;
% blank(~OuterPart)=blank(~OuterPart)+UxE1;
% blank(OuterPart)=blank(OuterPart)+UxE2;
% Ux=blank;
% blank(~OuterPart)=blank(~OuterPart)+UyE1;
% blank(OuterPart)=blank(OuterPart)+UyE2;
% Uy=blank;
% 
% %Reshaping to whole grid dimensions
% [Xx,Yy,Ux,Uy]=ReshapeData( rowcount,colcount,Xx,Yy,Ux,Uy );
% %drawing
% Scl=0.01;
% DrawDeformedGrid2d( Xx,Yy,Ux,Uy,Scl,cmap );


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Running analytical solution, converting the numerical results 
%from cartesian to polar stress tensors and drawing comparative figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

%Now appending the two sets of cartesian stress components
SxxE1=[SxxE1,SxxE2];
SyyE1=[SyyE1,SyyE2];
SxyE1=[SxyE1,SxyE2];

%Now appending the two sets X Y points, putting in the same orientation as
%the tensors (rows)
X=[XE1;XE2]';
Y=[YE1;YE2]';

%%%%%%%%%%%%%%%%%%
%Creating radial coordinates for each observation point and then converting
%tensors from cart to polar coordinates. Eq 6.96 Pollard and Fletcher 2005
[TH,R] = cart2pol(X,Y);
[ SRR,STT,SRT ] = StressTensorTransformation2d(SxxE1,SyyE1,SxyE1,cos(TH),cos((pi/2)-TH));

% bad=abs(SRR)>10;SRR(bad)=0;%removing spurious values
% bad=abs(STT)>10;STT(bad)=0;
% bad=abs(SRT)>10;SRT(bad)=0;


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
figure;plot(rOvb,SrrNorm,'LineWidth',2,'color','blue');
hold on
plot(rOvb,SttNorm,'LineWidth',2,'color','g');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Plotting numerical results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%Octave can't handle transparent objects
% isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
% if  isOctave==1;
% scatter(rOvb_a_r_b(:),SrrNorm_a_r_b(:),12,'k');
% scatter(rOvb_a_r_b(:),SttNorm_a_r_b(:),12,'k');
% elseif isOctave==0;
% scatter(rOvb_a_r_b(:),SrrNorm_a_r_b(:),12,'k');
% scatter(rOvb_a_r_b(:),SttNorm_a_r_b(:),12,'k');
% end
% title('Analytical - lines, Srr (blue) Stt (pink), Numerical result as points'), xlabel('Distance from Hole centre, Normalised to Ring 1 Radius')
% ylabel('Srr & Stt normalised to driving stress');



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
[TH,Rin] = cart2pol(XE1,YE1);[TH,ROut] = cart2pol(XE2,YE2);clear TH
%Calculating soultion for points same as points used in numerical calc
[rOvbPnts,SrrNormPnts,SttNormPnts]=AnnulusEq_Func_pnts(muE1,nuE1,muE2,nuE2,(Rin')/1e3,(ROut')/1e3);
%Calculating polar stress residuals.
ResidualSrr=SrrNormPnts'-SrrNorm_a_r_b;
ResidualStt=SttNormPnts'-SttNorm_a_r_b;
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
