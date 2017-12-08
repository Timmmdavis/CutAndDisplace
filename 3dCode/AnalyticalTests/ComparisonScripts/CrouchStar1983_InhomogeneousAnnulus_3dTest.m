%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you do not need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%STEP 1: Import the fault surface.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading Gocad ascii data.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  string='Tube_FixedInnnerTris.ts'; %Loading inner annlus, free bnd elastic 1
 [ PointsFBE1F,TrianglesFBE1F ] = GoCadAsciiReader( string );
 FDsz=numel(PointsFBE1F(:,1));
    
        %VeryHighSampleSetup
%  string='Tube_1000mRad_4000FacesONVHS.ts'; %Loading inner annlus, free bnd elastic 1
% [ PointsFBE1,TrianglesFBE1 ] = GoCadAsciiReader( string );		
%   string='Tube_2000mRad_4000FacesONVHS.ts'; %Loading interface elastic 1
% [ PointsIFE1E2,TrianglesIFE1E2 ] = GoCadAsciiReader( string );		
%   string='Tube_2000mRad_4000FacesINVHS.ts'; %Loading interface elastic 2
% [ PointsIFE2E1,TrianglesIFE2E1 ] = GoCadAsciiReader( string );		        
        
%  string='Tube_1000mRad_1500FacesONVHS.ts'; %Loading inner annlus, free bnd elastic 1
% [ PointsFBE1,TrianglesFBE1 ] = GoCadAsciiReader( string );		
%   string='Tube_2000mRad_1500FacesONVHS.ts'; %Loading interface elastic 1
% [ PointsIFE1E2,TrianglesIFE1E2 ]= GoCadAsciiReader( string );		
%   string='Tube_2000mRad_1500FacesINVHS.ts'; %Loading interface elastic 2
% [ PointsIFE2E1,TrianglesIFE2E1 ] = GoCadAsciiReader( string );		


%  %%%        HighSampleSetup 
 string='Tube_1000mRad_1500FacesONHS.ts'; %Loading inner annlus, free bnd elastic 1
 [ PointsFBE1,TrianglesFBE1 ] = GoCadAsciiReader( string );		
  string='Tube_2000mRad_1500FacesONHS.ts'; %Loading interface elastic 1
[ PointsIFE1E2,TrianglesIFE1E2 ]  = GoCadAsciiReader( string );		
  string='Tube_2000mRad_1500FacesINHS.ts'; %Loading interface elastic 2
[ PointsIFE2E1,TrianglesIFE2E1 ] = GoCadAsciiReader( string );	

 
 
% %      %MediumSampleSetup
%   string='Tube_1000mRad_1500FacesONMS.ts'; %Loading inner annlus, free bnd elastic 1
%  [ PointsFBE1,TrianglesFBE1 ] = GoCadAsciiReader( string );		
%    string='Tube_2000mRad_1500FacesONMS.ts'; %Loading interface elastic 1
%  [ PointsIFE1E2,TrianglesIFE1E2 ] = GoCadAsciiReader( string );		
%    string='Tube_2000mRad_1500FacesINMS.ts'; %Loading interface elastic 2
%  [ PointsIFE2E1,TrianglesIFE2E1 ] = GoCadAsciiReader( string );


% % %    %LowSampleSetup
%  string='Tube_1000mRad_1500FacesONLS.ts'; %Loading inner annlus, free bnd elastic 1
% [ PointsFBE1,TrianglesFBE1 ] = GoCadAsciiReader( string );		
%   string='Tube_2000mRad_1500FacesONLS.ts'; %Loading interface elastic 1
% [ PointsIFE1E2,TrianglesIFE1E2 ] = GoCadAsciiReader( string );		
%   string='Tube_2000mRad_1500FacesINLS.ts'; %Loading interface elastic 2
% [ PointsIFE2E1,TrianglesIFE2E1 ] = GoCadAsciiReader( string );	
 	
%Adding Fixed parts of E1 to E1 inner ann
[Points,Triangles,BoundaryFlag] = DataAppender3d( PointsFBE1,PointsFBE1F,TrianglesFBE1,TrianglesFBE1F );

%Adding interface E1E2, Points towards E2 
[Points,Triangles,BoundaryFlag] = DataAppender3d( Points,PointsIFE1E2,Triangles,TrianglesIFE1E2,BoundaryFlag,2 );

%Adding interface E2E1, Points towards E1
[Points,Triangles,BoundaryFlag] = DataAppender3d( Points,PointsIFE2E1,Triangles,TrianglesIFE2E1,BoundaryFlag,3 );

% %How to add Fixed parts of E2 to data
%[Points,Triangles,BoundaryFlag] = DataAppender3d( Points,PointsFBE1F,Triangles,TrianglesFBE1F,BoundaryFlag,5 );

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

% Defining elastic constants
mu=500;           	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
[ ~,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu );

% Getting surface ready for script. 
Points=[Points(:,1),Points(:,2),Points(:,3),(Points(:,4)-freesurface_height)];
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 

title('All Surfaces Loaded');

%Locked Els
Fdisp=(BoundaryFlag == 1 | BoundaryFlag == 5);



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

[Sxx,Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,Mu,Sf,strain ] = CreateBlankVars;
	%%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%

%Elastic 2 properties
nu2=0.25;
mu2=100; %0.5
lambda2=(2*mu2*nu2)/(1-(2*nu2));
%Appending these to the back of E1
nu=[nu;nu2];
mu=[mu;mu2];
lambda=[lambda;lambda2];
  
    
Sxx = 0.002;       			
Syy = 0.002; 
Szz = 0; 
Sxy = 0;        			
Sxz = 0;
Syz = 0; 
Option='F';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
[Dss,Dds,Dn,Sxx,Syy,Szz,Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,...
 Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
 FaceNormalVector,BoundaryFlag,strain,Mu,Sf,Option,Triangles,Points);


%Removing fixed disps from BoundaryFlag
BoundaryFlag(Fdisp)=[]; %removing

    % Removing any fixed elements
if any(Fdisp)==1
    [Triangles,FaceNormalVector,MidPoint,P1,P2,P3]...
    = RemovingFixedEls3d(Triangles,FaceNormalVector,MidPoint,P1,P2,P3,Fdisp);
end
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing figures of slip distribution on the crack.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,Dss,Dds,Dn)

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
    
% clear size
% [X,Y] = meshgrid(-3000:50:3000,-3000:50:3000);
% % % % [X,Y] = meshgrid(1200:50:3500,-3581:50:-1500);
% dimx = length(X(:,1));
% dimy = length(X(1,:));
% Y=Y-freesurface_height;

% a=1000;%Hole edge
% b=2000;%Interface
% c=3000;%Limit of data
% X= [linspace(a,b,50),linspace(b,c,50)]';

X=linspace(1000,3000,21)';
Y=zeros(size(X));

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

%Removing points in hole (Holewidth)
%Removing points near elements where BEM solution breaks down (tol width)
%Removing points really far away from the hole so the graph is nicer
Holewidth=1000; AnnWidth=2000; Toldis=100;
bad= (R<Holewidth+Toldis | R<AnnWidth+Toldis & R>AnnWidth-Toldis | R>Holewidth+AnnWidth)  ;
X(bad)=nan;Y(bad)=nan;Theta(bad)=nan;R(bad)=nan;
disp('Any points closer than 100m from dislocations removed') %toldis

%Flagging points in 2nd elastic
OuterPart=R>2000;
XE2=X(OuterPart);YE2=Y(OuterPart);
XE1=X(~OuterPart);YE1=Y(~OuterPart);
ZE1=zeros(size(YE1));
ZE2=zeros(size(YE2));

%Drawing the observation points on one of the previous figures
%Octave can't handle transparent objects
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
scatter(XE1(:),YE1(:),'b');
scatter(XE2(:),YE2(:),'r');
elseif isOctave==0
figure;scatter(XE1(:),YE1(:),'b','filled','MarkerFaceAlpha',1/8);hold on
scatter(XE2(:),YE2(:),'r','filled','MarkerFaceAlpha',1/8);hold off
end
title('Observation Points, Blue E1 Red E2'), xlabel('x'), ylabel('y')
        

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses/Disps at defined observation points in XYZ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

	    
%E1 Elastic 1
E1Bits =(BoundaryFlag == 0 | BoundaryFlag == 2);
DnE1=Dn(E1Bits );
DssE1=Dss(E1Bits );
DdsE1=Dds(E1Bits );
Part=1; %Which elastic we want to extract
[MidPointE1,P1E1,P2E1,P3E1,muE1,lambdaE1,FaceNormalVectorE1,nuE1,NUME1,FdispE1,FB_E1,IF_E1]...
= ExtractElasticParts3d( Part,BoundaryFlag,MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,Fdisp );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	
%Calculating the stress on the points that make up the annulus (elastic 1).
%0 input stress as the original analytical calc is internal pressure not
%boundary stress. its the same Tn anyway for this geometry
 [StressTTotalE1,StrainTTotalE1,StressTChgE1,StrainTChgE1,StressTRegE1,StrainTRegE1]...
   =CalculateStressOnSurroundingPoints3d(DssE1,DdsE1,DnE1,muE1,lambdaE1,XE1,YE1,ZE1...
   ,0,0,0,0,0,0,P1E1,P2E1,P3E1,halfspace,nuE1);
 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Calculate Displacements.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
[UxE1,UyE1,UzE1] = CalculateDisplacementOnSurroundingPoints3d...
(DssE1,DdsE1,DnE1,nuE1,XE1,YE1,ZE1,P1E1,P2E1,P3E1,halfspace);


%Extract stresses
SxxE1=StressTTotalE1(:,1)';SyyE1=StressTTotalE1(:,2)';SxyE1=StressTTotalE1(:,4)';
sz=numel(Sxx);
 
 
 %E2 Elastic 2
E2Bits =(BoundaryFlag == 3 | BoundaryFlag == 4);
DnE2=Dn(E2Bits );
DssE2=Dss(E2Bits );
DdsE2=Dds(E2Bits );
Part=2; %Which elastic we want to extract
[MidPointE2,P1E2,P2E2,P3E2,muE2,lambdaE2,FaceNormalVectorE2,nuE2,NUME2,FdispE2,FB_E2,IF_E2]...
= ExtractElasticParts3d( Part,BoundaryFlag,MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,Fdisp );
    
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%Calculating the stress on the points that make up the annulus (elastic 1).
%0 input stress as the original analytical calc is internal pressure not
%boundary stress. its the same Tn anyway for this geometry
[StressTTotalE2,StrainTTotalE2,StressTChgE2,StrainTChgE2,StressTRegE2,StrainTRegE2]...
     =CalculateStressOnSurroundingPoints3d(DssE2,DdsE2,DnE2,muE2,lambdaE2,XE2,YE2,ZE2...
	 ,0,0,0,0,0,0,P1E2,P2E2,P3E2,halfspace,nuE2);

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Calculate Displacements.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	 
	 
[UxE2,UyE2,UzE2] = CalculateDisplacementOnSurroundingPoints3d...
(DssE2,DdsE2,DnE2,nuE2,XE2,YE2,ZE2,P1E2,P2E2,P3E2,halfspace);

 %Extract stresses
SxxE2=StressTTotalE2(:,1)';SyyE2=StressTTotalE2(:,2)';SxyE2=StressTTotalE2(:,4)';
 
%Now appending the two sets of cartesian stress components
Sxx=[SxxE1,SxxE2];
Syy=[SyyE1,SyyE2];
Sxy=[SxyE1,SxyE2];

Ux=[UxE1;UxE2];
Uy=[UyE1;UyE2];

%Now appending the two sets X Y points, putting in the same orientation as
%the tensors (rows)

X=[XE1;XE2]';
Y=[YE1;YE2]';

hold on
quiver(X,Y,Ux',Uy'); axis('equal')
hold off

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating and Plotting analytical results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
a=1;  
b=2;       
r= [linspace(a,b,50),linspace(b,3,50)];
[rOvb,SrrNorm,SttNorm]=CrouchStar1983_InhomogeneousAnnulus(muE1,nuE1,muE2,a,b,r);
%Plotting analytical lines
figure;plot(rOvb,SrrNorm,'LineWidth',2,'color','blue');
hold on
plot(rOvb,SttNorm,'LineWidth',2,'color','g');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting numerical results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scatter(rOvb_a_r_b(:),SrrNorm_a_r_b(:),24,'k','filled');
scatter(rOvb_a_r_b(:),SttNorm_a_r_b(:),24,'k','filled');
%Cleaning the figure up:
title('Stress across elastic interface'), 
xlabel('Distance, interface at 1')
ylabel('Stress components')
grid on
plot([1 1],[-1 1.5],'k--') %interface shows as a dashed line
legend('show')
legend('Stt','Srr','Numerical results')
titlesz=25;
fntsz=21;
ChangeFontSizes(fntsz,titlesz);

%Calculating analytical solution stresses for points used in numerical calc
[~,Rin] = cart2pol(XE1,YE1);[~,ROut] = cart2pol(XE2,YE2);clear TH
%Calculating soultion for points same as points used in numerical calc
r= [(Rin'),(ROut')];
a=1000; %Same sizes as 3d sol:
b=2000;
[rOvbPnts,SrrNormPnts,SttNormPnts]=CrouchStar1983_InhomogeneousAnnulus(muE1,nuE1,muE2,a,b,r);

%Removing nans (these are removed in analytical func
SrrNorm_a_r_b(isnan(SrrNorm_a_r_b))=[];
SttNorm_a_r_b(isnan(SttNorm_a_r_b))=[];
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

