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
%STEP 1: Import the fault surface.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading Gocad ascii data.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 string='SphereFix.ts';
 [ PointsF,TrianglesF ] = GoCadAsciiReader( string );    
    
%string='HalfShell_750Faces.ts'; %Only half the surface
%string='SphereUniformDistributionON_46Faces.ts'; %~46 pnts
%string='SphereUniformDistributionON2_500Faces.ts';  %375 pnts
%string='SphereUniformDistributionON_96Faces.ts';
string='SphereUniformDistributionHS_ON_1500Faces.ts';  %1500 pnts
%string='SphereUniformDistribution5120FacesSide.ts';  %1500 pnts
%string='Sphere.ts';  %1500 pnts

[ Points,Triangles ] = GoCadAsciiReader( string );


 %ONLY FOR FIXING DISP
 %Getting some sizes
 IF_Sz=numel(TrianglesF(:,1));
 Sz1=numel(Points(:,1));    
 %Now adding the fixing interior points
 PointsF(:,1)= PointsF(:,1)+Sz1;
 TrianglesF=TrianglesF+Sz1; 
 Points=[Points;PointsF];
 Triangles=[Triangles;TrianglesF];

 %the half length of each axis of the ellipsoid
 a=1; 
 b=1;
 c=2;
 Points=[Points(:,1),Points(:,2)*a,Points(:,3)*b,(Points(:,4)*c)];


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
mu=1;           	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
[ ~,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu );

% Getting surface ready for script. 
Points=[Points(:,1),Points(:,2),Points(:,3),(Points(:,4)-freesurface_height)];
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 

%Locked triangles
Fdisp= zeros(size(Triangles(:,1))); 
Fdisp(end-IF_Sz+1:end,:)=Fdisp(end-IF_Sz+1:end,:)+1;


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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. Choose to define in stress or strain. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	%%%%%%%%%%%%%%
    %StrainInput       
	%%%%%%%%%%%%%%
    
% %Put to 1 to define the stresses defined in 'stress input' as strain values    
% strain=0; 

	%%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%
    

Sxx = -1;                 
Syy = 0; 
Szz = 0;
Sxy = 0;
Sxz = 0;
Syz = 0; 
Option='B'; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Dss,Dds,Dn,Sxx,Syy,Szz,Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,...
 Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
 FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points);


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
    

%Define the number of cells on the square grid
cells=6;   
%How many extra cells you want to add away from the faults, 
%reduce to 0 if having half space issues 
padding=1;  
[maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ,sz]=MinMaxDataExtents3d(Points(:,2:4),cells,padding);
[X,Y,Z] = meshgrid(mingriX:sz:maxgriX,mingriY:sz:maxgriY,mingriZ:sz:maxgriZ);    
    
x = mingriX:sz:maxgriX; 
y = mingriY:sz:maxgriY;
z = mingriZ:sz:maxgriZ;
 

%Drawing figure of surface and the observation points
figure;trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',(.2),'FaceColor', [0.5 0 0.9 ]);
hold on
scatter3(X(:),Y(:),Z(:),5,'k')  %Showing surface and obs points
xlabel('x'); ylabel('y');zlabel('z'); axis('equal'); title('Surface and Obs Points');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses/Disps at defined observation points in XYZ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
[StressTTotal,StrainTTotal,StressTChg,StrainTChg,StressTReg,StrainTReg]=...
CalculateStressOnSurroundingPoints3d(Dss,Dds,Dn,mu,lambda,X,Y,Z,Sxx,...
Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Calculate Displacements.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
[Ux,Uy,Uz] = CalculateDisplacementOnSurroundingPoints3d...
(Dss,Dds,Dn,nu,X,Y,Z,P1,P2,P3,halfspace);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP AA: Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Running only if in MATLAB
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1;
	disp('Analytical solution doesn t work in octave, mex files don t compile') 
elseif isOctave==0;
    Sxx=StressTReg(1,1);
    Syy=StressTReg(1,2);
    Szz=StressTReg(1,3);
    Sxy=StressTReg(1,4);
    Sxz=StressTReg(1,5);
    Syz=StressTReg(1,6);
	[X,Y,Z,sigx,sigy,sigz,sigxy,sigxz,sigyz,UxAn,UyAn,UzAn]=MengEshelbyFunc(a,b,c,x,y,z,Sxx,Syy,Szz,Sxy,Sxz,Syz,E,nu);
end

[Sxx,Syy,Szz,Sxy,Sxz,Syz ] = ExtractCols( StressTTotal );
[Exx,Eyy,Ezz,Exy,Exz,Eyz ] = ExtractCols( StrainTTotal );


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: Visualisation and analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%Loop version of Dave Healy's is point outside Func
tol_dist = 0.4; %Distance away from hole, numerical solution will have artefacts close to boundary
[External] = InsideEllipsoid(X, Y, Z, a, b, c,tol_dist);
ExternalL=logical(External);
ExternalL_CV=ExternalL(:);%col vec of flag

%Graphically showing the principal stress directions from the numerical
%model
X(~ExternalL)=nan;    
Y(~ExternalL)=nan;    
X(~ExternalL)=nan; 
StressTTotal(~ExternalL_CV,:)=nan;

%Draws stress ellipsoids
DrawStressEllipsoidsPrincipal(StressTTotal,StrainTTotal,X,Y,Z,'Triangles',Triangles,'Points',Points);
DrawStressEllipsoidsPrincipal(StressTChg,StrainTChg,X,Y,Z,'Triangles',Triangles,'Points',Points);


[Sxx,Syy,Szz,Sxy,Sxz,Syz ] = ExtractCols( StressTTotal );

%Calculating eigenvalues (no shear)
[S1,S2,S3,S1dir,S2dir,S3dir] = EigCalc3d(Sxx(:),Syy(:),Szz(:),Sxy(:),Sxz(:),Syz(:));
[Sig1,Sig2,Sig3] = EigCalc3d(sigx(:),sigy(:),sigz(:),sigxy(:),sigxz(:),sigyz(:));

S1=reshape(S1,(size(X)));
S2=reshape(S2,(size(X)));
S3=reshape(S3,(size(X)));
Sig1=reshape(Sig1,(size(X)));
Sig2=reshape(Sig2,(size(X)));
Sig3=reshape(Sig3,(size(X)));
Ux=reshape(Ux,(size(X)));
Uy=reshape(Uy,(size(X)));
Uz=reshape(Uz,(size(X)));

%Removing interior points
S1(~ExternalL)=nan;    Sig1(~ExternalL)=nan;
S2(~ExternalL)=nan;    Sig2(~ExternalL)=nan;
S3(~ExternalL)=nan;    Sig3(~ExternalL)=nan;
Ux(~ExternalL)=nan;    UxAn(~ExternalL)=nan;
Uy(~ExternalL)=nan;    UyAn(~ExternalL)=nan;
Uz(~ExternalL)=nan;    UzAn(~ExternalL)=nan;

%Drawing isosurfaces from the two solutions:
IsoContoursPrincipalStressPercentiles...
( S1,S2,S3,X,Y,Z,'Triangles',Triangles,'Points',Points,'Alpha',0.2);
title('Numerical solution');		

IsoContoursPrincipalStressPercentiles...
( Sig1,Sig2,Sig3,X,Y,Z,'Triangles',Triangles,'Points',Points,'Alpha',0.2);
title('Analytical solution');

%Quiver, coloured, scaled and shows positive/negative parts
X(~ExternalL)=nan;  
Y(~ExternalL)=nan;     
Z(~ExternalL)=nan;



%Running tests to check match
S1Res=S1-Sig1;
bad=isnan(S1Res);
fprintf('MaxResidualS1 %i.\n',max(S1Res(~bad))) %prints on one line unlike 'disp

S2Res=S2-Sig2;
bad=isnan(S2Res);
fprintf('MaxResidualS2 %i.\n',max(S2Res(~bad))) 

S3Res=S3-Sig3;
bad=isnan(S3Res);
fprintf('MaxResidualS3 %i.\n',max(S3Res(~bad)))

UxRes=Ux-UxAn;
bad=isnan(UxRes);
fprintf('MaxResidualUx %i.\n',max(UxRes(~bad))) 

UyRes=Uy-UyAn;
bad=isnan(UyAn);
fprintf('MaxResidualUy %i.\n',max(UyAn(~bad))) 


P=[max(S1Res(~bad)),max(S2Res(~bad)),max(S3Res(~bad))];
if any(P>0.1)
    error('Error in surrounding grid is too large. This no longer matches analytical solutions, check the comparative images')
else
	%do nothing for the moment, checking disp also is ok first
end

P=[max(UxRes(~bad)),max(UyRes(~bad))];
if any(P>0.1)
    error('Error in surrounding grid is too large. This no longer matches analytical solutions, check the comparative images')
else
    disp('Everything looks good, checked max stress and disp residuals and flagged errors above 0.01')
end


