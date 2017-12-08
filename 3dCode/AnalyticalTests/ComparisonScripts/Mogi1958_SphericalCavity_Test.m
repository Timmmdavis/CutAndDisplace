% Test: Spherical chamber with internal pressure causing ground
% displacement
% 
% Solution: Lisowski, M., 2007. Analytical volcano deformation source
% models. In'Volcano Deformation'(pp. 279-304). Springer Berlin Heidelberg.
% Equation: 8.15. This is a simplified analytical equation for the
% deformation induced by a point source in a half space. This requires a
% finite 'magma chamber' or sphere. An internal pressure is specified and
% this inflates and displaces the ground surface. The solution is relative
% to 0,0 and in radial coordinates. For simplicity the comparison with the
% script uses a magma chamber placed at 0,0 and a strip of observation
% points are placed along the positive X axis. The radius used is 1 and the
% elastic constants are poisons ratio 0.25 and mu=1. These can be changed
% if needed.
%  
% Proof: This shows the half space solution works correctly for meshed 3d
% objects subject to user defined tractions on each element.


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
 string='SphereFix4.ts';
 [ PointsF,TrianglesF ] = GoCadAsciiReader( string );  
 PointsF(:,2:4)=PointsF(:,2:4)*3;
 

%string='SphereUniformDistributionON_46Faces.ts'; %~46 pnts
string='SphereUniformDistributionON2_500Faces.ts';  %375 pnts
%string='SphereUniformDistribution5120FacesSides.ts';  %1500 pnts
[ Points,Triangles ] = GoCadAsciiReader( string );


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Appending and flagging the fixed data points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Appending the triangles and points to the places we are going to fix. 
[Points,Triangles,Fdisp] = DataAppender3d( Points,PointsF,Triangles,TrianglesF);
	

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define full/halfspace and elastic constants.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Defining halfspace. 1 is on, 0 is off. The code will run much faster with
% this off as the calculations are simpler in a full space.
halfspace = 1; 

% Defining the height of the freesurface relative to the value of 0 of
% your imported surfaces. 
freesurface_height = 50;

% Defining elastic constants
mu=1;           	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
[ ~,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu );

% Getting surface ready for script. 
Points=[Points(:,1),Points(:,2),Points(:,3),(Points(:,4)-freesurface_height)];
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 

% Drawing the fixed calculated trianglulation
hold on
trisurf(Triangles(logical(Fdisp),:),Points(:,2),Points(:,3),Points(:,4),'facecolor', 'red');
hold off

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


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option E = Defining boundary conditions as traction at the element
    %centres. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
Tss = 0;     
Tds = 0;        
Tn = 1;   

Option='E'; %Traction defined on elements (pressure)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(Option,'A')    
    [Dss,Dds,Dn,Sxx,Syy,Szz,Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,...
     Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
     FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points);
end

Sxx = 0;  
Syy = 0; 
Szz = 0; 
Sxy = 0;   
Sxz = 0;
Syz = 0;



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
    

%Obs points
X=linspace(0,100,100);
Y=zeros(1,numel(X));
Z=zeros(1,numel(X));

% %Grid if needed:
% [X,Y] = meshgrid(-50:1:50,-50:1:50);	
% Z=zeros(size(X));


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses/Disps at defined observation points in XYZ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Calculate Displacements.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ux,Uy,Uz] = CalculateDisplacementOnSurroundingPoints3d...
(Dss,Dds,Dn,nu,X,Y,Z,P1,P2,P3,halfspace);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: Visualisation and analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%Depth of chamber
D=50; %depth
Radius=1; %Radius of the magma chaber 
P=Tn; %Pressure

[u,v,w]=Mogi1958_SphericalCavity(D,X',Y',Radius,P,nu,mu);


figure;
hold on
plot(X,u,'b','LineWidth',2.5)
plot(X,w,'g','LineWidth',2.5)

XN = linspace(0,100,21);
UxN = interp1(X,Ux,XN,'pchip');
UzN = interp1(X,Uz,XN,'pchip');

scatter(XN,UxN,24,'k','filled')
scatter(XN,UzN,24,'k','filled')
xlabel('Radial distance from chamber centre'); ylabel('Ground displacement'); title ('Mogi point source displacement');
titlesz=25;
fntsz=21;
grid on
ChangeFontSizes(fntsz,titlesz);

ResUx=(u-Ux);
ResUy=(w-Uz);
fprintf('MaxResidualUx %i.\n',max(abs(ResUx(:))))
fprintf('MaxResidualUy %i.\n',max(abs(ResUy(:))))

%If Error is above certain value flip and error
P=(max(abs(ResUx(:)))+max(abs(ResUy(:))));
if any(P>3e-5)
    error('Error This deviates too far from the analytical solution, check the comparative disp figure')
else
    disp('Everything looks good, tolerance checks max disp residuals and flags errors above 3e-5')
end

%calculating error in %
Error_perc_x=((Ux./u).*100)-100; bad=Error_perc_x==inf;  Error_perc_x(bad)=0;  MaxError_perc_x=max(max(Error_perc_x));
Error_perc_z=((Uz./w).*100)-100; bad=Error_perc_z==inf;  Error_perc_z(bad)=0;  MaxError_perc_z=max(max(Error_perc_z));
