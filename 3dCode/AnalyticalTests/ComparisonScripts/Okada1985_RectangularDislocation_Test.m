% Test: Rectangular fault slip causing ground surface displacement
% 
% Solution: Okada, Y., 1985. Surface deformation due to shear and tensile
% faults in a half-space.�Bulletin of the seismological society of
% America,�75(4), pp.1135-1154. Code for this written in MATLAB by Fran�ois
% Beauducel :
% http://uk.mathworks.com/MATLABcentral/fileexchange/25982-okada--surface-deformation-due-to-a-finite-rectangular-source
% A fault that matches the one used in the Okada code is made from two
% triangles in the TDE code. A dipslip displacement 1 is used and the
% ground surface displaced. It checks that all ground surface displacements
% are within 1e-8 figures. The fault is at a different orientation to the
% Cartesian axes: a dip of 70 and azimuth of 150.
% 
% Proof: This shows the superposition is working correctly and the original
% M Nikko TDE functions have not changed. These were independently tested
% agasit Okada in the original TDE solution paper.


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


   string='OkadaFault_ON_2Faces.ts';
 [ Points,Triangles ] = GoCadAsciiReader( string );

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define full/halfspace and elastic constants.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


SecondSurface=0;
% Defining halfspace. 1 is on, 0 is off. The code will run much faster with
% this off as the calculations are simpler in a full space.
halfspace = 1; 

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

%Locked Els
Fdisp= zeros(size(Triangles(:,1))); 

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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Constant Slip across the fault surface.
	%Runs a constant slip on every element. Output 'stresses' are blank
	%arrays of 0's as this option uses displacement boundary conditions.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
sz=zeros(numel(P1(:,1)),1);
Dss =0;      Dss    = sz+Dss;   %Positive = LeftLatMovement 
Dds =1;      Dds    = sz+Dds;   %Positive = Reverse movement
Dn  =0;      Dn     = sz+Dn;    %Positive = Opening movement
clear sz 

Sxx = 0; 				
Syy = 0; 
Szz = 0; 
Sxy = 0; 					
Sxz = 0; 					
Syz = 0;
Option='A'; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XYZ observation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Simply define flat observation plane of points. XYZ with defined step size 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
[X,Y,Z] = meshgrid(-5:0.5:10,-5:0.5:10,0);
rowcount = length(X(:,1)); 
colcount = length(X(1,:));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option D = Load a secondary Fault Surface
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Secondary Surface
% SecondSurface=1; %Flag set
% string='GoCadExportHector.ts';
% [ PointsObs,TrianglesObs ] = GoCadAsciiReader( string );
% PointsObs(:,2:4)=PointsObs(:,2:4)/0.5E4;
% PointsObs=[PointsObs(:,1),PointsObs(:,2)-max(PointsObs(:,2)),PointsObs(:,3)-max(PointsObs(:,3)),(PointsObs(:,4)-freesurface_height)];
% hold on; trisurf(TrianglesObs,PointsObs(:,2),PointsObs(:,3),PointsObs(:,4),'FaceAlpha',(.8),'facecolor', 'cyan');
% [MidPointObs,FaceNormalVectorObs] = MidPointCreate(PointsObs,TrianglesObs,0);
% X=MidPointObs(:,1);
% Y=MidPointObs(:,2);
% Z=MidPointObs(:,3);
% rowcount = length(X(:,1)); 
% colcount = length(X(1,:));


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses/Disps at defined observation points in XYZ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%  [StressTTotal,StrainTTotal,StressTChg,StrainTChg,StressTReg,StrainTReg]=...
%      CalculateStressOnSurroundingPoints3d(Dss,Dds,Dn,mu,lambda,X,Y,Z,Sxx,Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu);

%Testing the obs influence building function
 [StrainTTotal,StrainTChg,StrainTReg]...
   =ObsPointsStressCalculator3d(Dn,Dss,Dds,mu,lambda...
   ,X,Y,Z,Sxx,Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu);

[ exx,eyy,ezz,exy,exz,eyz] = ExtractCols( StrainTTotal );
[ Sxx,Syy,Szz,Sxy,Sxz,Syz ] = HookesLaw3dStrain2Stress(exx,eyy,ezz,exy,exz,eyz,lambda,mu );


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Calculate Displacements.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ux,Uy,Uz] = CalculateDisplacementOnSurroundingPoints3d...
(Dss,Dds,Dn,nu,X,Y,Z,P1,P2,P3,halfspace);

if SecondSurface==0
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %STEP 7: VISUALISATION AND ANALYSIS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %Calculating the depth of the centre of the fault using Pythagoras theorem
    dip=70;
    width=5;
    TipDepth=1;
    MidDeptha=sind(dip)*width;
    MidDepth=MidDeptha/2+TipDepth;
    %Running the Okada function, see: http://uk.mathworks.com/MATLABcentral/fileexchange/25982-okada--surface-
    %deformation-due-to-a-finite-rectangular-source
    %For details
    [uE,uN,uZ,uNN,uNE,uEN,uEE]  = Okada1985_RectangularDislocation(X,Y,MidDepth,60,dip,5,width,90,1,0,'plot') ;
    hold on

    uxy=-(uEN(:)+uNE(:))/2;

    %Tri dislocations
    [E1,E2,E1dir,E2dir]=EigCalc2d(exx,eyy,exy);
    %Okada dislocations
    [E1ok,E2ok,E1dirok,E2dirok]=EigCalc2d(-uEE(:),-uNN(:),-uxy);

    [E1,E2,E1ok,E2ok ]=ReshapeData2d( rowcount,colcount,E1,E2,E1ok,E2ok );


    %Plotting the displacement induced in the TDE code and its fault
    figure;subplot(1,2,1),quiver3(X(:),Y(:),Z(:),Ux,Uy,Uz)
    xlabel('x'); ylabel('y'); axis('equal'); title('DisplacementAndFaultTDECode')
    hold on
    trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceColor', 'cyan', 'faceAlpha',0.8);
    hold off
    %Drawing the displacement vectors induced by the Okada function
    subplot(1,2,2),quiver3(X,Y,Z,uE,uN,uZ)
    xlabel('x'); ylabel('y'); axis('equal'); title('Displacement OkadaCode')
    hold off

    %Reshaping vectors before calculating residual
    Ux=reshape(Ux,(size(uE)));
    Uy=reshape(Uy,(size(uE)));
    Uz=reshape(Uz,(size(uE)));

    %Calculating residuals
    Resx=uE-Ux;
    Resy=uN-Uy;
    Resz=uZ-Uz;

    figure;subplot(1,3,1),contourf(X,Y,Resx);xlabel('x'); ylabel('y'); axis('equal'); title('Disp residual ux');colorbar;
    subplot(1,3,2),contourf(X,Y,Resy);xlabel('x'); ylabel('y'); axis('equal'); title('Disp residual uy');colorbar;
    subplot(1,3,3),contourf(X,Y,Resz);xlabel('x'); ylabel('y'); axis('equal'); title('Disp residual uz');colorbar;


    %Plotting Tri Dis displacements vs the Okada displacements
    DrawDeformedGrid2d( X,Y,Ux,Uy,cmap2,(E1+E2),'Scale',5);
    hax1=gca;no1=gcf;title('Triangular dislocations');
    DrawDeformedGrid2d( X,Y,uE,uN,cmap2,(E1+E2),'Scale',5);
    hax2=gca;no2=gcf;title('Rectangular dislocation');
    hf2=figure;
    s1=subplot(1,2,1);
    s2=subplot(1,2,2);
    pos1=get(s1,'Position');pos2=get(s2,'Position');
    delete(s1);delete(s2);
    hax_1=copyobj(hax1,hf2);
    hax_2=copyobj(hax2,hf2);
    set(hax_1, 'Position', pos1);colorbar;colormap(cmap2);
    xlabel('x'); ylabel('y');gca;% title('Triangles');
    set(hax_2, 'Position', pos2);colorbar;colormap(cmap2);
    xlabel('x'); ylabel('y');gca;%title('gggg');
    close(no1);close(no2);
    titlesz=25;
    fntsz=21;
    ChangeFontSizes(fntsz,titlesz);
    text(-25,13,'Ground displacements and dilatation, comparison of different analyical sources','Fontsize',25,'FontWeight','bold')

    bad=isnan(Resx);Resx=Resx(~bad);%removing nans to get proper mean
    fprintf('MaxResidualUx %i.\n',max(abs(Resx(:))))

    bad=isnan(Resy);Resy=Resy(~bad);%removing nans to get proper mean
    fprintf('MaxResidualUy %i.\n',max(abs(Resy(:))))

    bad=isnan(Resz);Resz=Resz(~bad);%removing nans to get proper mean
    fprintf('MaxResidualUz %i.\n',max(abs(Resz(:))))

    P=max(abs(Resx(:)))+max(abs(Resy(:)))+max(abs(Resz(:)));
    if any(P>1e-8) %Around the precision and errors introduced. Small values. 
        error('Calculated displacements are a long way from the analytical solutions')
    else
        disp('Everything looks good, tolerance checks max stress residuals and flags errors above 1e-8')
    end

    %Calculating residuals
    Resex=exx--uEE(:);
    Resey=exy-(-(uEN(:)+uNE(:))/2);

else
    
    Mu=ones(size(Sxx))*0.6; %Coeff Friction
    Cohesion=zeros(size(Sxx)); %Coeff Friction
    [CSS] = CalculateCoulombStressOnPlane(MidPointObs,FaceNormalVectorObs,...
    Sxx,Syy,Szz,Sxy,Sxz,Syz,Mu,Cohesion,PointsObs,TrianglesObs,cmap2 );
    
end
