% Test: Rectangular fault slip causing ground surface displacement
% 
% Solution: Okada, Y., 1985. Surface deformation due to shear and tensile
% faults in a half-space. Bulletin of the seismological society of
% America, 75(4), pp.1135-1154. Code for this written in MATLAB by François
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

%   Copyright 2017, Tim Davis, The University of Aberdeen



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you do not need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

clear;close all

% %===== add file paths ==========================
pathstring = pwd;                                   %Get the address of the current working directory
if ispc; parts = strsplit(mfilename('fullpath'), '\');       %Getting the address of the script working directory and splitting into cell array
else; parts = strsplit(mfilename('fullpath'), '/');  end
[~,n] = find(~cellfun(@isempty,strfind(parts,'BEM_DDM_MATLAB'))); %finding the scripts root directory and its location (n)
if ispc; addpath(genpath(strjoin(parts(1,1:n),'\'))); %adding all folders to the path that are in this dir 
else; addpath(genpath(strjoin(parts(1,1:n),'/'))); end;%mac/linux
cd(pathstring)                                       %jumping back to the current dir

	%Script to run dislocations on a fault surface

	% This follows the same unit convention as Poly3d so ''assumes dimensionally consistent units for physical quantities with
	% dimensions of length or stress. Thus if you use kilometers for any quantity with
	% dimensions of length (e.g. a coordinate position), you must use kilometers for all
	% quantities with dimensions of length (e.g. displacement discontinuity).''
	% Andrew Lyell Thomas, Masters thesis, Stanford University. 


cmap = colormap_cpt('Ccool-warm'); %Loads a colourmap to be used for figures. This is Kenneth Morelands diverging colourmap.
cmap2 = colormap_cpt('Ccool-warm2');

		
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: Import the fault surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading formatted ascii data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    %  Loads the triangles that define the fault surface 
    %  Splits this into two files that define the vertex's and which vertex's
    %  make which triangle. Use Gocad ascii format and rejig like the
    %  example files. Note column 1 for the points must start at 1.
    %  Reexport Gocad ascii file to get this. 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Loading Gocad ascii data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   string='OkadaFault_ON_2Faces.ts';
 [ Points,Triangles ] = GoCadAsciiReader( string,pathstring );

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define Full/Halfspace and Elastic constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %Defining halfspace. 1 is on, 0 is off. The code will run much faster with this off as the calculations are simpler in a full space. 
	%Use this if the algorithm is throwing errors for vertex's above 0m.
    
halfspace = 1; 

	%Defining the height of the freesurface relative to the value of 0 of
	%your imported surfaces. Used when deformation is not at
	%sealevel. See the output MaxZ height if your half space clashes with
	%your points and adjust until this outputs as a negative value. 
freesurface_height = 0;

    %Defining elastic constants, All three are required to run this whole script.  See equations below if you need to calculate these 
	%from other constants. 
	%Note nu is needed in the displacement calculation.  
	
mu=1;%0.4;%4000 	%Shear Mod, mu or G. Relates shear stress to shear strain. 
lambda=1;%4000;  %Lame's constant  0 perfectly compressible like cork, Infinate for incompressible material like rubber.
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;

	%Equations you can use to calculate elastic parameters that you need above. 
	
% E = 10000;                        %Young's Modulus
% K = E/3*(1-2*nu);                 %Bulk Modulus.    Equation 8.25 Pollard
% mu = E/(2*(1+nu));                %Shear Modulus.   Equation 8.26 Pollard  which is Geo convention 
% lambda= E*nu/((1+nu)*(1-2*nu));   %Lamé's  constant,  Equation 8.27 Pollard
% nu =lambda/(2*(mu+lamda);            %Poisson's ratio, Equation 8.28 Pollard


	%Function to extract the three vertex's for each triangle and 
	%triangles defined by a combination of vertex Id's. This is used later by the algorithms 
    %The max Z value is output as text in the command window giving an idea
    %of minimum free surface height. 
Points=[Points(:,1),Points(:,2),Points(:,3),(Points(:,4)-freesurface_height)];
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 
	%Fictitious disp flag, if set to one then this forces the triangles with this to 0 displacement. 
Fdisp= zeros(size(Triangles(:,1))); 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, Remote, Remote
%no opening, tractions on the elements or user imported slip value for each face. 

%Stress tensor outputs from this part of the toolbox i.e. Sxx Syy etc are simply the user defined remote stress
%components (boundary conditions). These are used later in the stress calculations. If strain is input these are converted to stress 
%using Hooke's law. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Constant Slip across the fault surface
	%Runs a constant slip on every triangle. A mixture of two slips will give an inclined slip.
	%Output 'stresses' are blank arrays of 0's as this option is not driven by boundary conditions. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
a=size(P1);
b=a(1,1);
c=zeros(b,1); 
    
StrikeSlipDisp  =0;      StrikeSlipDisp     = c+StrikeSlipDisp;     %Positive = LeftLatMovement 
DipSlipDisp     =1;      DipSlipDisp		= c+DipSlipDisp;        %Positive = Reverse movement
TensileSlipDisp =0;      TensileSlipDisp	= c+TensileSlipDisp;    %Positive = Opening movement
clear a b c             	%a b and c are used to create an array of slips the size of the number of triangles on the fault surface.

Sxx = 0; 					
Syy = 0; 
Szz = 0; 
Sxy = 0; 					
Sxz = 0; 					
Syz = 0;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting the calculated slips.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure('name','Step 3, option D'); colormap(cmap),trisurf(TR,StrikeSlipDisp);
% xlabel('x'); ylabel('y'); axis('equal'); title('SSDisp'),colorbar;
% figure('name','Step 3, option D'); colormap(cmap),trisurf(TR,DipSlipDisp);
% xlabel('x'); ylabel('y'); axis('equal'); title('DSDisp'),colorbar;
% figure('name','Step 3, option D'); colormap(cmap),trisurf(TR,TensileSlipDisp);
% xlabel('x'); ylabel('y'); axis('equal'); title('TSDisp'),colorbar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting the calculated slips.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XYZ observation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Simply define flat observation plane of points. XYZ with defined step size 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
[X,Y,Z] = meshgrid(-5:0.5:10,-5:0.5:10,0);
rowcount = length(X(:,1)); 
colcount = length(X(1,:));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses  on dispersed XYZ
%observation points. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%TotalStress,Stress and StressChange are Xx12 arrays of strain tensors and stress tensors. XX YY ZZ XY XZ YZ
	%Strain is the first 6 coloumns and stress the last 6. Stress is the regional stress 
	%StressChange is the stress on the points from the event.
	%TotalStress is the driving stress and stress change from the event added. 
	
 [TotalStrainStress,StressStrainChange,StressStrainRemoteSxx,Syy,Szz,Sxy,Sxz,Syz]=CalculateStressOnSurroundingPoints(StrikeSlipDisp,DipSlipDisp,TensileSlipDisp,mu,lambda,X,Y,Z,Sxx,Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu);
[exx,eyy,ezz,exy,exz,eyz,Sxx,Syy,Szz,Sxy,Sxz,Syz ] = ExtractCols( TotalStrainStress );


     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Calculate Displacements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  [Ux,Uy,Uz] = CalculateDisplacementOnSurroundingPoints( StrikeSlipDisp,DipSlipDisp,TensileSlipDisp, nu, X,Y,Z, P1, P2, P3,halfspace);


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
[uE,uN,uZ,uNN,uNE,uEN,uEE]  = okada85(X,Y,MidDepth,60,dip,5,width,90,1,0,'plot') ;
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
DrawDeformedGrid2d( X,Y,Ux,Uy,5,cmap2,(E1+E2));hax1=gca;no1=gcf;title('Triangular dislocations');
DrawDeformedGrid2d( X,Y,uE,uN,5,cmap2,(E1+E2));hax2=gca;no2=gcf;title('Rectangular dislocation');
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

