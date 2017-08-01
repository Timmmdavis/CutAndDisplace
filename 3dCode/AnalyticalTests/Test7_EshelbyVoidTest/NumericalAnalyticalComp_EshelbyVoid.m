%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you do not need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

%   Copyright 2017, Tim Davis, The University of Aberdeen
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
		
	%Precision:
		%Due to the fairly large nature of some models your machine may run out of memory. Simplicity is key when getting to know your model before 
		%increasing the complexity and sampling. 
		%In functions that are creating large arrays the precision has been reduced to 'single' in this script. See example below. Remove the word 'single' if you need a higher precision.
		%From the MATLAB help - ''Use double-precision to store values greater than approximately 3.4 x 10^38 or less than approximately -3.4 x 10^38. 
		%For numbers that lie between these two limits , you can use either double- or single-precision, but single requires less memory.''
		%Preferably you would change your units /models scale before increasing the precision and doubling the memory usage.  
		%SSinfmatrix = zeros(NUM*NUM,12,'single'); 
		
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

% S = load('TwoTriangle.dat');          %Loads the dat file -	GocadAsciiexport formatted with no text
% Points = S([1:4],[1,2,3,4]);          %Points XYZ that are the fault. First numbers define which rows, second which columns.
% Triangles = S([5:6],[1,2,3]);         %which points make a triangle
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Loading Gocad ascii data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 string='SphereFix.ts';
 [ PointsF,TrianglesF ] = GoCadAsciiReader( string,pathstring );    
    
%string='HalfShell_750Faces.ts'; %Only half the surface
%string='SphereUniformDistributionON_46Faces.ts'; %~46 pnts
%string='SphereUniformDistributionON2_500Faces.ts';  %375 pnts
%string='SphereUniformDistributionON_96Faces.ts';
string='SphereUniformDistributionHS_ON_1500Faces.ts';  %1500 pnts
%string='SphereUniformDistribution5120FacesSide.ts';  %1500 pnts
%string='Sphere.ts';  %1500 pnts

[ Points,Triangles ] = GoCadAsciiReader( string,pathstring );


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define Full/Halfspace and Elastic constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %Defining halfspace. 1 is on, 0 is off. The code will run much faster with this off as the calculations are simpler in a full space. 
	%Use this if the algorithm is throwing errors for vertex's above 0m.
    
halfspace = 0; 

	%Defining the height of the freesurface relative to the value of 0 of
	%your imported surfaces. Used when deformation is not at
	%sealevel. See the output MaxZ height if your half space clashes with
	%your points and adjust until this outputs as a negative value. 
freesurface_height =0;

    %Defining elastic constants, All three are required to run this whole script.  See equations below if you need to calculate these 
	%from other constants. 
	%Note nu is needed in the displacement calculation.  
	
mu=1;%0.4;%4000 	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
lambda=mu*(2*nu/(1-2*nu));%4000;  %Lame's constant  0 perfectly compressible like cork, Infinite for incompressible material like rubber.
E = mu*(2*(1+nu)) ; 
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
Fdisp(end-IF_Sz+1:end,:)=Fdisp(end-IF_Sz+1:end,:)+1;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, Remote, Remote
%no opening, tractions on the elements or user imported slip value for each face. 

%Stress tensor outputs from this part of the toolbox i.e. Sxx Syy etc are simply the user defined remote stress
%components (boundary conditions). These are used later in the stress calculations. If strain is input these are converted to stress 
%using Hooke's law. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
[Sxx,Syy,Szz,Sxy,Sxz,Syz,Tnn,Tss,Tds,Mu,Sf,strain ] = CreateBlankVars;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. Choose to define in stress or strain. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	%%%%%%%%%%%%%%
    %StrainInput        If halfspace = 1, Z stress components are removed
	%%%%%%%%%%%%%%
    
strain=1;                   %Put to 1 to define the stresses defined in 'stress input' as strain values
    
	%%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%
    

Sxx = -1;                    %Positive is tension
Syy = 0; 
Szz = 0;
Sxy = 0;        			%Positive = left lat when frac strikes NS and right lat when EW. See Pollard Fletcher Diagram - 6.13
Sxz = 0;% 0.05;
Syz = 0; 
Option='B'; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[StrikeSlipDisp,DipSlipDisp,TensileSlipDisp,Sxx,Syy,Szz,...
Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tnn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removing any fixed triangles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removing any fixed triangles
if any(Fdisp)==1
[Triangles,FaceNormalVector,MidPoint,P1,P2,P3]...
    = RemovingFixedEls3d(Triangles,FaceNormalVector,MidPoint,P1,P2,P3,Fdisp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting the calculated slips.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,StrikeSlipDisp,DipSlipDisp,TensileSlipDisp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Animate Fault movement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To do this just call 'Animate' in the cmd window. To change the
%parameters, save movie or the camera angle go into the 'Animate' function

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XYZ observation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Simply define flat observation plane of points. XYZ with defined step size 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cells=6;   %Define the number of observations points on the grid around the fault (not extended) %6
padding=1;  %How much extra bumph you want to add away from the grid. DISTANCE not cells.        %1
[maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ,sz]=MinMaxDataExtents3d(Points,cells,padding);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses and/or displacements on dispersed XYZ
%observation points. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%TotalStrainStress,StrainStressChange and StrainStressRemote are Xx12 arrays of strain tensors and stress tensors. XX YY ZZ XY XZ YZ
	%Strain is the first 6 coloumns and stress the last 6. StressStrainRemoteSxx is the regional stress 
	%StressStrainChange is the stress on the points from the event.
	%TotalStressStrain is the driving stress and stress change from the event added. 
	
 [TotalStrainStress,StrainStressChange,StrainStressRemote,Syy,Szz,Sxy,Sxz,Syz]=CalculateStressOnSurroundingPoints(StrikeSlipDisp,DipSlipDisp,TensileSlipDisp,mu,lambda,X,Y,Z,Sxx,Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu);

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Calculate Displacements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  [Ux,Uy,Uz] = CalculateDisplacementOnSurroundingPoints( StrikeSlipDisp,DipSlipDisp,TensileSlipDisp, nu, X,Y,Z, P1, P2, P3,halfspace);
  
  
  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP AA: Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sxx=StrainStressRemote(1,7);
Syy=StrainStressRemote(1,8);
Szz=StrainStressRemote(1,9);
Sxy=StrainStressRemote(1,10);
Sxz=StrainStressRemote(1,11);
Syz=StrainStressRemote(1,12);

%Running only if in MATLAB
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1;
disp('Analytical solution doesn t work in octave, mex files don t compile') 
elseif isOctave==0;
[X,Y,Z,sigx,sigy,sigz,sigxy,sigxz,sigyz,UxAn,UyAn,UzAn]=MengEshelbyFunc(a,b,c,x,y,z,Sxx,Syy,Szz,Sxy,Sxz,Syz,E,nu);
end

%CreatingStrainStressVec
[ exx,eyy,ezz,exy,exz,eyz ] = HookesLaw3dStress2Strain( sigx(:),sigy(:),sigz(:),sigxy(:),sigxz(:),sigyz(:),lambda,mu );
ZersAndStress=[exx(:),eyy(:),ezz(:),exy(:),exz(:),eyz(:),sigx(:),sigy(:),sigz(:),sigxy(:),sigxz(:),sigyz(:)];
sigxnr=sigx(:)-StrainStressRemote(1,7);
[ exx,eyy,ezz,exy,exz,eyz ] = HookesLaw3dStress2Strain( sigxnr(:),sigy(:),sigz(:),sigxy(:),sigxz(:),sigyz(:),lambda,mu );
ZersAndStressChng=[exx(:),eyy(:),ezz(:),exy(:),exz(:),eyz(:),sigxnr(:),sigy(:),sigz(:),sigxy(:),sigxz(:),sigyz(:)];




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: VISULISATION AND ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Loop version of Dave Healys is point outside Func
tol_dist = 0.4; %Distance away from hole, numerical solution will have artefacts close to boundary
[External] = IsPointOutsideL(X, Y, Z, a, b, c,tol_dist);
ExternalL=logical(External);
ExternalL_CV=ExternalL(:);%col vec of flag

%Graphically showing the principal stress directions from the numerical
%model
X(~ExternalL)=nan;    
Y(~ExternalL)=nan;    
X(~ExternalL)=nan; 
TotalStrainStress(~ExternalL_CV,:)=nan;

%Draws stress ellipsoids
DrawStressEllipsoidsPrincipal(StrainStressChange,X,Y,Z,Triangles,Points);
DrawStressEllipsoidsPrincipal(TotalStrainStress,X,Y,Z,Triangles,Points);

%Draws principal stress ellipsoids
DrawS1S2S3Directions(StrainStressChange,X,Y,Z,Triangles,Points);
DrawS1S2S3Directions(TotalStrainStress,X,Y,Z,Triangles,Points);

Sxx=TotalStrainStress(:,7);
Syy=TotalStrainStress(:,8);
Szz=TotalStrainStress(:,9);
Sxy=TotalStrainStress(:,10);
Sxz=TotalStrainStress(:,11);
Syz=TotalStrainStress(:,12);
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
Ux(~ExternalL)=nan;     UxAn(~ExternalL)=nan;
Uy(~ExternalL)=nan;     UyAn(~ExternalL)=nan;
Uz(~ExternalL)=nan;     UzAn(~ExternalL)=nan;

%finding percentile values for 3d plotting
Sig1Ord=sort(Sig1(:));
S1_75=Sig1Ord(round((numel(Sig1)/100)*75,0)); %75th percentile
S1_50=Sig1Ord(round((numel(Sig1)/100)*50,0)); %50th percentile
S1_25=Sig1Ord(round((numel(Sig1)/100)*25,0)); %25th percentile
Sig2Ord=sort(Sig2(:));
S2_75=Sig2Ord(round((numel(Sig2)/100)*75,0)); %75th percentile
S2_50=Sig2Ord(round((numel(Sig2)/100)*50,0)); %50th percentile
S2_25=Sig2Ord(round((numel(Sig2)/100)*25,0)); %25th percentile
Sig3Ord=sort(Sig3(:));
S3_75=Sig3Ord(round((numel(Sig3)/100)*75,0)); %75th percentile
S3_50=Sig3Ord(round((numel(Sig3)/100)*50,0)); %50th percentile
S3_25=Sig3Ord(round((numel(Sig3)/100)*25,0)); %25th percentile
% S1_75 = prctile(Sig1(:),75);%75th percentile
% S1_50 = prctile(Sig1(:),50);%50th percentile
% S1_25 = prctile(Sig1(:),25);%25th percentile
% S2_75 = prctile(Sig2(:),75);%75th percentile
% S2_50 = prctile(Sig2(:),50);%50th percentile
% S2_25 = prctile(Sig2(:),25);%25th percentile
% S3_75 = prctile(Sig3(:),75);%75th percentile
% S3_50 = prctile(Sig3(:),50);%50th percentile
% S3_25 = prctile(Sig3(:),25);%25th percentile


%Drawing isosurfaces of the principal stresses. 
figure_handle = figure;subplot(2,3,1), pneg=patch(isosurface(X,Y,Z,S1,S1_50)) ;
        set(pneg,'facecolor',RGB2Colour(255,48,48),'EdgeColor','none'); %red
        hold on
        p0 = patch(isosurface(X,Y,Z,S1,S1_25)) ; 
        set(p0,'FaceColor',RGB2Colour(255,165,0),'EdgeColor','none')    %orange 
        ppos = patch(isosurface(X,Y,Z,S1,S1_75)) ; 
        set(ppos,'FaceColor',RGB2Colour(255,255,0),'EdgeColor','none')  %yellow
grid on ; 
daspect([1 1 1]) ; 
view(3); 
axis equal ; 
camlight ; 
lighting gouraud ; 
title({'\fontsize{14}Sig1TDEScript','\fontsize{8}Percentiles,Red=50th Orange=25th Yellow=75th'})
hold off;
subplot(2,3,4), pneg=patch(isosurface(X,Y,Z,Sig1,S1_50)) ;
        set(pneg,'FaceColor',RGB2Colour(255,48,48),'EdgeColor','none');
        hold on
        p0 = patch(isosurface(X,Y,Z,Sig1,S1_25)) ; 
        set(p0,'FaceColor',RGB2Colour(255,165,0),'EdgeColor','none')    
        ppos = patch(isosurface(X,Y,Z,Sig1,S1_75)) ; 
        set(ppos,'FaceColor',RGB2Colour(255,255,0),'EdgeColor','none')    
grid on ; 
daspect([1 1 1]) ; 
view(3); 
axis equal ; 
camlight ; 
lighting gouraud ; 
title({'\fontsize{14}Sig1AnalyticalEshelby','\fontsize{8}Percentiles,Red=50th Orange=25th Yellow=75th'})
hold off;
subplot(2,3,2),pneg=patch(isosurface(X,Y,Z,S2,S2_50)) ; %subplot(2,3,2),
        set(pneg,'FaceColor',RGB2Colour(192,255,62),'EdgeColor','none');  %olive green
        hold on
        p0 = patch(isosurface(X,Y,Z,S2,S2_25)) ; 
        set(p0,'FaceColor',RGB2Colour(32,178,170),'EdgeColor','none');    %turquoise
        ppos = patch(isosurface(X,Y,Z,S2,S2_75)) ; 
        set(ppos,'FaceColor',RGB2Colour(224,102,255),'EdgeColor','none'); %purple
grid on ; 
daspect([1 1 1]) ; 
view(3); 
axis equal ; 
camlight ; 
lighting gouraud ; 
title({'\fontsize{14}Sig2TDEScript','\fontsize{8}Percentiles,Olive green=50th Turquoise=25th Purple=75th'})
hold off;
subplot(2,3,5),pneg=patch(isosurface(X,Y,Z,Sig2,S2_50)) ; %subplot(2,3,2),
        set(pneg,'FaceColor',RGB2Colour(192,255,62),'EdgeColor','none');
        hold on
        p0 = patch(isosurface(X,Y,Z,Sig2,S2_25)) ; 
        set(p0,'FaceColor',RGB2Colour(32,178,170),'EdgeColor','none');     
        ppos = patch(isosurface(X,Y,Z,Sig2,S2_75)) ; 
        set(ppos,'FaceColor',RGB2Colour(224,102,255),'EdgeColor','none');
grid on ; 
daspect([1 1 1]) ; 
view(3); 
axis equal ; 
camlight ; 
lighting gouraud ; 
title({'\fontsize{14}Sig2AnalyticalEshelby','\fontsize{8}Percentiles,Olive green=50th Turquoise=25th Purple=75th'})
hold off;
subplot(2,3,3),pneg=patch(isosurface(X,Y,Z,S3,S3_50)) ;
        set(pneg,'FaceColor',RGB2Colour(148,148,148),'EdgeColor','none'); %grey
        hold on
        p0 = patch(isosurface(X,Y,Z,S3,S3_25)) ; 
        set(p0,'FaceColor',RGB2Colour(152,245,255),'EdgeColor','none');   %blue    
        ppos = patch(isosurface(X,Y,Z,S3,S3_75)) ; 
        set(ppos,'FaceColor',RGB2Colour(148,0,211),'EdgeColor','none'); %dark purp  
grid on ; 
daspect([1 1 1]) ; 
view(3); 
axis equal ; 
camlight ; 
lighting gouraud ; 
title('Sig3TDEScript')
title({'\fontsize{14}Sig3TDEScript','\fontsize{8}Percentiles,Grey=50th Blue=25th Dark purple=75th'})
hold off;
subplot(2,3,6),pneg=patch(isosurface(X,Y,Z,Sig3,S3_50)) ; %subplot(2,3,2),
        set(pneg,'FaceColor',RGB2Colour(148,148,148),'EdgeColor','none'); %grey
        hold on
        p0 = patch(isosurface(X,Y,Z,Sig3,S3_25)) ; 
        set(p0,'FaceColor',RGB2Colour(152,245,255),'EdgeColor','none');   %blue       
        ppos = patch(isosurface(X,Y,Z,Sig3,S3_75)) ; 
        set(ppos,'FaceColor',RGB2Colour(148,0,211),'EdgeColor','none'); %dark purp 
grid on ; 
daspect([1 1 1]) ; 
view(3); 
axis equal ; 
camlight ; 
lighting gouraud ;
title({'\fontsize{14}Sig3AnalyticalEshelby','\fontsize{8}Percentiles,Grey=50th Blue=25th Dark purple=75th'})
hold off;
%Linking all axes together to spin at the same time
all_ha = findobj( figure_handle, 'type', 'axes', 'tag', '' );
hlink = linkprop(all_ha,{'CameraPosition','CameraUpVector'});
rotate3d on


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

%%%if you wanted to draw the isosurfaces in OCTAVE this is probably the best way (no lighting funcs exist)
% %Drawing isosurfaces of the principal stresses. 
% figure_handle = figure;subplot(1,3,1), pneg=patch(isosurface(X,Y,Z,S1,S1_50)) ;
%         set(pneg,'facecolor',normc([255;48;48])); %red
%         hold on
%         p0 = patch(isosurface(X,Y,Z,S1,S1_25)) ; 
%         set(p0,'FaceColor',normc([255;165;0]))    %orange 
%         ppos = patch(isosurface(X,Y,Z,S1,S1_75)) ; 
%         set(ppos,'FaceColor',normc([255;255;0]))  %yellow
% grid on ; 
% daspect([1 1 1]) ; 
% view(3); 
% axis equal ; 
% title({'\fontsize{14}Sig1TDEScript','\fontsize{8}Percentiles,Red=50th Orange=25th Yellow=75th'})
% hold off;
% subplot(1,3,3),pneg=patch(isosurface(X,Y,Z,S3,S3_50)) ;
%         set(pneg,'FaceColor',normc([148;148;148])); %grey
%         hold on
%         p0 = patch(isosurface(X,Y,Z,S3,S3_25)) ; 
%         set(p0,'FaceColor',normc([152;245;255]));   %blue    
%         ppos = patch(isosurface(X,Y,Z,S3,S3_75)) ; 
%         set(ppos,'FaceColor',normc([148;0;211])); %dark purp  
% grid on ; 
% daspect([1 1 1]) ; 
% view(3); 
% axis equal ; 
% title('Sig3TDEScript')
% title({'\fontsize{14}Sig3TDEScript','\fontsize{8}Percentiles,Grey=50th Blue=25th Dark purple=75th'})
% hold off;
% subplot(1,3,2),pneg=patch(isosurface(X,Y,Z,S2,S2_50)) ; %subplot(2,3,2),
%         set(pneg,'FaceColor',normc([192;255;62]));  %olive green
%         hold on
%         p0 = patch(isosurface(X,Y,Z,S2,S2_25)) ; 
%         set(p0,'FaceColor',normc([32;178;170]));    %turquoise
%         ppos = patch(isosurface(X,Y,Z,S2,S2_75)) ; 
%         set(ppos,'FaceColor',normc([224;102;255])); %purple
% grid on ; 
% daspect([1 1 1]) ; 
% view(3); 
% axis equal ; 
% hold off;
% %Linking all axes together to spin at the same time
% all_ha = findobj( figure_handle, 'type', 'axes', 'tag', '' );
% hlink = linkprop(all_ha,{'CameraPosition','CameraUpVector'});
% rotate3d on

