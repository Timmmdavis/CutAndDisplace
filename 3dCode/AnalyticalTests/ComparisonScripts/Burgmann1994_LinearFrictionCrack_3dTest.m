
%   Copyright 2017, Tim Davis, The University of Aberdeen%

%Test: 2d Planar fracture with increasing frictional profile: Comparison of
% slip distributions along a planar fracture with linearly increasing
% friction (sliding friction/frictional strength) from 0 at the centre to a
% set value at the tips.
% Also testing slip reduces with normal stress (mu). Havent yet found a
% analytical solution for this.
% 
% Solution: From Burgmann, R., Pollard, D.D. and Martel, S.J., 1994. Slip
% distributions on faults: effects of stress gradients, inelastic
% deformation, heterogeneous host-rock stiffness, and fault interaction.
% Journal of Structural Geology, 16(12), pp.1675-1690. This test calculates
% the analytical slip profile for a planar fracture with its friction
% increasing linearly to its tips. This fracture is subjected to a set
% shear stress. In the code a planar fracture lying along the X axis is
% loaded with a shear stress (Sxy) with 100 elements where each has its own
% set frictional strength. The shearing slip distribution is calculated.
% Checks the residual value is smaller than a threshold.  Res=Slip from
% solver ' Slip from solution.
% 
% Proof: This shows the frictional properties for each element in the code
% is working correctly. It shows with enough elements this results in a
% accurate slip profile for a fracture with friction. This does not need to
% test the resultant stress in the surrounding material as this is
% dependant on each slip contribution and is accurate as shown in the other
% tests.


% Setup: I have put this into 2 parts each with two subparts. Firstly I run
% a tall vertical fracture under Sxy and check it matches the analytical
% slip profile for its strike slip shearing. Secondly I check this fractures
% slip is reduced with a normal stress.
% Third and fourth I run the same two tests but for a flat lying fracture
% long in YY that is loaded under shearing Sxz. For this I check the dip
% slip shearing not strike slip. 








%%
%PART 1.1 Checking a vertical fracture matches the analytical frictional
%profile when loaded under Sxy.


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

 string='BladeFaultShort_500Tris_Vertical_NorthDip_0.5.ts'; %low sampling, still passes
 %string='BladeFaultShort_500Tris_Vertical_SouthDip_0.5.ts'; %low sampling, still passes
 %string='BladeFaultShort_2000Tris_Vertical_NorthDip_0.5.ts'; 
 %string='BladeFaultShort_2000Tris_Vertical_SouthDip_0.5.ts';
 [ Points,Triangles ] = GoCadAsciiReader( string );

 disp('Flattening surface in Y')
Points(:,3)=Points(:,3)*0; 

disp('Scaling on')
SclCrck=0.0001;
Points(:,2)=Points(:,2).*SclCrck;
Points(:,3)=Points(:,3).*SclCrck;
Points(:,4)=Points(:,4).*SclCrck;

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
	%Fictitious disp flag, if set to one then this forces the triangles with this to 0 displacement. 
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
[Sxx,Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,Mu,Sf,strain ] = CreateBlankVars;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option C = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. This option includes frictional contact
    %properties on the fault surface, elements cannot interpenetrate and
    %slip is reduced by the frictional parameters. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
Sxx = 0;       			
Syy = -0.05; 
Szz = 0; 
Sxy = 1;        			
Sxz = 0;
Syz = 0; 
% %Finding distance of points along X from 0,0. 
X=(MidPoint(:,1));
ne=size(MidPoint(:,1));
Mu  = 0.0;  Mu=repmat(Mu,ne);     %Coefficient of friction
%Creating a linear friction gradient across fractures X axis that increases
%with its X distance. As the fracture has a half length of 1 the max value
%will be 1 at the tips. 
Sf  = (abs(X))/SclCrck; 
%trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Sf);
Option='C'; %slip from uniform remote stress with friction

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Dss,Dds,Dn,~,~,~,~,~,~]=SlipCalculator3d(MidPoint,Sxx,...
 Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
 FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points);


%redefine driving stress
Sxy = -1;        			%Positive is right lateral movement  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[DssNegDr,~,~,~,~,~,~,~,~]=SlipCalculator3d(MidPoint,Sxx,...
 Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
 FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points);

   
	%%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%
	
Sxx = 0;       			
Syy = -0.05; 
Szz = 0;
Sxy = 1;        			
Sxz = 0;
Syz = 0; 
Option='B'; %slip No Fric


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
[DssNoFr,~,~,~,~,~,~,~,~]=SlipCalculator3d(MidPoint,Sxx,...
 Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
 FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing figures of slip distribution on the crack.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,Dss,Dds,Dn)


%Filtering values from the calculation as  fractures tips along Y have less
%slip so will not match the 2d analytical 'plane strain/stress'  profile 
FracMid=(MidPoint(:,3))<(1*SclCrck) & (MidPoint(:,3))>-(1*SclCrck); %Only finding slips in the central Y strip of fault
Profile=-Dss(FracMid); %Grabbing values of slip with friction
Profile2=-DssNoFr(FracMid); %Grabbing values of slip for no friction
Profile3=DssNegDr(FracMid); %Grabbing values of slip for neg friction
Srt=(MidPoint((FracMid),1)); %Grabbing X points within Y distance to sort on
C = [Srt,Profile,Profile2,Profile3];
jnk = sortrows(C); %Sorted along X axis for slips. 
%Grabbing values from filtered data
x=jnk(:,1); %Points along crack
PosFr=jnk(:,2);
NoFr=jnk(:,3);
NegFr=jnk(:,4);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP AA: Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%
%Linear increasing fricstion model of Burgmann Pollard and Martel 1994.
%Using Points from calculation as observation points. 
%%%%%%
Sg=1;
a = 1*SclCrck;  %Unit half length        
Sr =  abs(Sxy); %Driving stress (positive is extension)
G=mu; %ShearMod
nu=0.25;%Poissions ratio

%Pollard Burgmann_1994_eq5_UniformRemoteStress
[Slip_UniformRemote]=PollardSegall1987_FractureSlipProfile(G,nu,0,Sr,x,a);
%Pollard Burgmann_1994_eq14_LinearIncreasingFriction
[Slip_IncreasingFriction]=Burgmann1994_FractureLinearFrictionSlipProfile(G,nu,Sr,x,a);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual, plotting and checking for errors. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%For plots
%Grabbing the right hand side of the crack
rightcrack=x>0;
%Spacing and points we will interpolate on
RN = linspace(0,0.95*SclCrck,20);
%Making sure we have unique points or the interpolation is unhappy. 
[R2, ~] = uniquify(x(rightcrack));
%Interpolating the numerical results for plotting (Makes it easier to see
%errors)
ShearDispInt = interp1(R2,PosFr(rightcrack),RN,'pchip');
ShearDispNegFrInt = interp1(R2,NoFr(rightcrack),RN,'pchip');
ShearDispNoFrInt = interp1(R2,NegFr(rightcrack),RN,'pchip');
%Drawing the figure
figure;
hold on
%Plotting analytical profiles
plot(x(rightcrack),Slip_IncreasingFriction(rightcrack),'b','LineWidth',2.5);
plot(x(rightcrack),Slip_UniformRemote(rightcrack),'g','LineWidth',2.5);
%Plotting interpolated numerical disps every 20th halflength
scatter(RN,abs(ShearDispInt),24,'k','filled');
scatter(RN,abs(ShearDispNegFrInt),24,'k','filled');
scatter(RN,abs(ShearDispNoFrInt),24,'k','filled');
%Adding titles etc
title('Slip Distribution');
xlabel('Distance from crack centre')
ylabel({'Crack wall displacement'}) 
grid on
legend('show')
legend('Linear friction profile','No Friction','Numerical results')
%Making the plot nicer
titlesz=25;
fntsz=21;
ChangeFontSizes(fntsz,titlesz);


%Calculating residual
FrictionResidual=Slip_IncreasingFriction-abs(jnk(:,2));
FrictionResidualNegDirStress=Slip_IncreasingFriction-abs(jnk(:,4));

fprintf('MaxSlipDifference %i.\n',max(abs(FrictionResidual(:))))

%Printing error if the match is too far from an expected tolerance. 
if  max(abs(FrictionResidual)) > 0.12 || max(abs(FrictionResidualNegDirStress)) > 0.12 
    error('Not a good match to analytical solution')
else
    disp('Everything looks good, tolerance checks frictional slip residuals and flags differences to analytical profile above 0.11')
end

disp('Now checking mu changes strike slip disp when there is a compressional normal stress on the vertical fracture')   







%%
%PART 1.2 Checking a vertical fractures slip is reduced when I load this
%with a compressional normal stress along with shearing and give it a
%coefficient of friction.

%%


    %%%%%%%%%%%%%% 
    %StressInput - Friction 2, checking normal stress induces changes in
    %crack wall disp
	%%%%%%%%%%%%%%
    
Sxx = 0;       			%Positive is tension 
Syy = -2; 
Szz = 0;
Sxy = 1;        			%Positive is right lateral movement  
Sxz = 0;
Syz = 0; 
% %Finding distance of points along X from 0,0. 
ne=size(MidPoint(:,1));
Mu  = 0.2;  Mu=repmat(Mu,ne);     %Coefficient of friction
Sf  = zeros(size(Mu));  
%trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Sf);
Option='C'; %slip from uniform remote stress with friction

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Dss,Dds,Dn,~,~,~,~,~,~]=SlipCalculator3d(MidPoint,Sxx,...
 Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
 FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting the calculated slips.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,Dss,Dds,Dn)

MaxDisp=max(abs(Dss));

if  MaxDisp < 0.9*SclCrck || MaxDisp > 1*SclCrck
    error('Not a good match, the normal stress does not reduce the friction corrrectly')
else
    disp('Everything looks good. Mu and Normal stress relate')
end

clear a AdRs0 AdRs05 AdRs1 C DipSlipDisp DipSlipDispNegDr DipSlipDispNoFr DirPart FaceNormalVector Fdisp FracMid 
clear freesurface_height FrictionResidualNegDirStress G halfspace jnk lambda mu MidPoint Mu ne nu P1 P2 P3 parts Points pr Profile 
clear Profile2 Profile3 Sf Sg Slip_IncreasingFriction Slip_UniformRemote Sr Srt strain StrikeSlipDisp StrikeSlipDispNegDr 
clear StrikeSlipDispNoFr string Sxx Sxy Sxz Syy Syz Szz TensileSlipDisp TensileSlipDispNegDr TensileSlipDispNoFr Triangles x X














%%
%PART 2.1 Checking a flat fracture matches the analytical frictional
%profile when loaded under Sxz.
		
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: Import the fault surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Loading Gocad ascii data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %string='BladeFaultShort_500Tris_nodip.ts';
 string='BladeFaultShort_500Tris_Flat_WestDip_0.1.ts'; %low sampling, still passes
 %string='BladeFaultShort_500Tris_Flat_EastDip_0.1.ts'; %low sampling, still passes
 %string='BladeFaultShort_2000Tris_Flat_WestDip_0.1.ts'; 
 %string='BladeFaultShort_2000Tris_Flat_EastDip_0.1.ts'; 

 [ Points,Triangles ] = GoCadAsciiReader( string );

 
disp('Flattening surface in Z')
Points(:,4)=Points(:,4)*0; 

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
	%Fictitious disp flag, if set to one then this forces the triangles with this to 0 displacement. 
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
[Sxx,Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,Mu,Sf,strain ] = CreateBlankVars;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option C = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. This option includes frictional contact
    %properties on the fault surface, elements cannot interpenetrate and
    %slip is reduced by the frictional parameters. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
	%%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%
	
Sxx = 0;       			
Syy = 0; 
Szz = -0.05;
Sxy = 0;        			
Sxz = 1;
Syz = 0; 
Option='B'; 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[DssNoFr,DdsNoFr,DnNoFr,~,~,~,~,~,~]=SlipCalculator3d(MidPoint,Sxx,...
 Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
 FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points);

    %%%%%%%%%%%%%% 
    %StressInput - Friction
	%%%%%%%%%%%%%%
    
Sxx = 0;       			%Positive is compression 
Syy = 0; 
Szz = -0.05;
Sxy = 0;        			%Positive is right lateral movement  
Sxz = 1;
Syz = 0; 
% %Finding distance of points along X from 0,0. 
X=(MidPoint(:,1));
ne=size(MidPoint(:,1));
Mu  = 0.0;  Mu=repmat(Mu,ne);     %Coefficient of friction
%Creating a linear friction gradient across fractures X axis that increases
%with its X distance. As the fracture has a half length of 1 the max value
%will be 1 at the tips. 
Sf  = abs(X); 
%trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Sf);
Option='C'; %slip from uniform remote stress with no friction



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 [Dss,Dds,Dn,~,~,~,~,~,~]=SlipCalculator3d(MidPoint,Sxx,...
 Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
 FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points);

%redefine driving stress
Sxz = -1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [DssNegDr,DdsNegDr,DnNegDr,~,~,~,~,~,~]=SlipCalculator3d(MidPoint,Sxx,...
 Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
 FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing figures of slip distribution on the crack.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,Dss,Dds,Dn)


%Filtering values from the calculation as  fractures tips along Y have less
%slip so will not match the 2d analytical 'plane strain/stress'  profile 
FracMid=(MidPoint(:,2))<1 & (MidPoint(:,2))>-1; %Only finding slips in the central Y strip of fault
Profile=Dds(FracMid); %Grabbing values of slip with friction
Profile2=DdsNoFr(FracMid); %Grabbing values of slip for no friction
Profile3=DdsNegDr(FracMid); %Grabbing values of slip for no friction
Srt=(MidPoint((FracMid),1)); %Grabbing X points within Y distance to sort on
C = [Srt,Profile,Profile2,Profile3];
jnk = sortrows(C); %Sorted along X axis for slips. 
%Grabbing values from filtered data
x=jnk(:,1); %Points along crack
PosFr=jnk(:,2);
NoFr=jnk(:,3);
NegFr=jnk(:,4);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP AA: Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%
%Linear increasing fricstion model of Burgmann Pollard and Martel 1994.
%Using Points from calculation as observation points. 
%%%%%%
Sg=1;
a = 1;  %Unit half length        
Sr =  abs(Sxz); %Driving stress (positive is tension)
G=mu; %ShearMod
nu=0.25;%Poissions ratio

%Pollard Burgmann_1994_eq5_UniformRemoteStress
[Slip_UniformRemote]=PollardSegall1987_FractureSlipProfile(G,nu,0,Sr,x,a);
%Pollard Burgmann_1994_eq14_LinearIncreasingFriction
[Slip_IncreasingFriction]=Burgmann1994_FractureLinearFrictionSlipProfile(G,nu,Sr,x,a);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual, plotting and checking for errors. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Calculating residual
FrictionResidual=Slip_IncreasingFriction-abs(jnk(:,2));
FrictionResidualNegDirStress=Slip_IncreasingFriction-abs(jnk(:,4));

%%For plots
%Grabbing the right hand side of the crack
rightcrack=x>0;
%Spacing and points we will interpolate on
RN = linspace(0,0.95,20);
%Making sure we have unique points or the interpolation is unhappy. 
[R2, ~] = uniquify(x(rightcrack));
%Interpolating the numerical results for plotting (Makes it easier to see
%errors)
ShearDispInt = interp1(R2,PosFr(rightcrack),RN,'pchip');
ShearDispNegFrInt = interp1(R2,NoFr(rightcrack),RN,'pchip');
ShearDispNoFrInt = interp1(R2,NegFr(rightcrack),RN,'pchip');
%Drawing the figure
figure;
hold on
%Plotting analytical profiles
plot(x(rightcrack),Slip_IncreasingFriction(rightcrack),'b','LineWidth',2.5);
plot(x(rightcrack),Slip_UniformRemote(rightcrack),'g','LineWidth',2.5);
%Plotting interpolated numerical disps every 20th halflength
scatter(RN,abs(ShearDispInt),24,'k','filled');
scatter(RN,abs(ShearDispNegFrInt),24,'k','filled');
scatter(RN,abs(ShearDispNoFrInt),24,'k','filled');
%Adding titles etc
title('Slip Distribution');
xlabel('Distance from crack centre')
ylabel({'Crack wall displacement'}) 
grid on
legend('show')
legend('Linear friction profile','No Friction','Numerical results')
%Making the plot nicer
titlesz=25;
fntsz=21;
ChangeFontSizes(fntsz,titlesz);

fprintf('MaxSlipDifference %i.\n',max(abs(FrictionResidual(:))))

%Printing error if the match is too far from an expected tolerance. 
if  max(abs(FrictionResidual)) > 0.12 || max(abs(FrictionResidualNegDirStress)) > 0.12 
    error('Not a good match to analytical solution')
else
    disp('Everything looks good, tolerance checks frictional slip residuals and flags differences to analytical profile above 0.11')
end

disp('Now checking mu changes dip slip disp when there is a compressional normal stress on the flat lying fracture') 









%%
%PART 2.2 Checking a flat fractures slip is reduced when I load this
%with a compressional normal stress along with shearing and give it a
%coefficient of friction.


    %%%%%%%%%%%%%% 
    %StressInput - Friction 2, checking normal stress induces changes in
    %crack wall disp
	%%%%%%%%%%%%%%
    
Sxx = 0;       			%Positive is tension 
Syy = 0; 
Szz = -2;
Sxy = 0;        			%Positive is right lateral movement  
Sxz = 1;
Syz = 0; 
% %Finding distance of points along X from 0,0. 
ne=size(MidPoint(:,1));
Mu  = 0.2;  Mu=repmat(Mu,ne);     %Coefficient of friction
Sf  = zeros(size(Mu));  
%trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Sf);

Option='C'; %slip from uniform remote stress with no friction

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 [Dss,Dds,Dn,~,~,~,~,~,~]=SlipCalculator3d(MidPoint,Sxx,...
 Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
 FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing figures of slip distribution on the crack.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,Dss,Dds,Dn)

MaxDisp=max(abs(Dds));

if  MaxDisp < 0.9 || MaxDisp > 1
    error('Not a good match, the normal stress does not reduce the friction corrrectly')
else
    disp('Everything looks good. Mu and Normal stress relate')
end

