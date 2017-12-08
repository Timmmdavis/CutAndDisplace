% Test: Planar fracture with increasing frictional profile: Comparison of
% slip distributions along a planar fracture with linearly increasing
% friction (sliding friction/frictional strength) from 0 at the centre to a
% set value at the tips.
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

%   Copyright 2017, Tim Davis, The University of Aberdeen

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
     
%Single flat user defined fracture
x=linspace(-1,1,200);                        %List of x points
y=zeros(1,numel(x));
Pointsxy=[x;y]';
mystruct.line1=(1:(length(Pointsxy(:,1)))-1);
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

%Reshaping the list of XY into usable variables
[MidPoint,HalfLength,P1,P2,LineNormalVector,Lines,Points,Fdisp  ]...
 = CreateElements2d( Pointsxy,mystruct,Fdisp );

%Drawing fracture
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
    %Option C = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. This option includes frictional contact
    %properties on the fault surface, elements cannot interpenetrate and
    %slip is reduced by the frictional parameters. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%
       
%Input Stress
Sxx = 0; 					%Positive is tension
Syy = 0;
Sxy = 1;                   

%Increasing Friction Eq14 Burgmann, Fig 6 ratio
%This will be the friction at the tips
Sg=1;
Mu  = 0;     Mu=repmat(Mu,numel(HalfLength),1);  %Coefficient of friction

%Creating a linear friction profile
Sfa=linspace(0,Sg,numel(HalfLength)/2);
Sf  = [fliplr(Sfa),Sfa]';   
Option='C'; %slip from uniform remote stress with friction

%Loading fault that lies in XY and calculating slip
[Ds,Dn,Sxx,Syy,Sxy]=SlipCalculator2d(...
    MidPoint,HalfLength,Sxx,Syy,Sxy,Tn,Ts,nu,E,halfspace,LineNormalVector,strain,Fdisp,Mu,Sf,Option);

Sxy = -1;                   %Positive is left lateral movement
%checking negative stress gives same result
[DsNegFr,DnNegFr,Sxx,Syy,Sxy]=SlipCalculator2d(...
    MidPoint,HalfLength,Sxx,Syy,Sxy,Tn,Ts,nu,E,halfspace,LineNormalVector,strain,Fdisp,Mu,Sf,Option);

%No friction
Option='B'; %slip from uniform remote stress with friction
Sxy = 1;                   %Positive is left lateral movement
[DsNoFr,DnNoFr,Sxx,Syy,Sxy]=SlipCalculator2d(...
    MidPoint,HalfLength,Sxx,Syy,Sxy,Tn,Ts,nu,E,halfspace,LineNormalVector,strain,Fdisp,Mu,Sf,Option);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing figures of slip distribution on the crack.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

% Plot displacement discontinuity magnitude for all elements
GraphSlipVsElNum( Dn,Ds )

%Creating directions and magnitudes of slip for plotting
PlotOpeningVsShearOnEls( LineNormalVector,Dn,Ds,P1,P2,MidPoint,HalfLength )

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP AA: Analytical solution, Burgmann Pollard and Martel 1994.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%
%Linear increasing friction model of Burgmann Pollard and Martel 1994.
%%%%%%
   
Sr = abs(Sxy); %Driving stress (positive is tension)

%Pollard Burgmann_1994_eq5_UniformRemoteStress
[Slip_UniformRemote]=PollardSegall1987_FractureSlipProfile(mu,nu,0,Sr,MidPoint(:,1),1);
%Pollard Burgmann_1994_eq14_LinearIncreasingFriction
[Slip_IncreasingFriction]=Burgmann1994_FractureLinearFrictionSlipProfile(mu,nu,Sr,MidPoint(:,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual,drawing figures and checking for errors. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%For plots
%Grabbing the right hand side of the crack
rightcrack=MidPoint(:,1)>0;
%Spacing and points we will interpolate on
RN = linspace(0,0.95,20);
%Making sure we have unique points or the interpolation is unhappy. 
[R2, ~] = uniquify(MidPoint(rightcrack,1));
%Interpolating the numerical results for plotting (Makes it easier to see
%errors)
ShearDispInt      = interp1(R2,Ds(rightcrack),RN,'pchip');
ShearDispNegFrInt = interp1(R2,DsNegFr(rightcrack),RN,'pchip');
ShearDispNoFrInt  = interp1(R2,DsNoFr(rightcrack),RN,'pchip');
%Drawing the figure
figure;
hold on
%Plotting analytical profiles
plot(MidPoint(rightcrack,1),Slip_IncreasingFriction(rightcrack),'b','LineWidth',2.5);
plot(MidPoint(rightcrack,1),Slip_UniformRemote(rightcrack),'g','LineWidth',2.5);
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
Ds=abs(Ds);
FrictionResidual=(abs(Slip_IncreasingFriction)-abs(Ds));
FrictionResidualNegDirStress=abs(Slip_IncreasingFriction)-abs(DsNegFr);
MaxAnalyticalShear=max(abs(Slip_IncreasingFriction));
SlipResPerc=(100/MaxAnalyticalShear).*FrictionResidual; %slip residual % relative to max slip from analyical

fprintf('Maximum Residual as a percent of the maximum analytical shear value %i.\n',max(abs(SlipResPerc)))

if  max(abs(FrictionResidual)) > 0.12 || max(abs(FrictionResidualNegDirStress)) > 0.12
    error('Not a good match to analytical solution')
else
    disp('Everything looks good, checked the max % difference is below 0.5')
end


%%
%PART 1.2 Checking a fractures slip is reduced when I load this
%with a compressional normal stress along with shearing and give it a
%coefficient of friction.

%%
    %%%%%%%%%%%%%% 
    %StressInput - Friction 2, checking normal stress induces changes in
    %crack wall disp
	%%%%%%%%%%%%%%
    
Sxx = 0;       			%Positive is tension 
Syy = -2; 
Sxy = 1;        			


Mu  = 0.2;     Mu=repmat(Mu,numel(HalfLength),1);  %Coefficient of friction
Sf  = zeros(size(Mu));  
Option='C';
[Ds,Dn,Sxx,Syy,Sxy]=SlipCalculator2d(...
    MidPoint,HalfLength,Sxx,Syy,Sxy,Tn,Ts,nu,E,halfspace,LineNormalVector,strain,Fdisp,Mu,Sf,Option);

MaxDisp=max(abs(Ds));

if  MaxDisp < 0.89 || MaxDisp > 0.91
    error('Not a good match, the normal stress does not reduce the friction corrrectly')
else
    disp('Everything looks good. Mu and Normal stress relate')
end

%set(0,'DefaultFigureVisible','on') %Always on 2 start
