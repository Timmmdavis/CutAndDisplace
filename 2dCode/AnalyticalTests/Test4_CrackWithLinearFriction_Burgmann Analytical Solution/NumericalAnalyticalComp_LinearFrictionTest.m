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


clear;close all
cmap = colormap_cpt('Ccool-warm'); %Loads a colourmap to be used for figures. This is Kenneth Morelands diverging colourmap.
cmap2 = colormap_cpt('Ccool-warm2');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: import the fault surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading XY ascii data or manually creating fractures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Single flat user defined fracture
x=linspace(-1,1,500);                        %List of x points
y=zeros(1,numel(x));
Pointsxy=[x;y]';
mystruct.line1=(1:(length(Pointsxy(:,1)))-1);
Fnms=fieldnames(mystruct);
nf=numel(Fnms);clear Fnms

    %Checks and Creates Fdisp
chk=exist('Fdisp','var');
if chk==0 %exists in workspace
Fdisp=zeros(size(Pointsxy(:,1)));
end


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

%Elastic parameters
mu=1;            	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
E = mu*(2*(1+nu)) ;                     %Young's Modulus
lambda= E*nu/((1+nu)*(1-2*nu));   %Lamé's  constant,  Equation 8.27 Pollard

    %Reshaping the list of XY into usable variables
[ x,y,xe,ye,HalfLength,Beta,CosB,Points,NormAng,NUM,XBEG,XEND,YBEG,YEND]...
 = CreateElements2d( Pointsxy,mystruct );

    %Plotting figure of the the loaded fracture/fractures
hold on    
line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','r')
nx=cos(NormAng);
ny=cos((pi/2)-NormAng);
quiver(xe,ye,nx,ny)
title('fractures'), xlabel('x'), ylabel('y');axis equal
%Fixed elements shown as blue
line([Points(logical(Fdisp),1)';Points(logical(Fdisp),2)'],[Points(logical(Fdisp),3)';Points(logical(Fdisp),4)'],'color','b')
hold off

 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, Remote, Remote
%no opening, tractions on the elements or user imported slip value for each face. 

%Stress tensor outputs from this part of the toolbox i.e. Sxx Syy etc are simply the user defined remote stress
%components (boundary conditions). These are used later in the stress calculations. If strain is input these are converted to stress 
%using Hooke's law. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    
[Sxx,Syy,Sxy,Tn,Ts,Mu,Sf,strain ] = CreateBlankVars;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. Choose to define in stress or strain. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
    %%%%%%%%%%%%%% 
    %StressInput - Friction
	%%%%%%%%%%%%%%
       
%Input Stress
Sxx = 0; 					%Positive is tension
Syy = 0;
Sxy = 1;                   

%Increasing Friction Eq14 Burgmann, Fig 6 ratio
%This will be the friction at the tips
Sg=1;
Mu  = 0;     Mu=repmat(Mu,NUM,1);  %Coefficient of friction

%Creating a linear friction profile
Sfa=linspace(0,Sg,NUM/2);
Sf  = [fliplr(Sfa),Sfa]';   
Option='C'; %slip from uniform remote stress with friction

%Loading fault that lies in XY and calculating slip
[ShearDisp,TensileDisp,Sxx,Syy,Sxy]=SlipCalculator2d(x,y,xe,ye,HalfLength,Beta,Sxx,Syy,Sxy,Tn,Ts,nu,E,halfspace,NUM,NormAng,strain,Fdisp,Mu,Sf,Option);

Sxy = -1;                   %Positive is left lateral movement
%checking negative stress gives same result
[ShearDispNegFr,TensileDispNegFr,Sxx,Syy,Sxy]=SlipCalculator2d(x,y,xe,ye,HalfLength,Beta,Sxx,Syy,Sxy,Tn,Ts,nu,E,halfspace,NUM,NormAng,strain,Fdisp,Mu,Sf,Option);


%NO FRIC
Option='B'; %slip from uniform remote stress with friction
Sxy = 1;                   %Positive is left lateral movement
[ShearDispNoFr,TensileDispNoFr,Sxx,Syy,Sxy]=SlipCalculator2d(x,y,xe,ye,HalfLength,Beta,Sxx,Syy,Sxy,Tn,Ts,nu,E,halfspace,NUM,NormAng,strain,Fdisp,Mu,Sf,Option);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3 Drawing Figures: Drawing Figures of slip distribution on the crack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Plot displacement discontinuity magnitude for all elements
GraphSlipVsElNum( NUM,TensileDisp,ShearDisp )

%Creating directions and magnitudes of slip for plotting
PlotOpeningVsShearOnEls( NormAng,TensileDisp,ShearDisp,Points,x,y,HalfLength )

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP AA: Analytical solution, Burgmann Pollard and Martel 1994.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%
%Linear increasing friction model of Burgmann Pollard and Martel 1994.
%%%%%%

a = 1;  %Unit half length        
Sr = abs(Sxy); %Driving stress (positive is tension)
G=mu; %ShearMod
pr=0.25;%Poissions ratio

%Pollard Burgmann_1994_eq5_UniformRemoteStress
Slip_UniformRemote=(2*Sr).*((1-pr)/G).*(sqrt(a.^2-xe.^2));
%Pollard Burgmann_1994_eq14_LinearIncreasingFriction
Slip_IncreasingFriction=((1-pr)/G).*(2.*(((Sr.*(sqrt(a.^2-xe.^2))))-((Sg.*(1./pi)).*(((sqrt(a.^2-xe.^2))+(xe.^2./a).*acosh(a./xe))))));
Slip_IncreasingFriction=real(Slip_IncreasingFriction);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual,drawing figures and checking for errors. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%For plots
%Grabbing the right hand side of the crack
rightcrack=xe>0;
%Spacing and points we will interpolate on
RN = linspace(0,0.95,20);
%Making sure we have unique points or the interpolation is unhappy. 
[R2, ~] = uniquify(xe(rightcrack));
%Interpolating the numerical results for plotting (Makes it easier to see
%errors)
ShearDispInt = interp1(R2,ShearDisp(rightcrack),RN,'pchip');
ShearDispNegFrInt = interp1(R2,ShearDispNegFr(rightcrack),RN,'pchip');
ShearDispNoFrInt = interp1(R2,ShearDispNoFr(rightcrack),RN,'pchip');
%Drawing the figure
figure;
hold on
%Plotting analytical profiles
plot(xe(rightcrack),Slip_IncreasingFriction(rightcrack),'b','LineWidth',2.5);
plot(xe(rightcrack),Slip_UniformRemote(rightcrack),'g','LineWidth',2.5);
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
ShearDisp=abs(ShearDisp);
FrictionResidual=(abs(Slip_IncreasingFriction)-abs(ShearDisp));
FrictionResidualNegDirStress=abs(Slip_IncreasingFriction)-abs(ShearDispNegFr);
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


Mu  = 0.2;     Mu=repmat(Mu,NUM,1);  %Coefficient of friction
Sf  = zeros(size(Mu));  
Option='C';
[ShearDisp,TensileDisp,Sxx,Syy,Sxy]=SlipCalculator2d(x,y,xe,ye,HalfLength,Beta,Sxx,Syy,Sxy,Tn,Ts,nu,E,halfspace,NUM,NormAng,strain,Fdisp,Mu,Sf,Option);


MaxDisp=max(abs(ShearDisp));

if  MaxDisp < 0.89 || MaxDisp > 0.91
    error('Not a good match, the normal stress does not reduce the friction corrrectly')
else
    disp('Everything looks good. Mu and Normal stress relate')
end

%set(0,'DefaultFigureVisible','on') %Always on 2 start
