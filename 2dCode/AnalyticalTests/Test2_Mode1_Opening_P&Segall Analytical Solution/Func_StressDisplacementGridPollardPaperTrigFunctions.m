  

function [SxxChange,SyyChange,SxyChange,u,v]=Func_StressDisplacementGridPollardPaperTrigFunctions(G,PR,Syy,Sxy,spacing,minx,maxx,cmap)
% Returns Cartesian displacements and stresses at grid points
% for a loaded crack using the formulation from Pollard, D. D., and P. Segall. 
% "Theoretical displacements and stresses near fractures in rock: with 
% applications to faults, joints, veins, dikes, and solution surfaces." 
% Fracture mechanics of rock 277.349 (1987): 277-349.
% This uses the solutions based on  trigonometry rather than complex
% variables as it is easier to understand. 

%   Copyright 2017, Tim Davis, The University of Aberdeen

%Assuming the crack lies in the X axis
%G shear mod
%PR nu (Poisson's ratio)
%Syy driving stress normal to the crack
%Sxy driving stress shearing the crack
% minx spacing maxx make a square grid around the crack a is defined below
% so check this first then make your grid
%cmap is a colourmap file, you wont need this if you don’t make the figures at the
%base

%Setting up grid to analyse stresses and displacements on. It doesn't need
%to be a grid but this is a decent way to visulise this. 
% % spacing=spacing;
% % minx=minx; maxx=maxx;
[x,y] = meshgrid(minx:spacing:maxx); %large grid
x=x(:);
y=y(:);

% Unit half-length displacement discontinuity. This is going to have its 
% centre at 0,0. 
a = 1;  
%Defining stress. Extension positive "Engineering convention"
Syy=Syy;
Sxy=Sxy;
Syz=0; %Out of Plane stresses all done at the bottom of code
Sxz=0; %Out of Plane stresses all done at the bottom of code

% Now calculating the locations and angles of each point relative to the
% two ends (positive and negative) and the centre of the fracture. 
% See figure 8.3

% Centre of the crack
r = sqrt(x.^2 + y.^2);      %   Array giving each points radial distance from the centre. 
theta = atan2(y,-x);        %   Angle measured from the x axis. 
theta=abs(theta-pi);        %   making this start from x axis and finish at 2piR.


% Positive end of the crack
xp = x-a;
rp = sqrt(xp.^2 + y.^2);
tp = atan2(y,-xp);          %   theta here is not measured away from positve X
tpp=abs(tp-pi);             %   making this start from x axis and finish at 2piR as described in paper. 

xn = x+a;
rn = sqrt(xn.^2 + y.^2);
tn = atan2(y,-xn);          %   theta here is not measured away from positve X
tnn=abs(tn-pi);             %   making this start from x axis and finish at 2piR. 


% Creating the well used variables that are often used
R=sqrt(rp.*rn);     %as described in pollard 8.31a
THETA=(tpp+tnn)/2;  %as described in pollard 8.31a
R1=R.^-1;R3=R.^-3;  %Radius


% Calculate displacements eq 8.33a and 8.33b Pollard and Segall 1984
% These are seperated into seperate components based on related driving stress so the equation as a whole is easier to read 
vsyy=Syy.*(2*(1-PR)*(R.*sin(THETA)-r.*sin(theta))-r.*sin(theta).*(r.*R1.*cos(theta-THETA)-1));
usyy=Syy.*((1-(2*PR))*(R.*cos(THETA)-r.*cos(theta))-r.*sin(theta).*(r.*R1.*sin(theta-THETA)));
usxy=Sxy.*(2*(1-PR)*(R.*sin(THETA)-r.*sin(theta))+r.*sin(theta).*(r.*R1.*cos(theta-THETA)-1));
vsxy=(Sxy.*((1-(2*PR))*(R.*cos(THETA)-r.*cos(theta))+r.*sin(theta).*(r.*R1.*sin(theta-THETA)))).*-1; %Why do I need to do this?

v=(vsyy+vsxy)/2*G;%The 2/G corresponds to the left hand side of the eq 8.33 that i have pushed over
u=(usyy+usxy)/2*G;

% Calculate the stresses from eq 8.44 onwards from Pollard and Segall 1984
% Note this contains no elastic constants related to rigidity. The
% resultant stress related to a fracture and is only dependant on input stress magntitude
% on the fracture surface (this stress must already alude to the strength of the elastic solid) 
SyyChange_syy= Syy.*(r.*R1.*cos(theta-THETA)-1+a^2.*r.*R3.*sin(theta).*sin(3*THETA));   %8.44a
SyyChange_sxy= Sxy.*(a^2.*r.*R3.*sin(theta).*cos(3*THETA));
SyyChange= SyyChange_syy + SyyChange_sxy;%Can add remote here if want to find total stress, see eq

SxyChange_sxy= Sxy.*(r.*R1.*cos(theta-THETA)-1-a^2.*r.*R3.*sin(theta).*sin(3*THETA));   %8.44b
SxyChange_syy= Syy.*((a^2).*r.*R3.*sin(theta).*cos(3*THETA));                             
SxyChange=SxyChange_syy + SxyChange_sxy; %Can add remote here if want to find total stress, see eq

SxxChange_syy= Syy.*(r.*R1.*cos(theta-THETA)-1-a^2.*r.*R3.*sin(theta).*sin(3*THETA));   %8.44c 
SxxChange_sxy= Sxy.*(2.*r.*R1.*sin(theta-THETA)-a^2.*r.*R3.*sin(theta).*cos(3*THETA));
SxxChange= SxxChange_syy + SxxChange_sxy; %Can add remote here if want to find total stress, see eq


%Reshaping to Col Vec so this can be plotted
num=(sqrt(numel(x)));
x=reshape(x,num,num);
y=reshape(y,num,num);
u=reshape(u,num,num);
v=reshape(v,num,num);
SyyChange=reshape(SyyChange,num,num);
SxyChange=reshape(SxyChange,num,num);
SxxChange=reshape(SxxChange,num,num);

%Octave ContourF doesn't like 'nans', replacing with '-infs'
%MATLAB unsuprisingly handles this fine
SyyChange(isnan(SyyChange)) = -inf;    SxxChange(isnan(SxxChange)) = -inf;    SxyChange(isnan(SxyChange)) = -inf;   

if Sxy+Syy>0.5
    %Drawing Figures
    %figure;contourf(x,y,u);title('displacement ux');%caxis([-0.25 0.2]) match twodd
    %figure;contourf(x,y,v);title('displacement uy');%caxis([-0.7 0.6]);match twodd
    figure;colormap(cmap),contourf(x,y,SyyChange);title('Syy Change, P&S solution'),colorbar;
    figure;colormap(cmap),contourf(x,y,SxxChange);title('Sxx Change, P&S solution'),colorbar;
    figure;colormap(cmap),contourf(x,y,SxyChange);title('Sxy Change, P&S solution'),colorbar;
    figure;quiver(x,y,u,v);title('displacement, P&S solution');
end

if Syz>0 %Doing if user has chosen to run with out of plane stress
    % Calculate the stresses from eq 8.44 onwards from Pollard and Segall 1984
    % Note this contains no elastic constants related to rigidity. The
    % resultant stress related to a fracture and is only dependant on input stress magntitude
    % on the fracture surface (this stress must already alude to the strength of the elastic solid) 
    SyzChange_syz= Syz.*(r.*R1.*cos(theta-THETA)-1);                                        %8.44d
    SyzChange= SyzChange_syz;   %Can add remote here if want to find total stress, see eq
    SyzChange=reshape(SyzChange,num,num);
    figure;colormap(cmap),contourf(x,y,SyzChange);title('Syz Change');
end

if Sxz>0 %Doing if user has chosen to run with out of plane stress
    % Calculate the stresses from eq 8.44 onwards from Pollard and Segall 1984
    % Note this contains no elastic constants related to rigidity. The
    % resultant stress related to a fracture and is only dependant on input stress magntitude
    % on the fracture surface (this stress must already alude to the strength of the elastic solid) 
    SxzChange_sxz= Sxz.*(r.*R1.*sin(theta-THETA)-1);                                        %8.44d
    SxzChange= SxzChange_sxz;   %Can add remote here if want to find total stress, see eq
    SxzChange=reshape(SxzChange,num,num);
    figure;colormap(cmap),contourf(x,y,SxzChange);title('Sxz Change');
end

% SxxChange=SxxChange(:);
% SyyChange=SyyChange(:);
% SxyChange=SxyChange(:);
