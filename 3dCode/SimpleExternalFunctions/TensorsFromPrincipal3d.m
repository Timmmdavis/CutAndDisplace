%Two horizontal components of stress
%Calculate cartesian stress tensor components from given computed/observed
%principal stresses
%This tensor can then be input as your driving stress without any 'dodgy'
%coordinate transforms/rotations having to be performed on the fault surface. 
%The way this is setup we presume we have the three principal stresses and
%two are horizontal, the other is vertical. We define the az and magnitude
%and the tensors are given that we input into the script

%   Copyright 2017, Tim Davis, The University of Aberdeen

clear
close all


%% DefineMagnitude
%Horizontal stresses mag
S11= -50;     %Magnitude of component 1
S22= -233.8000;	%Magnitude of component 2
%Vertical stress mag
S33=-233.8000; 	%Magnitude of component 3

%% DefineDirection
Dir1=90-24.8; 	%direction (az) of component 1 (degrees)
Dir2=90+Dir1; 	%direction (az) of component 2 (degrees)
%Dir 3 we presume is vertical so do not define this

%% Checking Dir 1&2 are 90 degrees apart
Test=round((asind(abs(sind(Dir1)))+asind(abs(sind(Dir2)))),8)==90; %rounding to 8 digits as precision can be an issue with this test
if Test==0
    error('Input Azimuth angles are not 90 degrees apart')
end

%% Calculating direction cosines
Ax1=abs(sind(Dir1)); %Angle away from X axis for Dir 1
Ay1=abs(cosd(Dir1)); %Angle away from Y axis for Dir 1
Az1=0;
if Dir1>180
	Ax1=-Ax1;
end
%Correcting sign for Y
if 	Dir1>90 && Dir1<270
	Ay1=-Ay1;
end

Ax2=abs(sind(Dir2)); %Angle away from X axis for Dir 1
Ay2=abs(cosd(Dir2)); %Angle away from Y axis for Dir 1
Az2=0;
if Dir2>180
	Ax2=-Ax2;
end
%Correcting sign for Y
if 	Dir2>90 && Dir2<270
	Ay2=-Ay2;
end

%Dipping directly down cosines
Ax3=0;
Ay3=0;
Az3=1;

%accumulating new axes to direction cosines
xAXIS=[Ax1,Ax2,Ax3];
yAXIS=[Ay1,Ay2,Ay3];
zAXIS=[Az1,Az2,Az3];

[ CartSxx,CartSyy,CartSzz,CartSxy,CartSxz,CartSyz ] = StressTensorTransformation3d(S11,S22,S33,0,0,0,xAXIS,yAXIS,zAXIS);


%Printing to command window, you can just cut and paste this into the
%script to use as the driving stress
disp('Calculated cartesian tensors:')
fprintf('Sxx= %i;\n',CartSxx)
fprintf('Syy= %i;\n',CartSyy) 
fprintf('Szz= %i;\n',CartSzz) 
fprintf('Sxy= %i;\n',CartSxy) 
fprintf('Sxz= %i;\n',CartSxz) 
fprintf('Syz= %i;\n',CartSyz) 

  

