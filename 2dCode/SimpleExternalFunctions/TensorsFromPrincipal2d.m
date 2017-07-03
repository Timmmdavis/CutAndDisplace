%Two horizontal components of stress
%Calculate cartesian stress tensor components from given computed/observed
%principal stresses
%This tensor can then be input as your driving stress without any 'dodgy'
%coordinate transforms/rotations having to be performed on the fault surface. 
%The way this is setup we presume we have the two principal stresses in plane
%We define the 'azimuth' (angle away from 0 clockwise) and magnitude
%and the tensors are given that we input into the script

%   Copyright 2017, Tim Davis, The University of Aberdeen
clear
close all


%% DefineMagnitude
%Horizontal stresses mag
Sz1= -50;     %Magnitude of component 1
Sz2= -233.8000;	%Magnitude of component 2


%% DefineDirection
Dir1=90-24.8; 	%direction (az) of component 1 (degrees)
Dir2=90+Dir1; 	%direction (az) of component 2 (degrees)


%% Checking Dir 1&2 are 90 degrees apart
Test=round((asind(abs(sind(Dir1)))+asind(abs(sind(Dir2)))),8)==90; %rounding to 8 digits as precision can be an issue with this test
if Test==0
    error('Input Azimuth angles are not 90 degrees apart')
end

%% Calculating direction cosines
Ax1=abs(sind(Dir1)); %Angle away from X axis for Dir 1
Ay1=abs(cosd(Dir1)); %Angle away from Y axis for Dir 1
if Dir1>180
	Ax1=-Ax1;
end
%Correcting sign for Y
if 	Dir1>90 && Dir1<270
	Ay1=-Ay1;
end

Ax2=abs(sind(Dir2)); %Angle away from X axis for Dir 1
Ay2=abs(cosd(Dir2)); %Angle away from Y axis for Dir 1
if Dir2>180
	Ax2=-Ax2;
end
%Correcting sign for Y
if 	Dir2>90 && Dir2<270
	Ay2=-Ay2;
end

%checking all vectors are 90 to each other
%this should give a result of 0 if there is an angle of 90 between the vectors 
DotPro12=Ax1.*Ax2+Ay1.*Ay2;
if any(round(DotPro12,14)) %rounding as we end up with tiny values that ruin this check
    error('Direction cosines are at not 90 degrees to each other')
end    

CosAx=Ax1;
CosAy=Ay1;

[ CartSxx,CartSyy,CartSxy ] = StressTensorTransformation2d(Sz1,Sz2,0,CosAx,CosAy );


%Printing to command window, you can just cut and paste this into the
%script to use as the driving stress
disp('Calculated cartesian tensors:')
fprintf('Sxx= %i;\n',CartSxx)
fprintf('Syy= %i;\n',CartSyy) 
fprintf('Sxy= %i;\n',CartSxy) 

  

