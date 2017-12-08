function [CartSxx,CartSyy,CartSxy]=TensorFromPrincipal2d( S1,S2,DirS1 )
% TensorsFromPrincipal2d: Function that computes the Cartesian stress
%                   tensor from the magnitude and directions of defined
%                   principal stresses. Use to compute tensors for the BEM
%                   boundary conditions. This prints these to your console.
%
%               
% usage #1: Prints results to cmd window
% TensorFromPrincipal2d( S1,S2,DirS1 )
%
% usage #2: Outputs these for futher use
% [CartSxx,CartSyy,CartSxy]=TensorFromPrincipal2d( S1,S2,DirS1 )
%
% Arguments: (input)
%     S1,S2       - The magnitudes of the principal stresses. 
%
%     DirS1       - The azimuth of the S1 component (degrees). 
%                  (angle away from N clockwise)
%
% Arguments: (output)
% CartSxx,CartSyy,CartSxy - The magnitudes of the calculated Cartesian
%                          stresses.
%
% Example usage:
%
%  S1=0;
%  S2=-2;
%  DirS1=45;
%  TensorFromPrincipal2d( S1,S2,DirS1 );
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%% DefineMagnitude
%Horizontal stresses mag
Sz1= S1;     %Magnitude of component 1
Sz2= S2;	 %Magnitude of component 2


%% DefineDirection
Dir1=90-DirS1; 	%direction (az) of component 1 (degrees)
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

if nargout==0
    %Printing to command window, you can just cut and paste this into the
    %script to use as the driving stress
    disp('Calculated cartesian tensors:')
    fprintf('Sxx= %i;\n',CartSxx)
    fprintf('Syy= %i;\n',CartSyy) 
    fprintf('Sxy= %i;\n',CartSxy) 
else
    %do nothing, you get these out the script and probably dont need these
    %displaying
end
  

