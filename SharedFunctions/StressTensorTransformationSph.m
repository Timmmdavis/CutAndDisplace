function [ Srr,Spp,Stt,Srp,Srt,Spt  ] = StressTensorTransformationSph(Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,Az,El)
%StressTensorTransformation Converts 3d stress tensors into a spherical
%coordinate system. 

%Pxx,Pyy,Pzz,Pxy,Pxz,Pyz input cart tensors
%Az,El from MATLAB func : [Az,El,r] = cart2sph(x,y,z)
%Srr,Spp,Stt,Srp,Srt,Spt 

%Pxx,Pyy,Pxy etc are the tensors in the current coordinate system
%A11,A12,A13,A21,A22,A23,A31,A32,A33 are the direction cosines of the old
%coordinate system

%Az and el, see dir in:  https://uk.mathworks.com/help/MATLAB/ref/cart2sph.html 

%   Copyright 2017, Tim Davis, The University of Aberdeen
%Equation 6.95 & 6.96, Pollard and Fletcher.
%Function works perfectly fine on col vectors of stress. 

%A quick low down:
%Quat is a matrix containing directions (quaternary). The first row is the first
%new coordinate direction we are transforming too. If this was a vector
%pointing along the X axis the row would be: [1,0,0] or [A11,A12,A13] as it
%only points along X and is 90 to the other two directions
%the second row is the 2nd direction. For the example above this would be
%the y axis, ie [0,1,0]
%etc
%Tensor is simply a matrix of the tensors in their original input coordinates


%making sure everything is col vecs
[Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,Az,El] = RowVecToCol( Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,Az,El );

%http://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Polar_Coords/Polar_Coords.htm
Th=90-rad2deg(El);
Th=deg2rad(Th);
Phi=Az;



%Preallocating blank arrays of right size 
P11=zeros(size(Pxx));   %in cart xx
P22=P11;                %in cart yy
P33=P11;                %in cart zz
P12=P11;                %in cart xy
P13=P11;                %in cart xz
P23=P11;                %in cart yx


for i=1:numel(Pxx)

% Performing calculation for cartesian stresses
%Eq 2.23, Pollard, arranging cosines in table
Quat=[sin(Th(i,:)).*cos(Phi(i,:)),    sin(Th(i,:)).*sin(Phi(i,:)),    cos(Th(i,:)); %in cart x axis
      cos(Th(i,:)).*cos(Phi(i,:)),    cos(Th(i,:)).*sin(Phi(i,:)),    -sin(Th(i,:)); %in cart y axis
      -sin(Phi(i,:)),                 cos(Phi(i,:)),               0];%in cart z axis

%Chucking the tensor together, if you have somehow managed to interpret/compute shear components of
%stress from field observations you could put them in here too
%Tensor=[Pxx(:,i),Pxy(:,i),Pxz(:,i);Pxy(:,i),Pyy(:,i),Pyz(:,i);Pxz(:,i),Pyz(:,i),Pzz(:,i)];
Tensor=[Pxx(i,:),Pxy(i,:),Pxz(i,:);Pxy(i,:),Pyy(i,:),Pyz(i,:);Pxz(i,:),Pyz(i,:),Pzz(i,:)];
%http://continuummechanics.org/stressxforms.html
%Computing rotation
CartStress=Quat*Tensor*Quat';

% Extracting variables 
P11(i,:)=CartStress(1,1);
P22(i,:)=CartStress(2,2);
P33(i,:)=CartStress(3,3);
P12(i,:)=CartStress(1,2);
P13(i,:)=CartStress(1,3);
P23(i,:)=CartStress(2,3);
end

Srr=P11;
Spp=P22;
Stt=P33;
Srp=P12;
Srt=P13;
Spt=P23;

end

