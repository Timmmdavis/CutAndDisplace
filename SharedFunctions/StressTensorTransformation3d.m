function [ P11,P22,P33,P12,P13,P23 ] = StressTensorTransformation3d(Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,C1,C2,C3)
% StressTensorTransformation3d: Converts 3d stress/strain tensors into a new
%                   coordinate system. Can also be used to convert tensors
%                   between spherical and Cartesian coordinates.
%                   Equations 6.88 - 6.93 of Pollard and Fletcher.
%                   This could be sped up with vector multiplication, see:
%                   StressTensorTransformation2d.m as an example. This
%                   would also allow it to work on vector that isnt col
%                   vects. 
%               
% usage #1: Just changing the coordinate axis orientation. 
% [ P11,P22,P33,P12,P13,P23 ] = StressTensorTransformation3d(Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,C1,C2,C3)
%
% usage #2: Convert spherical to Cart tensors:
% %If "th" is azimuth and "phi" elevation as in 
% %Davis et al 2017. Stress concentrations around voids. 
% [th,phi,r] = cart2sph(x,y,z);
% C1=[sin(phi).*cos(th),cos(phi).*cos(th),-sin(th) ];
% C2=[sin(phi).*sin(th),cos(phi).*sin(th),cos(th)  ];
% C3=[cos(phi)         ,-sin(phi)         ,zer     ];
% [ Sxx,Syy,Szz,Sxy,Sxz,Syz ] = StressTensorTransformation3d(Srr,Spp,Stt,Srp,Srt,Spt,C1,C2,C3)
%
% usage #3: Convert Cart to spherical tensors:
% %If "th" is azimuth and "phi" elevation as in 
% %Davis et al 2017. Stress concentrations around voids. 
% [th,phi,r] = cart2sph(x,y,z);
% C1=[sin(phi).*cos(th),sin(phi).*sin(th),cos(phi) ];
% C2=[cos(phi).*cos(th),cos(phi).*sin(th),-sin(phi)];
% C3=[-sin(th)         ,cos(th)          ,zer    ];
%[ Srr,Spp,Stt,Srp,Srt,Spt ] = StressTensorTransformation3d(Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,C1,C2,C3)
%
% Arguments: (input)
% Pxx,Pyy,Pzz
% Pxy,Pxz,Pyz       - 3D stress tensors (or strain).
%
% C1                - New Coordinate axis X (n*3) or (1*3)
%
% C2                - New Coordinate axis Y (n*3) or (1*3)
%
% C3                - New Coordinate axis Z (n*3) or (1*3)
%
% Arguments: (output)
% P11,P22,P33
% P12,P13,P23       - The new tensors. 
% 
%
% Example usage:
%
% %Creating points randomly scattered on a sphere
% n=500;
% theta=pi*rand(n,1);
% phi=2*pi*rand(n,1);
% [x,y,z]=sph2cart(theta,phi,ones(n,1));
% Dat=[x,y,z];
% %Create some arrays
% one=ones(size(x));
% zer=zeros(size(x));
% %Create new coords
% [th,phi,r] = cart2sph(x,y,z);
% C1=[sin(phi).*cos(th),sin(phi).*sin(th),cos(phi) ];
% C2=[cos(phi).*cos(th),cos(phi).*sin(th),-sin(phi)];
% C3=[-sin(th)         ,cos(th)          ,zer    ];
% %No shear tensors
% Pxx=one;
% Pyy=one*1.5;
% Pzz=one*2;
% Pxy=zer;
% Pxz=zer;
% Pyz=zer;
% %Transforming to polar tensors
% [ P11,P22,P33,P12,P13,P23 ] = StressTensorTransformation3d(Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,C1,C2,C3);
% %Transforming back to cart tensors
% C1=[sin(phi).*cos(th),cos(phi).*cos(th),-sin(th) ];
% C2=[sin(phi).*sin(th),cos(phi).*sin(th),cos(th)  ];
% C3=[cos(phi)         ,-sin(phi)         ,zer     ];
% [ P11,P22,P33,P12,P13,P23 ] = StressTensorTransformation3d(P11,P22,P33,P12,P13,P23,C1,C2,C3);
% %Drawing: 
% DrawS1S2S3Directions([zer,zer,zer,zer,zer,zer,P11,P22,P33,P12,P13,P23],x,y,z,'Scale',0.1)
% 
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen



%A quick low down:
%Quat is a matrix containing directions (quaternary). The first row is the first
%new coordinate direction we are transforming too. If this was a vector
%pointing along the X axis the row would be: [1,0,0] or [A11,A12,A13] as it
%only points along X and is 90 to the other two directions
%the second row is the 2nd direction. For the example above this would be
%the y axis, ie [0,1,0]
%etc
%Tensor is simply a matrix of the tensors in their original input coordinates

%If the new coordiante system is a single value we repeat it. 
if isscalar(C1(:,1)) || isscalar(C2(:,1)) || isscalar(C3(:,1)) %checks if its just 1 row
    C1=repmat(C1,size(Pxx(:,1)));
    C2=repmat(C2,size(Pxx(:,1)));
    C3=repmat(C3,size(Pxx(:,1)));
end  

%Extracting parts of the new coordinates. 
A11=C1(:,1);A12=C1(:,2);A13=C1(:,3);
A21=C2(:,1);A22=C2(:,2);A23=C2(:,3);
A31=C3(:,1);A32=C3(:,2);A33=C3(:,3);
    
%making sure everything is col vecs
[ Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,A11,A12,A13,A21,A22,A23,A31,A32,A33 ] = RowVecToCol( Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,A11,A12,A13,A21,A22,A23,A31,A32,A33 );

for i = 1:numel(Pxx)
%checking all vectors are 90 to each other
%this should give a result of 0 if there is an angle of 90 between the vectors 
DotPro12=A11(i).*A12(i)+A21(i).*A22(i)+A31(i).*A32(i);
DotPro13=A11(i).*A13(i)+A21(i).*A23(i)+A31(i).*A33(i);
DotPro23=A12(i).*A13(i)+A22(i).*A23(i)+A32(i).*A33(i);
    if any(round(DotPro12,14)) || any(round(DotPro13,14)) || any(round(DotPro23,14)) %rounding so no precision issues close to eps
        error('Direction cosines are not at 90 degrees to each other')
    end    
end

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
    Quat=[A11(i,:),A12(i,:),A13(i,:); %in cart dir cosines of the x axis
          A21(i,:),A22(i,:),A23(i,:); %in cart y axis
          A31(i,:),A32(i,:),A33(i,:)];%in cart z axis

    %Chucking the tensor together, if you have somehow managed to interpret/compute shear components of
    %stress from field observations you could put them in here too
    Tensor=[Pxx(i,:),Pxy(i,:),Pxz(i,:);
            Pxy(i,:),Pyy(i,:),Pyz(i,:);
            Pxz(i,:),Pyz(i,:),Pzz(i,:)];

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

end

