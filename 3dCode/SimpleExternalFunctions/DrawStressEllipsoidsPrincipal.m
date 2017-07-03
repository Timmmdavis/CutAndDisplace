function  DrawStressEllipsoidsPrincipal(StressOutput,X,Y,Z,Triangles,Points)
%Draws the principal stress ellipsoids, pretty heavy function can be quite
%slow. I have added a catch if you are putting too much data into this.
% Input 'StressOutput' is from the stress calculation 12*n vector, strain then
% stress tensors. 
% Input 'XYZ' is the locations of the points we have the tensors for

% Things that can but don’t need to be supplied:
% Triangles is the triangles the code used
% Points the list of points 

%   Copyright 2017, Tim Davis, The University of Aberdeen

% Properties that need to be manually changed: 
%Manual scaling property, changes line and ellipsoid lengths
Scl=5000; %1
%Amount of points/faces that make up each ellipsoid in plot, the points are
%smpl^2. Reduce this no if figure is taking too long to draw
smpl=15;

figure;

%Loading external Octave packages, Normc etc are in here
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
pkg load general
pkg load miscellaneous
elseif isOctave==0
%do nothing
end

if numel(X)>1000
    disp('The "DrawPrincipalStresses" function needs to draw a lot of points, its slow.')
    disp('I recommend using only a few at first to get the parameters right.')
    disp('Leaving func, too many points')
    return
end


%Forcing XYZ to column vecs 
Xcv=X(:); %col vecs
Ycv=Y(:); %col vecs
Zcv=Z(:); %col vecs

%Extracting stress
Sxx=StressOutput(:,7);
Syy=StressOutput(:,8);
Szz=StressOutput(:,9);
Sxy=StressOutput(:,10);
Sxz=StressOutput(:,11);
Syz=StressOutput(:,12);
%Calculating principal stresses
[S1,S2,S3,S1dir,S2dir,S3dir] = EigCalc3d(Sxx(:),Syy(:),Szz(:),Sxy(:),Sxz(:),Syz(:));

%Calculating principal strains to calculate dilatation. (if vol change is
%positive or negative for drawing
Exx=StressOutput(:,1);
Eyy=StressOutput(:,2);
Ezz=StressOutput(:,3);
Exy=StressOutput(:,4);
Exz=StressOutput(:,5);
Eyz=StressOutput(:,6);
[E1,E2,E3] = EigCalc3d(Exx(:),Eyy(:),Ezz(:),Exy(:),Exz(:),Eyz(:));
Dilatation=E1+E2+E3; %(Change/Orig vol)
Expanded=Dilatation(:)>0; %flag, positive is expansion

%Now drawing first subplot
%This is a figure showing vectors that are the principal directions.
%Coloured for S1S2S3 and sign, see fig title for col details 
for i=1:numel(S1)
%Eq 2.23, Pollard, arranging cosines in table
Quat=[S1dir(i,1),S2dir(i,1),S3dir(i,1); %x
      S1dir(i,2),S2dir(i,2),S3dir(i,2); %y
      S1dir(i,3),S2dir(i,3),S3dir(i,3)];%z
%Creating the ellipsoid
[x,y,z] = ellipsoid(0,0,0,S1(i)*Scl,S2(i)*Scl,S3(i)*Scl,smpl); %Very low sampling
%Putting these in col vecs inside var to apply the transformation 
pnts=[x(:),y(:),z(:)];
%Preallocating array for loop
NwPnts=zeros(size(pnts));
%Running loop
for j=1:numel(x)
NwPnts(j,:)=(Quat*pnts(j,:)')'; %Eq 2.24 Pollard and Fletcher
end
%Reshaping for drawing surf plot
xn=reshape(NwPnts(:,1),size(x));
yn=reshape(NwPnts(:,2),size(x));
zn=reshape(NwPnts(:,3),size(x));
%Moving the points to the correct location in space
xn=xn+Xcv(i,:);
yn=yn+Ycv(i,:);
zn=zn+Zcv(i,:);

[ C ] = RGB2Colour(0,191,255);

%Drawing a single ellipsoid then setting the colour based on the total volume change sign at this
%infinitesimal point. 
if Expanded(i)==1
    hSurface=surf(xn, yn, zn);
    if  isOctave==1
    set(hSurface,'FaceColor','red','EdgeColor',([0.3;0.3;0.3]));%lighting not setup for surfaces in octave
    elseif isOctave==0
    set(hSurface,'FaceColor','red','FaceLighting','gouraud','EdgeColor',([0.3;0.3;0.3])); 
    end
else
    
    hSurface=surf(xn, yn, zn);
    if  isOctave==1
    set(hSurface,'FaceColor',C,'EdgeColor',([0.3;0.3;0.3]));%lighting not setup for surfaces in octave
    elseif isOctave==0
    set(hSurface,'FaceColor',C,'FaceLighting','gouraud','EdgeColor',([0.3;0.3;0.3])); 
    end
end    
hold on
end
xlabel('x'); ylabel('y');zlabel('z');axis equal ;
%Only setting tranparency if MATLAB
if  isOctave==1
elseif isOctave==0
alpha(.8);%transparency value
end
title({'\fontsize{14}Stress ellipsoids','\fontsize{8}Blue, net compressional volume change, Red, Net extensional volume change'})

if nargin>4
hold on
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',(.2),'FaceColor', [0.5 0 0.9 ]);
hold off
end

end

