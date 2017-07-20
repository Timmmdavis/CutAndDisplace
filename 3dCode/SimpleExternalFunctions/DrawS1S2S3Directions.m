function  DrawS1S2S3Directions(StressOutput,X,Y,Z,Triangles,Points )
%Draws the pricipal stress directions as vectors scaled 
%and coloured dependant on magnitude sign and direction:

% Input 'StressOutput' is from the stress calculation 12*n vector, strain then
% stress tensors. i.e.. exx,eyy,ezz,exy,exz,eyz,Sxx,Syy,Szz,Sxy,Sxz,Syz
% Input 'XYZ' is the locations of the points we have the tensors for

% Things that can but do not need to be supplied:
% Triangles is the triangles the code used
% Points the list of points 

%S1,Ext red,Comp yellow, S2,Ext green,Comp olivegreen, S3,Ext blue,Comp grey

% Properties that need to be manually changed: 
% Manual scaling property, changes line and ellipsoid lengths
Scl=5000; %1
% Line width property, MATLAB default is .5
Lw=2;   %5

%   Copyright 2017, Tim Davis, The University of Aberdeen

%Loading external Octave packages, Normc etc are in here
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
pkg load general
pkg load miscellaneous
elseif isOctave==0
%do nothing
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


%Now drawing first subplot
%This is a figure showing vectors that are the principal directions.
%Coloured for S1S2S3 and sign, see fig title for col details 
figure;
%Scaling directions by thier magnitudes
S1dirt=(bsxfun(@times,S1dir,abs(S1))); 
S2dirt=(bsxfun(@times,S2dir,abs(S2))); 
S3dirt=(bsxfun(@times,S3dir,abs(S3))); 
%Flagging negative and positive stresses, colour changes dependant on sign
S1Flag=S1>0; %Above 0 is extensional in this conv
S2Flag=S2>0; %Above 0 is extensional in this conv
S3Flag=S3>0; %Above 0 is extensional in this conv
%Extracting XYZ of the scaled vectors
S1x=S1dirt(:,1);S2x=S2dirt(:,1);S3x=S3dirt(:,1);
S1y=S1dirt(:,2);S2y=S2dirt(:,2);S3y=S3dirt(:,2);
S1z=S1dirt(:,3);S2z=S2dirt(:,3);S3z=S3dirt(:,3);
%Makes sure there is no scaling going on between holds, have to use manual scaling.
%MATLAB usually autoscales in a quiver otherwise
S=0;

hold on
%If loop, if its just pos values doing the neg one too ruins the axis lims
%etc
if sum(S1Flag) == numel(S1Flag)
%Drawing the positive part of S1 (extensional)
a = quiver3(Xcv(S1Flag) ,Ycv(S1Flag) ,Zcv(S1Flag) ,(S1x(S1Flag))*Scl, (S1y(S1Flag))*Scl, (S1z(S1Flag)*Scl) ,S);
%Drawing the opposite dir (making it full), these point in two directions
a2 = quiver3(Xcv(S1Flag) ,Ycv(S1Flag) ,Zcv(S1Flag) ,(-S1x(S1Flag))*Scl, (-S1y(S1Flag))*Scl, (-S1z(S1Flag)*Scl) ,S);
elseif sum(~S1Flag) == numel(S1Flag)
%Now drawing the negative part of S1 (compressional)
aneg= quiver3(Xcv(~S1Flag),Ycv(~S1Flag),Zcv(~S1Flag),(S1x(~S1Flag))*Scl,(S1y(~S1Flag))*Scl,(S1z(~S1Flag))*Scl,S);
%Drawing the opposite dir (making it full), these point in two directions
aneg2= quiver3(Xcv(~S1Flag),Ycv(~S1Flag),Zcv(~S1Flag),(-S1x(~S1Flag))*Scl,(-S1y(~S1Flag))*Scl,(-S1z(~S1Flag))*Scl,S);
else %do both
%Drawing the positive part of S1 (extensional)
a = quiver3(Xcv(S1Flag) ,Ycv(S1Flag) ,Zcv(S1Flag) ,(S1x(S1Flag))*Scl, (S1y(S1Flag))*Scl, (S1z(S1Flag)*Scl) ,S);
%Drawing the opposite dir (making it full), these point in two directions
a2 = quiver3(Xcv(S1Flag) ,Ycv(S1Flag) ,Zcv(S1Flag) ,(-S1x(S1Flag))*Scl, (-S1y(S1Flag))*Scl, (-S1z(S1Flag)*Scl) ,S);
%Now drawing the negative part of S1 (compressional)
aneg= quiver3(Xcv(~S1Flag),Ycv(~S1Flag),Zcv(~S1Flag),(S1x(~S1Flag))*Scl,(S1y(~S1Flag))*Scl,(S1z(~S1Flag))*Scl,S);
%Drawing the opposite dir (making it full), these point in two directions
aneg2= quiver3(Xcv(~S1Flag),Ycv(~S1Flag),Zcv(~S1Flag),(-S1x(~S1Flag))*Scl,(-S1y(~S1Flag))*Scl,(-S1z(~S1Flag))*Scl,S);
end

%Doing the same as above for S1 but for S2 and S3
if sum(S2Flag) == numel(S2Flag)
b = quiver3(Xcv(S2Flag) ,Ycv(S2Flag) ,Zcv(S2Flag) ,(S2x(S2Flag))*Scl, (S2y(S2Flag))*Scl, (S2z(S2Flag))*Scl ,S);
b2 = quiver3(Xcv(S2Flag) ,Ycv(S2Flag) ,Zcv(S2Flag) ,(-S2x(S2Flag))*Scl, (-S2y(S2Flag))*Scl, (-S2z(S2Flag))*Scl ,S);
elseif sum(~S2Flag) == numel(S2Flag)
bneg= quiver3(Xcv(~S2Flag),Ycv(~S2Flag),Zcv(~S2Flag),(S2x(~S2Flag))*Scl,(S2y(~S2Flag))*Scl,(S2z(~S2Flag))*Scl,S);
bneg2= quiver3(Xcv(~S2Flag),Ycv(~S2Flag),Zcv(~S2Flag),(-S2x(~S2Flag))*Scl,(-S2y(~S2Flag))*Scl,(-S2z(~S2Flag))*Scl,S);
else %do both
b = quiver3(Xcv(S2Flag) ,Ycv(S2Flag) ,Zcv(S2Flag) ,(S2x(S2Flag))*Scl, (S2y(S2Flag))*Scl, (S2z(S2Flag))*Scl ,S);
b2 = quiver3(Xcv(S2Flag) ,Ycv(S2Flag) ,Zcv(S2Flag) ,(-S2x(S2Flag))*Scl, (-S2y(S2Flag))*Scl, (-S2z(S2Flag))*Scl ,S);
bneg= quiver3(Xcv(~S2Flag),Ycv(~S2Flag),Zcv(~S2Flag),(S2x(~S2Flag))*Scl,(S2y(~S2Flag))*Scl,(S2z(~S2Flag))*Scl,S);
bneg2= quiver3(Xcv(~S2Flag),Ycv(~S2Flag),Zcv(~S2Flag),(-S2x(~S2Flag))*Scl,(-S2y(~S2Flag))*Scl,(-S2z(~S2Flag))*Scl,S);
end

if sum(S3Flag) == numel(S3Flag)
c = quiver3(Xcv(S3Flag) ,Ycv(S3Flag) ,Zcv(S3Flag) ,(S3x(S3Flag))*Scl, (S3y(S3Flag))*Scl, (S3z(S3Flag))*Scl ,S);
c2 = quiver3(Xcv(S3Flag) ,Ycv(S3Flag) ,Zcv(S3Flag) ,(-S3x(S3Flag))*Scl, (-S3y(S3Flag))*Scl, (-S3z(S3Flag))*Scl ,S);
elseif sum(~S3Flag) == numel(S3Flag)
cneg= quiver3(Xcv(~S3Flag),Ycv(~S3Flag),Zcv(~S3Flag),S3x(~S3Flag)*Scl,S3y(~S3Flag)*Scl,S3z(~S3Flag)*Scl,S);
cneg2= quiver3(Xcv(~S3Flag),Ycv(~S3Flag),Zcv(~S3Flag),-S3x(~S3Flag)*Scl,-S3y(~S3Flag)*Scl,-S3z(~S3Flag)*Scl,S);
else %do both
c = quiver3(Xcv(S3Flag) ,Ycv(S3Flag) ,Zcv(S3Flag) ,(S3x(S3Flag))*Scl, (S3y(S3Flag))*Scl, (S3z(S3Flag))*Scl ,S);
c2 = quiver3(Xcv(S3Flag) ,Ycv(S3Flag) ,Zcv(S3Flag) ,(-S3x(S3Flag))*Scl, (-S3y(S3Flag))*Scl, (-S3z(S3Flag))*Scl ,S);
cneg= quiver3(Xcv(~S3Flag),Ycv(~S3Flag),Zcv(~S3Flag),S3x(~S3Flag)*Scl,S3y(~S3Flag)*Scl,S3z(~S3Flag)*Scl,S);
cneg2= quiver3(Xcv(~S3Flag),Ycv(~S3Flag),Zcv(~S3Flag),-S3x(~S3Flag)*Scl,-S3y(~S3Flag)*Scl,-S3z(~S3Flag)*Scl,S);
end

title({'\fontsize{14}Principal directions','\fontsize{8}S1,Ext red,Comp yellow, S2,Ext green,Comp olivegreen, S3,Ext blue,Comp grey, scaled by magnitude'})
axis equal ;xlabel('x'); ylabel('y');zlabel('z');
%Setting styles for each set of lines

[ C1 ] = RGB2Colour(255,255,0);
[ C2 ] = RGB2Colour(192,255,62);
[ C3 ] = RGB2Colour(148,148,148);


if exist('a','var');set (a, 'color', 'red','showarrowhead', 'off','linewidth',Lw);end
if exist('a2','var');set (a2, 'color', 'red','showarrowhead', 'off','linewidth',Lw);end
if exist('aneg','var');set (aneg, 'color', C1,'showarrowhead', 'off','linewidth',Lw);end %yellow 
if exist('aneg2','var');set (aneg2, 'color', C1,'showarrowhead', 'off','linewidth',Lw);end%yellow
if exist('b','var');set (b, 'color', 'green','showarrowhead', 'off','linewidth',Lw);end %yellow 
if exist('b2','var');set (b2, 'color', 'green','showarrowhead', 'off','linewidth',Lw);end%yellow
if exist('bneg','var');set (bneg, 'color', C2,'showarrowhead', 'off','linewidth',Lw);end %olive green 
if exist('bneg2','var');set (bneg2, 'color', C2,'showarrowhead', 'off','linewidth',Lw);end%olive green
if exist('c','var');set (c, 'color', 'blue','showarrowhead', 'off','linewidth',Lw);end %yellow 
if exist('c2','var');set (c2, 'color', 'blue','showarrowhead', 'off','linewidth',Lw);end%yellow
if exist('cneg','var');set (cneg, 'color', C3,'showarrowhead', 'off','linewidth',Lw);end %olive green 
if exist('cneg2','var');set (cneg2, 'color', C3,'showarrowhead', 'off','linewidth',Lw);end%olive green
%set(c,'Visible','off');set(cneg,'Visible','off');set(c2,'Visible','off');set(cneg2,'Visible','off') %To turn of certain things do this
% ax.projection = 'perspective';

if nargin>4
hold on
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',(.2),'FaceColor', [0.5 0 0.9 ]);
hold off
end

end

