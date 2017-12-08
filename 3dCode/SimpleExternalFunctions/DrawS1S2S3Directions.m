function  DrawS1S2S3Directions(StressTensor,X,Y,Z,varargin)
% DrawS1S2S3Directions: Draws the pricipal stress directions as vectors 
%                   coloured dependant on magnitude sign and direction:
%                   S1,Ext red,     Comp yellow, 
%                   S2,Ext green,   Comp olivegreen,
%                   S3,Ext blue,    Comp grey.
%
%                   This also draws a triangulated surface if you supply
%                   this.
%               
% usage #1:
% DrawS1S2S3Directions(StressTensor,X,Y,Z,'Scale',5 )
%
% usage #2: with fracture surface
% DrawS1S2S3Directions(StressTensor,X,Y,Z,'Scale',5,'Triangles',Triangles,'Points',Points )
%
% Arguments: (input)
% StressTensor      - From the stress calculation 6*n vector, 
%                    [Sxx,Syy,Szz,Sxy,Sxz,Syz]
%
%
% XYZ               - The XYZ locations of the points where the tensors are.
%
% Additional arguments: (input)
%
% Scale,LineWidth   - Properties changing the length and the width of the
%                    resultant lines. Default to '1' and '2' if not set.
%                    MATLAB default width is 0.5.
%
% Points            - Columns 2 3 and 4 are the XYZ locations of one the
%                    corner points of a triangle. Column 1 is the index.
%                    Not needed unless you want to draw a surface too. 
%
% Triangles         -  Triangles is a list where each row contains 3 index
%                    locations in "Points" which contains the XYZ location
%                    of each corner of the triangle.
%                    Not needed unless you want to draw a surface too. 
%
% Arguments: (output)
% N/A - This just draws a figure. 
%
% Example usage:
%
% %Random Points on a sphere. 
% Theta = 2*pi*rand(1,1e3);
% Phi = asin(-1+2*rand(1,1e3));
% [X,Y,Z] = sph2cart(Theta',Phi',1);
% Sxx=ones(size(X)).*Theta';
% Syy=ones(size(X)).*-Theta';
% Szz=ones(size(X)).*Phi';
% Sxy=Syy;
% Sxz=Syy;
% Syz=Syy;
% StressTensor=[Sxx,Syy,Szz,Sxy,Sxz,Syz];
% DrawS1S2S3Directions(StressTensor,X,Y,Z,'Scale',0.01)
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Checking if additional arguments have been input into the function:
[ varargin,Scale ]     = AdditionalArgsInVaragin( 'Scale',    varargin,1 );
[ varargin,LineWidth ] = AdditionalArgsInVaragin( 'LineWidth',varargin,2 );
[ varargin,Points ]    = AdditionalArgsInVaragin( 'Points',   varargin,[] );
[ ~,Triangles ]        = AdditionalArgsInVaragin( 'Triangles',varargin,[] );


%Loading external Octave packages, Normc etc are in here
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
    pkg load general
    pkg load miscellaneous
    elseif isOctave==0
    %do nothing
end

%Creating colours
[ C1 ] = RGB2Colour(255,255,0);
[ C2 ] = RGB2Colour(192,255,62);
[ C3 ] = RGB2Colour(148,148,148);

%Forcing XYZ to column vecs 
Xcv=X(:); %col vecs
Ycv=Y(:); %col vecs
Zcv=Z(:); %col vecs

%Extracting strain & stress
[Sxx,Syy,Szz,Sxy,Sxz,Syz ] = ExtractCols( StressTensor );

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
[S1x,S1y,S1z] = ExtractCols( S1dirt );
[S2x,S2y,S2z] = ExtractCols( S2dirt );
[S3x,S3y,S3z] = ExtractCols( S3dirt );
%Makes sure there is no scaling going on between holds, have to use manual scaling.
%MATLAB usually autoscales in a quiver otherwise
S=0;

hold on

[Pdir,Pdirneg,Ndir,Ndirneg]=QuiverLines(Xcv,Ycv,Zcv,S1x,S1y,S1z,S1Flag,Scale,S);

if exist('Pdir','var');		set(Pdir,'color','red','showarrowhead','off','linewidth',LineWidth);	end		%red
if exist('Pdirneg','var');	set(Pdirneg,'color','red','showarrowhead','off','linewidth',LineWidth);	end		%red
if exist('Ndir','var');		set(Ndir,'color',C1,'showarrowhead','off','linewidth',LineWidth);		end		%yellow
if exist('Ndirneg','var');	set(Ndirneg,'color',C1,'showarrowhead','off','linewidth',LineWidth);	end		%yellow

[Pdir,Pdirneg,Ndir,Ndirneg]=QuiverLines(Xcv,Ycv,Zcv,S2x,S2y,S2z,S2Flag,Scale,S);

if exist('Pdir','var');		set(Pdir,'color','green','showarrowhead','off','linewidth',LineWidth);	end		%green
if exist('Pdirneg','var');   set(Pdirneg,'color','green','showarrowhead','off','linewidth',LineWidth);end	%green
if exist('Ndir','var');		set(Ndir,'color',C2,'showarrowhead','off','linewidth',LineWidth);		end		%olivegreen
if exist('Ndirneg','var');	set(Ndirneg,'color',C2,'showarrowhead','off','linewidth',LineWidth);	end		%olivegreen

[Pdir,Pdirneg,Ndir,Ndirneg]=QuiverLines(Xcv,Ycv,Zcv,S3x,S3y,S3z,S3Flag,Scale,S);

if exist('Pdir','var');		set(Pdir,'color','blue','showarrowhead','off','linewidth',LineWidth);	end		%blue
if exist('Pdirneg','var');	set(Pdirneg,'color','blue','showarrowhead','off','linewidth',LineWidth);end		%blue
if exist('Ndir','var');		set(Ndir,'color',C3,'showarrowhead','off','linewidth',LineWidth);		end		%olivegreen
if exist('Ndirneg','var');	set(Ndirneg,'color',C3,'showarrowhead','off','linewidth',LineWidth);	end		%olivegreen

%set(c,'Visible','off');set(cneg,'Visible','off');set(c2,'Visible','off');set(cneg2,'Visible','off') %To turn of certain things do this
% ax.projection = 'perspective';

title({'\fontsize{14}Principal directions','\fontsize{8}S1,Ext red,Comp yellow, S2,Ext green,Comp olivegreen, S3,Ext blue,Comp grey, scaled by magnitude'})
axis equal ;xlabel('x'); ylabel('y');zlabel('z');
WhiteFigure;
%Setting styles for each set of lines

chk=exist('Triangles','var');
if chk==1 
    hold on
    trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',(.2),'FaceColor', [0.5 0 0.9 ]);
    hold off
end

function [Pdir,Pdirneg,Ndir,Ndirneg]=QuiverLines(X,Y,Z,Snx,Sny,Snz,Flag,Scl,S)
    
    if sum(Flag) == numel(Flag)

        %Drawing the positive part of S1 (extensional)
        Pdir = quiver3(X(Flag) ,Y(Flag) ,Z(Flag) ,(Snx(Flag))*Scl, (Sny(Flag))*Scl, (Snz(Flag)*Scl) ,S);
        %Drawing the opposite dir (making it full), these point in two directions
        Pdirneg = quiver3(X(Flag) ,Y(Flag) ,Z(Flag) ,(-Snx(Flag))*Scl, (-Sny(Flag))*Scl, (-Snz(Flag)*Scl) ,S);
        Ndir=[];Ndirneg=[];

    elseif sum(~Flag) == numel(Flag)

        %Now drawing the negative part of S1 (compressional)
        Ndir= quiver3(X(~Flag),Y(~Flag),Z(~Flag),(Snx(~Flag))*Scl,(Sny(~Flag))*Scl,(Snz(~Flag))*Scl,S);
        %Drawing the opposite dir (making it full), these point in two directions
        Ndirneg= quiver3(X(~Flag),Y(~Flag),Z(~Flag),(-Snx(~Flag))*Scl,(-Sny(~Flag))*Scl,(-Snz(~Flag))*Scl,S);
        Pdir=[];Pdirneg=[];

    else %do both

        %Drawing the positive part of S1 (extensional)
        Pdir = quiver3(X(Flag) ,Y(Flag) ,Z(Flag) ,(Snx(Flag))*Scl, (Sny(Flag))*Scl, (Snz(Flag)*Scl) ,S);
        %Drawing the opposite dir (making it full), these point in two directions
        Pdirneg = quiver3(X(Flag) ,Y(Flag) ,Z(Flag) ,(-Snx(Flag))*Scl, (-Sny(Flag))*Scl, (-Snz(Flag)*Scl) ,S);
        %Now drawing the negative part of S1 (compressional)
        Ndir= quiver3(X(~Flag),Y(~Flag),Z(~Flag),(Snx(~Flag))*Scl,(Sny(~Flag))*Scl,(Snz(~Flag))*Scl,S);
        %Drawing the opposite dir (making it full), these point in two directions
        Ndirneg= quiver3(X(~Flag),Y(~Flag),Z(~Flag),(-Snx(~Flag))*Scl,(-Sny(~Flag))*Scl,(-Snz(~Flag))*Scl,S);

    end
    
end  %end internal func


end

