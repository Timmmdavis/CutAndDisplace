function [ hlink ] = IsoContoursPrincipalStressPercentiles( S1,S2,S3,X,Y,Z,varargin)
% IsoContoursPrincipalStressPercentiles:  Draws Isocontours of the 3 principal
%                   stresses for gridded data in XYZ. The resultant diagram
%                   is 3 linked subplots for S1, S2 and S3. The isocontours
%                   are drawn at the percentiles of the input data chosen
%                   by the user. 
%
%                   Mesh surfaces can also be input (last arguments) to
%                   draw these in the diagrams too. 
%               
% usage #1: simple, drawing the 75 50 and 25th percentile
% IsoContoursPrincipalStressPercentiles( S1,S2,S3,X,Y,Z )
%
% usage #2: drawing the lowest percentile levels at the 15th percentile
% IsoContoursPrincipalStressPercentiles( S1,S2,S3,X,Y,Z,'P_Low',15 )
%
% usage #3: As usage 1 but also drawing a triangulated mesh. 
% IsoContoursPrincipalStressPercentiles( S1,S2,S3,X,Y,Z,'Triangles',Triangles,'Points',Points)
%
% usage #4: As usage 1 but making the contours slightly transparent. 
% IsoContoursPrincipalStressPercentiles( S1,S2,S3,X,Y,Z,'Alpha',0.2)
%
%
% Arguments: (input)
% S1,S2 & S3        - Gridded principal stress magnitudes at grid points
%                    that correpond to those in XYZ, Assuming with the
%                    labelling that S1 is the most extensional stress.
%
% XYZ               - The XYZ locations of the points where the principal
%                    stresses are
%
% Additional arguments: (input)
%
% P_Low,P_Mid,P_High- The percentiles of each principal stress that the 
%                    stresses are drawn at. 
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
% Alpha             -  The transparency of the isocontours. 0 being
%                     see-through and 1 opaque. 
%
% Arguments: (output)
% hlink             - Property linking of the figures, only needed if you are
%                    modifying this property later. 
%
% Example usage:
%
% %Just drawing principal stresses that scale with coordinate. 
% [X,Y,Z]=meshgrid(-1:0.1:1,-1:0.1:1,-1:0.1:1);
% S1=X;S2=Y;S3=Z;
% IsoContoursPrincipalStressPercentiles( S1,S2,S3,X,Y,Z )
% %Now changing S1 to show the difference. 
% S1=abs(S1.*Y.^2)./(Z);
% IsoContoursPrincipalStressPercentiles( S1,S2,S3,X,Y,Z )
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Finding if additional arguments have been input
[ varargin,Alpha ]  = AdditionalArgsInVaragin( 'Alpha', varargin,1 );
[ varargin,P_Low ]  = AdditionalArgsInVaragin( 'P_Low', varargin,25 );
[ varargin,P_Mid ]  = AdditionalArgsInVaragin( 'P_Mid', varargin,50 );
[ varargin,P_High ] = AdditionalArgsInVaragin( 'P_High',varargin,75 );
[ varargin,Points ] = AdditionalArgsInVaragin( 'Points',varargin,[] );
[ ~,Triangles ] = AdditionalArgsInVaragin( 'Triangles', varargin,[] );

%finding percentile values for 3d plotting
%                                   If defaults:
S1_High = Percentile(S1(:),P_High); %75th percentile
S1_Mid  = Percentile(S1(:),P_Mid);  %50th percentile
S1_Low  = Percentile(S1(:),P_Low);  %25th percentile
S2_High = Percentile(S2(:),P_High); %75th percentile
S2_Mid  = Percentile(S2(:),P_Mid);  %50th percentile
S2_Low  = Percentile(S2(:),P_Low);  %25th percentile
S3_High = Percentile(S3(:),P_High); %75th percentile
S3_Mid  = Percentile(S3(:),P_Mid);  %50th percentile
S3_Low  = Percentile(S3(:),P_Low);  %25th percentile


%Drawing isosurfaces of the principal stresses. 
figure_handle = figure;
%Getting figure and inflating this to full screen size
set(figure_handle, 'Position', get(0,'Screensize'));

%Drawing S1
subplot(1,3,1), 
hold on
C_Low= normr([255,165,0]);%orange 
C_Mid= normr([255,48,48]);%red
C_High=normr([255,255,0]);%yellow
DrawIsosurface(X,Y,Z,S1,S1_Low,S1_Mid,S1_High,C_Low,C_Mid,C_High,Alpha)
title({'\fontsize{14}Sig1...MostExtensional',...
    ['\fontsize{8}Percentiles,Red=',num2str(P_Mid),...
     'th Orange=',num2str(P_Low),...
     'th Yellow=',num2str(P_High),'th']})
if ~isempty(Triangles)
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',0.8,'FaceColor', [0.5 0 0.9 ]);
end
hold off;

%Drawing S2
subplot(1,3,2), 
hold on
C_Low= normr([32,178,170]);%turquoise
C_Mid= normr([192,255,62]);%olive green
C_High=normr([255,255,0]); %purple
DrawIsosurface(X,Y,Z,S2,S2_Low,S2_Mid,S2_High,C_Low,C_Mid,C_High,Alpha)
title({'\fontsize{14}Sig2...2ndPrincipal',...
    ['\fontsize{8}Percentiles,OliveGreen=',num2str(P_Mid),...
     'th Turquoise=',num2str(P_Low),...
     'th Purple=',num2str(P_High),'th']})

if ~isempty(Triangles)
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',0.8,'FaceColor', [0.5 0 0.9 ]);
end
hold off;

%Drawing S3
subplot(1,3,3), 
hold on
C_Low= normr([152,245,255]);%blue   
C_Mid= normr([148,148,148]);%grey
C_High=normr([148,0,211]);  %dark purple
DrawIsosurface(X,Y,Z,S3,S3_Low,S3_Mid,S3_High,C_Low,C_Mid,C_High,Alpha)
title({'\fontsize{14}Sig3...MostCompressional',...
    ['\fontsize{8}Percentiles,Grey=',num2str(P_Mid),...
     'th Blue=',num2str(P_Low),...
     'th DarkPurple=',num2str(P_High),'th']})
if ~isempty(Triangles)
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',0.8,'FaceColor', [0.5 0 0.9 ]);
end
hold off;

WhiteFigure;
%Linking all axes together to spin axes at the same time
all_ha = findobj( figure_handle, 'type', 'axes', 'tag', '' );
hlink = linkprop(all_ha,{'CameraPosition','CameraUpVector'});
rotate3d on

%Internal function to find the Percentiles. 
%https://de.mathworks.com/matlabcentral/answers/127944-how-to-pick-the-j-th-percentile-of-a-vector
function val = Percentile(arr, pct)
    len = length(arr);
    ind = floor(pct/100*len);
    newarr = sort(arr);
    val = newarr(ind);
end

%Internal function to draw the Isosurfaces at chosen locations. 
function DrawIsosurface(X,Y,Z,Sx,Sx_Low,Sx_Mid,Sx_High,C_Low,C_Mid,C_High,Alpha)
    
    P1 = patch(isosurface(X,Y,Z,Sx,Sx_Low)) ; 
    set(P1,'FaceColor',C_Low,'EdgeColor','none','FaceAlpha',Alpha)   

    P2=patch(isosurface(X,Y,Z,Sx,Sx_Mid)) ;
    set(P2,'FaceColor',C_Mid,'EdgeColor','none','FaceAlpha',Alpha);

    P3 = patch(isosurface(X,Y,Z,Sx,Sx_High)) ; 
    set(P3,'FaceColor',C_High,'EdgeColor','none','FaceAlpha',Alpha) 

    %Drawing things:
    grid on ; 
    daspect([1 1 1]) ; 
    view(3); 
    axis equal ; 
    camlight ; 
    lighting gouraud ;
    
end

end

