function [ hlink ] = IsoContoursPrincipalStressPercentiles( S1,S2,S3,X,Y,Z,Triangles,Points )
%IsoContoursPrincipalStressPercentiles Draw Isocontours of the 3 principal
%stresses. Big linked diagram
%Does this at percentiles of the stress. This can be changed if you want 
%to pick an actual value~line 19

%S1,S2,S3 = input stresses as vectors
%X,Y,Z = Gridded observation points
%Triangles,Points = the boundary surface (not needed)

%Colour descriptions are described in plot titles

%Example usage:
%[ hlink ] =IsoContoursPrincipalStressPercentiles( S1,S2,S3,X,Y,Z,Triangles,Points );

%   Copyright 2017, Tim Davis, The University of Aberdeen

High=20; %% 75 would be the 75th percentile
Mid=15;  %% 50 would be the 50th percentile
Low=10;  %% 25 would be the 50th percentile

%finding percentile values for 3d plotting
S1_High = prctile(S1(:),High);%75th percentile
S1_Mid = prctile(S1(:),Mid);%50th percentile
S1_Low = prctile(S1(:),Low);%25th percentile
S2_High = prctile(S2(:),High);%75th percentile
S2_Mid = prctile(S2(:),Mid);%50th percentile
S2_Low = prctile(S2(:),Low);%25th percentile
S3_High = prctile(S3(:),High);%75th percentile
S3_Mid = prctile(S3(:),Mid);%50th percentile
S3_Low = prctile(S3(:),Low);%25th percentile


%Drawing isosurfaces of the principal stresses. 
figure_handle = figure;subplot(1,3,1), pneg=patch(isosurface(X,Y,Z,S1,S1_Mid)) ;
        set(pneg,'facecolor',normc([255;48;48]),'EdgeColor','none'); %red
        hold on
        p0 = patch(isosurface(X,Y,Z,S1,S1_Low)) ; 
        set(p0,'FaceColor',normc([255;165;0]),'EdgeColor','none')    %orange 
        ppos = patch(isosurface(X,Y,Z,S1,S1_High)) ; 
        set(ppos,'FaceColor',normc([255;255;0]),'EdgeColor','none')  %yellow
grid on ; 
daspect([1 1 1]) ; 
view(3); 
axis equal ; 
camlight ; 
lighting gouraud ; 
title({'\fontsize{14}Sig1...MostExtensional','\fontsize{8}Percentiles,Red=50th Orange=25th Yellow=75th'})
if nargin>6
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',(.2),'FaceColor', [0.5 0 0.9 ]);
end
hold off;

subplot(1,3,2),pneg=patch(isosurface(X,Y,Z,S2,S2_Mid)) ; %subplot(2,3,2),
        set(pneg,'FaceColor',normc([192;255;62]),'EdgeColor','none');  %olive green
        hold on
        p0 = patch(isosurface(X,Y,Z,S2,S2_Low)) ; 
        set(p0,'FaceColor',normc([32;178;170]),'EdgeColor','none');    %turquoise
        ppos = patch(isosurface(X,Y,Z,S2,S2_High)) ; 
        set(ppos,'FaceColor',normc([224;102;255]),'EdgeColor','none'); %purple
grid on ; 
daspect([1 1 1]) ; 
view(3); 
axis equal ; 
camlight ; 
lighting gouraud ; 
title({'\fontsize{14}Sig2...2ndPrincipal','\fontsize{8}Percentiles,Olive green=50th Turquoise=25th Purple=75th'})
if nargin>6
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',(.2),'FaceColor', [0.5 0 0.9 ]);
end
hold off;

subplot(1,3,3),pneg=patch(isosurface(X,Y,Z,S3,S3_Mid)) ;
        set(pneg,'FaceColor',normc([148;148;148]),'EdgeColor','none'); %grey
        hold on
        p0 = patch(isosurface(X,Y,Z,S3,S3_Low)) ; 
        set(p0,'FaceColor',normc([152;245;255]),'EdgeColor','none');   %blue    
        ppos = patch(isosurface(X,Y,Z,S3,S3_High)) ; 
        set(ppos,'FaceColor',normc([148;0;211]),'EdgeColor','none'); %dark purp  
grid on ; 
daspect([1 1 1]) ; 
view(3); 
axis equal ; 
camlight ; 
lighting gouraud ; 
title('Sig3TDEScript')
title({'\fontsize{14}Sig3...MostCompressional','\fontsize{8}Percentiles,Grey=50th Blue=25th Dark purple=75th'})
if nargin>6
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',(.2),'FaceColor', [0.5 0 0.9 ]);
end
hold off;

%Linking all axes together to spin at the same time
all_ha = findobj( figure_handle, 'type', 'axes', 'tag', '' );
hlink = linkprop(all_ha,{'CameraPosition','CameraUpVector'});
rotate3d on


end

