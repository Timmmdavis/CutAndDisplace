function Animate3d(Save,fps,Scale,Seperation,AxisOff,AxisLimPad,FaceNormalVector,Dn,Dds,Dss,P1,P2,P3)
% Animate2d: Drawing an animation of the displacement of elements and
%                    exporting as a MP4 (film). %This was used to create
%                    the animation shown in
%                    "Path\Images\NashPointArray.gif"
%               
% usage #1:
% Animate2d(Save,fps,Scale,Seperation,AxisOff,AxisLimPad,LineNormalVector,Dn,Ds,Points)
%
%
% Arguments: (input)
% Save              - Flag, when "1" a movie is saved to the current
%                     directory called 'Animated_Displacement.Mp4".
%                     Slightly different in Octave where a load of images
%                     are printed and you need to make a gif with these. 
%
%
% fps               - The frames per second of the MP4. (24 default)
%
% Scale             - Scales the displacement, if you are not seeing
%                     displaecment in the animation increase this. If you
%                     see too little decrease this value.
%
% Seperation        - The Start seperation of the elements. If 0 these will
%                     overlap and its hard to tell which side is the
%                     element normal. When a postive number the red lines
%                     are the positive side of the element.
%
% AxisOff           - Turns the axis ticks off and titles. 
%
% AxisLimPad        - Additional padding for the final axis limits of the
%                     animation. If 0 the elements will reach the very edge
%                     of the axis limits on the last frame
%
% FaceNormalVector  - The Direction cosines of each element. 
%                    [CosAx,CosAy,CosAz];
%
% Dn,Dds,Dss        - The displacement value at each element (normal,
%                     dipslip and strikeslip)
%
% P1,P2,P3          - n*3 Column vector where each 'P' represents the
%                    different corner points of one of the triangles (XYZ).
%                    Not as efficient in terms of storage but easier to
%                    understand.
%
%
% Arguments: (output)
% N/A
%
% Example usage:
%
% %Run "MainFrame.m" past "SlipCalculator2d.m" then call the function: 
% Save=1; 
% fps=15; 
% Scale=10; 
% Seperation=0.01; 
% AxisOff=0;
% AxisLimPad=0.1; %padding the axis lims
% Animate3d(Save,fps,Scale,Seperation,AxisOff,AxisLimPad,FaceNormalVector,Dn,Dds,Dss,P1,P2,P3)
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%%
%Calculating new points and tris, less compact form for drawing
%Triangles will go 123,456,567 and point rows will correspond.
PointsLst=[P1,P2,P3];
TrianglesLst=1:numel(FaceNormalVector);
TrianglesLst=reshape(TrianglesLst,3,[])';

%% Setting up the the movement vectors  

[ StrikeSlipCosine,DipSlipCosine ] = CalculateDSandSSDirs( FaceNormalVector );

%Calculate magnitude of displacement vectors (Cart coords)
SS_XYZ=(bsxfun(@times,Dss,StrikeSlipCosine));
DS_XYZ=(bsxfun(@times,Dds,DipSlipCosine));
TS_XYZ=(bsxfun(@times,Dn,FaceNormalVector));

%Final displacement
FinD_XYZ=SS_XYZ+DS_XYZ+TS_XYZ;
FinD_XYZ=FinD_XYZ*Scale;

%% First we want the axis limits of our data (set limits of animation to this)

%Creating fig
figure;

Seperation=[FaceNormalVector(:,1).*Seperation,FaceNormalVector(:,2).*Seperation,FaceNormalVector(:,3).*Seperation];

%Start axis lims
%The triangles after being seperated by chonsen amount
[PointsNew2P]=DisplaceTris(PointsLst,Seperation);
[PointsNew2N]=DisplaceTris(PointsLst,-Seperation);
trisurf(TrianglesLst,PointsNew2P(:,1),PointsNew2P(:,2),PointsNew2P(:,3));
hold on
trisurf(TrianglesLst,PointsNew2N(:,1),PointsNew2N(:,2),PointsNew2N(:,3),'FaceColor','m');
xlabel('x'); ylabel('y');axis('equal');
xlS=xlim;%Got start lims
ylS=ylim;
zlS=zlim;



%Final axis lims
%Getting final state
[PointsNew3]=DisplaceTris(PointsLst,FinD_XYZ+Seperation);
[PointsNew3N]=DisplaceTris(PointsLst,-FinD_XYZ-Seperation);
clf('reset')
%Draw this
trisurf(TrianglesLst,PointsNew3(:,1),PointsNew3(:,2),PointsNew3(:,3));
hold on
trisurf(TrianglesLst,PointsNew3N(:,1),PointsNew3N(:,2),PointsNew3N(:,3),'FaceColor','m');
xlabel('x'); ylabel('y');axis('equal');
xlF=xlim;%final axis limits 
ylF=ylim;
zlF=zlim;

%Getting the limits that will be used
xlB=[xlS;xlF];	ylB=[ylS;ylF];	zlB=[zlS;zlF];
xl=[min(xlB(:,1))-AxisLimPad,max(xlB(:,2))+AxisLimPad];
yl=[min(ylB(:,1))-AxisLimPad,max(ylB(:,2))+AxisLimPad];
zl=[min(zlB(:,1))-AxisLimPad,max(zlB(:,2))+AxisLimPad];
xlim(xl);		ylim(yl);		zlim(zl)
%Turning of the axes if desired
if AxisOff==1
    axis off; %Turn of the axes
    set(gcf,'NumberTitle','off')%hide title
end

%If you want to change the camera angle do so here: (add breakpoint, change
%angle, continue run)
test=1;

%# save all Camera* properties
g.CameraPosition     = get(gca, 'CameraPosition');  
g.CameraTarget       = get(gca, 'CameraTarget');    
g.CameraUpVector     = get(gca, 'CameraUpVector');
g.CameraUpVectorMode = get(gca, 'CameraUpVectorMode');
g.CameraViewAngle    = get(gca, 'CameraViewAngle');

%% Now the loop starts where the animation takes place.

%Setting up some properties
N=3;        %for logical indexing the XYZ
sz=50;      %the amount of figures that will be drawn. This times the pause time will give the run time of the figure
inv = linspace(0,1,sz); %linear vector that slowly forces the triangles to move to the max disp
Frame=getframe(gcf); %Start array 'Frame'
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB

%Running animation
for i=1:sz

    %Calculating the points of the triangles as these move, this is along the
    %normal direction
    [PointsNew3]=DisplaceTris(PointsLst,[FinD_XYZ(:,1).*inv(i),FinD_XYZ(:,2).*inv(i),FinD_XYZ(:,3).*inv(i)]+Seperation);
    %Neg Dir
    [PointsNew3N]=DisplaceTris(PointsLst,[FinD_XYZ(:,1).*-inv(i),FinD_XYZ(:,2).*-inv(i),FinD_XYZ(:,3).*-inv(i)]-Seperation);
    
    %Restting the figure on each loop
    clf('reset')
    %Drawing the movement in the first direction, default colourmap
    colormap default
    trisurf(TrianglesLst,PointsNew3(:,1),PointsNew3(:,2),PointsNew3(:,3));
    hold on     %Drawing the movement in the other direction (Pink)
    trisurf(TrianglesLst,PointsNew3N(:,1),PointsNew3N(:,2),PointsNew3N(:,3),'FaceColor','m');
    
    %Setting the axis limits and creating the title
    xlabel('x'); ylabel('y');axis('equal');  
    xlim(xl)
    ylim(yl)
    zlim(zl)
    title('Animation of triangle displacement, Cols=ElNormSde'),
    titlesz=14;
    fntsz=12;
    ChangeFontSizes(fntsz,titlesz);

    %# set previous Camera* properties
    if ~isempty(g)
        set(gca, g); 
    end
    if AxisOff==1
        axis off; %Turn of the axes
        set(gcf,'NumberTitle','off') %hide title
    end

    %# save all Camera* properties
    g.CameraPosition     = get(gca, 'CameraPosition');  
    g.CameraTarget       = get(gca, 'CameraTarget');    
    g.CameraUpVector     = get(gca, 'CameraUpVector');
    g.CameraUpVectorMode = get(gca, 'CameraUpVectorMode');
    g.CameraViewAngle    = get(gca, 'CameraViewAngle');
    
    %Forcing this to draw
    drawnow
    %Pausing so the timing is not dependant on the machine
    pause(0.005)

    if Save==1
        [ Frame ] = GrabFramesAndSaveFilm( i,Frame,isOctave );
    end  

end %end lp

%Images/frames are now created

%Images/frames are now created
if Save==1

    if  isOctave==1
        
        disp('Octave has issues Animating, use the created images and create a gif/movie online')

    elseif isOctave==0
        
        % Exporting video if wanted
        %movie(Frame,loops,fps); %movie(Frame,n(loops),fps(default12))
        v=VideoWriter('Animated Displacement','MPEG-4');
        v.FrameRate=fps;
        open(v)
        writeVideo(v,Frame)
        disp('duration will be this many seconds, change framerate to increase')
        disp(v.Duration)
        close(v)
        
    end
end

end %end func
 
%Internal func
function [PointsNewLoc]=DisplaceTris(Points,Movement)
    %Takes inputs and moves each triangles points by the amount specified in the 3*3 vector 'Movement'
    N=3;
    PointsNewMvX=bsxfun(@plus,Points(:,1:N:end),(Movement(:,1))); 
    PointsNewMvY=bsxfun(@plus,Points(:,2:N:end),(Movement(:,2))); 
    PointsNewMvZ=bsxfun(@plus,Points(:,3:N:end),(Movement(:,3))); 
    PointsNewLoc=[reshape(PointsNewMvX.',1,[])',reshape(PointsNewMvY.',1,[])',reshape(PointsNewMvZ.',1,[])'];
end