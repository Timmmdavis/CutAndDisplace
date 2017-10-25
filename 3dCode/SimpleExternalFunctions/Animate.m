%function  Animate
%Draws the displacement of the triangles
%Drawing animation of slip and exporting a film
%This was used to create the animation shown in
%"BEM_DDM_MATLAB\SharedFunctions\NashPointArray.gif"

%   Copyright 2017, Tim Davis, The University of Aberdeen

%Example usage: run slip calc to find your faults displacement then call 
%Animate
%in your cmd window,
%If you have the save options on it saves a MP4 in your current dir
%Variables defined below are important:


%Setting up some variables that you can use to change the animation
Save=0; %logical flag, 1 for saving the film and 0 for not. Saves in the current directory

%Movie properties if you save this
fps=15; %frames per second (24 default)
CreateVid=1; %1 for video, 0 for images (can be used to create a gif and works in octave). 

%Scaling displacement for animation
Scale=500; %Scale the displacement so its bigger or smaller. 1 is the same
mo=0.1; %the original seperation between the front and back faces, not
%needed if there is tensile displacement but for a shear fault it makes 
%a better animation. Put to 0 if not wanted
pad=0.1;%Adding some paddding to your axis limits if needed
Axisoff=1; %turning of the axes

%%%%
%Starting func
%%%%

%Creating fig
figure;

[ StrikeSlipCosine,DipSlipCosine ] = CalculateDSandSSDirs( FaceNormalVector );

%Calculate magnitude of displacement vectors (Cart coords)
SS_XYZ=(bsxfun(@times,StrikeSlipDisp,StrikeSlipCosine));
DS_XYZ=(bsxfun(@times,DipSlipDisp,DipSlipCosine));
TS_XYZ=(bsxfun(@times,TensileSlipDisp,FaceNormalVector));

TotalMovementXYZ=SS_XYZ+DS_XYZ+TS_XYZ;
TotalMovementXYZ=TotalMovementXYZ*Scale;

%Calculating new points and tris, less compact form for drawing
%Triangles will go 123,456,567 and point rows will correspond.
PointsNew=[];
for i=1:numel(Triangles(:,1))
    %very lazy way of doing this, seems fast enough..
    PointsNew=[PointsNew,Points(Triangles(i,1),2:4),Points(Triangles(i,2),2:4),Points(Triangles(i,3),2:4)];
end
sz=numel(PointsNew)/3;
PointsNew=reshape(PointsNew.',[],sz)';
TrianglesNew=1:numel(PointsNew(:,1));
sz=numel(TrianglesNew)/3;
TrianglesNew=reshape(TrianglesNew.',[],sz)';
PointsNew2=reshape(PointsNew.',[],sz)'; %each row now contains t1XYZ t2XYZ t2XYZ etc

%Calculating the points of the triangles as these move, this is along the
%normal direction
[PointsNew2P]=DisplaceTris(PointsNew2,[FaceNormalVector(:,1).*mo,FaceNormalVector(:,2).*mo,FaceNormalVector(:,3).*mo]);
PointsNew2P=reshape(PointsNew2P.',[],sz)'; %each row now contains t1XYZ t2XYZ t2XYZ etc
%opposite to the normal direction
[PointsNew2N]=DisplaceTris(PointsNew2,[-FaceNormalVector(:,1).*mo,-FaceNormalVector(:,2).*mo,-FaceNormalVector(:,3).*mo]);
PointsNew2N=reshape(PointsNew2N.',[],sz)'; %each row now contains t1XYZ t2XYZ t2XYZ etc

trisurf(TrianglesNew,PointsNew2(:,1),PointsNew2(:,2),PointsNew2(:,3),'FaceColor','m');
xlabel('x'); ylabel('y');axis('equal');
xlS=xlim;%final axis limits 
ylS=ylim;
zlS=zlim;

clf('reset')
%Getting final state, want to set the axes to this
[PointsNew3]=DisplaceTris(PointsNew2P,TotalMovementXYZ);
[PointsNew3N]=DisplaceTris(PointsNew2P,-TotalMovementXYZ);

%Draw this
trisurf(TrianglesNew,PointsNew3(:,1),PointsNew3(:,2),PointsNew3(:,3));
hold on
trisurf(TrianglesNew,PointsNew3N(:,1),PointsNew3N(:,2),PointsNew3N(:,3),'FaceColor','m');
xlabel('x'); ylabel('y');axis('equal');
%final axis limits 
xlF=xlim;		ylF=ylim;		zlF=zlim;
xlB=[xlS;xlF];	ylB=[ylS;ylF];	zlB=[zlS;zlF];
xl=[min(xlB(:,1))-pad,max(xlB(:,2))+pad];
yl=[min(ylB(:,1))-pad,max(ylB(:,2))+pad];
zl=[min(zlB(:,1))-pad,max(zlB(:,2))+pad];
xlim(xl);		ylim(yl);		zlim(zl)
if Axisoff==1
axis off; %Turn of the axes
set(gcf,'NumberTitle','off')%hide title
end

%PAUSE AND CHANGE TO DESIRED CAMERA ANGLE HERE:
test=1;

%# save all Camera* properties
g.CameraPosition     = get(gca, 'CameraPosition');  
g.CameraTarget       = get(gca, 'CameraTarget');    
g.CameraUpVector     = get(gca, 'CameraUpVector');
g.CameraUpVectorMode = get(gca, 'CameraUpVectorMode');
g.CameraViewAngle    = get(gca, 'CameraViewAngle');


%Setting up some properties
N=3;        %for logical indexing the XYZ
sz=50;      %the amount of figures that will be drawn. This times the pause time will give the run time of the figure
inv = linspace(0,1,sz); %linear vector that slowly forces the triangles to move to the max disp



if Save==1
if  CreateVid==0
f = warndlg('this func is about to spit loads of images in your current directory, do you really want to continue?, ctrl+C in cmd window to quit, if you do continue make sure your in a directory that you do not mind filling with images');
drawnow     % Necessary to print the message
%pause
waitfor(f);
elseif CreateVid==1
%do nothing
end
end

%Running animation
for i=1:sz

%Calculating the points of the triangles as these move, this is along the
%normal direction
[PointsNew3]=DisplaceTris(PointsNew2P,[TotalMovementXYZ(:,1).*inv(i),TotalMovementXYZ(:,2).*inv(i),TotalMovementXYZ(:,3).*inv(i)]);
%Neg Dir
[PointsNew3N]=DisplaceTris(PointsNew2N,[TotalMovementXYZ(:,1).*-inv(i),TotalMovementXYZ(:,2).*-inv(i),TotalMovementXYZ(:,3).*-inv(i)]);


%Restting the figure on each loop
clf('reset')
%Drawing the movement in the first direction, default colourmap
colormap default
trisurf(TrianglesNew,PointsNew3(:,1),PointsNew3(:,2),PointsNew3(:,3));
hold on
%Drawing the movement in the other direction, setting to pink so we know
%these are the other faces
trisurf(TrianglesNew,PointsNew3N(:,1),PointsNew3N(:,2),PointsNew3N(:,3),'FaceColor','m');
%Setting the axis limits and creating the title
xlabel('x'); ylabel('y');axis('equal');  
xlim(xl)
ylim(yl)
zlim(zl)
title('Animation of triangle displacement'),

titlesz=25;fntsz=21;
ChangeFontSizes(fntsz,titlesz);

%# set previous Camera* properties
if ~isempty(g)
    set(gca, g); 
end
if Axisoff==1
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

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  CreateVid==0
print ('Image','-djpeg'); %creating image
oldname='Image.jpg'; %its current name 
name=(['Image' num2str(i)]); %adding number
movefile(oldname, [name '.jpg']); %renaming

elseif CreateVid==1
% Store Each frame (Turn on when creating Movie)
M(i)=getframe(gcf); % leaving gcf out crops the frame in the movie.
end
  
else
%no saving
end    

end

%Images/frames are now created


if Save==1

if  CreateVid==0
% % % disp('Octave has issues Animating, use the created images and create a gif/movie online')
% % %if octave imwrite didnt crash I would use this script:
% ImageType = '*.jpg'; %we are looking for pngs
% Pngs=dir(ImageType); %searching current dir for these		
% L = {Pngs.name};	 %grabbing these as a string
% L=L';
% % %L needs sorting in numerical order
% Lnum = regexprep(L,'Image','');			%removing everything but numbers
% Lnum = regexprep(Lnum,'.jpg','');
% Lnum=str2double(Lnum); %create a double
% LL = [L, num2cell(Lnum)]; %append this to the cell (need to convert num to double or wont let you sort)
% [trash idx] = sort([LL{:,2}], 'ascend'); %create a sort index
% LL=LL(idx,:); %sort on index
% file_name = LL(:,1);	 %grabbing the sorted names
% file_name_NoExt = regexprep(file_name,'.jpg',''); %
% filename = 'AnimatedGif.gif';
% for i=1:length(file_name)
% clear RGB map1 map A 
% [RGB,map1]  = imread(file_name{i});
% % RGB = ind2rgb(RGB,map1);
% [A,map] = rgb2ind(RGB); %[A,map] = rgb2ind(RGB,256); in MATLAB we define quality
% 	if i == 1
% 	imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',fps); %crashes on my build of octave
%     else
%     imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',fps);
%     end
% end    

elseif CreateVid==1
% Exporting video if wanted
%movie(M,loops,fps); %movie(M,n(loops),fps(default12))
v=VideoWriter('Animated Displacement','MPEG-4');
v.FrameRate=fps;
open(v)
writeVideo(v,M)
disp('duration will be this many seconds, change framerate to increase')
disp(v.Duration)
close(v)
end
else
end    

%Clearing all vars so there is no clutter in the workspace
clear Save fps Scale mo StrikeSlipCosine DipSlipCosine SS_XYZ DS_XYZ TS_XYZ TotalMovementXYZ TotalMovementXYZ
clear PointsNew sz TrianglesNew PointsNew2 N PointsNewMvX PointsNewMvY PointsNewMvZ PointsNewP PointsNew2P PointsNewN PointsNew2N
clear test g v M hAxes inv PointsNew3 PointsNew3N

function [PointsNewLoc]=DisplaceTris(Points,Movement)
%Takes inputs and moves each triangles points by the amount specified in the 3*3 vector 'Movement'
N=3;
PointsNewMvX=bsxfun(@plus,Points(:,1:N:end),(Movement(:,1))); 
PointsNewMvY=bsxfun(@plus,Points(:,2:N:end),(Movement(:,2))); 
PointsNewMvZ=bsxfun(@plus,Points(:,3:N:end),(Movement(:,3))); 
PointsNewLoc=[reshape(PointsNewMvX.',1,[])',reshape(PointsNewMvY.',1,[])',reshape(PointsNewMvZ.',1,[])'];
end