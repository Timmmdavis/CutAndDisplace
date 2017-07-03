%function  Animate
%Draws the displacement of the triangles
%Drawing animation of slip

%   Copyright 2017, Tim Davis, The University of Aberdeen

%Example usage: run slip calc to find your faults displacement then call 
%Animate2D
%in your cmd window,
%If you have the save options on it saves a MP4 in your current dir
%Variables defined below are important:

%Setting up some variables that you can use to change the animation

%Movie properties if you save this
Save=0; %logical flag, 1 for saving the film and 0 for not. Saves in the current directory
%Movie properties if you save this
fps=15; %frames per second (24 default)

%Scaling displacement for animation
Scale=0.1; %Scale the displacement so its bigger or smaller. 1 is the same
mo=0.1; %the original seperation between the dislocation faces, not
%needed if there is tensile displacement but for a shear fault it makes 
%a better animation. Put to 0 if not wanted


%%%%
%Starting func
%%%%

%Creating fig
figure;

TensileCosine=[-cos(NormAng),cos(Beta)];
ShearCosine=[cos(Beta),-cos(NormAng)];

TS_XY=(bsxfun(@times,TensileDisp,TensileCosine));
SS_XY=(bsxfun(@times,ShearDisp,ShearCosine));


TotalMovementXY=TS_XY+SS_XY;
TotalMovementXY=TotalMovementXY*Scale;


%Quickly moving the inner and out triangles away a miniscule amount for
%drawing
%Calculating the points of the triangles as these move, this is along the
%normal direction
PointsNewMvXP=(bsxfun(@plus,Points(:,1:2),-cos(NormAng).*-mo)); 
PointsNewMvYP=(bsxfun(@plus,Points(:,3:4),cos(Beta).*-mo)); 
PointsNewP=[PointsNewMvXP,PointsNewMvYP];

%Calculating the points of the triangles as these move, this is opposite to
%the normal direction (and the DS and SS dir cosines), if needed you can draw these
%in the func that creates them.
PointsNewMvXN=(bsxfun(@minus,Points(:,1:2),-cos(NormAng).*-mo)); 
PointsNewMvYN=(bsxfun(@minus,Points(:,3:4),cos(Beta).*-mo)); 
PointsNewN=[PointsNewMvXN,PointsNewMvYN];

%Getting limits if you need to fix these
% figure; trisurf(TrianglesNew,PointsNew(:,1),PointsNew(:,2),PointsNew(:,3));xlabel('x'); ylabel('y'); axis('equal');
% xl = xlim;yl = ylim;zl = zlim;
% axislimits =[xl,yl,zl];
figure; h=axes; hAxes.DrawMode = 'fast';hold on
for i=1:numel(PointsNewP(:,1));
plot(PointsNewP(i,1:2),PointsNewP(i,3:4),'r');
end
for i=1:numel(PointsNewN(:,1));
plot(PointsNewN(i,1:2),PointsNewN(i,3:4),'b');
end
title('Animation of displacement'), xlabel('x'), ylabel('y'),axis equal;


%PAUSE AND CHANGE zoom, size etc
test=1;

%Setting up some properties
N=3;        %for logical indexing the XYZ
sz=50;      %the amount of figures that will be drawn. This times the pause time will give the run time of the figure
inv = linspace(0,1,sz); %linear vector that slowly forces the triangles to move to the max disp



if Save==1
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1;
f = warndlg('Octave is about to spit loads of images in your current directory, do you really want to continue?, ctrl+C in cmd window to quit, if you do continue make sure your in a directory that you don’t mind filling with images');
drawnow     % Necessary to print the message
%waitfor(f);
elseif isOctave==0;
%do nothing
end
end

%Running animation
for i=1:sz;

%Calculating the points of the elements as these move, this is along the
%normal direction
PointsNewMvXP=(bsxfun(@plus,PointsNewP(:,1:2),TotalMovementXY(:,1).*inv(i))); 
PointsNewMvYP=(bsxfun(@plus,PointsNewP(:,3:4),TotalMovementXY(:,2).*inv(i))); 

PointsNewMvXN=(bsxfun(@minus,PointsNewN(:,1:2),TotalMovementXY(:,1).*inv(i))); 
PointsNewMvYN=(bsxfun(@minus,PointsNewN(:,3:4),TotalMovementXY(:,2).*inv(i))); 

%Restting the figure on each loop
clf('reset')
hold on
line([PointsNewMvXP(:,1)';PointsNewMvXP(:,2)'],[PointsNewMvYP(:,1)';PointsNewMvYP(:,2)'],'color','r')
line([PointsNewMvXN(:,1)';PointsNewMvXN(:,2)'],[PointsNewMvYN(:,1)';PointsNewMvYN(:,2)'],'color','b')
scatter([PointsNewMvXP(:,1);PointsNewMvXP(:,2)],[PointsNewMvYP(:,1);PointsNewMvYP(:,2)],5,'r')
scatter([PointsNewMvXN(:,1);PointsNewMvXN(:,2)],[PointsNewMvYN(:,1);PointsNewMvYN(:,2)],5,'b')
%Setting the axis limits and creating the title
xlabel('x'); ylabel('y');axis('equal');  
title('Animation of displacement'),

% % %Drawing the movement in the first direction, default colourmap
% % hold on
% % for ii=1:numel(PointsNewP(:,1));
% % plot(PointsNewP(ii,1:2),PointsNewP(ii,3:4),'r');
% % end
% % for ii=1:numel(PointsNewP(:,1));
% % plot(PointsNewP(ii,1:2),PointsNewP(ii,3:4),'o', 'markersize', 3,'markeredgecolor','r');
% % end
% % %Drawing the movement in the other direction, setting to pink so we know
% % %these are the other faces
% % for ii=1:numel(PointsNewN(:,1));
% % plot(PointsNewN(ii,1:2),PointsNewN(ii,3:4),'b');
% % end
% % for ii=1:numel(PointsNewN(:,1));
% % plot(PointsNewN(ii,1:2),PointsNewN(ii,3:4),'o', 'markersize', 3,'markeredgecolor','b');
% % end
% % %Setting the axis limits and creating the title
% % xlabel('x'); ylabel('y');axis('equal');  
% % title('Animation of displacement'),

%Forcing this to draw
drawnow
%Pausing so the timing is not dependant on the machine
pause(0.00005)



if Save==1

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1;
print 'Image'.jpg; %creating image
oldname='Image.jpg'; %its current name 
name=(['Image' num2str(i)]); %adding number
movefile(oldname, [name '.jpg']); %renaming

elseif isOctave==0;
% Store Each frame (Turn on when creating Movie)
M(i)=getframe(gcf); % leaving gcf out crops the frame in the movie.
end
  
else
%no saving
end    

end

%Images/frames are now created


if Save==1

if  isOctave==1;
disp('Octave has issues Animating, use the created images and create a gif/movie online')

%if octave imwrite didnt crash I would use this script:
% % ImageType = '*.jpg'; %we are looking for pngs
% % Pngs=dir(ImageType); %searching current dir for these		
% % L = {Pngs.name};	 %grabbing these as a string
% % L=L';
% % %L needs sorting in numerical order
% % Lnum = regexprep(L,'Image','');			%removing everything but numbers
% % Lnum = regexprep(Lnum,'.jpg','');
% % Lnum=str2double(Lnum); %create a double
% % LL = [L, num2cell(Lnum)]; %append this to the cell (need to convert num to double or wont let you sort)
% % [trash idx] = sort([LL{:,2}], 'ascend'); %create a sort index
% % LL=LL(idx,:); %sort on index
% % file_name = LL(:,1);	 %grabbing the sorted names
% % file_name_NoExt = regexprep(file_name,'.jpg',''); %
% % filename = 'AnimatedGif.gif';
% % for i=1:length(file_name)
% % clear RGB map1 map A 
% % [RGB,map1]  = imread(file_name{i});
% % % RGB = ind2rgb(RGB,map1);
% % [A,map] = rgb2ind(RGB); %[A,map] = rgb2ind(RGB,256); in MATLAB we define quality
% % 	if i == 1;
% % 	imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',fps); %crashes on my build of octave
% %     else
% %     imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',fps);
% %     end
% % end    



elseif isOctave==0;
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

