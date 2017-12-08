function Animate2d(Save,fps,Scale,Seperation,AxisOff,AxisLimPad,LineNormalVector,Dn,Ds,Points)
% Animate2d: Drawing an animation of the displacement of elements and
%                    exporting as a MP4 (film)
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
% LineNormalVector  - The Direction cosines of each element. [CosAx,CosAy];
%
% Dn,Ds             - The displacement value at each element (normal and
%                     shear)
%
% Points            - The Start X and Y and End X and Y end points of each
%                     element, n*4 vector
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
% Animate2d(Save,fps,Scale,Seperation,AxisOff,AxisLimPad,LineNormalVector,Dn,Ds,Points)
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%% Setting up the the movement vectors  

%Grab direction cosines of each element
CosAx=LineNormalVector(:,1);
CosAy=LineNormalVector(:,2);

%Compute the directions of the slip components.
DnCosine=LineNormalVector;
DsCosine=[-CosAy,CosAx];
%Cartesian directions and magnitudes of slip
Dn_xy=(bsxfun(@times,Dn,DnCosine));
Ds_xy=(bsxfun(@times,Ds,DsCosine));

%Final displacement
FinD_XY=Dn_xy+Ds_xy;
FinD_XY=FinD_XY*Scale;

%% First we want the axis limits of our data (set limits of animation to this)

%Creating fig
figure;
%Start axis lims
%Moving by the desired seperation
PointsNewMvXP=(bsxfun(@plus,Points(:,1:2),CosAx.*Seperation)); 
PointsNewMvYP=(bsxfun(@plus,Points(:,3:4),CosAy.*Seperation)); 
PointsNewP=[PointsNewMvXP,PointsNewMvYP];
%Neg side
PointsNewMvXN=(bsxfun(@minus,Points(:,1:2),CosAx.*Seperation)); 
PointsNewMvYN=(bsxfun(@minus,Points(:,3:4),CosAy.*Seperation)); 
PointsNewN=[PointsNewMvXN,PointsNewMvYN];
line([PointsNewMvXP(:,1)';PointsNewMvXP(:,2)'],[PointsNewMvYP(:,1)';PointsNewMvYP(:,2)'],'color','r')
line([PointsNewMvXN(:,1)';PointsNewMvXN(:,2)'],[PointsNewMvYN(:,1)';PointsNewMvYN(:,2)'],'color','b')
title('Animation of displacement, Red=ElNormSde'), xlabel('x'), ylabel('y'),axis equal;
xlS=xlim; %Got start lims
ylS=ylim;


%Final axis lims
%Getting the final state (want axis limits as this)
PointsNewMvXP=(bsxfun(@plus,PointsNewP(:,1:2),FinD_XY(:,1))); 
PointsNewMvYP=(bsxfun(@plus,PointsNewP(:,3:4),FinD_XY(:,2))); 
%Neg Dir
PointsNewMvXN=(bsxfun(@minus,PointsNewN(:,1:2),FinD_XY(:,1))); 
PointsNewMvYN=(bsxfun(@minus,PointsNewN(:,3:4),FinD_XY(:,2))); 
%Draw this
clf %Clear old fig
line([PointsNewMvXP(:,1)';PointsNewMvXP(:,2)'],[PointsNewMvYP(:,1)';PointsNewMvYP(:,2)'],'color','r')
line([PointsNewMvXN(:,1)';PointsNewMvXN(:,2)'],[PointsNewMvYN(:,1)';PointsNewMvYN(:,2)'],'color','b')
title('Animation of displacement, Red=ElNormSde'), xlabel('x'), ylabel('y'),axis equal;
xlF=xlim; %Got end lims	
ylF=ylim;	

%Put a breakpoint here to pause and change zoom if wanted
test=1;

%Getting the limits that will be used
xlB=[xlS;xlF];	ylB=[ylS;ylF];
xl=[min(xlB(:,1))-AxisLimPad,max(xlB(:,2))+AxisLimPad];
yl=[min(ylB(:,1))-AxisLimPad,max(ylB(:,2))+AxisLimPad];
xlim(xl);		
ylim(yl);	
%Turning of the axes if desired
if AxisOff==1
    axis off; %Turn of the axes
    set(gcf,'NumberTitle','off')%hide title
end



%% Now the loop starts where the animation takes place.

%Setting up some properties
sz=50;      %the amount of figures that will be drawn. This times the pause time will give the run time of the figure
inv = linspace(0,1,sz); %linear vector that slowly forces elements
Frame=getframe(gcf); %Start array 'Frame'
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB

%Running animation
for i=1:sz

    %Calculating the points of the elements as these move, this is along the
    %normal direction
    PointsNewMvXP=(bsxfun(@plus,PointsNewP(:,1:2),FinD_XY(:,1).*inv(i))); 
    PointsNewMvYP=(bsxfun(@plus,PointsNewP(:,3:4),FinD_XY(:,2).*inv(i))); 
    %Neg Dir
    PointsNewMvXN=(bsxfun(@minus,PointsNewN(:,1:2),FinD_XY(:,1).*inv(i))); 
    PointsNewMvYN=(bsxfun(@minus,PointsNewN(:,3:4),FinD_XY(:,2).*inv(i))); 

    %Restting the figure on each loop
    clf('reset')
    hold on
    line([PointsNewMvXP(:,1)';PointsNewMvXP(:,2)'],[PointsNewMvYP(:,1)';PointsNewMvYP(:,2)'],'color','r')
    line([PointsNewMvXN(:,1)';PointsNewMvXN(:,2)'],[PointsNewMvYN(:,1)';PointsNewMvYN(:,2)'],'color','b')
    scatter([PointsNewMvXP(:,1);PointsNewMvXP(:,2)],[PointsNewMvYP(:,1);PointsNewMvYP(:,2)],5,'r')
    scatter([PointsNewMvXN(:,1);PointsNewMvXN(:,2)],[PointsNewMvYN(:,1);PointsNewMvYN(:,2)],5,'b')
    
    %Setting the axis limits and creating the title
    xlabel('x'); ylabel('y');axis('equal');  
    xlim(xl)
    ylim(yl)
    title('Animation of displacement, Red=ElNormSde'),
    titlesz=14;
    fntsz=12;
    ChangeFontSizes(fntsz,titlesz);

    %Forcing this to draw
    drawnow
    %Pausing so the timing is not dependant on the machine
    pause(0.00005)

    if Save==1
        [ Frame ] = GrabFramesAndSaveFilm( i,Frame,isOctave );
    end    

end %end lp

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



