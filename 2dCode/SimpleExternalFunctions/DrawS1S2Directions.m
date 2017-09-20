function DrawS1S2Directions( X,Y,S1dir,S2dir,Points )
%DrawS1S2 If S1 and S2 principal directions have been calculated this
%function draw these as a quiver plot Red and Blue

%   Copyright 2017, Tim Davis, The University of Aberdeen
%Change Scl to make the lines longer or shorter
scl=0.25;
%Change Thick
thk=1;

%Extracting parts for visability
S1dirx=S1dir(:,1);
S1diry=S1dir(:,2);
S2dirx=S2dir(:,1);
S2diry=S2dir(:,2);


figure,
hold on
q1 = quiver(X,Y,S1dirx,S1diry,thk,'.'); 
q15 = quiver(X,Y,-S1dirx,-S1diry,thk,'.');
q2 = quiver(X,Y,S2dirx,S2diry,thk,'.');
q25 = quiver(X,Y,-S2dirx,-S2diry,thk,'.');
q1.Color = 'red';q15.Color = 'red';q2.Color = 'b';q25.Color = 'b';
q1.AutoScaleFactor=scl;q15.AutoScaleFactor=scl;q2.AutoScaleFactor=scl;q25.AutoScaleFactor=scl;

if nargin==5
line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','r')
end
axis equal,
xlabel('x'); ylabel('y'); axis('equal'); title('S1 (red) S2 (blue) directions'); 
hold off

end

