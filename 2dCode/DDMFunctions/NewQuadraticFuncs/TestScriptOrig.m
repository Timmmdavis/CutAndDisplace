%Script to test the func 'Quad coeff func'
c

%For this lying at the point Ux is only ~0 when nu (pr) is 0.5
% x=linspace(-1,1,50);
% y=zeros(size(x))-0.0001;

%Using some thing like this will allow you to draw the disps correctly.
% x=linspace(-1,1,50);
% y=x;%zeros(size(x))-0;

 X = linspace(-3,3,60); %grid size
 Y = linspace(-3,3,60); 
 [x,y] = meshgrid(X,Y);  
 x=x(:)';
 y=y(:)';

xe=0;%[-0.928571428571429];
ye=0;
a=1;%0.0714; 
Beta=degtorad(145);

pr=0.5;
G=1;

Sd=[1,1,1];
Nd=[0,0,0];


[StressDisp] = ...
Quad_coeff_func(x,y,xe,ye,a,Beta,Sd,Nd,pr,G);

Sxy=StressDisp(:,3)';
Syy=StressDisp(:,2)';
Sxx=StressDisp(:,1)';

Ux=StressDisp(:,4)';
Uy=StressDisp(:,5)';
quiver(x,y,Ux,Uy);axis('equal')

%Calclating 2d EigenValues
[S1,S2,S1dir,S2dir]=EigCalc2d(Sxx,Syy,Sxy);

figure
scatter(x,y,15,S2);axis('equal');
DrawS1S2Directions(x(:),y(:),S1dir,S2dir )