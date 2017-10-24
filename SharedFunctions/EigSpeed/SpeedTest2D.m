n=1000000; %1000000;

Sxx=rand(n,1);
Syy=randn(n,1)*2;
Sxy=ones(n,1);


%My current func
[S1,S2,S1dir,S2dir] = EigCalc2d(Sxx,Syy,Sxy);
disp('Done part 1')
%New func
[E1,E2,E1dir,E2dir] = EigCalc2dSpeed(Sxx,Syy,Sxy);
disp('Done part 2')