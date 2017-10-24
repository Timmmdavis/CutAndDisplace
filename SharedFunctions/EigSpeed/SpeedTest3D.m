n=10000; %1000000;

Sxx=rand(n,1);
Syy=randn(n,1)*2;
Szz=ones(n,1)*3;
Sxy=ones(n,1);
Sxz=ones(n,1);
Syz=ones(n,1);

%My current func
[S1,S2,S3,S1dir,S2dir,S3dir] = EigCalc3d(Sxx,Syy,Szz,Sxy,Sxz,Syz);
disp('Done part 1')
%New func
[E1,E2,E3,E1dir,E2dir,E3dir] = EigCalc3dSpeed(Sxx,Syy,Szz,Sxy,Sxz,Syz);
disp('Done part 2')