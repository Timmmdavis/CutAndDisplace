function [SXX,SYY,SXY,SRR,STT,SRT]=KirschSolutionFunction(X,Y,ri,pm,Sxx,Syy)

% Script from Pollard/Fletcher Fundamentals of structural geology
% Don't reproduce without the above statement ^^^
% fig_06_33
% plot stress components for circular hole in infinite plate
% Internal Pressure
% Jaeger and Cook (1979)
% equations (6.108) - (6.110)

if nargin<5
    Sxx=0;
    Syy=0;
end    

% Sxx=-Sxx;
% Syy=-Syy;

%ri = 1;%radius
%Sxx = 0;%Max remote
%Syy = 0;%Min remote

pm=pm(1,1);%InternalPressure

% x = linspace(-4,4,161)+eps; y = linspace(-4,4,161);
%[X,Y] = meshgrid(x,y);
[TH,R] = cart2pol(X,Y);

ST = sin(TH); S2T = sin(2*TH); ST2 = ST.^2; 
CT = cos(TH); C2T = cos(2*TH); CT2 = CT.^2;
R2 = (ri./R).^2; R4 = R2.^2;
% Polar stress components
SRR = -(0.5*(-Sxx-Syy)*(1-R2))+(pm*R2)-(0.5*(-Sxx+Syy)*((1-4*R2+3*R4).*C2T));
STT = -(0.5*(-Sxx-Syy)*(1+R2))-(pm*R2)+(0.5*(-Sxx+Syy)*((1+3*R4).*C2T));
SRT = 0.5*(-Sxx+Syy)*((1+2*R2-3*R4).*S2T);
Tol=0;%1e-4
SRR(R<ri-Tol) = nan; STT(R<ri-Tol) = nan; SRT(R<ri-Tol) = nan;

rowcount = size(SRR,1);
colcount = size(SRR,2);
[ SXX,SYY,SXY  ] = StressTensorTransformation2d(SRR(:),STT(:),SRT(:),cos(TH(:)),cos((pi/2)-TH(:)));
[SXX,SYY,SXY ]=ReshapeData2d( rowcount,colcount,SXX,SYY,SXY  );
% % Spherical to Cartesian stress components
% SXX = SRR.*CT2+STT.*ST2-2*SRT.*CT.*ST;
% SYY = SRR.*ST2+STT.*CT2+2*SRT.*CT.*ST;
% SXY = (SRR-STT).*CT.*ST+SRT.*(CT2-ST2);

% % Cartesian to Spherical stress components
% SRR2 = SXX.*CT2+SYY.*ST2+2*SXY.*CT.*ST;
% STT2 = SXX.*ST2+SYY.*CT2-2*SXY.*CT.*ST;
% SRT2 = ST.*CT.*(SYY-SXX)+SXY.*(CT2-ST2);

