function [Sv,I1,I2] = StressInvariants2D(Sxx,Syy,Sxy)
%StressInvariants2D Calculate the vonmises stress and stress invariants
%from the stress tensors in 2D

%Von Mises stress: https://en.wikipedia.org/wiki/Von_Mises_yield_criterion:
%21 April 2020, at 13:01 (UTC).
Sv=sqrt(Sxx.^2-(Sxx.*Syy)+Syy.^2+(3.*Sxy.^2));
%https://www.continuummechanics.org/principalstress.html : 2D invariants
I1=(Sxx+Syy);
I2=(Sxx.*Syy)-Sxy.^2;

end

