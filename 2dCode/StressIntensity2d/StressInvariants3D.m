function [Sv,I1,I2,I3] = StressInvariants3D(Sxx,Syy,Szz,Sxy,Sxz,Syz)
%StressInvariants3D Calculate the vonmises stress and stress invariants
%from the stress tensors in 3D

%Von Mises stress: https://en.wikipedia.org/wiki/Von_Mises_yield_criterion:
%21 April 2020, at 13:01 (UTC).
Sv=sqrt(0.5*((Sxx+Syy)^2+(Syy+Szz)^2+(Szz+Szz)^2+6*(Sxy^2+Sxz^2+Syz^2)));
%Pollard and fletcher Eq.9.15
I1=(Sxx+Syy+Szz);
I2=-(Sxx.*Syy+Syy*Szz+Szz*Sxx)-Sxy.^2+Sxz^2+Syz^2;
I3=Sxx*Syy*Szz+2*(Sxy*Sxz*Syz)-Sxx*Syz^2-Syy*Sxz^2-Szz*Sxy^2;

end

