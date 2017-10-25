function [ CSS ] = CalculateCoulombShearStress( Tn,Ts,Mu,C )
%CoulombShearStress Calculates the coulomb shear stress on a plane
%       We assume tension positive convention and if the plane is not in
%       contact neglect the resisting coefficient of friction force 
%       Equation from Pollard and Fletcher Book Eq 9.40 
%Tn is the normal traction on the plane
%Ts is the shear traction abs(maximal) on the plane
%Mu is the coefficient of friction
%C is the cohesion

Sz=numel(Tn);
CSS=zeros(Sz,1);
for i=1:Sz
    if Tn(i)>=C(i)  %The plane is not in contact
    CSS(i)=abs(Ts(i))-C(i);
    else            %Normal Eq applies
    CSS(i)=abs(Ts(i))+(Mu(i).*Tn(i))-C(i);
    end
end

end

