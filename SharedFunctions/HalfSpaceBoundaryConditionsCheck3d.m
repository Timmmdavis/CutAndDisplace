function  HalfSpaceBoundaryConditionsCheck3d(Z,Pzz,Pyz,Pxz)
%Checks that at the halfspace surface the equilibrium condition for no z stress
%components is met. The components must be 0 at this surface. 
% In this func I assume that stresses increases linearly with depth, this may not be the case. See typical text books like
% Jaeger and cook. This function pulling an hissy fit doesn't necessarily mean your inputs are wrong. 
% its just a decent first check that you have throught about your inputs. 

%   Copyright 2017, Tim Davis, The University of Aberdeen


if any(Pzz~=0)
    if numel(Pzz)==1 %no gradient so it can't reach 0 at the halfspace
    HSError %pass to internal func
    end
    %finding if the stress reduces to 0 at the free surface. 
    [InterceptPzz]=GradCalc(Pzz,Z);
    InterceptPzz = round(InterceptPzz,9); %round as we can have tiny values
    if InterceptPzz~=0
        HSError
    end
end
    
if any(Pxz~=0)
    if numel(Pxz)==1 %no gradient so it can't reach 0 at the halfspace
    HSError %pass to internal func
    end
    %finding if the stress reduces to 0 at the free surface. 
    [InterceptPxz]=GradCalc(Pxz,Z);
    if InterceptPxz~=0
        HSError
    end
end

if any(Pyz~=0)
    if numel(Pyz)==1 %no gradient so it can't reach 0 at the halfspace
    HSError %pass to internal func
    end
    %finding if the stress reduces to 0 at the free surface. 
    [InterceptPyz]=GradCalc(Pyz,Z);
    if InterceptPyz~=0
        HSError
    end
end
    
end

function [YIntercept]=GradCalc(Y,X)
%subfunction to find the intercept of a straight line from two known points. C in Y=Mx+C  
    %getting the min and max Y values and the X values. 
    [Y1,Y1Index]=max(Y);
    [Y2,Y2Index]=min(Y);
    X1=X(Y1Index);
    X2=X(Y2Index);
    %calculating the gradient of the straight line
    grad=(Y2-Y1)/(X2-X1);
    %check if this meets the Y axis at 0
    YIntercept=Y1-(X1*grad);
end

function HSError
            disp ('your boundary conditions have tractions that cannot exist on the half space surface')  
            disp('The stress change from your solution will satisfy the condition that no Z tensors exist on the halfspace')
            disp('But the total stress if you are defining your remote as a far field stress will violate the equilibrium condition')
            disp('See McTigue, D.F. and Segall, P., 1988. Geophysical Research Letters, 15(6), pp.601-604')
            disp('Use remote stress conditions where z stress tensor components are 0 or grade to 0 where they meet the halfspace surface') 
            error('redefine boundary conds')
end

