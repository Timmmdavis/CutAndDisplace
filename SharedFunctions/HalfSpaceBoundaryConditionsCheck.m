function  HalfSpaceBoundaryConditionsCheck(Zcoord,varargin)
% HalfSpaceBoundaryConditionsCheck: Checks that at the halfspace surface
%                   the condition for no normal/shear in 3D stress
%                   components is met. Throws an error if there is an
%                   issue. Note this assumes stresses vary linearly in the
%                   z-axis which may not be the case. This function pulling
%                   an hissy fit doesn't necessarily mean your inputs are
%                   wrong its just a decent first check that you have
%                   throught about your inputs.
%   
% usage #1:
% HalfSpaceBoundaryConditionsCheck3d(Z,Pzz,Pyz,Pxz)
%
% Arguments: (input)
% Zcoord            - The 'z' coordinate of the data. (typically a vector
%                     with element midpoints).
%
% varargin          - Stresses at these coordiantes. 
%
% Arguments: (output)
% N/A               - Throws error if there is an issue. No outputs.
%
% Example usage 1: No error:
%
%  Zcoord=-4:0.1:-3;
%  Pzz=-8:0.2:-6;
%  Pxz=zeros(size(Ycoord));
%  Pyz=zeros(size(Ycoord));
%  HalfSpaceBoundaryConditionsCheck(Zcoord,Pzz,Pyz,Pxz)
%
% Example usage 2: Normal/Shear tractions at free surface (2D). 
%
%  Ycoord=-4:0.1:-3;
%  Pyy=-4:0.2:-2;
%  Pxy=zeros(size(Ycoord));
%  HalfSpaceBoundaryConditionsCheck( Ycoord,Pyy,Pxy )
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Flips for every input
for i=1:numel(varargin)
    strVarName = inputname(i+1);
    Stress=varargin{i};
    if any(Stress~=0)
        if numel(Stress)==1 %no gradient so it can't reach 0 at the halfspace
            HSError %pass to internal func
        end
        % % %finding if the stress reduces to 0 at the free surface. 
        % % [InterceptStress]=GradCalc(Stress,Zcoord);
        % % InterceptStress = round(InterceptStress,9); %round as we can have tiny values
        % % if InterceptStress~=0
        % %     %Drawing stress with depth
        % %     scatter(Stress,Zcoord,'b'); hold on
        % %     scatter(InterceptStress,0,'r')
        % %     xlabel(strVarName); ylabel('depth');title('Red dot should be 0 stress') 
        % %     %Throwing error            
        % %     HSError
        % % end
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
disp ('Your boundary conditions have tractions that cannot exist on the half-space surface')  
%error('Redefine boundary conditions or remove function call')
end

