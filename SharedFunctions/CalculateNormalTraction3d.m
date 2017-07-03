function [ NormalTraction ] = CalculateNormalTraction3d( Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,CosAx,CosAy,CosAz )
%CalculateNormalTraction Calculates normal traction on plane
    %Stresses are col vectors
    %CosAx CosAy and CosAz are the direction cosines of the planes NORMAL

%   Copyright 2017, Tim Davis, The University of Aberdeen

    % Normal traction on the planes. 
    % Equation 6.49 Pollard and Fletcher Fundamentals
    AA=(bsxfun(@times,Pxx,CosAx.^2));
    BB=(bsxfun(@times,Pyy,CosAy.^2));
    CC=(bsxfun(@times,Pzz,CosAz.^2));
    DD=(bsxfun(@times,2*Pxy,CosAx.*CosAy));
    EE=(bsxfun(@times,2*Pyz,CosAy.*CosAz));
    FF=(bsxfun(@times,2*Pxz,CosAz.*CosAx));
    NormalTraction=AA+BB+CC+DD+EE+FF;

end

