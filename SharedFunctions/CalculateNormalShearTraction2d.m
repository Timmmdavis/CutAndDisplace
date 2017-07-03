function [ Tn,Ts ] = CalculateNormalShearTraction2d( Tx,Ty,CosAx,CosAy)
%CalculateTractionInChosenDirection Calculates traction on a plane in any
%direction the user chooses. The user just needs to supply the direction
%cosines of the direction they want this in. 
%   Tx Ty are cartesian traction components for the plane. I presume these are
%   column vectors. 
%   CosAx CosAy the direction cosines of the plane. Again I presume these are
%   column vectors. 

%   Copyright 2017, Tim Davis, The University of Aberdeen

    %Converts the traction XY to normal and shear traction components.
    %Eq 6.54 in Pollard, D.D. and Fletcher, R.C., 2005. Fundamentals of
    %structural geology. Cambridge University Press.
    Tn = (bsxfun(@times,Tx,CosAx))+(bsxfun(@times,Ty,CosAy));
    Ts = (bsxfun(@times,-Tx,CosAy))+(bsxfun(@times,Ty,CosAx));

end

