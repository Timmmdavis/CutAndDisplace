function [ TractionInDirection ] = CalculateTractionInChosenDirection3d( Tx,Ty,Tz,CosAx,CosAy,CosAz,ChosenDirectionCos )
%CalculateTractionInChosenDirection Calculates traction on a plane in any
%direction the user chooses. The user just needs to supply the direction
%cosines of the direction they want this in. 
%   Tx Ty and Tz are cartesian traction components for the plane. I presume these are
%   column vectors. 
%   CosAx CosAy and CosAz are the direction cosines of the plane. Again I presume these are
%   column vectors. 
%   ChosenDirectionCos is a 3*n matrix of the direction you want to know
%   the traction in. the first col is the CosAx component etc

%   Copyright 2017, Tim Davis, The University of Aberdeen

    %Pollard eq 6.52
    one=   (bsxfun(@times,(bsxfun(@times,Tx,(1-(CosAx.^2))))- (bsxfun(@times,Ty,CosAx.*CosAy))  -  (bsxfun(@times,Tz,CosAx.*CosAz))  ,ChosenDirectionCos(:,1)));
    two=   (bsxfun(@times,(bsxfun(@times,Tx,-CosAx.*CosAy)) + (bsxfun(@times,Ty,(1-(CosAy.^2))))-  (bsxfun(@times,Tz,CosAy.*CosAz))  ,ChosenDirectionCos(:,2)));
    three= (bsxfun(@times,(bsxfun(@times,Tx,-CosAz.*CosAx)) - (bsxfun(@times,Ty,CosAz.*CosAy))  +  (bsxfun(@times,Tz,(1-(CosAz.^2)))),ChosenDirectionCos(:,3)));
    TractionInDirection=one+two+three; 

end

