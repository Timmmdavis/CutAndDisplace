function [ Tx,Ty,Tz ] = TractionVectorCartesianComponents3d(  Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,CosAx,CosAy,CosAz )
%TractionVectorCartesianComponents Calculates the cartesian components of the
%traction vector
    %Turning these stress col vecs to traction (vecs)
    %Stresses are col vectors
    %CosAx CosAy and CosAz are the direction cosines of the planes NORMAL
	%Tx Ty and Tz are cartesian components of the 3D traction vec t(n)

%   Copyright 2017, Tim Davis, The University of Aberdeen
	
    %Equation 6.40 and 6.41 in David Pollards Book. 
    %T1 [sxx sxy sxz] * [n]
    %T2 [sxy syy syz] * [n]
    %T3 [sxz syz szz] * [n]
    Tx = (bsxfun(@times,Pxx,CosAx))+(bsxfun(@times,Pxy,CosAy))+(bsxfun(@times,Pxz,CosAz)); 
    Ty = (bsxfun(@times,Pxy,CosAx))+(bsxfun(@times,Pyy,CosAy))+(bsxfun(@times,Pyz,CosAz));
    Tz = (bsxfun(@times,Pxz,CosAx))+(bsxfun(@times,Pyz,CosAy))+(bsxfun(@times,Pzz,CosAz));

end

