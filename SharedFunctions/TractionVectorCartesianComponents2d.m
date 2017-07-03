function [ Tx,Ty ] = TractionVectorCartesianComponents2d( Pxx,Pyy,Pxy,CosAx,CosAy )
%TractionVectorCartesianComponents2d Calculates the cartesian components of the
%traction vector
    %Turning these stress col vecs to traction (vecs)
    %Stresses are col vectors
    %CosAx CosAy are the direction cosines of the planes NORMAL
    %Equation 6.40 and 6.41 in David Pollards Book. 
    %T1 [sxx sxy] * [n]
    %T2 [sxy syy] * [n]
    %Stresses are col vectors
    %CosAx CosAy are the direction cosines of the planes NORMAL
	%Tx and Ty are cartesian components of the 2D traction vec

%   Copyright 2017, Tim Davis, The University of Aberdeen

    %Converts the input stress to a traction XY.
    %Equation 6.55 in Pollard, D.D. and Fletcher, R.C., 2005. Fundamentals of
    %structural geology. Cambridge University Press.
    Tx = (bsxfun(@times,Pxx,CosAx))+(bsxfun(@times,Pxy,CosAy));
    Ty = (bsxfun(@times,Pxy,CosAx))+(bsxfun(@times,Pyy,CosAy));

end

