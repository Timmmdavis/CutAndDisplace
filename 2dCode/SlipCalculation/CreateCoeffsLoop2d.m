function [infmatrix,Dinfmatrix]=CreateCoeffsLoop2d(infmatrix,...
    Dinfmatrix,NUM,x,y,MidPoint,HalfLength,LineNormalVector,Ds,Dn,nu,E,halfspace,FD)
% CreateCoeffsLoop2d: Loop that calls the DDM functions of C&S
%             and fills a large column mat with influence coefficients.
%
% usage #1:
%[infmatrix]=CreateCoeffsLoop2d(infmatrix,...
%    NUM,x,y,MidPoint,a,LineNormalVector,Ds,Dn,nu,E,halfspace,StringHS,StringFS)
%
% Arguments: (input)
%  infmatrix  - 5*n 'vector' of 0's to be filled with coefficients. 
%
%     NUM     - The number of elements
%
%     x & y   - The element midpoints (moved slightly down normal direction)
%
%    MidPoint - The element midpoints 
%
% HalfLength  - An array of each each elements half length
%
%       nu    - The Poisson's ratio
%
%       E     - The Young's modulus
%
%  halfspace  - Defines if we work out the coefficients in a half or whole
%              space
%
% LineNormalVector - The direction cosines, CosAx (Nx) and CosAy in a list
%                   for each element. 
%
%      Ds,Dn  - Either one will be one or 0, depending if you are creating
%              a ds or dn inf matrix.
%
% Arguments: (output)
% DsTn,DsTs,DnTn,DnTs - Square influence matricies of much a displacement
%                   of one element (first part of name) effects the traction
%                   on another element (last part of name).
%
% DsUx,DsUy,DnUx,DnUy - Square influence matricies of much a displacement
%                   of one element (first part of name) effects the
%                   displacement at a point at the midpoint of another
%                   element (not the element displacement itself). 
%
% FD          - Flag telling the user if any elements are going to be fixed
%              (if this is the case displacement influence matricies are
%              filled).
%
%
% Example usage:
%
% [ DsTn,DsTs,DnTn,DnTs,DsUx,DsUy,DnUx,DnUy]...
% = CalculateInfluenceMatrices2d(halfspace,MidPoint,a,nu,E,LineNormalVector,Fdisp )
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

NUM2=size(MidPoint(:,1),1);

%Direction cosines 
CosAx=LineNormalVector(:,1);
CosAy=LineNormalVector(:,2); 
Beta=atan2(CosAx,-CosAy);
%ElementMidPoints
xe=MidPoint(:,1);
ye=MidPoint(:,2);

%Strings for progress
if Ds==1
    StringHS='1/2 CalculatingShearDispInfMatrixHS';
    StringFS='1/2 CalculatingShearDispInfMatrixFS';
else
    StringHS='2/2 CalculatingNormalDispInfMatrixHS';
    StringFS='2/2 CalculatingNormalDispInfMatrixFS';
end
    
%Starting progress bar function
if halfspace==1
    %Getting progress bar string
    progressbar(StringHS)  
else
    %Getting progress bar string
    progressbar(StringFS)          
end        

for i=1:NUM2
    %Creating size and space that will be filled in the influence
    %matrix in each loop
    first = (i-1)*NUM+1;
    last = i*NUM;
    if halfspace==1         %Calling modified TWODD function
        infmatrix(first:last,:) = LDstressHS(x,y,xe(i),ye(i),HalfLength(i),Beta(i),Ds,Dn,nu,E);
        if FD==1
            Dinfmatrix(first:last,:) = LDdispHS(x,y,xe(i),ye(i),HalfLength(i),Beta(i),Ds,Dn,nu);
        end
    else %Fullspace
        infmatrix(first:last,:) = LDstressFS(x,y,xe(i),ye(i),HalfLength(i),Beta(i),Ds,Dn,nu,E);
        if FD==1
            Dinfmatrix(first:last,:) = LDdispFS(x,y,xe(i),ye(i),HalfLength(i),Beta(i),Ds,Dn,nu);
        end
    end       
    %Updating progress bar
    progressbar(i/NUM2) 
end

end