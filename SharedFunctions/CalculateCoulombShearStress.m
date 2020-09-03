function [ CSS ] = CalculateCoulombShearStress( Tn,Ts,Mu,C )
% CoulombShearStress: Calculates the Coulomb shear stress on a plane
%                   We assume tension positive convention and if the plane
%                   is not in contact neglect the resisting coefficient of
%                   friction force. Can be used for column vectors of
%                   inputs. 
%                   This script assumes tension positive convention. 
%                   For problems where you have stress tensors and a plane
%                   of known geometry use the script:
%                   "CalculateCoulombStressOnPlane".
%
%                   Equation from Pollard and Fletcher, 2005, Eq. 9.40. 
%
%               
% usage #1:
% [ CSS ] = CalculateCoulombShearStress( Tn,Ts,Mu,C )
%
% Arguments: (input)
% Tn                - The normal traction on the plane. (can be column vec
%                    of values)
%
% Ts                - The shear traction abs(maximal) on the plane. (can be
%                    column vec of values)
%
% Mu                - The coefficient of friction. (single value or col
%                     vec).
%
% C                 - The cohesive strength on the plane (in some texts
%                    called "sliding friction").
%
%
% Arguments: (output)
% CSS              - The value of Coulomb stress change. 
%
% Example usage 1:
%
% [ CSS ] = CalculateCoulombShearStress( -1,0.5,0.6,0.2 )
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Getting some sizes for the loop and initiating arrays to fill.
Sz=numel(Tn);
CSS=zeros(Sz,1);

%If C and Mu are single values we just repeat these here before the loop.
if isequal(size(C),([1,1]))
    Ones=ones(Sz,1);
    C=C.*Ones;
    Mu=Mu.*Ones;
end

%Looping and calculating for each input value of Tn and Ts etc. 
for i=1:Sz
    if Tn(i)>=C(i)  %The plane is not in contact.
        CSS(i)=Ts(i)-C(i);
    else            %Normal equation applies.
        CSS(i)=Ts(i)+(Mu(i).*Tn(i))-C(i);
    end
end

end

