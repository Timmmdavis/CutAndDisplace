function [varargout]=FlipColVecs(varargin)
% FlipColVecs: Pass as many arguments as you want in and get the same
%                   amount out, just flipped (up down). The inputs can be
%                   different sizes and it can deal will cell arrays. 
%   
% usage #1:
% [varargout]=FlipColVecs(varargin)
%
% Arguments: (input)
% varargin          - Vector data arguments.
%
% Arguments: (output)
% varargout       	- Vector data flipped so base value (end) is now the
%                     top value (1).
%
% Example usage 1:
%
% %Grabbing each column from a matrix 'a':
% VectorA=(1:1000)';
% VectorB=(1:2:1000)';
% figure;hold on; plot(VectorA); plot(VectorB);
% [VectorA,VectorB]=FlipColVecs(VectorA,VectorB)
% plot(VectorA); plot(VectorB);
% xlabel('Row number in vector'); ylabel('Value'); 
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Gets the inputs size
InputsSize=nargin;

vararginExtra=struct;
WholeArraySz=1; %max(size(varargin));
StrucSzs=[]; %Data on how big each part of the varargin was
Names=[];
for i=1:InputsSize
    %Grab the current varargin
    data=varargin{i};  
    %If the current is a structure we are going to do this to every value
    %in here:
    if isstruct(data)
        fnames = fieldnames(data);
        StrucSz=numel(fnames);
        WholeArraySz=WholeArraySz+StrucSz;
        for j=1:StrucSz
            ArrayNm=fnames{j};
            vararginExtra.(['X',num2str(WholeArraySz-j)])=data.(ArrayNm);
        end
        fnames=flipud(fnames); %So the correct names come out
    else
        
        StrucSz=1;
        fnames=[];
        vararginExtra.(['X',num2str(WholeArraySz)])=data;
        WholeArraySz=WholeArraySz+1;
    end
    StrucSzs=[StrucSzs,StrucSz];
    Names=[Names,fnames];
       
end

%Flips for every input
for i=1:numel(fieldnames(vararginExtra))
    %Grab the current varargin
    data=vararginExtra.(['X',num2str(i)]); 
    vararginExtra.(['X',num2str(i)])=flipud(data);
end

StructureArray=struct;
for i=1:InputsSize
    
    %The amount of args up to this pnt. 
    DataIndx=sum(StrucSzs(1:i-1));
  
    if StrucSzs(i)==1
    varargin{i}=vararginExtra.(['X',num2str((DataIndx+1))]);
    else
        for j=1:StrucSzs(i)
            StructureArray.(Names{DataIndx+j})=vararginExtra.(['X',num2str((DataIndx)+j)]);
        end
    varargin{i}=StructureArray;
    StructureArray=struct; %Resetting
    end
    
end
 
%Assign data to outputs
varargout=varargin;

end

