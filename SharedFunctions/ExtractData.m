function [ varargout ] = ExtractData( Flag,type,varargin )
% ExtractData: Grabs columns or rows specified by a flag, works on 
%                   structure arrays too. 
%   
% usage #1:
% [ varargout ] = ExtractCols( Data )
%
% Arguments: (input)
% Flag              - A Flag that corresponds to the row or column size
%                    used for extraction or deletion.
%
% Type              - 1 = Extract rows specified by flag.
%                     2 = Remove  rows specified by flag.
%
% varargin          - Any number of arrays you are going to change using
%                    the flag. Can be arrays or structures. 
%
% Additional arguments: (input)
%
% 'Columns'         - Call "'Columns',1" in the varargin and you will extract
%                    columns with the flag not rows
%
% Scale             - Property changing the length of the 
%                    resultant ellipsoid. (Single value)  
%
% Arguments: (output)
% varargout       	- Edited arrays from the varargin. Best if the number
%                    of outputs matches the number of inputs.
%
% Example usage 1:
%
% %Grabbing each column from a matrix 'a':
% n=1000; 
% a=[rand(n,1),zeros(n,1),ones(n,1)];
% %Rows in first col of a that are below 0.5;
% Flag=a(:,1)<0.5;
% %Grab rows that meet the condition in the flag
% [ aCol1BelowPnt5 ] = ExtractData( Flag,1,a );
% MaxAndMinBelowPnt5=[min(aCol1BelowPnt5(:,1)),max(aCol1BelowPnt5(:,1))]
% %Grab rows that don't meet the condition in the flag
% [ aCol1AbovePnt5 ] = ExtractData( Flag,2,a );
% MaxAndMinAbovePnt5=[min(aCol1AbovePnt5(:,1)),max(aCol1AbovePnt5(:,1))]
% %Now we flip a to show how to act on cols
% a=a';
% %Rows in first col of a that are below 0.25;
% Flag=a(1,:)<0.5;
% [ aRow1BelowPnt5 ] = ExtractData( Flag,1,a,'Columns',1 );
% MaxAndMinBelowPnt5Rw=[min(aRow1BelowPnt5(1,:)),max(aRow1BelowPnt5(1,:))]
% 
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Checking if additional arguments have been input into the function:
[ varargin,Columns ]  = AdditionalArgsInVaragin( 'Columns', varargin,0 );

InputsSize=numel(varargin);

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


for i=1:numel(fieldnames(vararginExtra))
    
    %Grab the current varargin
    data=vararginExtra.(['X',num2str(i)]); 

        if type==1
        %We extract the data defined by the flag

        if Columns==1
            if numel(Flag)~=numel(data(1,:))
                error('Input column sizes must match flag size')
            end
        %Extract the columns based on the flag
        data=data(:,Flag); 
        else
            if numel(Flag)~=numel(data(:,1))
                error('Input row sizes must match flag size')
            end
        %Get the rows (default)
        data=data(Flag,:);
        end

    elseif type==2
    %We remove the data defined by the flag        

        if Columns==1
            if numel(Flag)~=numel(data(1,:))
                error('Input column sizes must match flag size')
            end            
        %Remove the columns based on the flag            
        data=data(:,~Flag); 
        else
            if numel(Flag)~=numel(data(:,1))
                error('Input row sizes must match flag size')
            end            
        %Remove the rows (default)            
        data=data(~Flag,:);
        end

    end
    %Set the current varargin to go out. 
    vararginExtra.(['X',num2str(i)])=data;     
    
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

