function [ Points,Triangles ] = GoCadAsciiReaderOld( Fid )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Open the text file.
fileID = fopen(Fid,'r');

if fileID == -1
error('The chosen Gocad file doesn't exist in the current working directories')
end



S = textscan(fileID,'%s');%'%*5c');
fclose(fileID);
G=S{1,1};             %Table created

% C = textscan(fid, '%s %n %n %n %n', 'delimiter', ',', ...
%                  'treatAsEmpty', {'NA', 'na'}, ...
%                  'commentStyle', '//');


Size=numel(G);
IndexG = strfind(G, 'VRTX');
Index = find(not(cellfun('isempty', IndexG)));
FirstIndex = Index(1,1);
G=G(FirstIndex:end,1);                 %Removing Header data

Size=numel(G);
IndexG = strfind(G, 'TRGL');
Index = find(not(cellfun('isempty', IndexG)));
FirstIndex = Index(1,1);                %Finding first instance of triangles

Pnts=G(1:(FirstIndex-1),1);             %Splitting into points and triangles
Trgls=G(FirstIndex:end,1); 
Size=numel(Trgls);
Trgls=Trgls(1:(end-1),1);               %Removes cell called 'end'

Pnts = reshape(Pnts,5,[]);
Pnts = Pnts(2:5,:);
Pnts = Pnts.';
Points = str2double(Pnts);

Trgls = reshape(Trgls,4,[]);
Trgls = Trgls(2:4,:);
Trgls = Trgls.';
Triangles = str2double(Trgls);

%Fixing a bug, this codebase relies on the fact the cols are the same as
%numbers the triangles are pointing to in the first row of 'Points'.
%Some program exports do not get this right. 
for i=1:numel(Points(:,1));
    if Points(i,1) ~= i          %if the row is not equal to the first column on that row. 
        
    %any numbers correponding to this - 'Points(i,1)' need to be changed
    %to this 'i' in triangles
    NumToChange=Points(i,1);
    ToChangeTo=i;
    Bad = Triangles == NumToChange;   
    Triangles(Bad)=ToChangeTo;
    Points(i,1)=ToChangeTo;
    end
end   

end

