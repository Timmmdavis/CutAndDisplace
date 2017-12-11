%Script to add file paths. 

%If you change the name of the folder containing the scripts then change
%the variable 'FolderName' to match. 

%Change this is you change the top level directory name
FolderName='CutAndDisplace';

%Get the address of the current working directory
pathstring = pwd;      

%Splitting this into a cell array
if ispc
    parts = strsplit(mfilename('fullpath'), '\'); %Windows      
else
    parts = strsplit(mfilename('fullpath'), '/'); %Mac/Linux
end

%Finding the scripts root directory and its location (n)
if verLessThan('matlab', '9.1') %(2016b)
    [~,n] = find(~cellfun(@isempty,strfind(parts,FolderName)));
else %Cleaner version:
    [~,n] = find(contains(parts,FolderName));     
end

%Adding all folders to the path that are in this dir 
if ispc
    addpath(genpath(strjoin(parts(1,1:n),'\'))); %Windows  
else
    addpath(genpath(strjoin(parts(1,1:n),'/'))); %Mac/Linux
end

%Jumping back to the current directory
cd(pathstring)                                       

%Clear all the variables just created
clear