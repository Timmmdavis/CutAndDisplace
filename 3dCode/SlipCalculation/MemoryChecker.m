
function MemoryChecker(NUM,DispInfCreate,infNum) 
%Checks the arrays are not so big the computer is going to start pushing in
%and out of memory and freeze. A matrix the size of the inf matrix size is
%created and its checked for the amount of space this will take up. 
%NUM*NUM * 6 is the size of each inf matrix. 6 as 6 tensors
%DispInfCreate is a flag, if its true then we are also creating XYZ disp
%inf matrices so each inf matrix is NUM*NUM * 9
%infNum, default is 3 inf matrices, instances where there are more and
%less exist. 

%   Copyright 2017, Tim Davis, The University of Aberdeen
if ~isunix %Need to add, otherwise fails in linux+mac

A = exist('infNum');
if A==0
    %infnum not brought in, we create the default
    infNum=3;
end    

%Calculating free memory before any matrices are created
[user,sys] = memory;
Orig=sys.PhysicalMemory.Available;   %Availible RAM

if DispInfCreate==0
%Creating one inf matrix of corrrect size
Dssinfmatrix = ones(NUM*NUM,6);
else
Dssinfmatrix = ones(NUM*NUM,9);
end

%Checking the sizes off all matrices, removing all values but the size of
%the matrix and name. 
a=whos;
clear Dssinfmatrix %free up memory again asap
field = 'size';
s = rmfield(a,field);
field = 'class'; s = rmfield(s,field);
field = 'global'; s = rmfield(s,field);
field = 'sparse'; s = rmfield(s,field);
field = 'complex'; s = rmfield(s,field);
field = 'nesting'; s = rmfield(s,field);
field = 'persistent'; two = rmfield(s,field);
field = 'bytes'; s = rmfield(two,field);
c = struct2cell(s);
loc=strfind(c,'Dssinfmatrix');
field = 'name'; bytes = rmfield(two,field);
Index = find(not(cellfun('isempty', loc)));
%extracting the size of the inf matrix and converting from cell to num
MatrixSize=bytes(Index); MatrixSize2=struct2cell(MatrixSize);MatrixSize=cell2mat(MatrixSize2);  %disp ('matrixsize'); disp (MatrixSize);

%The code uses 3 influence matrices but during allocation and calculations probably closer to 4. 
SmallestProbableSize=MatrixSize*infNum;
ProbableSize=MatrixSize*(infNum+1);
MemLeft=Orig-ProbableSize; %Calculating system memory left
MemLeftSmall=Orig-SmallestProbableSize; %Calculating system memory left

end

function BarChart(Orig,ProbableSize,SmallestProbableSize) 
%Setting up fig properties
barchrt=[Orig ProbableSize SmallestProbableSize]; 
barchrt=barchrt./10^9; %converting from bytes to GB 
x=[1 2 3];  
Line=[Orig Orig Orig];Line=Line./10^9;
plot(x,Line,'--','Color','red');hold on
Names=['YourRAM___________';'UpperInfMatrixSize';'LowerInfMatrixSize'];
bar(x,barchrt);
set(gca,'XTick',1:3,'XTickLabel',Names)
ylabel('RAM,GB'); title('FreeRAM Vs InfluenceMatrixSize');
hold off
end

if ~isunix %Need to add, otherwise fails in linux+mac

%Now doing if statement and pausing if the matrices are too big.     
if  Orig <= ProbableSize && Orig >= SmallestProbableSize
    disp 'Dangerously close to using all available RAM, if you are patient you MAY be ok';
    disp 'Approximation of the memory you will have, between:';
    disp(MemLeftSmall)
    disp 'and';
    disp(MemLeft)
    BarChart(Orig,ProbableSize,SmallestProbableSize) 
    f = warndlg('continue?, if you want to quit don’t close this, use (ctrl+c) in the cmd window');
    drawnow     % Necessary to print the message
    waitfor(f) %waiting until figure is closed
    disp('Carrying on');
    
elseif ProbableSize>Orig
    disp 'Influence matrcies are going to exceed RAM, use smaller mesh';
    disp 'Approximation of the memory you will have:';
    disp(MemLeft) %Printing size of memory left
    disp 'quit? (ctrl+c)';
    figure;BarChart(Orig,ProbableSize,SmallestProbableSize) 
    pause
    
else    
    %continuing uninterrupted
end
end

clear all

end
