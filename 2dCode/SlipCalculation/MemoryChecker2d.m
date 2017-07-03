
function MemoryChecker2d(NUM) 
%Checks the arrays are not so big the computer is going to start pushing in
%and out of memory and freeze. A matrix the size of the inf matrix size is
%created and its checked for the amount of space this will take up. 

%   Copyright 2017, Tim Davis, The University of Aberdeen
if ~isunix %Need to add, otherwise fails in linux+mac

%Calculating free memory before any matrices are created
[user,sys] = memory;
Orig=sys.PhysicalMemory.Available;   %Availible RAM

%Creating one inf matrix of corrrect size
SdInfMatrix = ones(NUM*NUM,5);

%Checking the sizes off all matrices, removing all values but the size of
%the matrix and name. 
a=whos;
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
loc=strfind(c,'SdInfMatrix');
field = 'name'; bytes = rmfield(two,field);
Index = find(not(cellfun('isempty', loc)));
%extracting the size of the inf matrix and converting from cell to num
MatrixSize=bytes(Index); MatrixSize2=struct2cell(MatrixSize);MatrixSize=cell2mat(MatrixSize2);  %disp ('matrixsize'); disp (MatrixSize);

%The code uses 3 influence matrices but during allocation and calculations probably closer to 4. 
SmallestProbableSize=MatrixSize*2;
ProbableSize=MatrixSize*3;
MemLeft=Orig-ProbableSize; %Calculating system memory left
MemLeftSmall=Orig-SmallestProbableSize; %Calculating system memory left

end %linux if loop

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

if ~isunix

%Now doing if statement and pausing if the matrices are too big.     
if  Orig <= ProbableSize && Orig >= SmallestProbableSize
    disp 'Dangerously close to using all available RAM, if you are patient you MAY be ok';
    disp 'Approximation of the memory you will have, between:';
    disp(MemLeftSmall)
    disp 'and';
    disp(MemLeft)
    disp 'continue? (press any key), quit? (ctrl+c)';
    BarChart(Orig,ProbableSize,SmallestProbableSize) 
    pause
    
elseif ProbableSize>Orig
    disp 'Influence matrcies are going to exceed RAM, use less elements';
    disp 'Approximation of the memory you will have:';
    disp(MemLeft) %Printing size of memory left
    disp 'quit? (ctrl+c)';
    BarChart(Orig,ProbableSize,SmallestProbableSize) 
    pause
    
else    
    %continuing uninterrupted
end

end

clear all
%OLD CODE
%Now creating a second large inf matrix to take up more ram. MATLAB has
% % %smart memory allocation so I need to do something to the matrix during 
% % %allocation to a new variable. 
% % DSinfmatrix=SSinfmatrix*3;
% % 
% % %Calculating memory used for 2 inf matrices
% % [user,sys] = memory;
% % Avail2=sys.PhysicalMemory.Available;     %With a second matrix size used
% % 
% % %Now calculating the size of the new matrix
% % MatrixSize=Avail-Avail2;

end
