function MemoryChecker(Num,Reps) 
%MemoryChecker
%Checks the arrays are not so big the computer is going to start pushing in
%and out of memory and freeze. A matrix the size of the inf matrix size is
%created and its checked for the amount of space this will take up. 

%   Copyright 2017, Tim Davis, The University of Aberdeen

if ~isunix %Need to add, otherwise fails in linux+mac

%Calculating free memory before any matrices are created
[~,sys] = memory;
Orig=sys.PhysicalMemory.Available;   %Availible RAM

%Creating one inf matrix of corrrect size
InfMatrix = ones(Num^2,1);

%Getting the list of variables in the func
a=whos;
%Assuming that our inf matrix is the largest mat
MatrixSize = max(getfield(a, 'bytes'));

%The code uses the amount of influence matrices defined in Reps but during
%allocation and calculations probably closer to Reps+1.
SmallestProbableSize=MatrixSize*Reps;
ProbableSize=(MatrixSize*Reps)+1; 
MemLeft=Orig-ProbableSize; %Calculating system memory left
MemLeftSmall=Orig-SmallestProbableSize; %Calculating system memory left

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
    disp 'Influence matrices are going to exceed RAM, use less elements';
    disp 'Approximation of the memory you will have:';
    disp(MemLeft) %Printing size of memory left
    disp 'quit? (ctrl+c)';
    BarChart(Orig,ProbableSize,SmallestProbableSize) 
    pause
    
else    
    %continuing uninterrupted
end

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

%Clears mats within the func
clear InfMatrix

end
