function MemoryCheckerOctave(NUM,Reps) 
%Checks the arrays are not so big the computer is going to start pushing in
%and out of memory and freeze. I am not sure how Octave Ram Allocation
%works and where vars are stored. This may not be as crutial as MATLAB on Windows

%   Copyright 2017, Tim Davis, The University of Aberdeen
%Creating one inf matrix of corrrect size
infmatrix = ones(NUM^2,1);

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
loc=strfind(c,'infmatrix');
field = 'name'; bytes = rmfield(two,field);
Index = find(not(cellfun('isempty', loc)));
MatrixSize=bytes(Index); 
MatrixSize2=struct2cell(MatrixSize);
matrixsize=cell2mat(MatrixSize2);  disp ('matrixsize'); disp (matrixsize);
estMem=matrixsize*Reps;
jmem = javamem ();
FreeMem=jmem(3,1);
Freemem=cell2mat(FreeMem);
MemLeft=Freemem-matrixsize;

if estMem>Freemem
    disp 'Size of influence matrices will be bigger than Octaves free memory, may or may not be an issue';
    disp 'Approximation of the memory you will have:';
    disp(MemLeft); %Printing size of memory left
    disp('https://www.gnu.org/software/octave/doc/v4.0.0/How-can-I-handle-memory-limitations_003f.html');
    disp 'continue (type 'continue' into the cmd window), quit? (press: ctrl+c)';
    pause

else    

end

end
