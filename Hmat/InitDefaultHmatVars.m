function [EPS,MaxRank,minn,n1,n2,n3,n4] = InitDefaultHmatVars( A )

if size(A) == size(A')
%If its a prime you will put in 1 and the prime as n1 and n2
SQRENUMCHK =sqrt(size(A,1));
PRIMECHK = isprime(size(A,1));
else
PRIMECHK=0;
SQRENUMCHK =0;
end


if size(A) == size(A')
	if any(PRIMECHK)
	n1=numel(A);
	n2=1;n3=n1;n4=n2;
	elseif isreal(SQRENUMCHK) && rem(SQRENUMCHK,1)==0
	n1=sqrt(size(A,1));n2=n1;n3=n1;n4=n1;
	else
	%if not
	D=[1; unique(cumprod(perms(factor(size(A,1))),2))];
	n1= D(ceil(end/2), :);n3=n1;
	n2= D(ceil(end/2)+1, :);n4=n2;
	end
else
	D=[1; unique(cumprod(perms(factor(size(A,1))),2))];
	n1= D(ceil(end/2), :);
	n2= D(ceil(end/2)+1, :);
	D=[1; unique(cumprod(perms(factor(size(A,2))),2))];
	n3= D(ceil(end/2), :);
	n4= D(ceil(end/2)+1, :);
end

EPS = 1e-6;

%%lower pres:
%MaxRank = 2;
%minn = 4;

%%Higher pres:
MaxRank=6;
minn=8;

%then just use: AA = HMatrix(A, [n1,n2],[n3,n4], 'S', [0,0], [0,0],  EPS, MaxRank, minn);

end
