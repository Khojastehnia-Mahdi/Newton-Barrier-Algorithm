% duplication matrix 
% This function is needed for the Newton_barrier algorithm.
% vec(A)=dup_n vech(A)
% input = n ( A is n*n)

function duplicationmatrix=dup_n(n)
    duplicationmatrix=zeros(n*n,n*(n+1)/2);
	for i=1:n
        for j=1:n
            if i>=j
                a=zeros(n*n,1);
                a((j-1)*n+i)=1;
                a((i-1)*n+j)=1;
                duplicationmatrix(:,((j-1)*(n-(j/2))+i))=a;
            end
        end
	end

end

