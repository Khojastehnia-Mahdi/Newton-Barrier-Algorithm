% vec function
% This function is needed for the Newton_barrier algorithm.
% input: matrix A

function vec_output=veC(A)
	[m,n]=size(A);
	if m==n
        for i=1:n
            for j=1:m
                vec_outputt(j+(i-1)*n)=A(j,i);
            end
        end
        vec_output=conj(vec_outputt');
       
	end
end

