% vech function
% This function is needed for the Newton_barrier algorithm.
% input : hermitian matrix A

function vech_output=vech(A)
    [mq,nq]=size(A);
	if mq==nq
        vech_outputt=[];
        for i=1:nq
            clear a
            for j=i:mq
                a(j-i+1)=A(j,i);
            end
            vech_outputt=[vech_outputt,a];
        end
        vech_output=conj(vech_outputt');
    else
        error('Input must be a square matrix.')
	end
        
end

