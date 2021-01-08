% inverse of vech function
% This function is needed for the Newton_barrier algorithm.
% note : vech(A)=xq and A is m*m
% inputs : xq and m

function  invvech_output=invvech(xq,m)
    invvech_output=zeros(m,m);
    for i=1:m
    	invvech_output((i-1)+1:end,i)=xq((i-1)*m-((i-2)*(i-1)/2)+1:i*m-((i-1)*i/2));
    end
	for i=1:m
        for j=1:m
            if i<j
                invvech_output(i,j)=conj(invvech_output(j,i));
            end
        end
	end
end
