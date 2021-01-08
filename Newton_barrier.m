% Newton_barrier algorithm
% capacity under the joint TPC+PAC constraints for a MIMO channel
% inputs: 
% W : Channel Gram matrix
% PT : total transmit power constraint
% P1 : per-antenna power constraints


function Newton_barrier_TPC_PAC(W,PT,P1)

% t00 : the initial value for t
t00=100;
t=t00;
t_barrier(1)=t;

% the value of t increases by the factor u
u=5;

diag_P1=diag(P1);
[m,~]=size(W);

% initial Tx covariance matrix
for i=1:m
    initial_entry_R(i)=(1/2)*min(P1(i),PT/m);
end
R_Newton{1}=diag(initial_entry_R);
x_Newton(:,1)=vech(R_Newton{1});

% the value of t_max
tmax=1e7;

zzz=1;
ii=1;

%the barrier method
while(t<tmax)

    R_before_newton{ii}=R_Newton{1};
    
    % gradient
    Z{1}=(((eye(m,m)+W*R_Newton{1}))^(-1))*W;
    diag_R{1}=diag(diag(R_Newton{1}));
    barrier_P1{1}=(diag_R{1}-diag_P1)^(-1);
    gradient_R_Ft{1}=Z{1}+(1/t)*((R_Newton{1})^(-1))+(1/t)*barrier_P1{1};
    gradient_x_Ft{1}=(dup_n(m))'*(veC(gradient_R_Ft{1})+(1/t)*...
    veC(eye(m,m))*(1/(trace(R_Newton{1})-PT)));
    
    x_Newton(:,1)=vech(R_Newton{1});
    r_w{1}=(gradient_x_Ft{1})';

    k=1;
    Norm_r_w_k=1;
    zzz=5*zzz*u;
    
    % Newton algorithm
    condition_newton=1;
    
    % accuracy of the Newton algorithm
    epsilon=10^(-8);
    while (condition_newton>epsilon)

        R_Newton{k};
        R_P1=R_Newton{k};
        
        % Hessian
        for i=1:m
            err=zeros(m,m);
            brr=zeros(m,m);
            err(i,i)=1/((R_P1(i,i)-P1(i))^2);
            brr(i,i)=1;
            er{i}=kron(brr,err);
        end
        [mk,nk]=size(er{i});
        b_P1=zeros(mk,nk);
        for i=1:m
            b_P1= b_P1+er{i};
        end
        hessian_x_Ft{k}=-(dup_n(m))'*(kron(Z{k},Z{k})+(1/t)...
        *kron((R_Newton{k})^(-1),(R_Newton{k})^(-1))+(1/t)*b_P1+...
        +(1/t)*(1/((PT-trace(R_Newton{k}))^(2)))*veC(eye(m,m))...
        *(veC(eye(m,m)))')*(dup_n(m));

        KKT_mat{k}=hessian_x_Ft{k};
        delta_w{k}=((-r_w{k})/KKT_mat{k});
        delta_x{k}=(delta_w{k})';
        
        
        % backtracking line search
        
        % residual norm is decreased by (a)
        a=0.3;
        
        % reduction in s is controlled by (B)
        B=0.5;
        
        s=1;
        condition_backtrack=10;
        Norm_r_new=100;
        aqq=0;
        while ((Norm_r_new>condition_backtrack)||(aqq>1))
    
            x_Newton(:,k+1)=x_Newton(:,k)+s*delta_x{k};
            R_Newton{k+1}=invvech(x_Newton(:,k+1),m);
            R_Newton{k+1};
            Z{k+1}=((eye(m,m)+W*R_Newton{k+1})^(-1))*W;
            diag_R{k+1}=diag(diag(R_Newton{k+1}));
            barrier_P1{k+1}=(diag_R{k+1}-diag_P1)^(-1);
            gradient_R_Ft{k+1}=Z{k+1}+(1/t)*((R_Newton{k+1})^(-1))+...
            (1/t)*(barrier_P1{k+1});
            gradient_x_Ft{k+1}=(dup_n(m))'*(veC(gradient_R_Ft{k+1})+...
            (1/t)*veC(eye(m,m))*(1/(trace(R_Newton{k+1})-PT)));
            r_w{k+1}=(gradient_x_Ft{k+1})';
            Norm_r_new=norm(r_w{k+1});
            condition_backtrack=(1-a*s)*norm(r_w{k});
            s=B*s;
            yuyu=eig(R_Newton{k+1});
            aqq=0;
            
            % checking the feasibility of covariance matrix
            if ((R_Newton{k+1})')==R_Newton{k+1}
                for jljl=1:length(yuyu)
                    if (yuyu(jljl))>0
                        aqq=0+aqq;
                    else
                         aqq=2+aqq;
                    end
                end
            else
                aqq=2.5;
            end

            new_R=R_Newton{k+1};
            
            for yryr=1:m
                if new_R(yryr,yryr)<=P1(yryr)
                    aqq=0+aqq;
                else
                    aqq=2+aqq;
                    R_Newton{k+1};
                end
            end

            if trace(new_R)<=PT
                aqq=0+aqq;
            else
                aqq=2+aqq;
                R_Newton{k+1};
            end


        end

        condition_newton=norm(r_w{k+1});
        residual_norm(k+1)=norm(r_w{k+1});
        norm(r_w{k+1});
        tyty(k)=norm(r_w{k+1});
        c_capacity(k)=log(det(eye(m,m)+W*R_Newton{k+1}));
        log(det(eye(m,m)+W*R_Newton{k+1}));
        k=k+1;
    end

cc_capacity{ii}=c_capacity;
clear c_capacity
ttt{ii}=tyty;
clear tyty
R_Newton{1}=R_Newton{k};
x_Newton(:,1)=vech(R_Newton{1});
t=u*t;
C(ii)=log(det(eye(m,m)+W*R_Newton{1}));
t_barrier(ii+1)=t;
ii=ii+1;
end
save('PAC_TPC_NB.mat')
end
