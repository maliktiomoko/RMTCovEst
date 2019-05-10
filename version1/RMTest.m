function [ C_est,cost ] = RMTest(x,C0,C,check_gradient,plot_cost,distance)
% This is a function to estimate the covariance matrix
% with the RMT estimator proposed in [CTZM18].
% To do so, we use a Riemannian gradient steepest-descent method on the
% symmetric positive definite manifold via the manopt toolbox.
% The manopt toolbox can be downloaded at: www.manopt.org
p=size(x,1);n=size(x,2);
switch distance
    case 'log'
        f=@(z) log(z);
    case 'Fisher'
        f=@(z) log(z).^2;
    case 'log1st'
        s=1;
        f=@(z) log(1+s*z);
    case 't'
        f=@(z) z;
    case 'Battacharrya'
        f=@(z) -(1/4)*log(z)+(1/2)*log(1+z)-(1/2)*log(2);
    case 'KL'
        f=@(z) -(1/2)*log(z)+(1/2)*z-1/2;
    case 'Wassertein'
    case 'Inverse_log'
        f=@(z) log(z);
    case 'Inverse_Fisher'
        f=@(z) log(z).^2;
    case 'Inverse_log1st'
        s=1;
        f=@(z) log(1+s*z);
    case 'Inverse_t'
        f=@(z) z;
    case 'Inverse_Battacharrya'
        f=@(z) -(1/4)*log(z)+(1/2)*log(1+z)-(1/2)*log(2);
    case 'Inverse_KL'
        f=@(z) -(1/2)*log(z)+(1/2)*z-1/2;
end
% SPD manifold
man = sympositivedefinitefactory(p);

%  compute the sample covaraince matrix
[U,Cs]  = eig(x*x'/n);
[cost_fun] = @(S) RMT_estim(S,Cs,n,distance);
rgrad    = @(S) RMT_estim_rgrad(S,Cs,n,distance);

% prepare problem structure for manopt
problem      = [];
problem.M    = man;
problem.cost = cost_fun;
problem.grad = rgrad;

% numerically check the gradient if asked (read explanations in the console)
if check_gradient==1
    checkgradient(problem);
end


options  = struct('verbosity',0,'statsfun',@statsfun,'maxiter',100);

%Steepest descent:We can try other(conjuguate gradient,pso)
[L_est,cost,info,opt_out] = steepestdescent(problem,C0,options);
C_est=U*L_est*U';

if plot_cost==1
   Matrix=[info.x];
   Matrix_reshape=reshape(Matrix,p,p,length([info.iter]));
if strcmp(distance,'Fisher') || strcmp(distance,'log') || strcmp(distance,'log1st') || strcmp(distance,'t') || strcmp(distance,'Battacharrya') ||...
   strcmp(distance,'KL')
   metric=@(D) mean(f(eig(D\C)));
else
  metric=@(D) mean(f(eig(D\(C^(-1)))));
end
Matrix_int=zeros(length([info.iter]),1);
for i=1:size(Matrix_reshape,3)
        Matrix_int(i)=metric(squeeze(Matrix_reshape(:,:,i)));
end
%Matrix_int(length([info.iter])+1)=metric(Cs^(-1));
 figure(2)
 plot(Matrix_int,'r*-','LineWidth',2,'MarkerSize',4)
 hold on
 plot((sqrt([info.cost])),'go-','LineWidth',2,'MarkerSize',4)
legend('real distance','cost function')
end
end