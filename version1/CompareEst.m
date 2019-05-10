%%%%%%%%%%%%%%%%%%% Script comparing the different method%%%%%%%%%%%%%%
clear all
clc
close all
addpath('Manopt')
addpath('Othermethods')
addpath('Utilities')
addpath('Dataset')


%Parameters of simulation (For the figure in the article we run over 100 simulations but it is slow)
%But it is possible to run over 1 to see the trend.
n_simulation=1;
%n_vec=512; % Choose this for figure 1
n_vec=floor(linspace(200,500,10));
%n_vec=[512]; c=1.5;
p=190;
% Which distance(Fisher,Battacharrya,KL,Wassertein,Inverse,log,log1st,t)?
    distance='Fisher';
% Which target covariance (manopt_rand,Wishart,dirac,toeplitz,manual)
    covariance='Wishart';
    param.toeplitz=0.9;param.dirac=[0.10;1;3;4];param.wishart=2;
% Which initialisation(shrinkage,ledoit-wolf,alternative shrinkage,manual)?
    initialisation='shrinkage';
    %initialisation='manual';
% Do you want to check the gradient (1) or not(0)
    gradient_check=0;
% Do you want to plot the cost Vs the real distance(1) or not(0)
    plot_cost=0;
    
    
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
for b=1:length(n_vec)
    %p=p_vec(b);n=p*c;
     n=n_vec(b);
     c=p/n;
    for fg=1:n_simulation
    fg
    %Choose the target matrix to reach
    switch covariance
        case 'manopt_rand'
            %Either a random matrix from Manopt
            man = sympositivedefinitefactory(p);
            C=man.rand();
        case 'dirac'
            v=[];
            %Either a dirac
            for i=1:length(param.dirac)
                if floor(i/2)==i/2
                    v=[v;param.dirac(i)*ones(floor(p/length(param.dirac)),1)];
                else
                    v=[v;param.dirac(i)*ones(floor(p/length(param.dirac)),1)];
                end
            end
            %p-length(v)
            if length(v)<p
                v=[v;param.dirac(length(param.dirac))*ones(p-length(v),1)];
            end
            C=diag(v);
        case 'toeplitz'
            %Either a toeplitz
            C=toeplitz(param.toeplitz.^(0:p-1));
        case 'Wishart'
            %Either Wishart matrix
             rd=zeros(p,param.wishart*p);
             for k=1:param.wishart*p
                 rd(:,k) = mvnrnd(zeros(1,p),toeplitz(0.99.^(0:p-1)));
             end
             % compute the sample covaraince matrix
              C  = rd*rd'/(param.wishart*p);
        case 'manual'
            %Choose one SPD matrix
         
    end
    
    %Simulate sample from C
    x=zeros(p,n);df=1;
    for k=1:n
        x(:,k) = mvnrnd(zeros(1,p),C);
    end
    
    %Define the fisher metric distance
    if strcmp(distance,'Fisher') || strcmp(distance,'log') || strcmp(distance,'log1st') || strcmp(distance,'t') || ...
            strcmp(distance,'Battacharrya') || strcmp(distance,'KL')
        metric=@(D) mean(f(eig(D\C)));
    else
        metric=@(D) mean(f(eig(D\(C^(-1)))));
    end
    % QUEST method
    now2=tic;
    [sigmahat,dhat,tauhat,speed,sigmahat2,dhat2,lambdah, ...
   lambdahat,exitflag,numiter,x0]=QuESTimate(x',0);
    time_quest_vec(fg)=toc(now2);

   % Initialisation for our method
   switch initialisation
       case 'shrinkage'
            % compute linear shrinkage coefficient
            y=x';
            Z=y.^2;k=0;Cs=x*x'/n;
            phi=sum(sum(Z'*Z/(n-k)-2*(y'*y).*Cs/(n-k)+Cs.^2));
            gama=norm(Cs-mean(diag(Cs))*eye(p),'fro')^2;
            shrinkage=max(0,min(1,phi/gama/(n-k)));

            %select starting point by linearly shrinking sample eigenvalues
            mlambda=mean(lambdah);ts=100;
            tauscaling=ts/mlambda;
            x0=mlambda+sqrt(1-shrinkage).*(lambdah-mlambda);
            C0=diag(x0);
            if strcmp(distance,'Fisher') || strcmp(distance,'log') || strcmp(distance,'log1st') || strcmp(distance,'t') || ...
            strcmp(distance,'Battacharrya') || strcmp(distance,'KL')
                C0=C0;
            else
                C0=C0^(-1);
            end
       case 'ledoit-wolf'
            C0=sigmahat;
       case 'alternative shrinkage'
           if strcmp(distance,'Fisher') || strcmp(distance,'log') || strcmp(distance,'log1st') || strcmp(distance,'t') || ...
            strcmp(distance,'Battacharrya') || strcmp(distance,'KL')
           [C0,a]=shrinkage_cov(x','oas');
           else
               [C0,a]=shrinkage_cov(x','rblw');
               C0=C0^(-1);
           end
       case 'manual'
           C0=eye(p);
   end
       %Compute the theorical and empirical SCM distance D(C,\hat{C})=int fm_{MP}
    if strcmp(distance,'Fisher') || strcmp(distance,'log') || strcmp(distance,'log1st') || strcmp(distance,'t') || ...
            strcmp(distance,'Battacharrya') || strcmp(distance,'KL')
    [distance_mod_vec(fg),distance_real_vec(fg)] = SCM_theo(p,n,x,C,f);
    else
        Cs=x*x'/n;
        distance_real_vec(fg)=metric(Cs^(-1));
        distance_mod_vec(fg)=distance_real_vec(fg);
    end
   
   %Make call of our method
   now1=tic();
    [C_est,cost]=RMTest(x,C0,C,gradient_check,plot_cost,distance);
    time_our_vec(fg)=toc(now1);
    
   % Compute the D(C_LW,C)
   if strcmp(distance,'Fisher')||strcmp(distance,'Battacharrya')||strcmp(distance,'KL')||...
           strcmp(distance,'log')||strcmp(distance,'log1st')||strcmp(distance,'t')
    distance_ledoi_vec(fg)=metric(sigmahat);
    distance_ledoi2_vec(fg)=metric(sigmahat2);
   else
       distance_ledoi_vec(fg)=metric(sigmahat^(-1));
       distance_ledoi2_vec(fg)=metric(sigmahat2^(-1));
   end
    
   %COmput the distance D(C_our,C)
    distance_our_vec(fg)=metric(C_est);
    
    %Cramer Rao Bound
    crame_rao_unbiased_vec(fg)=(1/p)*((p*(p+1))/((n)));
    end
    distance_real(b)=mean(distance_real_vec);
    erreur_real(b)=max((distance_real_vec-distance_real(b)));
    distance_mod(b)=mean(distance_mod_vec);
    erreur_mod(b)=max((distance_mod_vec-distance_mod(b)));
    distance_ledoi(b)=mean(distance_ledoi_vec);
    sigma_ledoi(b)=sqrt(var(distance_ledoi_vec))./distance_ledoi(b);
    distance_ledoi2(b)=mean(distance_ledoi2_vec);
    sigma_ledoi2(b)=sqrt(var(distance_ledoi2_vec))./distance_ledoi2(b);
    distance_our(b)=mean(distance_our_vec);
    sigma_our(b)=sqrt(var(distance_our_vec))./distance_our(b);
    crame_rao_unbiased(b)=mean(crame_rao_unbiased_vec);
    time_our(b)=mean(time_our_vec);
    time_quest(b)=mean(time_quest_vec);
end
% Plot all the figure on loglog axis
figure(3)
hold on
loglog(n_vec./p,distance_mod,'g*-','LineWidth',2,'MarkerSize',4)
loglog(n_vec./p,distance_real,'ro-','LineWidth',2,'MarkerSize',4)
loglog(n_vec./p,distance_our,'k+-','LineWidth',2,'MarkerSize',8)
loglog(n_vec./p,crame_rao_unbiased,'o--','LineWidth',2,'MarkerSize',4)
loglog(n_vec./p,distance_ledoi,'m*-','LineWidth',2,'MarkerSize',4)
loglog(n_vec./p,distance_ledoi2,'b.-','LineWidth',2,'MarkerSize',4)
legend('SCM\_theo','SCM\_calc','our','Crame\_rao\_unbiased','Ledoit\_Peche1','Ledoit\_Peche2')
xlabel('log(n/p)','FontSize',15)
ylabel('error','FontSize',15)
vec1=zeros(2*length(n_vec),1);
vec1(1:2:end)=n_vec./p;
vec1(2:2:end)=distance_mod;
sprintf('(%d,%d)',vec1)

vec2=zeros(2*length(n_vec),1);
vec2(1:2:end)=n_vec./p;
vec2(2:2:end)=distance_real;
sprintf('(%d,%d)',vec2)


vec3=zeros(2*length(n_vec),1);
vec3(1:2:end)=n_vec./p;
vec3(2:2:end)=distance_our;
sprintf('(%d,%d)',vec3)

vec4=zeros(2*length(n_vec),1);
vec4(1:2:end)=n_vec./p;
vec4(2:2:end)=crame_rao_unbiased;
sprintf('(%d,%d)',vec4)

vec5=zeros(2*length(n_vec),1);
vec5(1:2:end)=n_vec./p;
vec5(2:2:end)=distance_ledoi;
sprintf('(%d,%d)',vec5)

vec6=zeros(2*length(n_vec),1);
vec6(1:2:end)=n_vec./p;
vec6(2:2:end)=distance_ledoi2;
sprintf('(%d,%d)',vec6)