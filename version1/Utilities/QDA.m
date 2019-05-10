function [ error_c,error_LWa,error_LWb,error_our] = QDA( X,Y,x_test,y_test,sigmas )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
sigmas0=squeeze(sigmas(:,:,1));
sigmas1=squeeze(sigmas(:,:,2));
n=length(Y);%nombre d'�chantillons
No=n-sum(Y);%nombre de Y==0
p1=sum(Y)./n;
p=size(X,2);
mu_0=sum(X(Y==0,:))./(No);
mu_1=sum(X(Y==1,:))./(sum(Y));
%mu=zeros(size(X));
%mu(Y==0,:)=repmat(mu_0,No,1);
%mu(Y==1,:)=repmat(mu_1,sum(Y),1);
sigma_0=(X(Y==0,:)-repmat(mu_0,No,1))'*(X(Y==0,:)-repmat(mu_0,No,1))/(No);
sigma_1=(X(Y==1,:)-repmat(mu_1,sum(Y),1))'*(X(Y==1,:)-repmat(mu_1,sum(Y),1))/(sum(Y));
y1=X(Y==0,:)-mean(X(Y==0,:));
Z1=y1.^2;k=0;
phi1=sum(sum(Z1'*Z1/(n-k)-2*(y1'*y1).*sigma_0/(n-k)+sigma_0.^2));
gama1=norm(sigma_0-mean(diag(sigma_0))*eye(p),'fro')^2;
shrinkage1=max(0,min(1,phi1/gama1/(n-k)));

% select starting point by linearly shrinking sample eigenvalues
mlambda1=mean(eig(sigma_0));ts=100;
x01=mlambda1+sqrt(1-shrinkage1).*(eig(sigma_0)-mlambda1);
tauscaling=ts/mlambda1;
x_01=x01.*tauscaling;
%x_01=x01;
C01=diag(x_01);
y2=X(Y==1,:)-mean(X(Y==1,:));
Z2=y2.^2;k=0;
phi2=sum(sum(Z2'*Z2/(n-k)-2*(y2'*y2).*sigma_1/(n-k)+sigma_1.^2));
gama2=norm(sigma_1-mean(diag(sigma_1))*eye(p),'fro')^2;
shrinkage2=max(0,min(1,phi2/gama2/(n-k)));

% select starting point by linearly shrinking sample eigenvalues
mlambda2=mean(eig(sigma_1));ts=100;
x02=mlambda2+sqrt(1-shrinkage2).*(eig(sigma_1)-mlambda2);
tauscaling=ts/mlambda2;
x_02=x02.*tauscaling;
%x_02=x02;
[U0,V0]=eig(sigma_0);[U1,V1]=eig(sigma_1);
C02=diag(x_02);
%%%%Ledoit and Wolf
    [sigma_LW0,dhat,tauhat,speed,sigmahat20,dhat2,lambdah, ...
   lambdahat,exitflag,numiter,x0]=QuESTimate(X(Y==0,:),0,sigma_0);
[sigma_LW1,dhat,tauhat,speed,sigmahat21,dhat2,lambdah, ...
   lambdahat,exitflag,numiter,x0]=QuESTimate(X(Y==1,:),0,sigma_1);
%%%%%% Methode classique
%%%%%% Linear shrinkage as initialization
[C01,a]=shrinkage_cov(X(Y==0,:),'oas');
[C02,a]=shrinkage_cov(X(Y==1,:),'oas');
C01=C01^(-1);C02=C02^(-1);
    [ sigma_our0,cost ] = RMTest(X(Y==0,:)',C01,sigmas0,0,0,'Inverse_Fisher');
    [ sigma_our1,cost ] = RMTest(X(Y==1,:)',C02,sigmas1,0,0,'Inverse_Fisher');

%%%%%%%%%%%%% Calcul des w et b

alpha_c=0.5*(sigma_0^(-1)-sigma_1^(-1));
alpha_LWa=0.5*(sigma_LW0^(-1)-sigma_LW1^(-1));
alpha_LWb=0.5*(sigmahat20^(-1)-sigmahat21^(-1));
alpha_our=0.5*(sigma_our0-sigma_our1);
beta_c=mu_1*sigma_1^(-1)-mu_0*sigma_0^(-1);
beta_LWa=mu_1*sigma_LW1^(-1)-mu_0*sigma_LW0^(-1);
beta_LWb=mu_1*sigmahat21^(-1)-mu_0*sigmahat20^(-1);
beta_our=mu_1*sigma_our1-mu_0*sigma_our0;
gamma_c=0.5*mu_0*sigma_0^(-1)*mu_0'-0.5*mu_1*sigma_1^(-1)*mu_1'-0.5*log(det(sigma_0\sigma_1))-log((1-p1)/p1);
gamma_LWa=0.5*mu_0*sigma_LW0^(-1)*mu_0'-0.5*mu_1*sigma_LW1^(-1)*mu_1'-0.5*log(det(sigma_LW0\sigma_LW1))-log((1-p1)/p1);
gamma_LWb=0.5*mu_0*sigmahat20^(-1)*mu_0'-0.5*mu_1*sigmahat21^(-1)*mu_1'-0.5*log(det(sigmahat20\sigmahat21))-log((1-p1)/p1);

gamma_our=0.5*mu_0*sigma_our0*mu_0'-0.5*mu_1*sigma_our1*mu_1'+0.5*log(det(sigma_our0\sigma_our1))-log((1-p1)/p1);
c_test=(diag(x_test*alpha_c*x_test')+(beta_c*x_test')'+gamma_c);
y_c_test(c_test>0)=1;y_c_test(c_test<0)=0;
ntest=size(x_test,1);
error_c=sum(y_c_test==y_test)/ntest;
LW_testa=(diag(x_test*alpha_LWa*x_test')+(beta_LWa*x_test')'+gamma_LWa);
y_LW_testa(LW_testa>0)=1;y_LW_testa(LW_testa<0)=0;
error_LWa=sum(y_LW_testa==y_test)/ntest;
LW_testb=(diag(x_test*alpha_LWb*x_test')+(beta_LWb*x_test')'+gamma_LWb);
y_LW_testb(LW_testb>0)=1;y_LW_testb(LW_testb<0)=0;
error_LWb=sum(y_LW_testb==y_test)/ntest;
our_test=(diag(x_test*alpha_our*x_test')+(beta_our*x_test')'+gamma_our);
y_our_test(our_test>0)=1;y_our_test(our_test<0)=0;
error_our=sum(y_our_test==y_test)/ntest;
end

