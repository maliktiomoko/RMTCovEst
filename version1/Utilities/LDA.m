function [ omega_c,sigma_c,b_c,p,mu_0,mu_1,erreur_c,erreur_LWa,erreur_LWb,erreur_our ] = LDA( X,Y,x_test,y_test )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%Estimer les param�tres
n=length(Y);%nombre d'�chantillons
No=n-sum(Y);%nombre de Y==0
p1=sum(Y)./n;
p=size(X,2);
mu_0=sum(X(Y==0,:))./(No);
mu_1=sum(X(Y==1,:))./(sum(Y));
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
tauscaling=ts/mlambda1;
x_02=x02.*tauscaling;
%x_02=x02;
C02=diag(x_02);
%%%%Ledoit and Wolf
    [sigma_LW0,dhat,tauhat,speed,sigmahat20,dhat2,lambdah, ...
   lambdahat,exitflag,numiter,x0]=QuESTimate(X(Y==0,:),0,sigma_0);
[sigma_LW1,dhat,tauhat,speed,sigmahat21,dhat2,lambdah, ...
   lambdahat,exitflag,numiter,x0]=QuESTimate(X(Y==1,:),0,sigma_1);
%%%%%% Our method
[C01,a]=shrinkage_cov(X(Y==0,:),'rblw');
[C02,a]=shrinkage_cov(X(Y==1,:),'rblw');
    [ sigma_our0,cost ] = RMTest(X(Y==0,:)',C01,sigma_0,0,0,'Fisher');
    [ sigma_our1,cost ] = RMTest(X(Y==1,:)',C02,sigma_1,0,0,'Fisher');
%%%%%%%%%%%%%Calcul des  matrices de covariances
sigma_c=(No*sigma_0+sum(Y)*sigma_1)/n;
sigma_our=(No*sigma_our0+sum(Y)*sigma_our1)/n;
sigma_LWa=(No*sigma_LW0+sum(Y)*sigma_LW1)/n;
sigma_LWb=(No*sigmahat20+sum(Y)*sigmahat21)/n;

%%%%%%%%%%%%% Computation of w and b

omega_c=(mu_0-mu_1)*sigma_c^(-1);
omega_our=(mu_0-mu_1)*sigma_our^(-1);
 omega_LWa=(mu_0-mu_1)*sigma_LWa^(-1);
 omega_LWb=(mu_0-mu_1)*sigma_LWb^(-1);
b_c=0.5*(mu_1*sigma_c^(-1)*mu_1'-mu_0*sigma_c^(-1)*mu_0')+log((1-p1)/(p1));
b_our=0.5*(mu_1*sigma_our^(-1)*mu_1'-mu_0*sigma_our^(-1)*mu_0')+log((1-p1)/(p1));
 b_LWa=0.5*(mu_1*sigma_LWa^(-1)*mu_1'-mu_0*sigma_LWa^(-1)*mu_0')+log((1-p1)/(p1));
 b_LWb=0.5*(mu_1*sigma_LWb^(-1)*mu_1'-mu_0*sigma_LWb^(-1)*mu_0')+log((1-p1)/(p1));
pt=1000;
x=linspace(min(X(:,1)),max(X(:,1)),pt);
c_test=(omega_c*x_test'+b_c);y_c_test=zeros(size(c_test));
y_c_test(c_test<0)=1;y_c_test(c_test>0)=0;
ntest=size(x_test,1);
erreur_c=sum(y_c_test==y_test)./ntest;
LW_testa=(omega_LWa*x_test'+b_LWa);y_LW_testa=zeros(size(LW_testa));
y_LW_testa(LW_testa<0)=1;y_LW_testa(LW_testa>0)=0;
erreur_LWa=sum(y_LW_testa==y_test)./ntest;
LW_testb=(omega_LWb*x_test'+b_LWb);y_LW_testb=zeros(size(LW_testb));
y_LW_testb(LW_testb<0)=1;y_LW_testb(LW_testb>0)=0;
erreur_LWb=sum(y_LW_testb==y_test)./ntest;
our_test=(omega_our*x_test'+b_our);y_our_test(our_test>0)=0;
y_our_test(our_test<0)=1;y_our_test(our_test>0)=0;
erreur_our=sum(y_our_test==y_test)./ntest;
end

