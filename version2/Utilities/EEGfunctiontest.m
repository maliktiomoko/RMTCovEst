% Question 1c Implementation of the MLE;
% Lire les donn�es
clear all
close all
clc
load EEG_data.mat
init_data = EEG_data;
init_labels = EEG_labels;
n_train=110;n_test=100;p=100;plot_a=0;
vec_err_c=[];vec_err_LWa=[];vec_err_LWb=[];
vec_err_our=[];
for i=1:5
    i
    if i==1 || i==2
        jvec=[3 4 5];
    else
        jvec=[1 2];
    end
    for j=jvec
X1=init_data(:,init_labels==i);X1=X1(:,1:n_train);
X2=init_data(:,init_labels==j);X2=X2(:,1:n_train);
X=[X1 X2]';
Xtest1=init_data(:,init_labels==i);Xtest1=Xtest1(:,n_train+1:n_train+n_test);
Xtest2=init_data(:,init_labels==j);Xtest2=Xtest2(:,n_train+1:n_train+n_test);
X_test=[Xtest1 Xtest2]';
y=[zeros(n_train,1);ones(n_train,1)]';
y_test=[zeros(n_test,1);ones(n_test,1)]';
man = sympositivedefinitefactory(p);
mu1=zeros(1,p);
mu2=mu1+80/p;
mu=[mu1;mu2];
po=[0.5;0.5];
fg1=zeros(p,2*p);
fg2=zeros(p,2*p);
 sigmas(:,:,1)=eye(p);sigmas(:,:,2)=eye(p);
erreur_cqv=[];erreur_LWqv=[];erreur_ourqv=[];
error_cqv=[];error_LWqv=[];error_ourqv=[];
%[ omega_LDA_A,sigma_LDA_A,b_LDA_A,p_A,mu_0A,mu_1A,erreur_cv,erreur_LWav,erreur_LWbv,erreur_ourv ] = LDA( X,y,plot_a,X_test,y_test,sigmas);
[erreur_cv,erreur_LWav,erreur_LWbv,erreur_ourv ] = QDA( X,y,X_test,y_test,sigmas);

vec_err_c=[vec_err_c;erreur_cv];
vec_err_LWa=[vec_err_LWa;erreur_LWav];
vec_err_LWb=[vec_err_LWb;erreur_LWbv];
vec_err_our=[vec_err_our;erreur_ourv];
    end
end
%[ error_cqv(s),error_LWqv(s),error_ourqv(s)] = QDA( X,y,X_test,y_test,sigmas );

%error_c(i)=mean(error_cqv);error_LW(i)=mean(error_LWqv);error_our(i)=mean(error_ourqv);

%[ error_cq(i),error_LWq(i),error_ourq(i)] = QDA( X,y,X_test,y_test,sigmas );
[new_c,ind]=sort(vec_err_c);
new_LWa=vec_err_LWa(ind);
new_LWb=vec_err_LWb(ind);
new_our=vec_err_our(ind);
figure
hold on
plot(new_c(1:2:end),'r.-')
plot(new_LWa(1:2:end),'b*-')
plot(new_LWb(1:2:end),'k*-')
plot(new_our(1:2:end),'go-')
legend('SCM','LW1','LW2','our')
% figure
% hold on
% plot(n./p_vec,error_c,'r.-')
% plot(n./p_vec,error_LW,'b*-')
% plot(n./p_vec,error_our,'go-')
% legend('SCM','LW','our')
vec=zeros(2*6,1);
vec(1:2:end)=1:6;
vec(2:2:end)=new_c(1:2:end);
sprintf('(%d,%d)',vec)

vec2a=zeros(2*6,1);
vec2a(1:2:end)=1:6;
vec2a(2:2:end)=new_LWa(1:2:end);
sprintf('(%d,%d)',vec2a)

vec2b=zeros(2*6,1);
vec2b(1:2:end)=1:6;
vec2b(2:2:end)=new_LWb(1:2:end);
sprintf('(%d,%d)',vec2b)

vec3=zeros(2*6,1);
vec3(1:2:end)=1:6;
vec3(2:2:end)=new_our(1:2:end);
sprintf('(%d,%d)',vec3)

% vec4=zeros(2*length(p_vec),1);
% vec4(1:2:end)=n./p_vec;
% vec4(2:2:end)=error_cq;
% sprintf('(%d,%d)',vec4)
% 
% vec5=zeros(2*length(p_vec),1);
% vec5(1:2:end)=n./p_vec;
% vec5(2:2:end)=error_LWq;
% sprintf('(%d,%d)',vec5)
% 
% vec6=zeros(2*length(p_vec),1);
% vec6(1:2:end)=n./p_vec;
% vec6(2:2:end)=error_ourq;
% sprintf('(%d,%d)',vec6)