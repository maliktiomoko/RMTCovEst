% Question 1c Implementation of the MLE;
% Lire les donn�es
clear all
close all
clc
addpath('Manopt')
addpath('Othermethods')
addpath('Utilities')
addpath('Dataset')
%WHich applications? The synthetic one(synthetic) or the EEG dataset(eeg)
dataset='synthetic';
%Choose between LDA and QDA
application='QDA';
%Choose between Wishart and toeplitz for covariance 1 and 2
covariance1='Wishart';covariance2='Wishart';
%Number of simulations
% In the article we run over 10 simulations, but to be fast you can try
% single realization
n_simu=1;
switch dataset
    case 'synthetic'
        
        p_vec=floor(linspace(500,200,10));
        n=512;
        for i=1:length(p_vec)
            i
            M=2;p=p_vec(i);
            man = sympositivedefinitefactory(p);
            mu1=rand(1,p);
            if strcmp(application,'LDA')
                mu2=mu1+80/p;
            elseif strcmp(application,'QDA')
                if strcmp(covariance1,'toeplitz') && strcmp(covariance2,'toeplitz')
                    mu2=mu1+80/p;
                else
                    mu2=mu1+1/p;
                end
            end
            mu=[mu1;mu2];
            po=[0.5;0.5];
            fg1=zeros(p,2*p);
            fg2=zeros(p,2*p);
            erreur_cqv=[];erreur_LWqv=[];erreur_ourqv=[];
            error_cqv=[];error_LWqv=[];error_ourqv=[];
            for s=1:n_simu
                sigmas=zeros(p,p,2);
                for k=1:2*p
                    fg1(:,k) = mvnrnd(zeros(1,p),toeplitz(0.1.^(0:p-1)));
                    fg2(:,k) = mvnrnd(zeros(1,p),toeplitz(0.1.^(0:p-1)));
                end
                if strcmp(covariance1,'Wishart')
                    sigmas(:,:,1)=(fg1*fg1'/(2*p));
                elseif strcmp(covariance1,'toeplitz')
                    sigmas(:,:,1)=toeplitz(0.4.^(0:p-1));
                end
                if strcmp(covariance2,'Wishart')
                    sigmas(:,:,2)=(fg2*fg2'/(2*p));
                elseif strcmp(covariance2,'toeplitz')
                    sigmas(:,:,2)=toeplitz(0.4.^(0:p-1));
                end
                [ X,y ] = sim_misture_gaussian(n,mu,sigmas,M,po );plot_a=1;
                [ X_test,y_test ] = sim_misture_gaussian(n,mu,sigmas,M,po );
                y_test(y_test==2)=0;
                y(y==2)=0;
                if strcmp(application,'LDA')
                    [ omega_LDA_A,sigma_LDA_A,b_LDA_A,p_A,mu_0A,mu_1A,erreur_cv(s),erreur_LWav(s),erreur_LWbv(s),erreur_ourv(s) ] = LDA( X,y,X_test,y_test );
                elseif strcmp(application,'QDA')
                    [ erreur_cv(s),erreur_LWav(s),erreur_LWbv(s),erreur_ourv(s)] = QDA( X,y,X_test,y_test,sigmas );
                end
                erreur_cq(i)=mean(erreur_cv);erreur_LWa(i)=mean(erreur_LWav);erreur_LWb(i)=mean(erreur_LWbv);erreur_our(i)=mean(erreur_ourv);
            end
        end
            figure
            hold on
            plot(n./p_vec,erreur_cq,'r.-')
            plot(n./p_vec,erreur_LWa,'b*-')
            plot(n./p_vec,erreur_LWb,'k*-')
            plot(n./p_vec,erreur_our,'go-')
            legend('SCM','LW1','LW2','our')
    case 'eeg'
            load EEG_data.mat
            init_data = EEG_data;
            init_labels = EEG_labels;
            if strcmp(application,'LDA')
            n_train=150;n_test=1000;p=100;
            elseif strcmp(application,'QDA')
                n_train=500;n_test=1000;p=100;
            end
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
                    sigmas(:,:,1)=eye(p);sigmas(:,:,2)=eye(p);
                    erreur_cqv=[];erreur_LWqv=[];erreur_ourqv=[];
                    error_cqv=[];error_LWqv=[];error_ourqv=[];
                    if strcmp(application,'LDA')
                        [ omega_LDA_A,sigma_LDA_A,b_LDA_A,p_A,mu_0A,mu_1A,erreur_cv,erreur_LWav,erreur_LWbv,erreur_ourv ] = LDA( X,y,X_test,y_test);
                    elseif strcmp(application,'QDA')
                        [erreur_cv,erreur_LWav,erreur_LWbv,erreur_ourv ] = QDA( X,y,X_test,y_test,sigmas);
                    end

                    vec_err_c=[vec_err_c;erreur_cv];
                    vec_err_LWa=[vec_err_LWa;erreur_LWav];
                    vec_err_LWb=[vec_err_LWb;erreur_LWbv];
                    vec_err_our=[vec_err_our;erreur_ourv];
                end
            end
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
end