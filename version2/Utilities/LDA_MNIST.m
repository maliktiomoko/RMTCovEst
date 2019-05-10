% Question 1c Implementation of the MLE;
% Lire les donn�es
clear all
close all
clc
clear; clc;
addpath(genpath('utils'));
addpath(genpath('mnist'));

% Prepare data file
train_img_filename = 'mnist/train-images-idx3-ubyte';
train_lbl_filename = 'mnist/train-labels-idx1-ubyte';
test_img_filename = 'mnist/t10k-images-idx3-ubyte';
test_lbl_filename = 'mnist/t10k-labels-idx1-ubyte';
[train_image, train_label] = read_data(train_img_filename, train_lbl_filename, 10000, 0);
[test_image,  test_label] =  read_data(test_img_filename, test_lbl_filename, 2000, 0);
p=784;
for i=1
    i
X1=train_image(:,train_label==1);
X2=train_image(:,train_label==5);
X=[X1 X2]';n=size(X,2);
[features] = HOG_features(X);
y=[zeros(size(X1,2),1);ones(size(X2,2),1)];
X_test1=train_image(:,test_label==1);
X_test2=train_image(:,test_label==5);
X_test=[X_test1 X_test2]';
[features_test] = HOG_features(X_test);
plot_a=0;
y_test=[zeros(size(X_test1,2),1);ones(size(X_test2,2),1)];
[ omega_LDA_A,sigma_LDA_A,b_LDA_A,p_A,mu_0A,mu_1A,erreur_c(i),erreur_LW(i),erreur_our(i) ] = LDA( features,y,plot_a,features_test,y_test );
%[ error_cq(i),error_LWq(i),error_ourq(i)] = QDA( features,y',features_test,y_test',sigmas );
end
hold on
plot(n./p_vec,erreur_c,'r.-')
plot(n./p_vec,erreur_LW,'b*-')
plot(n./p_vec,erreur_our,'go-')
legend('SCM','LW','our')
figure
hold on
plot(n./p_vec,error_cq,'r.-')
plot(n./p_vec,error_LWq,'b*-')
plot(n./p_vec,error_ourq,'go-')
legend('SCM','LW','our')

vec=zeros(2*length(n./p_vec),1);
vec(1:2:end)=n./p_vec;
vec(2:2:end)=erreur_c;
sprintf('(%d,%d)',vec)
vec2=zeros(2*length(n./p_vec),1);
vec2(1:2:end)=n./p_vec;
vec2(2:2:end)=erreur_LW;
sprintf('(%d,%d)',vec2)
vec3=zeros(2*length(n./p_vec),1);
vec3(1:2:end)=n./p_vec;
vec3(2:2:end)=erreur_our;
sprintf('(%d,%d)',vec3)
vec4=zeros(2*length(n./p_vec),1);
vec4(1:2:end)=n./p_vec;
vec4(2:2:end)=error_cq;
sprintf('(%d,%d)',vec4)
vec5=zeros(2*length(n./p_vec),1);
vec5(1:2:end)=n./p_vec;
vec5(2:2:end)=error_LWq;
sprintf('(%d,%d)',vec5)
vec6=zeros(2*length(n./p_vec),1);
vec6(1:2:end)=n./p_vec;
vec6(2:2:end)=error_ourq;
sprintf('(%d,%d)',vec6)