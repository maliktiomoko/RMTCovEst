error_lda_training=[0.67;3;5.75];
error_lda_test=[1.8;4.15;4.33];
error_log_reg_training=[0;2;4];
error_log_reg_test=[3.4;4.3;2.27];
error_lr_training=[1.33;3;5.5];
error_lr_test=[2.07;4.15;4.23];
figure
hold on
plot(1:3,error_lda_training,'ro--')
plot(1:3,error_lr_training,'bo--')
plot(1:3,error_log_reg_training,'go--')
plot(1:3,error_lda_test,'ro-')
plot(1:3,error_lr_test,'bo-')
plot(1:3,error_log_reg_test,'go-')
legend('training\_LDA','training\_linear','training\_log\_reg','test\_LDA','test\_linear','test\_log\_reg')
xlabel('Numéro du set utilisé (A=1,B=2,C=3)','Fontsize',15)
ylabel('Erreur de classification (en %)','Fontsize',15)
title('Erreur de classification','Fontsize',15)
error_qda_training=[0.67;1.33;5.25];
error_qda_test=[2;2;3.83];
figure
plot(1:3,error_qda_training,'bo-')
hold on
plot(1:3,error_qda_test,'go-')
h=legend('training\_QDA','test\_QDA');
set(h,'Fontsize',15)
xlabel('data (A=1,B=2,C=3)','Fontsize',15)
ylabel('erreur de classification (en %)','Fontsize',15)
title('Erreur de classification pour les données d''entrainement et de test','Fontsize',15)
%%
figure
hold on
plot(1:3,error_lda_training,'ro--')
plot(1:3,error_lr_training,'bo--')
plot(1:3,error_log_reg_training,'go--')
plot(1:3,error_lda_test,'ro-')
plot(1:3,error_lr_test,'bo-')
plot(1:3,error_log_reg_test,'go-')
error_qda_training=[0.67;1.33;5.25];
error_qda_test=[2;2;2.93];
plot(1:3,error_qda_training,'ko--')
hold on
plot(1:3,error_qda_test,'ko-')
legend('training\_LDA','training\_linear','training\_log\_reg','test\_LDA','test\_linear','test\_log\_reg','training\_QDA','test\_QDA')
xlabel('Numéro du set utilisé (A=1,B=2,C=3)','Fontsize',15)
ylabel('Erreur de classification (en %)','Fontsize',15)
title('Erreur de classification','Fontsize',15)