function [ out_,out_1] = RMT_estim_rgrad(S,Cs,n2,distance)
%UNTITLED2 Summary of this function goes here

%   Gradient see papers Tiomoko, Couillet, Florent for more details
switch distance
    case 'Fisher'
        F=S\Cs;taille=size(S);p=taille(1);c2=p/n2;
        if sum(sum(isnan(F)))==0
        [K,V]=eig((F));
        lambda=sort(diag(V));
        zeta=sort((eig(diag(lambda)-(1/n2)*sqrt(lambda)*sqrt(lambda)')));
        out_prim=zeros(p,p);
        for i=1:length(lambda)
            lambdac=lambda;
            lambdac(i)=[];
        conn3=-1./lambda(i)+sum(1./(lambda(i)-zeta(1:end)))-sum(1./(lambda(i)-lambdac));
        c=-log(1-c2);
        out_prim(i,i)=(-2*conn3)*(log(lambda(i)))+(2+2*c)/lambda(i)+...
                    2*sum((1./(lambda(i)-zeta)).*log(zeta))-2*sum((1./(lambda(i)-lambdac)).*log(lambdac));
        end
        out_p=-(1/p)*out_prim;
        [diago,indi]=sort(diag(V)); diag_out=diag(out_p);
        for i=1:length(indi)
            pos=find(indi==i);
            vec(i)=diag_out(pos);
        end
        my_out=Cs*K*diag(vec)/K;
        out_1=real((my_out+transpose(my_out))/2);
        [a,r]=RMT_estim(S,Cs,n2,'Fisher');
        out_=2*r*sqrt(RMT_estim(S,Cs,n2,'Fisher'))*out_1;
        else
            out_=10*eye(p);
            out_1=out_;
        end
    case 'log'
        F=S\Cs;taille=size(S);p=taille(1);c2=p/n2;
        [K,V]=eig((F));
        lambda=sort(diag(V));
for r=1:length(lambda)
    out_prim(r,r)=1./lambda(r);
end
out_p=-(1/p)*out_prim;
        [diago,indi]=sort(diag(V)); diag_out=diag(out_p);
        for i=1:length(indi)
            pos=find(indi==i);
            vec(i)=diag_out(pos);
        end
        my_out=Cs*K*diag(vec)/K;
        out_1=real((my_out+transpose(my_out))/2);
        [a,r]=RMT_estim(S,Cs,n2,'log');
        out_=2*r*sqrt(RMT_estim(S,Cs,n2,'log'))*out_1;
    case 'log1st'
        F=S\Cs;taille=size(S);p=taille(1);c2=p/n2;
        [K,V]=eig((F));
        lambda=sort(diag(V));
m=@(z) c2*mean(1./(lambda-z))-(1-c2)./z;
    kappa_p=0;
    kappa_m=-10;
    while abs(kappa_p-kappa_m)>1e-7*abs(lambda(p)-lambda(1))
                        kappa_=(kappa_p+kappa_m)/2;
                        if 1-m(kappa_)<0
                            kappa_p=kappa_;
                        else
                            kappa_m=kappa_;
                        end
    end
kappa_0=(kappa_p+kappa_m)/2;
for r=1:length(lambda)
    out_prim(r,r)=1./(lambda(r)-kappa_0);
end
out_p=-(1/p)*out_prim;
        [diago,indi]=sort(diag(V)); diag_out=diag(out_p);
        for i=1:length(indi)
            pos=find(indi==i);
            vec(i)=diag_out(pos);
        end
        my_out=Cs*K*diag(vec)/K;
        out_1=real((my_out+transpose(my_out))/2);
        [a,r]=RMT_estim(S,Cs,n2,'log1st');
        out_=2*r*sqrt(RMT_estim(S,Cs,n2,'log1st'))*out_1;
    case 't'
        F=S\Cs;taille=size(S);p=taille(1);c2=p/n2;
        [K,V]=eig((F));
        lambda=sort(diag(V));
        zeta=sort((eig(diag(lambda)-(1/n2)*sqrt(lambda)*sqrt(lambda)')));
        mp=@(z) c2*mean(1./(lambda-z).^2)+(1-c2)./(z.^2);
        for r=1:length(lambda)
            lambdac=lambda;lambdac(r)=[];
            ro=0;
            for g=1:length(zeta)
                ro=ro-1./(mp(zeta(g))*((zeta(g)-lambda(r)).^2));
            end
            out_prim(r,r)=ro+p/c2;
        end
        out_p=-(1/p)*out_prim;
        [diago,indi]=sort(diag(V)); diag_out=diag(out_p);
        for i=1:length(indi)
            pos=find(indi==i);
            vec(i)=diag_out(pos);
        end
        my_out=Cs*K*diag(vec)/K;
        out_1=real((my_out+transpose(my_out))/2);
        [a,r]=RMT_estim(S,Cs,n2,'t');
        out_=2*r*sqrt(RMT_estim(S,Cs,n2,'t'))*out_1;
     case 'KL'
        F=S\Cs;taille=size(S);p=taille(1);c2=p/n2;
        [K,V]=eig((F));
        lambda=sort(diag(V));
        zeta=sort((eig(diag(lambda)-(1/n2)*sqrt(lambda)*sqrt(lambda)')));
mp=@(z) c2*mean(1./(lambda-z).^2)+(1-c2)./(z.^2);
        for r=1:length(lambda)
            lambdac=lambda;lambdac(r)=[];
            ro=0;
            for g=1:length(zeta)
                ro=ro-1./(mp(zeta(g))*((zeta(g)-lambda(r)).^2));
            end
            out_pr(r,r)=ro+p/c2;
        end
for r=1:length(lambda)
    out_prim(r,r)=(-1/2)*(1./lambda(r))+(1/2)*out_pr(r,r);
end
out_p=-(1/p)*out_prim;
        [diago,indi]=sort(diag(V)); diag_out=diag(out_p);
        for i=1:length(indi)
            pos=find(indi==i);
            vec(i)=diag_out(pos);
        end
        my_out=Cs*K*diag(vec)/K;
        out_1=real((my_out+transpose(my_out))/2);
        [a,r]=RMT_estim(S,Cs,n2,'KL');
        out_=2*r*sqrt(RMT_estim(S,Cs,n2,'KL'))*out_1; 
         case 'Battacharrya'
        s=1;
        F=S\Cs;taille=size(S);p=taille(1);c2=p/n2;
        [K,V]=eig((F));
        lambda=sort(diag(V));
m=@(z) c2*mean(1./(lambda-z))-(1-c2)./z;
    kappa_p=0;
    kappa_m=-10;
    while abs(kappa_p-kappa_m)>1e-7*abs(lambda(p)-lambda(1))
                        kappa_=(kappa_p+kappa_m)/2;
                        if 1-m(kappa_)<0
                            kappa_p=kappa_;
                        else
                            kappa_m=kappa_;
                        end
    end
kappa_0=(kappa_p+kappa_m)/2;
for r=1:length(lambda)
    out_prim(r,r)=(1/2)*(1./(lambda(r)-kappa_0))-(1/4*(1./lambda(r)));
end
out_p=-(1/p)*out_prim;
        [diago,indi]=sort(diag(V)); diag_out=diag(out_p);
        for i=1:length(indi)
            pos=find(indi==i);
            vec(i)=diag_out(pos);
        end
        my_out=Cs*K*diag(vec)/K;
        out_1=real((my_out+transpose(my_out))/2);
        [a,r]=RMT_estim(S,Cs,n2,'Battacharrya');
        out_=2*r*sqrt(RMT_estim(S,Cs,n2,'Battacharrya'))*out_1;
     case 'Inverse_Fisher'
        F=Cs*S;taille=size(S);p=taille(1);c2=p/n2;
        if sum(sum(isnan(F)))==0
        [K,V]=eig((F));
        lambda=sort(diag(V));
        zeta=sort((eig(diag(lambda)-(1/n2)*sqrt(lambda)*sqrt(lambda)')));
        
        
        out_prim=zeros(p,p);
        for i=1:length(lambda)
            lambdac=lambda;
            lambdac(i)=[];
        conn3=-1./lambda(i)+sum(1./(lambda(i)-zeta(1:end)))-sum(1./(lambda(i)-lambdac));
        c=-log(1-c2);
        out_prim(i,i)=(-2*conn3)*(log(lambda(i)))+(2+2*c)/lambda(i)+...
                    2*sum((1./(lambda(i)-zeta)).*log(zeta))-2*sum((1./(lambda(i)-lambdac)).*log(lambdac));
        end
        out_p=-(1/p)*out_prim;
        [diago,indi]=sort(diag(V)); diag_out=diag(out_p);
        for i=1:length(indi)
            pos=find(indi==i);
            vec(i)=diag_out(pos);
        end
        my_out=-S*(K*diag(vec)/K)*Cs*S;
        out_1=real((my_out+transpose(my_out))/2);
        [a,r]=RMT_estim(S,Cs,n2,'Inverse_Fisher');
        out_=2*r*sqrt(RMT_estim(S,Cs,n2,'Inverse_Fisher'))*out_1;
        else
            out_=10*eye(p);
            out_1=out_;
        end
    case 'Inverse_log'
        F=Cs*S;taille=size(S);p=taille(1);c2=p/n2;
        [K,V]=eig((F));
        lambda=sort(diag(V));
for r=1:length(lambda)
    out_prim(r,r)=1./lambda(r);
end
out_p=(1/p)*out_prim;
        [diago,indi]=sort(diag(V)); diag_out=diag(out_p);
        for i=1:length(indi)
            pos=find(indi==i);
            vec(i)=diag_out(pos);
        end
        my_out=-S*(K*diag(vec)/K)*Cs*S;
        out_1=real((my_out+transpose(my_out))/2);
        [a,r]=RMT_estim(S,Cs,n2,'Inverse_log');
        out_=2*r*sqrt(RMT_estim(S,Cs,n2,'Inverse_log'))*out_1;
    case 'Inverse_log1st'
        s=1;
        %altitude=1;
        F=Cs*S;taille=size(S);p=taille(1);c2=p/n2;
        [K,V]=eig((F));
        lambda=sort(diag(V));
m=@(z) c2*mean(1./(lambda-z))-(1-c2)./z;
    kappa_p=0;
    kappa_m=-10;
    while abs(kappa_p-kappa_m)>1e-7*abs(lambda(p)-lambda(1))
                        kappa_=(kappa_p+kappa_m)/2;
                        if 1-m(kappa_)<0
                            kappa_p=kappa_;
                        else
                            kappa_m=kappa_;
                        end
    end
kappa_0=(kappa_p+kappa_m)/2;
for r=1:length(lambda)
    out_prim(r,r)=-1./(lambda(r)-kappa_0)+1./lambda(r);
end
out_p=(1/p)*out_prim;
        [diago,indi]=sort(diag(V)); diag_out=diag(out_p);
        for i=1:length(indi)
            pos=find(indi==i);
            vec(i)=diag_out(pos);
        end
        my_out=-S*(K*diag(vec)/K)*Cs*S;
        out_1=real((my_out+transpose(my_out))/2);
        [a,r]=RMT_estim(S,Cs,n2,'Inverse_log1st');
        out_=2*r*sqrt(RMT_estim(S,Cs,n2,'Inverse_log1st'))*out_1;
    case 'Inverse_t'
        F=Cs*S;taille=size(S);p=taille(1);c2=p/n2;
        [K,V]=eig((F));
        lambda=sort(diag(V));
        for r=1:length(lambda)
            out_prim(r,r)=(1-c2)./(lambda(r).^2);
        end
        out_p=(1/p)*out_prim;
        [diago,indi]=sort(diag(V)); diag_out=diag(out_p);
        for i=1:length(indi)
            pos=find(indi==i);
            vec(i)=diag_out(pos);
        end
        my_out=-S*(K*diag(vec)/K)*Cs*S;
        out_1=real((my_out+transpose(my_out))/2);
        [a,r]=RMT_estim(S,Cs,n2,'Inverse_t');
        out_=2*r*sqrt(RMT_estim(S,Cs,n2,'Inverse_t'))*out_1;
    case 'Inverse_KL'
        F=Cs*S;taille=size(S);p=taille(1);c2=p/n2;
        [K,V]=eig((F));
        lambda=sort(diag(V));
for r=1:length(lambda)
    out_prim(r,r)=(-1/2)*(1./lambda(r))+(1/2)*((1-c2)./(lambda(r).^2));
end
out_p=(1/p)*out_prim;
        [diago,indi]=sort(diag(V)); diag_out=diag(out_p);
        for i=1:length(indi)
            pos=find(indi==i);
            vec(i)=diag_out(pos);
        end
        my_out=-S*(K*diag(vec)/K)*Cs*S;
        out_1=real((my_out+transpose(my_out))/2);
        [a,r]=RMT_estim(S,Cs,n2,'Inverse_KL');
        out_=2*r*sqrt(RMT_estim(S,Cs,n2,'Inverse_KL'))*out_1;
    case 'Inverse_Battacharrya'
        F=Cs*S;taille=size(S);p=taille(1);c2=p/n2;
        [K,V]=eig((F));
        lambda=sort(diag(V));
m=@(z) c2*mean(1./(lambda-z))-(1-c2)./z;
    kappa_p=0;
    kappa_m=-10;
    while abs(kappa_p-kappa_m)>1e-7*abs(lambda(p)-lambda(1))
                        kappa_=(kappa_p+kappa_m)/2;
                        if 1-m(kappa_)<0
                            kappa_p=kappa_;
                        else
                            kappa_m=kappa_;
                        end
    end
kappa_0=(kappa_p+kappa_m)/2;
for r=1:length(lambda)
    out_prim(r,r)=(1/2)*(-1./(lambda(r)-kappa_0)+1./lambda(r))-...
        (1/4)*(1./lambda(r));
end
out_p=(1/p)*out_prim;
        [diago,indi]=sort(diag(V)); diag_out=diag(out_p);
        for i=1:length(indi)
            pos=find(indi==i);
            vec(i)=diag_out(pos);
        end
        my_out=-S*(K*diag(vec)/K)*Cs*S;
        out_1=real((my_out+transpose(my_out))/2);
        [a,r]=RMT_estim(S,Cs,n2,'Inverse_Battacharrya');
        out_=2*r*sqrt(RMT_estim(S,Cs,n2,'Inverse_Battacharrya'))*out_1;
end
end