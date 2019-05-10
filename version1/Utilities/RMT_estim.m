function [ out,r ] = RMT_estim( S,Cs,n2,distance )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

switch distance
    case 'Fisher'
        p=size(S,1);
        F=Cs./S;
%         if sum(sum(isnan(F)))==0
        c2=p/n2;
        lambda=sort(diag(F));
        %zeta=sort((eig(diag(lambda)-(1/n2)*sqrt(lambda)*sqrt(lambda)')));
        lambdat=[0;lambda];
        zeta=zeros(p,1);
        for i=1:p
            kappa_p=lambdat(i+1);
            kappa_m=lambdat(i);
                while abs(kappa_p-kappa_m)>1e-6*abs(lambdat(1)-lambdat(end))
                    zeta_=(kappa_p+kappa_m)/2;
                    if (1/n2)*sum(lambda./(lambda-zeta_))<1
                        kappa_m=zeta_;
                    else
                        kappa_p=zeta_;
                    end
                end
                zeta(i)=(kappa_p+kappa_m)/2;  
        end
        ker_vec=(2/p)*((lambda./lambda').*log(abs(lambda./lambda'))./(1-lambda./lambda'));
        ker_vec(isnan(ker_vec))=0;
        ker=sum(sum((ker_vec)))-2;
        out=real((2/p)*sum(sum(((zeta./lambda').*log(abs(zeta./lambda'))./(1-zeta./lambda'))))-ker+(2/p)*sum(log(lambda))...
            -(1-c2)/c2*(log(1-c2)^2+sum(log(lambda).^2-log(zeta).^2))...
            -1/p*(2*sum(sum(dilog2(1-(zeta*ones(1,p))./(ones(p,1)*lambda'))-dilog2(1-(lambda*ones(1,p))./(ones(p,1)*lambda'))))-sum(log(lambda).^2)));
        r=sign(out);
        out=out.^2;
%         else
%             out=1000;
%         end
    case 'log'
        p=size(S,1);
        F=S\Cs;
        c2=p/n2;
        lambda_hatC2 = sort(eig(F));
        out=1/p*sum(log(lambda_hatC2))+(1-c2)/c2*log(1-c2)+1;
        r=sign(out);
        out=out.^2;
    case 't'
        F=S\Cs;
        lambda_hatC2 = sort(eig(F));
        out1=mean(lambda_hatC2);
        r=sign(out1);
        out=out1.^2;
    case 'log1st'
        s=1;
        p=size(S,1);
        F=S\Cs;c2=p/n2;
        lambda=sort(real(eig(F)));
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
        kappa_0_hatC2=(kappa_p+kappa_m)/2;
        out=(1+s*kappa_0_hatC2+log(abs(-s*kappa_0_hatC2)))/c2+1/p*sum(log(abs(1-lambda/kappa_0_hatC2)));
        r=sign(out);
        out=out.^2;
    case 'KL'
        p=size(S,1);
        F=S\Cs;
        c2=p/n2;
        lambda_hatC2 = sort(eig(F));
        out=-1/2*(1/p*sum(log(lambda_hatC2))+(1-c2)/c2*log(1-c2)+1)-(1/2)+(1/2)*(mean(lambda_hatC2));
        r=sign(out);
        out=out.^2;
    case 'Battacharrya'
        s=1;
        p=size(S,1);
        F=S\Cs;c2=p/n2;
        lambda=sort(real(eig(F)));
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
        kappa_0_hatC2=(kappa_p+kappa_m)/2;
        out=(1/2)*((1+s*kappa_0_hatC2+log(abs(-s*kappa_0_hatC2)))/c2+1/p*sum(log(abs(1-lambda/kappa_0_hatC2))))-...
            (1/4)*(1/p*sum(log(lambda))+(1-c2)/c2*log(1-c2)+1)-(1/2)*log(2);
        r=sign(out);
        out=out.^2;
    case 'Inverse_Fisher'
        p=size(S,1);
        F=Cs*S;
        if sum(sum(isnan(F)))==0
        c2=p/n2;
        lambda=sort(real(eig(F)));
        zeta=sort((eig(diag(lambda)-(1/n2)*sqrt(lambda)*sqrt(lambda)')));
        c1=0;
        ker_vec=(2/p)*((lambda./lambda').*log(abs(lambda./lambda'))./(1-lambda./lambda'));
        ker_vec(isnan(ker_vec))=0;
        ker=sum(sum((ker_vec)))-2;
        out=real((2/p)*sum(sum(((zeta./lambda').*log(abs(zeta./lambda'))./(1-zeta./lambda'))))-ker+(2/p)*sum(log(lambda))...
            -(1-c2)/c2*(log(1-c2)^2-log(1-c1)^2+sum(log(lambda).^2-log(zeta).^2))...
            -1/p*(2*sum(sum(dilog2(1-(zeta*ones(1,p))./(ones(p,1)*lambda'))-dilog2(1-(lambda*ones(1,p))./(ones(p,1)*lambda'))))-sum(log((1-c1)*lambda).^2)));
        r=sign(out);
        out=out.^2;
        else
            out=1000;
        end
    case 'Inverse_log'
        p=size(S,1);
        F=Cs*S;
        c2=p/n2;
        lambda_hatC2 = real(sort(eig(F)));
        out=-1/p*sum(log(lambda_hatC2))-(1-c2)/c2*log(1-c2)-1;
        r=sign(out);
        out=out.^2;
    case 'Inverse_log1st'
        p=size(S,1);s=1;
        F=Cs*S;c2=p/n2;
        lambda=sort(real(eig(F)));
        m=@(z) c2*mean(1./(lambda-z))-(1-c2)./z;
        kappa=zeros(p+1,1);lambda1=[0;lambda];
        kappa_p=0;
        kappa_m=-10;
        while abs(kappa_p-kappa_m)>1e-7*abs(lambda(p)-lambda(1))
            kappa_=(kappa_p+kappa_m)/2;
            if m(kappa_)>1/s
                kappa_p=kappa_;
            else
                kappa_m=kappa_;
            end
        end
        kappa(1)=(kappa_p+kappa_m)/2;
        for i=2:length(lambda1)
            kappa_p=lambda1(i);
            kappa_m=lambda1(i-1);
            while abs(kappa_p-kappa_m)>1e-7*abs(lambda(p)-lambda(1))
                kappa_=(kappa_p+kappa_m)/2;
                if m(kappa_)>1/s
                    kappa_p=kappa_;
                else
                    kappa_m=kappa_;
                end
            end
            kappa(i)=(kappa_p+kappa_m)/2;
        end
        out=(1/(c2))*sum(lambda-kappa(2:end))-(1/p)*sum(log(lambda))+(1/p)*sum(log(lambda-kappa(1)))+...
            (1/c2-1)*sum(log(lambda./kappa(2:end)))-1;
        r=sign(out);
        out=out.^2;
    case 'Inverse_t'
        p=size(S,1);
        F=Cs*S;c2=p/n2;
        lambda=sort(real(eig(F)));
        out=0;
        for i=1:length(lambda)
            lambdac=lambda;lambdac(i)=[];
            out=out-(1/p)*((c2/p)*sum(1./(lambda(i)-lambdac))-(1-c2)./lambda(i));
        end
        r=sign(out);
        out=out.^2;
    case 'Inverse_Battacharrya'
        p=size(S,1);s=1;
        F=Cs*S;c2=p/n2;
        lambda=sort(real(eig(F)));
        m=@(z) c2*mean(1./(lambda-z))-(1-c2)./z;
        kappa=zeros(p+1,1);lambda1=[0;lambda];
        kappa_p=0;
        kappa_m=-10;
        while abs(kappa_p-kappa_m)>1e-7*abs(lambda(p)-lambda(1))
            kappa_=(kappa_p+kappa_m)/2;
            if m(kappa_)>1/s
                kappa_p=kappa_;
            else
                kappa_m=kappa_;
            end
        end
        kappa(1)=(kappa_p+kappa_m)/2;
        for i=2:length(lambda1)
            kappa_p=lambda1(i);
            kappa_m=lambda1(i-1);
            while abs(kappa_p-kappa_m)>1e-7*abs(lambda(p)-lambda(1))
                kappa_=(kappa_p+kappa_m)/2;
                if m(kappa_)>1/s
                    kappa_p=kappa_;
                else
                    kappa_m=kappa_;
                end
            end
            kappa(i)=(kappa_p+kappa_m)/2;
        end
        out=(1/2)*((1/(c2))*sum(lambda-kappa(2:end))-(1/p)*sum(log(lambda))+(1/p)*sum(log(lambda-kappa(1)))+...
            (1/c2-1)*sum(log(lambda./kappa(2:end)))-1)-(1/4)*(-1/p*sum(log(lambda))-(1-c2)/c2*log(1-c2)-1)-(1/2)*log(2);
        r=sign(out);
        out=out.^2;
    case 'Inverse_KL'
        p=size(S,1);
        F=Cs*S;
        c2=p/n2;
        lambda_hatC2 = real(sort(eig(F)));
        out_t=0;
        for i=1:length(lambda_hatC2)
            lambdac=lambda_hatC2;lambdac(i)=[];
            out_t=out_t-(1/p)*((c2/p)*sum(1./(lambda_hatC2(i)-lambdac))-(1-c2)./lambda_hatC2(i));
        end
        out=-(1/2)*(-1/p*sum(log(lambda_hatC2))-(1-c2)/c2*log(1-c2)-1)+...
            (1/2)*out_t-(1/2);
        r=sign(out);
        out=out.^2;
end
end

