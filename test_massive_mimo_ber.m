clear;clc;
MM=[400 200 100];
KK=50;
p_db=-20:.5:10;
m=16;    % mQAM modulation order

max_err=100000;max_iter=1000000;
hMod = comm.RectangularQAMModulator('ModulationOrder',m,'BitInput',true,'NormalizationMethod','Average power');
hDemod = comm.RectangularQAMDemodulator('ModulationOrder',m,'BitOutput',true,'NormalizationMethod','Average power');
ber_sim=zeros(length(p_db),length(MM),length(KK));
ber_app=zeros(length(p_db),length(MM),length(KK));
ber_the=zeros(length(p_db),length(MM),length(KK));
counter=0;
for im=1:length(MM)
    M=MM(im);
    for ik=1:length(KK)
        K=KK(ik);
        for ip=1:length(p_db)
            counter=counter+1;
            pu=10^(p_db(ip)/10);
            if mod(p_db(ip),2)==0
                tic;
                fprintf('%2d/%d (M=%3d, K=%3d, Pu=%6.2f dB) ====>  ',counter,length(p_db)*length(MM)*length(KK),M,K,p_db(ip));
                nErr=0;
                nIter=0;
                while (nErr<max_err)&&(nIter<max_iter)
                    nIter=nIter+1;
                    u=randi([0 1],K*log2(m),1);
                    s=step(hMod,u);
                    H=1/sqrt(2)*(randn(M,K)+1j*randn(M,K));
                    n=1/sqrt(2)*(randn(M,1)+1j*randn(M,1));
                    y=sqrt(pu)*H*s+n;
                    A=H;%*inv(H'*H); %MRC
                    r=(A'*y)/sqrt(pu)/M;
                    w=step(hDemod,r);
                    old_number_of_dots=floor(nErr*10/max_err);
                    nErr=nErr+sum(w~=u);
                    new_number_of_dots=min(floor(nErr*10/max_err),10);
                    if (new_number_of_dots>old_number_of_dots)
                        for dots=1:new_number_of_dots-old_number_of_dots
                            fprintf('.');
                        end
                    end
                end
                ber_sim(ip,im,ik)=nErr/(nIter*K*log2(m));
                if (nIter==max_iter)
                    fprintf('ITER_MAX');
                end
                if (nErr==0)
                    fprintf(' -- No errors detected, there is no point to continue.\n');
                    break;
                end
                time=toc;
                fprintf(' BER=%d/%d= %e (%.3f sec)\n',nErr,nIter*K*log2(m),ber_sim(ip,im,ik),time);
            end
            MEAN=M*(exp(1/pu)*double(vpa(expint(K-1,sym(1/pu)),40)));
            VAR=(M*(M+1))*((exp(1/pu)*(K-2+1/pu)*double(vpa(expint(K-2,sym(1/pu)),40))-1)/(K-2))-(MEAN)^2;
            b=VAR/MEAN;
            a=MEAN^2/VAR;
            ber_app(ip,im,ik)=.2/(1.5*b/(m-1)+1)^(a);
            SINR=gamrnd(M,1,100000,1)./(gamrnd(K-1,1,100000,1)+1/pu);
            ber_the(ip,im,ik)=(4*(sqrt(m)-1))/(log2(m)*(sqrt(m)))*mean(qfunc(sqrt((3*SINR)/(m-1))));
        end
        semilogy(p_db(mod(p_db,2)==0),ber_sim((mod(p_db,2)==0),im,ik),'^');
        grid
        hold on
        semilogy(p_db,ber_app(:,im,ik));
        semilogy(p_db,ber_the(:,im,ik));
        legend('sim','app','the')
    end
end
