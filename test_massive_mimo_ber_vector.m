clear;clc;
MM=200;
KK=20;
p_db=-10:10;
m=4;    % mQAM modulation order
N=500;
max_err=1000;max_iter=100000/N;
hMod = comm.RectangularQAMModulator('ModulationOrder',m,'BitInput',true,'NormalizationMethod','Average power');
hDemod = comm.RectangularQAMDemodulator('ModulationOrder',m,'BitOutput',true,'NormalizationMethod','Average power');
ber_sim=zeros(length(p_db),length(MM),length(KK));

counter=0;
for im=1:length(MM)
    M=MM(im);
    for ik=1:length(KK)
        K=KK(ik);
        l=logical(kron(eye(N),ones(K,1)));
        for ip=1:length(p_db)
            counter=counter+1;
            pu=10^(p_db(ip));
            tic;
            fprintf('%2d/%d (M=%3d, K=%3d, Pu=%6.2f dB) ====>  ',counter,length(p_db)*length(MM)*length(KK),M,K,p_db(ip));
            nErr=0;
            nIter=0;
            while (nErr<max_err)&&(nIter<max_iter)
                nIter=nIter+1;
                u=randi([0 1],K*log2(m)*N,1);
                s=step(hMod,u);
                h=1/sqrt(2)*(randn(1,M*K*N)+1j*randn(1,M*K*N));
                H=sparse(repmat(1:M*N,1,K),kron(repmat(1:K:N*K,1,K)+kron(0:K-1,ones(1,N)),ones(1,M)),h);
                n=1/sqrt(2)*(randn(M*N,1)+1j*randn(M*N,1));
                y=sqrt(pu)*H*s+n;
                Y=sparse(1:M*N,kron(1:N,ones(1,M)),y);
                A=H;%*inv(H'*H); %MRC
                r=(A'*Y)/sqrt(pu)/M;
                r=full(r(l));
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
            ber_sim(ip,im,ik)=nErr/(nIter*K*log2(m)*N);
            if (nIter==max_iter)
                fprintf('ITER_MAX');
            end
            if (nErr==0)
                fprintf(' -- No errors detected, there is no point to continue.\n');
                break;
            end
            time=toc;
            fprintf(' BER=%d/%d= %e (%.3f sec)\n',nErr,nIter*K*log2(m)*N,ber_sim(ip,im,ik),time);
        end
    end
end
semilogy(p_db,ber_sim(:,1,1));
grid