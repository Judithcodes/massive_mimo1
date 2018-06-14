clear;clc;
channel='awgn';
M=16;
hMod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',true,'NormalizationMethod','Average power');
hDemod = comm.RectangularQAMDemodulator('ModulationOrder',M,'BitOutput',true,'NormalizationMethod','Average power');
eb_n0_db=-10:10;
eb_n0=10.^(eb_n0_db/10);
n0=1;
eb=eb_n0*n0;
es=log2(M)*eb;
N=10000;

ber_sim=zeros(size(eb_n0_db));
max_err=10000;max_iter=1000;
for i=1:length(eb_n0_db)
    tic;
    release(hMod);release(hDemod)
    set(hMod,'AveragePower',es(i));
    set(hDemod,'AveragePower',es(i));
    fprintf('%2d/%d (%6.2f dB) ====> ',i,length(eb_n0_db),eb_n0_db(i));
    nErr=0;
    nIter=0;
    while (nErr<max_err)&&(nIter<max_iter)
        nIter=nIter+1;
        u=randi([0 1],N*log2(M),1);
        x=step(hMod,u);
        n=sqrt(n0)/sqrt(2)*(randn(N,1)+1j*randn(N,1));
        if strcmp(channel,'rayleigh')
            h=1/sqrt(2)*(randn(N,1)+1j*randn(N,1));           
        else
            h=ones(N,1);           
        end
        y=(h.*x+n).*conj(h)./abs(h).^2;
        w=step(hDemod,y);
        old_number_of_dots=floor(nErr*10/max_err);
        nErr=nErr+sum(w~=u);
        new_number_of_dots=min(floor(nErr*10/max_err),10);
        if (new_number_of_dots>old_number_of_dots)
            for dots=1:new_number_of_dots-old_number_of_dots
                fprintf('.');
            end      
        end
    end
    ber_sim(i)=nErr/(nIter*N*log2(M));
    if (nIter==max_iter)
        fprintf('ITER_MAX');
    end
    if (nErr==0)
        fprintf(' -- No errors detected, there is no point to continue.\n');
        ber_sim(i+1:end)=0;
        break;
    end  
    time=toc;
    fprintf(' BER=%d/%d= %e (%.3f sec)\n',nErr,nIter*N*log2(M),ber_sim(i),time);    
end
semilogy(eb_n0_db,ber_sim,'-^','linewidth',2);
if strcmp(channel,'awgn')
    hold on
    ber_the=4*(sqrt(M)-1)/sqrt(M)/log2(M)*qfunc(sqrt(3*es/(M-1)));
    ber_app=.2*exp(-1.5*es/(M-1));
    semilogy(eb_n0_db,ber_app,'--g','linewidth',2);   
    semilogy(eb_n0_db,ber_the,'--k','linewidth',2);
end
grid