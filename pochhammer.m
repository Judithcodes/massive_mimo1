function p=pochhammer(n,a)
    if (floor(a)~=a)||(a>0)
        p=gamma(a+n)/gamma(a);
    else
        m=-a;
        if n<=m
            p=(-1)^n*factorial(m)/factorial(m-n);
        else
            p=0;
        end;
    end;