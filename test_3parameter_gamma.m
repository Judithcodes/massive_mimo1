a=3;
b=2;
z=1;
N=10000;
V2=gamrnd(a,b,1,N);
V3=V2+z;
[fs2,x2]=hist(V2,100);
fs2=fs2/N/(x2(2)-x2(1));
ft2=1/b^a/gamma(a)*x2.^(a-1).*exp(-x2/b);
[fs3,x3]=hist(V3,100);
fs3=fs3/N/(x3(2)-x3(1));
ft3=(1/b)^a/gamma(a)*(x3-z).^(a-1).*exp(-(x3-z)/b);

bar(x2,fs2);
hold on
plot(x2,ft2,'r');

bar(x3,fs3);
plot(x3,ft3,'g');
