m=16;
g=0:.01:10;
p1=4*(sqrt(m)-1)/sqrt(m)/log2(m)*qfunc(sqrt(3*g/(m-1)));
p2=4/log2(m)*qfunc(sqrt(3*g*log2(m)/(m-1)));
plot(g,p1)
hold on
plot(g,p2)
