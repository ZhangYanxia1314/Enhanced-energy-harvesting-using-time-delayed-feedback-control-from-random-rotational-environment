function Th1=tt(a,b)
global delta1 delta3 k alpha H w
syms x1
f=@(x1)1./sqrt(abs(2*H-(-delta1*x1.^2+delta3*x1.^4/2+k*w^2*x1.^2/(alpha^2+w^2))));
Th=0;tn=500;td=(b-a)/tn;
for ii=1:tn-1
    x1=a+td*ii;
    Th=Th+f(x1)*td;
end
Th1=roundn(Th,-4);
end
