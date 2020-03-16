function HH=Hh(w,x1,x2)
global delta1 delta3 k alpha
HH=x2^2/2-delta1*x1^2/2+delta3*x1^4/4+0.5*k*w^2*x1^2/(alpha^2+w^2);
end
