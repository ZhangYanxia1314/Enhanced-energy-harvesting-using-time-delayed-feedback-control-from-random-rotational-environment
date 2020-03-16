function XA=xa(H,w)
global delta1 delta3 k alpha
XA1=sqrt((delta1-(k*w^2)/(alpha^2+w^2)-sqrt((delta1-(k*w^2)/(alpha^2+w^2))^2+4*delta3*H))/delta3);
XA=roundn(XA1,-4);
end
