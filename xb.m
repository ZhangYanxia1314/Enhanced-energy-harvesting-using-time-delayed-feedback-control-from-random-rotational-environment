function XB=xb(H,w)
global delta1 delta3 k alpha
XB1=sqrt((delta1-(k*w^2)/(alpha^2+w^2)+sqrt((delta1-(k*w^2)/(alpha^2+w^2))^2+4*delta3*H))/delta3);
XB=roundn(XB1,-4);
end
