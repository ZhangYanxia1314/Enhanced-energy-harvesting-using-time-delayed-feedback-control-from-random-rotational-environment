%基于能量依赖的频率计算平均输出功率E[P]
clc;
clear all;
load('H_delta1=3_delta3=3_k=0.3_alpha=0.05.mat')
load('WH_delta1=3_delta3=3_k=0.3_alpha=0.05.mat')
w=Win;
beta=0.02;
delta1=3;
delta3=3;
kk=0.3;
alpha=0.05;
D=0.005;
c=0.3;
%u=0.2;
%v=0.03;
tau1=0.5;
tau2=0.5;
hx=81;
hy=81;
hd=0.049;
x=-2+(0:hx-1)*hd;
y=-2+(0:hy-1)*hd;
n=length(x);
Nu=[-0.02:0.002:0.02];
nd=length(Nu);
Nv=[-0.01:0.002:0.02];
ndv=length(Nv);
H=zeros(n,n);
P=zeros(n,n);
Ep=zeros(nd,ndv);
for n1=1:nd
    u=Nu(n1);
    for n2=1:ndv
        v=Nv(n2);
        for i=1:n
            for j=1:n
                w1=w(i,j);
                H(i,j)=1/2*y(j)^2-1/2*delta1*x(i)^2+1/4*delta3*x(i)^4+1/2*(kk*w1^2/(alpha^2+w1^2)-u*cos(w1*tau1)-v*w1*sin(w1*tau2))*x(i)^2;
                beta1=beta+kk*alpha/(alpha^2+w1^2)+u*sin(w1*tau1)/w1-v*cos(w1*tau2);
                P(i,j)=(1+c^2*w1^2)/D*exp(-beta1*(1+c^2*w1^2)/D*H(i,j));
            end
        end
        ss=sum(sum(P));
        pp=P/(ss*hd^2);
        Ev1=0;
        for i=1:n
            for j=1:n
                w1=w(i,j);
                Hv=NHi(i,j);
                if Hv>=0
                    Hv1=0;
                else
                    Hv1=sign(x(i));
                end
                Ev1=Ev1+(w1^2/(alpha^2+w1^2)*(x(i)-Hv1*sqrt(delta1/delta3))+alpha/(alpha^2+w1^2)*y(j))^2*pp(i,j)*hd*hd;
            end
        end
        Ep(n1,n2)=kk*alpha*Ev1;
    end
end
figure;
mesh(Nu,Nv,Ep')
xlabel('$\mu$','Interpreter','latex','fontsize',14);
ylabel('$\nu$','Interpreter','latex','fontsize',14);
zlabel('$E[P]$','Interpreter','latex','fontsize',14);
