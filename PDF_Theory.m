%基于能量依赖的频率求概率密度函数PDF-半解析解
clc;clear all;
load('WH_delta1=3_delta3=3_k=0.3_alpha=0.05.mat')
w=Win;
beta=0.02;
delta1=3;
delta3=3;
kk=0.3;
alpha=0.05;
D=0.005;
c=0.3;
tau1=1.3;
tau2=0.5;
u=0.01;
v=0.01;
hx=81;
hy=81;
hd=0.049;
x=-2+(0:hx-1)*hd;
y=-2+(0:hy-1)*hd;
n=length(x);
H=zeros(n,n);
P=zeros(n,n);
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
figure;
mesh(x,y,pp')
xlabel('$X$','Interpreter','latex');
ylabel('$\dot X$','Interpreter','latex');
zlabel('$p(X,\dot X)$','Interpreter','latex');
xlim([-2,2]);
ylim([-2,2]);
