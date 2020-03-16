%信噪比SNR，以SND-D-v为例
clc;
clear all;
load('WH_delta1=3_delta3=3_k=0.3_alpha=0.05.mat')
w=Win;
beta=0.02;
delta1=3;
delta3=3;
kk=0.3;
alpha=0.05;
c=0.3;
u=0;
tau1=0;
tau2=0.5;
F=0.01;
omega=0.02;
hx=81;
hy=81;
hd=0.049;
x=-2+(0:hx-1)*hd;
y=-2+(0:hy-1)*hd;
n=length(x);
w1=w(61,42); %稳态点处的频率
ND=0:0.001:0.03;
nd=length(ND);
Nv=-0.03:0.002:0.03;
ndv=length(Nv);
SNR=zeros(nd,ndv);
for b=1:nd
    D=ND(b);
    for b1=1:ndv
        v=Nv(b1);
        beta1=beta+kk*alpha/(alpha^2+w1^2)+u*sin(w1*tau1)/w1-v*cos(w1*tau2);
        delta=kk*w1^2/(alpha^2+w1^2)-u*cos(w1*tau1)-v*w1*sin(w1*tau2);
        xs2=(delta1-kk*w1^2/(alpha^2+w1^2))/delta3;
        R0=sqrt(2)/4/pi*(-beta-kk*alpha/(alpha^2+w1^2)+sqrt((beta+kk*alpha/(alpha^2+w1^2))^2+4*(delta1-kk*w1^2/(alpha^2+w1^2))))*...
            exp(beta1*(1+c^2*w1^2)/D*(-1/2*delta1*xs2+1/4*delta3*xs2^2+1/2*delta*xs2));
        R1=R0*beta1*(1+c^2*w1^2)/D*sqrt(xs2);
        SNR(b,b1)=pi*R1^2*F^2/(4*R0)/(1-R1^2*F^2/(2*R1^2+2*omega^2));
    end
end

figure;
mesh(ND,Nv,SNR');
xlim([0,0.03]);
ylim([-0.03,0.03]);
xlabel('$D$','Interpreter','latex');
ylabel('$\nu$','Interpreter','latex');
zlabel('$SNR$','Interpreter','latex');
hold on
plot3(ND(6)*ones(1,31),Nv,SNR(6,:)','r','linewidth',2)