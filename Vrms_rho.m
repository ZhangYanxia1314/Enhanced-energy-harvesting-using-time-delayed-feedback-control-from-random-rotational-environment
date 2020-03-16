%数值计算有效输出电压Vrms和能量转换率rho%
clc;
clear all;

beta=0.02;
delta1=3;
delta3=3;
kk=0.3;
alpha=0.05;
c=0.3;
u=-0.01;
vv=0.01;
tau1=0.6;
tau2=2.5;
f=0.01;
omega=0.05;
omega1=omega/2/pi;
n=100000;
dt=0.05;
t=0:dt:n*dt;
n2=50;%序列个数
x=zeros(1,n+1);
y=zeros(1,n+1);
v=zeros(1,n+1);
z=zeros(1,n+1);
%分布边界确定
n1=51; % 网格稀疏程度
a1=-2;
a2=2;
dh=(a2-a1)/(n1-1);
a_x=linspace(a1,a2,n1);
a_y=linspace(a1,a2,n1);
ax_mid=linspace(a1+dh,a2-dh,n1-1);
ay_mid=linspace(a1+dh,a2-dh,n1-1);
%对每个点进行循环分类，循环步数为n1
num=zeros(n1-1,n1-1);
tn1=round(tau1/dt);
tn2=round(tau2/dt);
tn=max(tn1,tn2);
ND=[0.001:0.001:0.03];
Nn=length(ND);
Evv=zeros(n2,Nn);
PM=zeros(n2,Nn);
PE=zeros(n2,Nn);
eta=zeros(n2,Nn);
for qi=1:Nn
    tic
    D=ND(qi);
    for m=1:n2
        N1=randn(1,n);
        for j=tn+1:n
            A1=y(j)*dt;
            B1=dt*(-1)*(beta*y(j)-delta1*x(j)+delta3*x(j)^3+kk*v(j))+dt*(z(j)+f*sin(omega*t(j))+u*x(j-tn1)+vv*y(j-tn2));
            C1=dt*(-alpha*v(j)+y(j));
            D1=dt*(-z(j)/c)+1/c*(sqrt(2*D*dt)*N1(1,j));
            
            A2=(y(j)+0.5*B1)*dt;
            B2=dt*(-1)*(beta*(y(j)+0.5*B1)-delta1*(x(j)+0.5*A1)+delta3*(x(j)+0.5*A1)^3+kk*(v(j)+0.5*C1))+dt*((z(j)+0.5*D1)+f*sin(omega*(t(j)+0.5*dt))+u*(x(j-tn1)+0.5*A1)+vv*(y(j-tn2)+0.5*B1));
            C2=dt*(-alpha*(v(j)+0.5*C1)+(y(j)+0.5*B1));
            D2=dt*(-(z(j)+0.5*D1)/c)+1/c*(sqrt(2*D*dt)*N1(1,j));
            
            A3=(y(j)+0.5*B2)*dt;
            B3=dt*(-1)*(beta*(y(j)+0.5*B2)-delta1*(x(j)+0.5*A2)+delta3*(x(j)+0.5*A2)^3+kk*(v(j)+0.5*C2))+dt*((z(j)+0.5*D2)+f*sin(omega*(t(j)+0.5*dt))+u*(x(j-tn1)+0.5*A2)+vv*(y(j-tn2)+0.5*B2));
            C3=dt*(-alpha*(v(j)+0.5*C2)+(y(j)+0.5*B2));
            D3=dt*(-(z(j)+0.5*D2)/c)+1/c*(sqrt(2*D*dt)*N1(1,j));
            
            A4=(y(j)+B3)*dt;
            B4=dt*(-1)*(beta*(y(j)+B3)-delta1*(x(j)+A3)+delta3*(x(j)+A3)^3+kk*(v(j)+C3))+dt*((z(j)+D3)+f*sin(omega*(t(j)+dt))+u*(x(j-tn1)+A3)+vv*(y(j-tn2)+B3));
            C4=dt*(-alpha*(v(j)+C3)+(y(j)+B3));
            D4=dt*(-(z(j)+D3)/c)+1/c*(sqrt(2*D*dt)*N1(1,j));
            
            x(j+1)=x(j)+(A1+2*A2+2*A3+A4)/6;
            y(j+1)=y(j)+(B1+2*B2+2*B3+B4)/6;
            v(j+1)=v(j)+(C1+2*C2+2*C3+C4)/6;
            z(j+1)=z(j)+(D1+2*D2+2*D3+D4)/6;
            pm(j)=D/(1+c^2*omega1^2)*abs(y(j+1))+abs(y(j+1)*f*sin(omega*t(j+1)));
            pe(j)=alpha*kk*v(j+1).^2;
        end
        v1=v(2*n/5:end).^2;
        Evv(m,qi)=sqrt(mean(v1));
        PM(m,qi)=sqrt(mean(pm(2*n/5:end).^2));
        PE(m,qi)=sqrt(mean(pe(2*n/5:end).^2));
        eta(m,qi)=PE(m,qi)./PM(m,qi);
    end
    qi
    toc
end

Evvm=sum(Evv)/n2; % the output RMS voltage
Eta=sum(eta)/n2*100;% the power conversion efficiency

%%%%%%%%%%%%%%%%%%%%%%%%%% 拟合曲线
Q1=polyfit(ND,Evvm,11);
Evn=zeros(1,Nn);
for i=1:Nn
    Evn(i)=polyval(Q1,ND(i));
end

figure;
plot(ND,Evvm,ND,Evn,'r-o');
xlabel('$D$','Interpreter','latex','fontsize',14);
ylabel('$Vrms$','Interpreter','latex','fontsize',14);


%%%%%%%%%%%%%%%%%%%%%%%%%% 拟合曲线
Q2=polyfit(ND,Eta,11);
Etn=zeros(1,Nn);
for i=1:Nn
    Etn(i)=polyval(Q2,ND(i));
end

figure;
plot(ND,Eta,ND,Etn,'b-o');
xlabel('$D$','Interpreter','latex','fontsize',14);
ylabel('\rho%','fontsize',14);
