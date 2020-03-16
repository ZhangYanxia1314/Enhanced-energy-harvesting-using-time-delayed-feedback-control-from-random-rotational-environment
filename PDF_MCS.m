%采用MCS计算概率密度函数PDF的数值解
clc;
clear all;
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
vv=0.01;
f=0;
omega=0.05;
n=200000;
dt=0.05;
t=0:dt:n*dt;
n2=150;%序列个数
x=zeros(1,n+1);
y=zeros(1,n+1);
v=zeros(1,n+1);
z=zeros(1,n+1);
%分布边界确定
n1=51; %网格稀疏程度
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
for m=1:n2
    N1=randn(1,n);
    tic
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
    end
    for p=tn+1:n
        j = sum(x(p)>a_x);
        k = sum(y(p)>a_y);
        if x(p)>min(a_x)&&x(p)<max(a_x)&&y(p)>min(a_y)&&y(p)<max(a_y)
            num(j,k) = num(j,k)+1;
        end
    end
    toc
    m
end

num=num/n2;
num_density=num/(dh^2*(n-tn));
figure;
mesh(ax_mid,ay_mid,num_density');
xlabel('$X$','Interpreter','latex');
ylabel('$\dot X$','Interpreter','latex');
zlabel('$p(X,\dot X)$','Interpreter','latex');
xlim([-2,2]);
ylim([-2,2]);
