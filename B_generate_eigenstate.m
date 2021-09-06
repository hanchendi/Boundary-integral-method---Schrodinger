clear
clc

load('eigv_statical.mat');
mode=5;
k=eigv_statical(mode);
N=300;

r=linspace(0,1,N);
delta_theta=2*pi/N;
theta=delta_theta:delta_theta:2*pi;
r=r';
b=0.1;
c(1)=0.1;
c(2)=0.2;
delta(1)=pi/2;
delta(2)=pi/2;
X1=r*cos(theta);
Y1=r*sin(theta);
Z=X1+sqrt(-1)*Y1;
W=(Z+b*Z.^2+c(1)*exp(sqrt(-1)*delta(1))*Z.^3+c(2)*exp(sqrt(-1)*delta(2))*Z.^4)/sqrt(1+2*b^2+3*c(1)^2+4*c(2)^2);
X=real(W);
Y=imag(W);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%µÚN²ã
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0=X(N,:);
y0=Y(N,:);
for i=1:N-1
    ds(i)=sqrt((x0(i+1)-x0(i))^2+(y0(i+1)-y0(i))^2);
end
ds(N)=sqrt((x0(1)-x0(N))^2+(y0(1)-y0(N))^2);

for i=1:N-1
    t1(i)=(y0(i+1)-y0(i))/(x0(i+1)-x0(i));
    t2(i)=x0(i+1)-x0(i);
end
t1(N)=(y0(1)-y0(N))/(x0(1)-x0(N));
t2(N)=x0(1)-x0(N);
for i=1:N
    if t1(i)>=0 && t2(i)>=0
        derivation(i)=atan(t1(i));
    elseif t1(i)<0 && t2(i)<0
        derivation(i)=atan(t1(i))+pi;
    elseif t1(i)>=0 && t2(i)<0
        derivation(i)=atan(t1(i))+pi;
    elseif t1(i)<0 && t2(i)>=0
        derivation(i)=atan(t1(i))+2*pi;
    end
end  

for i=1:N
    for j=1:N
        if i~=j
            s=(y0(j)-y0(i))/(x0(j)-x0(i));
            if s>=0 && (x0(j)-x0(i))>=0
                angle(i,j)=atan(s);
            elseif s<0 && (x0(j)-x0(i))<0
                angle(i,j)=atan(s)+pi;
            elseif s>=0 && (x0(j)-x0(i))<0
                angle(i,j)=atan(s)+pi;
            elseif s<0 && (x0(j)-x0(i))>=0
                angle(i,j)=atan(s)+2*pi;
            end
            rho(i,j)=sqrt((y0(j)-y0(i))^2+(x0(j)-x0(i))^2);
        end
    end
end

for i=1:N
    alpha(i,:)=angle(i,:)-derivation(i);
    B(i,:)=ds;
end

A=zeros(N,N);
A=sqrt(-1)/2.*k*sin(alpha).*besselh(1,1,k.*rho).*B;
for i=1:N
    A(i,i)=0;
end
[s1 s2]=eig(A);
s3=diag(abs(s2+eye(N)));
t1=find(s3==min(s3));
t1=t1(1);
u0=s1(:,t1)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%µÚN²ã
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi=zeros(N,N);
for t=1:N-1
    x=X(t,:);
    y=Y(t,:);
    
    for i=1:N
        for j=1:N
            rho(i,j)=sqrt((y0(j)-y(i))^2+(x0(j)-x(i))^2);
        end
    end
    A=-sqrt(-1)/2.*k.*besselh(0,1,k.*rho);
    phi(t,:)=-A*(u0'.*ds');
    disp(t)
end

X(:,N+1)=X(:,1);
Y(:,N+1)=Y(:,1);
phi(:,N+1)=phi(:,1);
%minus=5;
figure()
mesh(X,Y,abs(phi).^2);
axis([-0.9 1.5 -1.2 1.2])
axis off
phi_B=abs(phi).^2;
save([pwd,'/Boundary_',num2str(mode),'.mat'],'phi_B');