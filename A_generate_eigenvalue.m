clear
clc

% boundary integral method, find zero point start from kmin to kmax
% N boundary discretizations
kmin=15;
kmax=16;
dk=0.0001;
N=500;
delta_theta=2*pi/N;
theta=delta_theta:delta_theta:2*pi;
b=0.1;
c(1)=0.1;
c(2)=0.2;
delta(1)=pi/2;
delta(2)=pi/2;
x1=cos(theta);
y1=sin(theta);
z=x1+sqrt(-1)*y1;
w=(z+b*z.^2+c(1)*exp(sqrt(-1)*delta(1))*z.^3+c(2)*exp(sqrt(-1)*delta(2))*z.^4)/sqrt(1+2*b^2+3*c(1)^2+4*c(2)^2);
x=real(w);
y=imag(w);

for i=1:N-1
    ds(i)=sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
end
ds(N)=sqrt((x(1)-x(N))^2+(y(1)-y(N))^2);

for i=1:N-1
    t1(i)=(y(i+1)-y(i))/(x(i+1)-x(i));
    t2(i)=x(i+1)-x(i);
end
t1(N)=(y(1)-y(N))/(x(1)-x(N));
t2(N)=x(1)-x(N);
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
            s=(y(j)-y(i))/(x(j)-x(i));
            S(i,j)=s;
            if s>=0 && (x(j)-x(i))>=0
                angle(i,j)=atan(s);
            elseif s<0 && (x(j)-x(i))<0
                angle(i,j)=atan(s)+pi;
            elseif s>=0 && (x(j)-x(i))<0
                angle(i,j)=atan(s)+pi;
            elseif s<0 && (x(j)-x(i))>=0
                angle(i,j)=atan(s)+2*pi;
            end
            rho(i,j)=sqrt((y(j)-y(i))^2+(x(j)-x(i))^2);
        end
    end
end

for i=1:N
    alpha(i,:)=angle(i,:)-derivation(i);
    B(i,:)=ds;
end

t=1;
for k=kmin:dk:kmax
    A=zeros(N,N);
    A=sqrt(-1)/2.*k*sin(alpha).*besselh(1,1,k.*rho).*B;
    for i=1:N
        A(i,i)=0;
    end
    A=A+eye(N);
%     [U,S,V] = svd(A);
%     svd_1(t)=S(N,N);
%     svd_2(t)=S(N-1,N-1);
    result(t,1)=k;
    result(t,2)=abs(det(A));
    disp(k)
    t=t+1;
end

load('eigv_statical.mat');
figure()
plot(result(:,1),result(:,2),'b');hold on;
plot(eigv_statical(1:10),0,'*r');
set(gca,'FontSize',15)
axis([0 8 0 1])