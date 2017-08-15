% author :  Jiachi.Zou  @TUe
% function: simulating the process of EKF
% date:     08/25/2017
% 
R = 6371393;    % average radius of earth
dT = 0.1;       % time interval
Tx = 180*1e7/(pi*R);   
N = 100;        % #point
n=6;            % #state
q=[100;100;0.1;0.05*dT;0.73*dT;0.14*dT];    %���̱�׼��
r=[2791;7056;0.18;0.33;0.31;0.21];          %������׼��
Q=zeros(n);         %���̷���
R=zeros(n);         %����ֵ�ķ��� 
for i = 1:n
    Q(i,i) = q(i)^2;
    R(i,i) = r(i)^2;
end
f=@(x)[x(1)+x(4)*Tx*dT; x(2)+x(5)*Tx/cos(x(1)*1e-7*pi/180)*dT; x(3)-x(6)*dT; x(4); x(5); x(6)];  %״̬����
h=@(x)[x(1);x(2);x(3);x(4);x(5);x(6)];                   %��������
s=[0;0;0;0;0;0];                                %��ʼ״̬
%��ʼ��״̬
x=s+r.*(2*(rand(n,1)-0.5));                         
P = R;                               
xV = zeros(n,N);          
sV = zeros(n,N);         
zV = zeros(n,N);
u = [10;10;-10];
for k=1:N
  Ty = Tx/cos(s(1)*1e-7*pi/180);
  z = h(s) + r.*(2*(rand(n,1)-0.5));                     
  sV(:,k) = s;                             %ʵ��״̬
  zV(:,k) = z;                           %״̬����ֵ
  [x1,A]=jaccsd(f,x); %����f���ſɱȾ���
  x1(4:6) = x1(4:6) + u*dT;
  P=A*P*A'+Q;         %���̷���Ԥ��
  P1 = P
  [z1,H]=jaccsd(h,x1); %����h���ſɱȾ���
  K=P*H'*pinv(H*P*H'+R) %���������� 
  x=x1+K*(z-z1);        %״̬EKF����ֵ
  P=P-K*H*P;            %EKF����
  P2 = P
  xV(:,k) = x;          %save
  
  % update s
%   ds = s(4:6)*dT+0.5*u*dT^2;
%   s(1) = s(1)+ds(1)*Ty;
%   s(2) = s(2)+ds(2)*Tx;
%   s(3) = s(3)-ds(3);
%   s(4:6) = s(4:6)+u*dT;
  
  % random s
%   ds = 0.5*(2*(rand(3,1)-0.5));
%   s(1) = s(1)+ds(1)*Ty;
%   s(2) = s(2)+ds(2)*Tx;
%   s(3) = s(3)-ds(3);
%   u = 2*(ds-s(4:6)*dT)/dT^2;
%   s(4:6) = s(4:6)+u*dT;
%   u = u + q(4:6).*(2*(rand(3,1)-0.5));
  
  %random u
  u = 10*(2*(rand(3,1)-0.5));
  ds = s(4:6)*dT+0.5*u*dT^2; 
  s(1) = s(1)+ds(1)*Ty;
  s(2) = s(2)+ds(2)*Tx;
  s(3) = s(3)-ds(3);
  s(4:6) = s(4:6)+u*dT;
  u = u + q(4:6).*(2*(rand(3,1)-0.5));
end

for k=1:3
  subplot(3,1,k);  
  plot(1:N, sV(k,:), 'r-', 1:N, xV(k,:), 'b--', 1:N, zV(k,:), 'black+');
  title(sprintf('state %d', k));
end
figure;
for k=1:3
  subplot(3,1,k);  
  plot(1:N, sV(k+3,:), 'r-', 1:N, xV(k+3,:), 'b--', 1:N, zV(k+3,:), 'black+');
  title(sprintf('state %d', k+3));
end