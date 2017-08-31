% author :  Jiachi.Zou  @TUe
% function: simulating the process of EKF
% date:     08/25/2017
% 

%load data [h a_down]
data = load('ekf2_2.dat');

% temp_data = load('ekf2_test5.dat');
% data = [temp_data(:,2),temp_data(:,5)];

% qh = 0.35;
% qv = 0.15;
% rh = 4;
% rv = 0.9;
%qh = 0.35;
qh = 0.3;
qv = 0.1;
rh = 0.175;
rv = 0.12;

dT = 0.05;
N = length(data);   %data point number
n = 2;              %states number

F = [1 -dT;0 1];
B = [0;dT];
H = eye(2,2);

Q=zeros(n);
R=zeros(n);
q = [qh;qv];
%Q = B*B'*qv*qv
% Q = q*q';
Q(1,1) = qh*qh*dT*dT;
Q(2,2) = qv*qv*dT*dT;
R(1,1) = rh*rh;
R(2,2) = rv*rv;
%init P=R
P = R;
%P(2,2) = 1;
%init states
x=[data(1,1);0];  

u = data(:,2);
z = x;

xV = zeros(n,N);                   
zV = zeros(n,N);
vel = 0;

% init FIR filter
fir_order = 10;
b = myfir(100, 5, fir_order);
% do FIR
out = zeros(length(u), 1);
fir_len = fir_order+1;
fir_index = 0;
circle_buffer = zeros(1, fir_len);
for i = 1 : length(u)
    circle_buffer(fir_index+1) = u(i);
    fir_index = fir_index + 1;
    if fir_index >= fir_len
        fir_index = 0;
    end
    for j = 1 : fir_len
        out(i) = out(i) + circle_buffer(fir_index+1)*b(j);
        if fir_index ~= 0
            fir_index = fir_index - 1;
        else
            fir_index = fir_len - 1;
        end
    end
end

for k=1:N
    z(1) = data(k,1);
    pre_alt = x(1);
    
    x = F*x+B*out(k);
    P = F*P*F'+Q;  

    y = z - H*x;
    S = H*P*H' + R;
    K = P*H'/S; 
    x = x + K*y;
    P = P - K*H*P;
    
%     vel = vel*0.8 + (x(1)-pre_alt)/dT*0.2;
    vel = vel*0.8 + (pre_alt-x(1))/dT*0.2;
    z(2) = vel;
    xV(:,k) = x;
    zV(:,k) = z;
end
figure
plot(1:length(data), xV(1,:), 'r*-', 1:length(data), zV(1,:), 'b*-');
legend('卡尔曼高度',  '观测高度')
figure
plot(1:length(data), xV(2,:), 'r*-', 1:length(data), zV(2,:), 'b*-');
legend('卡尔曼速度',  '观测速度')
figure
plot(1:length(data), u, 'r*-', 1:length(data), out, 'b*-');
legend('加速度', 'FIR加速度')