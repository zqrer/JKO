addpath /home/qrzhang/software/Ipopt/Ipopt/contrib/MatlabInterface/
tic
clc;clear;close all;
order_n = 6;
T = 0.55;
t0 = 0.05;
total_base = round((T-t0)/0.25);
g = 1.0;
Error =  zeros(order_n,1);
Error1 = zeros(order_n,1);
N = 5;
total = total_base;
dt = (T-t0)/total;
for i = 1:order_n
   fprintf('\nStep = %d\n', i);
    [Error(i) Error1(i)] = analytic(N, dt, total, t0, g);
    N = N*2;
    total = total*2;
    dt = dt/2;
 %   fprintf('Error %d', i);
 %   disp(Error(i));
 %   fprintf('Error1 %d', i);
 %   disp(Error1(i))
 end
 
    fprintf('Error');
    disp(Error);
    fprintf('Error1 ' );
    disp(Error1);
 Order0 = log2(Error(1:end-1)./Error(2:end));
 Order1 = log2(Error1(1:end-1)./Error1(2:end));
 fprintf('order =');
 disp(Order0);
 fprintf('order1 =');
 disp(Order1);



