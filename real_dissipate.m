addpath /home/qrzhang/software/Ipopt/Ipopt/contrib/MatlabInterface/
tic
clc;clear;close all;
format long
N = 50;
x_max = 1; %[0, x_max]
dx = x_max/N;
dt = 0.0001;
total = 5000;
%dtbeta = 0.0;
dtbeta = zeros(N-1,1);
%dtbeta = sqrt(0.5)*dt*ones(N-1,1);
x_mid = zeros(N+1,1);
for i = 1:N+1
   x_mid(i) = (i-1)*dx;
end
x = x_mid(1:N)+dx/2;
Dx = x(2:N)-x(1:N-1);
Dx_mid  = x_mid(2:N+1)-x_mid(1:N);
g = 10.0;


D1 = 1;
z1 = 1;
c1 = zeros(N,1);
m1 = zeros(N-1,1);
M1 = zeros(N,1);
t0 = 0.0;

%V = -1.0.*g.*sin(pi.*x);
%V_mid = -1.0.*g.*sin(pi.*x_mid(2:N));
V = -1.0.*g.*cos(g*pi.*x);
V_mid = -1.0.*g.*cos(g*pi.*x_mid(2:N));
G = zeros(N-1,1);
G1 = zeros(N-1, 1);
for i = 1:N-1
   if V(i)~=V(i+1)
      G(i) = (1-exp(z1*(V(i)-V(i+1))*0.5))/(1-exp(z1*(V(i)-V(i+1))));
      G1(i) = -exp(z1*(V(i)-V(i+1))*0.5)/(1-exp(z1*(V(i)-V(i+1))))*z1*(V(i)-V(i+1))/Dx(i);
   else
      G(i) = 0.5;
      G1(i) = 1/Dx(i);
   end
end
%c1 = exp(g.*sin(pi.*x))+cos(2.*pi.*x)+g.*sin(pi.*x);
c1 = 2*pi.*x-sin(2*pi.*x);
M1 = exp(-1.0 .* V);
M = sum(c1.*Dx_mid)/sum(M1.*Dx_mid);
E_inf1 = sum(M.*M1.*log(M).*Dx_mid);
%m1 = -2.*pi.*sin(2.*pi.*x_mid(2:N))+pi.*g.*cos(pi.*x_mid(2:N))-pi.*g.*cos(pi.*x_mid(2:N)).*cos(2.*pi.*x_mid(2:N))-pi.*g.*g.*sin(pi.*x_mid(2:N)).*cos(pi.*x_mid(2:N)); 
m1 = 2*pi-2*pi.*cos(2*pi.*x_mid(2:N))+g^2*pi.*sin(g*pi.*x_mid(2:N)).*(2*pi.*x_mid(2:N)-sin(2*pi.*x_mid(2:N)));
m1 = -1 * m1 * D1;

result_c1 = zeros(total+1,N);
result_m1 = zeros(total+1,N-1);
result_c1(1,:) = c1;
result_m1(1,:) = m1;
result_E = zeros(total,1);


k = 0;
flag = 0;
j = 1;
while j <= total
   T = dt*j+t0;
%   [u1, info1] = opt_func(c1, m1, N, dx, dt, D1, z1, dtbeta, V, flag);
   [u1, info1] = opt_func1(c1, m1, N, Dx, Dx_mid, dt, D1, z1, dtbeta, V, flag, G, G1);
%    [u1, info1] = opt_func2(c1, m1, N, Dx, Dx_mid, dt, D1, z1, dtbeta, V, flag, g);
   if info1.status < 0 
      if info1.status == -2 && flag == 0 
         fprintf('\nRelaxation!! Result %d\n', j);
         flag = 1;
      elseif info1.status == -2 && flag == 1
         fprintf('\nRelaxation 2 !! Result %d\n', j);
         %[c1 m1 status] = real_relaxation(c1, m1, N, Dx, Dx_mid, dt, D1, z1, dtbeta, j, t0, g, x_mid, G, G1);
         [c1 m1 status] = relaxation1(c1, m1, N, Dx, Dx_mid, dt, D1, z1, dtbeta, V, j, G, G1);
          if status == 1
            fprintf('Result %d\n', j);
            fprintf('sum c1:\n');
            disp(sum(c1));
            fprintf('c1:\n');
            disp(c1);
            fprintf('m1: \n');
            disp(m1);
            if(j~=0)
               k = k+1;
               result_c1(k+1,:) = c1;
               result_m1(k+1,:) = m1;
               for i = 1:N-1
                 result_E(k) = result_E(k) + D1* (c1(i)*log(c1(i))+z1*V(i)*c1(i))*Dx_mid(i); 
%                     c_t = c1(i)*(1-G(i))+c1(i+1)*G(i);
 %                  result_E(k) = result_E(k) + D1*(c_t*log(c_t)+z1*V_mid(i)*c_t)*Dx(i);
               end
               fprintf('\n Count:%d',k);
            end
            flag = 0;
            j = j + 1;
         else
            fprintf('Failed Step = %d', j);
            break;
         end%if status == 1
      else
         break;
      end% if info1.status == -2 && flag == 0 
   else
       fprintf('Result %d\n', j);
       [c1 m1] = deal(u1{:});
       fprintf('sum c1: \n');
       disp(sum(c1));
       fprintf('c1: \n');
       disp(c1);
       fprintf('m1: \n');
       disp(m1);
       if(j~=0)
         k = k+1;
         result_c1(k+1,:) = c1;
         result_m1(k+1,:) = m1;
         for i = 1:N-1
            result_E(k) = result_E(k) + D1* (c1(i)*log(c1(i))+z1*V(i)*c1(i))*Dx_mid(i); 
%            c_t = c1(i)*(1-G(i))+c1(i+1)*G(i);
%            result_E(k) = result_E(k) + D1*(c_t*log(c_t)+z1*V_mid(i)*c_t)*Dx(i);
         end
           fprintf('\n Count:%d',k);
       end
       flag = 0;
       j = j + 1;
   end %if info1.status < 0
   
      disp(info1);

end
    fprintf('energy :');
    disp(result_E-E_inf1);
    
    time = 0:dt:0.5;
    save('real_dissipate10_50_new.mat','result_c1','x','result_m1','x_mid','result_E','time','E_inf1');

