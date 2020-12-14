addpath /home/qrzhang/software/Ipopt/Ipopt/contrib/MatlabInterface/
tic
clc;clear;close all;
format long
N = 40;
x_max = 1; %[0, x_max]
y_max = 1; %[0, x_max]
dx = x_max/N;
dy = y_max/N;
dt = 0.0001;
total = 1000;
dtbeta = dt;
x_mid = zeros(N+1,1);
y_mid = zeros(N+1,1);
for i = 1:N+1
   x_mid(i) = (i-1)*dx;
   y_mid(i) = (i-1)*dy;
end
x = x_mid(1:N)+dx/2;
y = y_mid(1:N)+dy/2;
dx_mid = dx;
dy_mid = dy;
g = 10.0;
D1 = 1;
z1 = 1;
c1 = zeros(N*N,1);
m1x = zeros(N*(N-1),1);
m1y = zeros(N*(N-1),1);
t0 = 0.0;
V = zeros(N,N);
for i = 1:N
   for l = 1:N
      V(i, l) = -1.0*g*cos(pi*x(i))-g*cos(pi*y(l));
   end
end

Gx = zeros(N-1, N);
Gy = zeros(N, N-1);
G1x = zeros(N-1, N);
G1y = zeros(N, N-1);
for i = 1:N-1
   for l = 1:N
   if(V(i,l)==V(i+1,l))
         Gx(i,l) = 0.5;
         G1x(i,l) = 1/dx;
   else
      Gx(i,l) = (1-exp(z1*(V(i,l)-V(i+1,l))*0.5))/(1-exp(z1*(V(i,l)-V(i+1,l))));
      G1x(i,l) = -exp(z1*(V(i,l)-V(i+1,l))*0.5)/(1-exp(z1*(V(i,l)-V(i+1,l))))*z1*(V(i,l)-V(i+1,l))/dx;
   end
end
end
for i = 1:N
   for l = 1:N-1
      if(V(i,l)==V(i,l+1))
         Gy(i,l) = 0.5;
         G1y(i,l) = 1/dy;
      else
      Gy(i,l) = (1-exp(z1*(V(i,l)-V(i,l+1))*0.5))/(1-exp(z1*(V(i,l)-V(i,l+1))));
      G1y(i,l) = -exp(z1*(V(i,l)-V(i,l+1))*0.5)/(1-exp(z1*(V(i,l)-V(i,l+1))))*z1*(V(i,l)-V(i,l+1))/dy;
   end
end
end

for i = 1:N
   for l =1:N
      c1((i-1)*N+l) = 2*pi*x(i)-sin(2*pi*x(i))+2*pi*y(l)-sin(2*pi*y(l));
   end
end
 for i = 1:N-1
    for l = 1:N
       m1x((i-1)*N+l) = 2*pi-2*pi*cos(2*pi*x(i))+g*pi*sin(pi*x(i))*(2*pi*x(i)-sin(2*pi*x(i))+2*pi*y(l)-sin(2*pi*y(l)));
    end
 end
 for i = 1:N
    for l = 1:N-1
       m1y((i-1)*(N-1)+l) = 2*pi-2*pi*cos(2*pi*y(l))+g*pi*sin(pi*y(l))*(2*pi*x(i)-sin(2*pi*x(i))+2*pi*y(l)-sin(2*pi*y(l)));
    end
 end
m1x = -1 * m1x * D1;
m1y = -1 * m1y * D1;

M1 = exp(-1.0*z1 .* V);
M = sum(c1*dx*dx)/sum(sum(M1*dx*dx));
E_inf1 = sum(sum(M*log(M)*dx*dx*M1));

result_c1 = zeros(total+1,N*N);
result_c1(1,:) = c1;
result_E = zeros(total,1);


k = 0;
flag = 0;
i = 1;
while i <= total
   [u1, info1] = opt_func1(c1, m1x, m1y, N, dx, dx_mid, dt, D1, z1, dtbeta, V, flag, Gx, Gy, G1x, G1y);
   if info1.status < 0 
      if info1.status == -2 && flag == 0 
         fprintf('\nRelaxation!! Result %d\n', i);
         flag = 1;
      elseif info1.status == -2 && flag == 1
         fprintf('\nRelaxation 2 !! Result %d\n', i);
         [c1 m1x m1y status] = relaxation1(c1, m1x, m1y, N, dx, dx_mid, dt, D1, z1, dtbeta, V, i, Gx, Gy, G1x, G1y);
          if status == 1
            fprintf('Result %d\n', i);
            fprintf('sum c1:\n');
            disp(sum(c1));
            fprintf('c1:\n');
            disp(c1);
            %fprintf('m1: \n');
            %disp(m1);
            if(i~=0)
               k = k+1;
               result_c1(k+1,:) = c1;
               for j = 1:N
                  for l = 1:N
                  result_E(k) = result_E(k) + D1* (c1((j-1)*N+l)*log(c1((j-1)*N+l))+z1*V(j,l)*c1((j-1)*N+l))*dx_mid*dx_mid; 
               end
            end
               fprintf('\n Count:%d',k);
            end
            flag = 0;
            i = i + 1;
         else
            fprintf('Failed Step = %d', j);
            break;
         end%if status == 1
      else
         break;
      end% if info1.status == -2 && flag == 0 
   else
       fprintf('Result %d\n', i);
       [c1 m1x m1y] = deal(u1{:});
       fprintf('sum c1: \n');
       disp(sum(c1));
       fprintf('c1: \n');
       disp(c1);
       if(i~=0)
         k = k+1;
         result_c1(k+1,:) = c1;
          for j = 1:N
            for l = 1:N
               result_E(k) = result_E(k) + D1* (c1((j-1)*N+l)*log(c1((j-1)*N+l))+z1*V(j,l)*c1((j-1)*N+l))*dx_mid*dx_mid; 
            end
          end
           fprintf('\n Count:%d',k);
       end
       flag = 0;
       i = i + 1;
   end %if info1.status < 0
   
      disp(info1);

end
    fprintf('energy :');
    disp(result_E-E_inf1);
    
    time = 0:dt:0.1;
    save('real_dissipate10_40.mat','result_c1','x','y','x_mid','y_mid','result_E','time','E_inf1');

