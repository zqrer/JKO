%addpath /home/qrzhang/software/Ipopt/Ipopt/contrib/MatlabInterface/
%tic
%clc;clear;close all;
function [Error Error1] = analytic(N, dt, total, t0, g)
%N = 10;
x_max = 1; %[-x_max, x_max]
y_max = 1; %[-x_max, x_max]
dx = x_max/N;
dy = y_max/N;
%dt = 0.1;
x_mid  = zeros(N+1,1);
y_mid  = zeros(N+1,1);
for i = 1:N+1
   x_mid(i) = (i-1)*dx;
   y_mid(i) = (i-1)*dy;
end
x = x_mid(1:N)+dx/2;
y = y_mid(1:N)+dy/2;
dx_mid = dx;
dy_mid = dy;
%g = 1.0;
D1 = 1;
z1 = 1;
dtbeta  = 0.0;
c1 = zeros(N*N,1);
c1_exa = zeros(N*N,1);
c1_err = zeros(N*N,1);
m1x = zeros((N-1)*N,1);
m1y = zeros((N-1)*N,1);
V = zeros(N,N);
for i = 1:N
   for l = 1:N
      V(i, l) = -1.0*g*x(i);
   end
end
%t0 = 0.1;
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
alpha = pi^2 + g^2/4;
for i = 1:N
   for l = 1:N
   c1((i-1)*N+l) = exp(-alpha*t0+g/2*x(i))*(pi*cos(pi*x(i))+g/2*sin(pi*x(i)))+pi * exp(g*(x(i)-0.5));
end
end

for i = 1:N-1
   for l = 1:N
      m1x((i-1)*N+l) = -(g^2/4+pi^2)*exp(-alpha*t0+g/2*x_mid(i+1))*sin(pi*x_mid(i+1));
   end
end

m1y = zeros((N-1)*N,1);
m1x = -1 * m1x * D1;
m1y = -1 * m1y * D1;

result_c1 = zeros(total+2,N*N);
result_c1(1,:) = c1;

k = 0;
flag = 0;
i = 1;
Error = 0;
Error1 = 0;
while i <= total

   %[u1, info1] = opt_func(c1, m1, N, dx, dt, D1, z1, dtbeta, V, flag);
   [u1, info1] = opt_func1(c1, m1x, m1y, N, dx, dx_mid, dt, D1, z1, dtbeta, V, flag, Gx, Gy, G1x, G1y);
   if info1.status < 0 
      if info1.status == -2 && flag == 0 
         fprintf('\nRelaxation!! Result %d\n', i);
         flag = 1;
      elseif info1.status == -2 && flag == 1
         fprintf('\nRelaxation 2 !! Result %d\n', i);
         [c1 m1x m1y status] = relaxation1(c1, m1x, m1y, N, dx, dx_mid, dt, D1, z1, dtbeta, V, i, Gx, Gy, G1x, G1y);

%         dt1 = dt/10;
%         flag1 = 0;
%         j1 = 1;
%         while j1 <= 10
%            [u1, info1] = opt_func(c1, m1, N, dx, dt1, D1, z1, beta, V, flag1);
%            if info1.status < 0
%               if info1.status == -2 && flag1 == 0
%                  fprintf('\nRelaxation!! Relaxation!! Result %d\n', j1);
%                  flag1 = 1;
%               else
%                  break;
%               end
%            else
%                  fprintf('\nRelaxation!! Small step!! Result %d\n', j1);
%                  [c1 m1] = deal(u1{:});
%                  flag1 = 0;
%                  j1 = j1 + 1;
%            end
%         end% while
%         if j1 > 10
          if status == 1
            fprintf('Result %d\n', i);
            fprintf('c1 %d:\n',sum(c1));
            disp(c1);
            %fprintf('m1:\n');
            %disp(m1 /dt);
            %if(j==10^k || j == 50 || j == 100  || j == 150 || j == 200 || j == 250)
            if(i~=0)
               T = dt*i+t0;
               for j = 1:N
                  for l = 1:N
                  c1_exa((j-1)*N+l) = exp(-alpha*T+g/2*x(j))*(pi*cos(pi*x(j))+g/2*sin(pi*x(j)))+pi*exp(g*(x(j)-0.5));
                  end
               end
               c1_err = abs(c1-c1_exa);
               fprintf('c1_err:\n');
               disp(c1_err);
               c1_err = abs((c1-c1_exa)*dx_mid*dx_mid);
               if Error < sum(c1_err)
                  Error = sum(c1_err);
               end
               Error1 = Error1 + sum(c1_err);
               k = k+1;
               result_c1(k+1,:) = c1;
               fprintf('\n Count:%d',k);
            end
            flag = 0;
            i = i + 1;
         else
            fprintf('Failed Step = %d', i);
            break;
         end
      else
         break;
      end% if
   else
       fprintf('Result %d\n', i);
       [c1 m1x m1y] = deal(u1{:});
       fprintf('c1 %d:\n',sum(c1));
       disp(c1);
       %if(j==10^k || j == 50 || j == 100  || j == 150 || j == 200 || j == 250)
       if(i~=0)
         T = dt*i + t0;
         for j = 1:N
            for l = 1:N
           c1_exa((j-1)*N+l) = exp(-alpha*T+g/2*x(j))*(pi*cos(pi*x(j))+g/2*sin(pi*x(j)))+pi*exp(g*(x(j)-0.5));
         end
      end
          c1_err = abs(c1-c1_exa);
          fprintf('c1_err:\n');
          disp(c1_err);
          c1_err = abs((c1-c1_exa)*dx_mid*dx_mid);
          if Error < sum(c1_err)
            Error = sum(c1_err);
          end
           Error1 = Error1 + sum(c1_err);
           k = k+1;
           result_c1(k+1,:) = c1;
           fprintf('\n Count:%d',k);
       end
       flag = 0;
       i = i + 1;
   end
   
      disp(info1);

end
    Error1 = Error1 *dt;
    fprintf('Error = ');
    disp(Error);
    fprintf('Error1 = ');
    disp(Error1);
    save('analytic1.mat','result_c1','x','y','x_mid','y_mid');
