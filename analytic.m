%addpath /home/qrzhang/software/Ipopt/Ipopt/contrib/MatlabInterface/
%tic
%clc;clear;close all;
function [Error Error1] = analytic(N, dt, total, t0, g)
%N = 10;
x_max = 1; %[-x_max, x_max]
dx = x_max/N;
%dt = 0.1;
x_mid  = zeros(N+1,1);
for i = 1:N+1
   x_mid(i) = (i-1)*dx;
end
x = x_mid(1:N)+dx/2;
Dx_mid = x_mid(2:N+1)-x_mid(1:N); %N vector
Dx = x(2:N) - x(1:N-1);
%g = 1.0;
Voltage = @(t) -1.0 * g * t;% unit = 1 
D1 = 1;
z1 = 1;
%dtbeta = (dt*dx)^2;
%dtbeta  = 0.0;
dtbeta = zeros(N-1,1);
%dtbeta = 8*dt*dx*ones(N-1,1);
c1 = zeros(N,1);
c1_exa = zeros(N,1);
c1_err = zeros(N,1);
m1 = zeros(N-1,1);
V = zeros(N,1);
%t0 = 0.1;

V = Voltage(x);
G = zeros(N-1,1);
G1 = zeros(N-1, 1);
for i = 1:N-1
   G(i) = (1-exp(z1*(V(i)-V(i+1))*0.5))/(1-exp(z1*(V(i)-V(i+1))));
   G1(i) = -exp(z1*(V(i)-V(i+1))*0.5)/(1-exp(z1*(V(i)-V(i+1))))*z1*(V(i)-V(i+1))/Dx(i);
end
alpha = pi^2 + g^2/4;
%for i = 1:N
%   c1(i) = exp(-alpha*t0+g/2*x(i))*(pi*cos(pi*x(i))+g/2*sin(pi*x(i)))+pi * exp(g*(x(i)-0.5));
%end
c1 = exp(-alpha*t0+g/2.*x).*(pi.*cos(pi.*x)+g/2.*sin(pi.*x))+pi.*exp(g.*(x-0.5));

%for i = 1:N-1
%   m1(i) = g/2*exp(-alpha*t0+g/2*x_mid(i))*(pi*cos(pi*x_mid(i))+g/2*sin(pi*x_mid(i)))+exp(-alpha*t0+g/2*x_mid(i))*(-1*pi*pi*sin(pi*x_mid(i))+g*pi/2*cos(pi*x_mid(i)))+g*pi*exp(g*(x_mid(i)-0.5))-g*(exp(-alpha*t0+g/2*x_mid(i))*(pi*cos(pi*x_mid(i))+g/2*sin(pi*x_mid(i)))+pi*exp(g*(x_mid(i)-0.5)));
%end

   %m1 = -g/2.*exp(-alpha*t0+g/2.*x_mid(2:N)).*(pi.*cos(pi.*x_mid(2:N))+g/2.*sin(pi.*x_mid(2:N)))+exp(-alpha*t0+g/2.*x_mid(2:N)).*(-1*pi*pi.*sin(pi.*x_mid(2:N))+g*pi/2.*cos(pi.*x_mid(2:N)));
m1 = -(g^2/4+pi^2).* exp(-alpha*t0+g/2.*x_mid(2:N)).*sin(pi.*x_mid(2:N));
m1 = -1 * m1 * D1;
result_c1 = zeros(total+2,N);
result_m1 = zeros(total+2,N-1);
result_c1(1,:) = c1;
result_m1(1,:) = m1;

k = 0;
flag = 0;
j = 1;
Error = 0;
Error1 = 0;
%while j <= 1e6
while j <= total

   %[u1, info1] = opt_func(c1, m1, N, dx, dt, D1, z1, dtbeta, V, flag);
   [u1, info1] = opt_func1(c1, m1, N, Dx, Dx_mid, dt, D1, z1, dtbeta, V, flag, G, G1);
   if info1.status < 0 
      if info1.status == -2 && flag == 0 
         fprintf('\nRelaxation!! Result %d\n', j);
         flag = 1;
      elseif info1.status == -2 && flag == 1
         fprintf('\nRelaxation 2 !! Result %d\n', j);
         %[c1 m1 status] = relaxation(c1, m1, N, dx, dt, D1, z1, dtbeta, V, j);
         [c1 m1 status] = relaxation1(c1, m1, N, Dx, Dx_mid, dt, D1, z1, dtbeta, V, j, G, G1);

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
            fprintf('Result %d\n', j);
            fprintf('c1 %d:\n',sum(c1));
            disp(c1);
            fprintf('m1:\n');
            disp(m1 /dt);
            %if(j==10^k || j == 50 || j == 100  || j == 150 || j == 200 || j == 250)
            if(j~=0)
               T = dt*j+t0;
               for i = 1:N
                  c1_exa(i) = exp(-alpha*T+g/2*x(i))*(pi*cos(pi*x(i))+g/2*sin(pi*x(i)))+pi*exp(g*(x(i)-0.5));
               end
               c1_err = abs(c1-c1_exa);
               fprintf('c1_err:\n');
               disp(c1_err);
               c1_err = abs((c1-c1_exa).*Dx_mid);
               if Error < sum(c1_err)
                  Error = sum(c1_err);
               end
               Error1 = Error1 + sum(c1_err);
               k = k+1;
               result_c1(k+1,:) = c1;
               result_m1(k+1,:) = m1;
               fprintf('\n Count:%d',k);
            end
            flag = 0;
            j = j + 1;
         else
            fprintf('Failed Step = %d', j);
            break;
         end
      else
         break;
      end% if
   else
       fprintf('Result %d\n', j);
       [c1 m1] = deal(u1{:});
       fprintf('c1 %d:\n',sum(c1));
       disp(c1);
       fprintf('m1:\n');
       disp(m1/dt);
       %if(j==10^k || j == 50 || j == 100  || j == 150 || j == 200 || j == 250)
       if(j~=0)
         T = dt*j + t0;
         for i = 1:N
           c1_exa(i) = exp(-alpha*T+g/2*x(i))*(pi*cos(pi*x(i))+g/2*sin(pi*x(i)))+pi*exp(g*(x(i)-0.5));
         end
          c1_err = abs(c1-c1_exa);
          fprintf('c1_err:\n');
          disp(c1_err);
          c1_err = abs((c1-c1_exa).*Dx_mid);
          if Error < sum(c1_err)
            Error = sum(c1_err);
          end
           Error1 = Error1 + sum(c1_err);
           k = k+1;
           result_c1(k+1,:) = c1;
           result_m1(k+1,:) = m1;
           fprintf('\n Count:%d',k);
       end
       flag = 0;
       j = j + 1;
   end
   
      disp(info1);

end
    Error1 = Error1 *dt;
    fprintf('Error = ');
    disp(Error);
    fprintf('Error1 = ');
    disp(Error1);
    save('analytic1.mat','result_c1','x','result_m1','x_mid');
