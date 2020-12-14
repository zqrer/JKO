addpath /home/qrzhang/software/Ipopt/Ipopt/contrib/MatlabInterface/
tic
clc;clear;close all;
format long
N = 4;
x_max = 1; %[0, x_max]
dx = x_max/N;
dt = 0.001;
total = 3000;
%dtbeta = 0.0;
dtbeta = zeros(N-1,1);
%dtbeta = sqrt(2)*dt*ones(N-1,1);
x_mid = zeros(N+1,1);
for i = 1:N+1
   x_mid(i) = (i-1)*dx;
end
x = x_mid(1:N)+dx/2;
Dx = x(2:N)-x(1:N-1);
Dx_mid  = x_mid(2:N+1)-x_mid(1:N);
g = 10.0;
Voltage = @(t) -1.0.* g.* t;% unit = 1 
D1 = 1;
z1 = 1;
c1 = zeros(N,1);
c1_exa = zeros(N,1);
c1_err = zeros(N,1);
m1 = zeros(N-1,1);
M1 = zeros(N,1);
t0 = 0.0;

%for i = 1:N-1
%   V(i) = Voltage(x(i)+dx/2);
%end
V = Voltage(x);
G = zeros(N-1,1);
G1 = zeros(N-1, 1);
for i = 1:N-1
   G(i) = (1-exp(z1*(V(i)-V(i+1))*0.5))/(1-exp(z1*(V(i)-V(i+1))));
   G1(i) = -exp(z1*(V(i)-V(i+1))*0.5)/(1-exp(z1*(V(i)-V(i+1))))*z1*(V(i)-V(i+1))/Dx(i);
end

alpha = pi^2 + g^2/4;

%for i = 1:N
   %c1(i) = exp(-alpha*t0+g/2*x(i))*(pi*cos(pi*x(i))+g/2*sin(pi*x(i)))+pi * exp(g*(x(i)-0.5));
   %M1(i) = exp(-1.0 * Voltage(x(i)));
%end

c1 = exp(-alpha*t0+g/2.*x).*(pi.*cos(pi.*x)+g/2.*sin(pi.*x))+pi.*exp(g.*(x-0.5));
M1 = exp(-1.0 .* V);
M = sum(c1.*Dx_mid)/sum(M1.*Dx_mid);
E_inf1 = sum(M.*M1.*log(M).*Dx_mid);

%   m1 = g/2.*exp(-alpha*t0+g/2.*x_mid).*(pi.*cos(pi.*x_mid)+g/2.*sin(pi.*x_mid))+exp(-alpha*t0+g/2.*x_mid).*(-1*pi*pi.*sin(pi.*x_mid)+g*pi/2.*cos(pi.*x_mid))+g*pi.*exp(g.*(x_mid-0.5))-g.*(exp(-alpha*t0+g/2.*x_mid).*(pi.*cos(pi.*x_mid)+g/2.*sin(pi.*x_mid))+pi.*exp(g.*(x_mid-0.5)));

m1 = -(g^2/4+pi^2).* exp(-alpha*t0+g/2.*x_mid(2:N)).*sin(pi.*x_mid(2:N));

m1 = -1 * m1 * D1;

result_c1 = zeros(total+1,N);
result_c1_exa = zeros(total+1,N);
result_m1 = zeros(total+1,N-1);
result_c1(1,:) = c1;
result_c1_exa(1,:) = c1;
result_m1(1,:) = m1;
result_E = zeros(total,1);
%result_E_exa = zeros(total,1);


k = 0;
flag = 0;
j = 1;
Error = 0;
Error1 = 0;
%while j <= 1e6
while j <= total

%   [u1, info1] = opt_func(c1, m1, N, dx, dt, D1, z1, dtbeta, V, flag);
   [u1, info1] = opt_func1(c1, m1, N, Dx, Dx_mid, dt, D1, z1, dtbeta, V, flag, G, G1);
%    [u1, info1] = opt_func2(c1, m1, N, Dx, Dx_mid, dt, D1, z1, dtbeta, V, flag, g);
   if info1.status < 0 
      if info1.status == -2 && flag == 0 
         fprintf('\nRelaxation!! Result %d\n', j);
         flag = 1;
      elseif info1.status == -2 && flag == 1
         fprintf('\nRelaxation 2 !! Result %d\n', j);
   %      [c1 m1 status] = relaxation(c1, m1, N, dx, dt, D1, z1, dtbeta, V, j);
         [c1 m1 status] = relaxation1(c1, m1, N, Dx, Dx_mid, dt, D1, z1, dtbeta, V, j, G, G1);
%         dt1 = dt/10;
%         flag1 = 0;
%         j1 = 1;
%         while j1 <= 10
%            [u1, info1] = opt_func(c1, m1, N, dx, dt1, D1, z1, dtbeta, V, flag1);
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
            fprintf('c1:\n');
            disp(c1);
           % fprintf('m1:\n');
           % disp(m1 /dt);
            if(j~=0)
               T = dt*j+t0;
               %for i = 1:N
               %   c1_exa(i) = exp(-alpha*T+g/2*x(i))*(pi*cos(pi*x(i))+g/2*sin(pi*x(i)))+pi*exp(g*(x(i)-0.5));
               %end
               %c1_exa = exp(-alpha*T+g/2.*x).*(pi.*cos(pi.*x)+g/2.*sin(pi.*x))+pi.*exp(g.*(x-0.5));
               c1_exa = exp(-alpha*T+g/2.*x).*(pi.*cos(pi.*x)+g/2.*sin(pi.*x))+pi.*exp(g.*(x-0.5));
               c1_err = abs(c1-c1_exa);
               fprintf('c1_exa  c1\n');
               disp(sum(c1_exa));
               disp(sum(c1));
               fprintf('c1_err:\n');
               disp(c1_err);
               c1_err = abs((c1-c1_exa).*Dx_mid);
               if Error < sum(c1_err)
                  Error = sum(c1_err);
               end
               Error1 = Error1 + sum(c1_err);
               k = k+1;
               result_c1_exa(k+1,:) = c1_exa;
               result_c1(k+1,:) = c1;
               result_m1(k+1,:) = m1;
               for i = 1:N
                  result_E(k) = result_E(k) + D1* (c1(i)*log(c1(i))+z1*V(i)*c1(i))*Dx_mid(i); 
               end
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
       fprintf('c1: \n');
       disp(c1);
      % fprintf('m1:\n');
      % disp(m1/dt);
       %if(j==10^k || j == 50 || j == 100  || j == 150 || j == 200 || j == 250)
       if(j~=0)
         T = dt*j + t0;
         %for i = 1:N
         %  c1_exa(i) = exp(-alpha*T+g/2*x(i))*(pi*cos(pi*x(i))+g/2*sin(pi*x(i)))+pi*exp(g*(x(i)-0.5));
         %end
          c1_exa = exp(-alpha*T+g/2.*x).*(pi.*cos(pi.*x)+g/2.*sin(pi.*x))+pi.*exp(g.*(x-0.5));
          c1_err = abs(c1-c1_exa);
          fprintf('c1_exa  c1\n');
          disp(sum(c1_exa));
          disp(sum(c1));
          fprintf('c1_err:\n');
          disp(c1_err);
          c1_err = abs((c1-c1_exa).*Dx_mid);
          if Error < sum(c1_err)
            Error = sum(c1_err);
          end
            Error1 = Error1 + sum(c1_err);
           k = k+1;
           result_c1_exa(k+1,:) = c1_exa;
           result_c1(k+1,:) = c1;
           result_m1(k+1,:) = m1;
          for i = 1:N-1
            result_E(k) = result_E(k) + D1* (c1(i)*log(c1(i))+z1*V(i)*c1(i))*Dx_mid(i); 
          end
           fprintf('\n Count:%d',k);
       end
       flag = 0;
       j = j + 1;
   end
   
      disp(info1);

end
    fprintf('Error = ');
    disp(Error);
    fprintf('Error1 = ');
    disp(Error1*dt);
    fprintf('energy1 :');
    disp(result_E-E_inf1);
    
    time = 0:dt:3.0;
    save('dissipate10_4_10.mat','result_c1','x','result_m1','x_mid','result_E','time','E_inf1','result_c1_exa');
