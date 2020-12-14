addpath /home/qrzhang/software/Ipopt/Ipopt/contrib/MatlabInterface/
tic
clc;clear;close all;
N = 21;
x_max = 10.0; %[-x_max, x_max]
dx = 2 * x_max/(N-1);
dt = 1e-2;
x  = zeros(N,1);
for i = 1:N
   x(i) = -x_max + (i-1)*dx;
end
x_mid = x(1:N-1)+dx/2;
Dx = x(2:N)-x(1:N-1);
Dx_mid  = zeros(N,1);
Dx_mid(1) = dx/2;
Dx_mid(N) = dx/2;
Dx_mid(2:N-1) = x_mid(2:N-1)-x_mid(1:N-2);
ec = 1.602189e-19;
Kb = 1.380650524e-23;
T = 298;
NA = 6.0220450e23;
epsilon_0 = 8.854187817e-12;
epsilon_s = 80; %solvent = water
lambda_D = sqrt((epsilon_0 * epsilon_s * Kb * T)/(ec * ec * 10 * NA))* 1e9; % Debye length for one ion species
mu = ec/(Kb *T );
phi = 1.0;
Voltage = @(t) phi .*mu .* (t+x_max)./(2*x_max);% unit = 1 
D1 = 0.133; % unit = 1
%D2 = 0.203;
z1 = 1;
%z2 = -1;
%dtbeta = zeros(N-1,1);
dtbeta = 0.0;
c1 = zeros(N,1);
%c2 = zeros(N,1);
m1 = zeros(N-1,1);
%m2 = zeros(N-1,1);
V = zeros(N-1,1);

%for i = 1:N-1
%   V(i) = Voltage(x(i)+dx/2);
%end
V = Voltage(x_mid);

c1_l = 1.0; % unit = 1
c1_r = 0.1; % unit = 1
K = (c1_r - c1_l)/(2*x_max);
%for i = 1:N
%   c1(i) = K*(x(i)+x_max)+c1_l;
%end
c1 = K.*(x+x_max)+c1_l;

%for i = 1:N-1
%   m1(i) = -1*(K + ((phi * mu)/(2*x_max))*(K*(x(i)+dx/2+x_max)+c1_l));
%end

  m1 = K + ((phi .* mu)./(2.*x_max)).*(K.*(x_mid+x_max)+c1_l);

m1 = -1.* m1 .* D1;


result_c1 = zeros(30,N);
result_m1 = zeros(30,N-1);
result_c1(1,:) = c1;
result_m1(1,:) = m1;

k = 0;
flag = 0;
j = 1;
while j <= 1e3

   [u1, info1] = opt_func(c1, m1, N, dx, dt, D1, z1, dtbeta, V, flag);
   %[u1, info1] = opt_func1(c1, m1, N, Dx, Dx_mid, dt, D1, z1, dtbeta, V, flag);

   if info1.status < 0 
      if info1.status == -2 && flag == 0 
         fprintf('\nRelaxation!! Result %d\n', j);
         flag = 1;
      elseif info1.status == -2 && flag == 1
         fprintf('\nRelaxation 2 !! Result %d\n', j);
         %[c1 m1 status] = relaxation(c1, m1, N, dx, dt, D1, z1, dtbeta, V, j);
         [c1 m1 status] = relaxation1(c1, m1, N, Dx, Dx_mid, dt, D1, z1, dtbeta, V, j);

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
            disp(m1 /lambda_D/ dt);
            if(j==10^k || j == 5e4 || j == 1e5  || j == 3e5 || j == 5e5 || j == 7e5 || j == 1e6)
               k = k+1;
               result_c1(k+1,:) = c1;
               result_m1(k+1,:) = m1 /lambda_D /dt;
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
       disp(m1 /lambda_D /dt);
       if(j==10^k || j == 5e4 || j == 1e5  || j == 3e5 || j == 5e5 || j == 7e5 || j == 1e6)
           k = k+1;
           result_c1(k+1,:) = c1;
           result_m1(k+1,:) = m1 / lambda_D /dt;
           fprintf('\n Count:%d',k);
       end
       flag = 0;
       j = j + 1;
   end
   
      disp(info1);

end
   save('result1.mat','result_c1','x','result_m1','x_mid');
