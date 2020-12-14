function [u, info] = opt_func1(cc, mx, my, N, dx, dx_mid, dt, D, z, dtbeta, V, flag, Gx, Gy, G1x, G1y)

   b = cc;
   if flag ~= 0
      temp = diag(-1 * ones(1, N-1),-1) + diag(ones(1,N),0);
      My_temp = temp(:,1:N-1);
      My = My_temp;
      for j = 1:N-1
             My = blkdiag(My, My_temp);
      end
      My = 1/dx_mid .* My;
      Mx_temp = eye(N*N) + diag(-1.*ones(1,N*(N-1)),-N);
      Mx = Mx_temp(:,1:(N-1)*N);
      Mx = 1/dx_mid .* Mx;
      Mx = sparse(Mx);
      My = sparse(My);
      c0 = b - Mx * mx - My * my;
   else
      c0 = cc;
   end
   u0 = {c0
         mx
         my};

   options.cl = [b; zeros(N*N,1)];
   options.cu = [b; Inf * ones(N*N,1)];

   options.auxdata = {N dx dx_mid dt D z dtbeta V Gx Gy G1x G1y};



    options.ipopt.jac_d_constant   = 'yes';
    options.ipopt.jac_c_constant   = 'yes';
    options.ipopt.hessian_approximation = 'limited-memory';
    %options.ipopt.derivative_test       = 'first-order';
    options.ipopt.mu_strategy      = 'adaptive';
    options.ipopt.max_iter         = 100000;
    options.ipopt.tol              = 1e-8;


    funcs.objective         = @objective;
    funcs.constraints       = @constraints;
    funcs.gradient          = @gradient;
    funcs.jacobian          = @jacobian;
    funcs.jacobianstructure = @jacobianstructure;

    [u info] = ipopt_auxdata(u0,funcs,options);



% ------------------------------------------------------------------
function f = objective (u, auxdata)
   [N dx dx_mid dt D z dtbeta V Gx Gy G1x G1y] = deal(auxdata{:});
   [c0 mx my] = deal(u{:});
   c = zeros(N*N,1);
   for j = 1:N*N
      if(c0(j)<1e-6)
         c(j) = 1e-6;
      else
         c(j) = c0(j);
      end
   end
      f = 0.0;
      for j = 1:N-1
         for l = 1:N
            f = f + mx((j-1)*N+l)^2/(c((j-1)*N+l)+Gx(j,l)*(c(j*N+l)-c((j-1)*N+l)))*dx*dx + dtbeta^2*((c(j*N+l)-c((j-1)*N+l))*G1x(j,l)*(log(c(j*N+l))-log(c((j-1)*N+l))))*dx; 
         end
      end
      for j = 1:N
         for l = 1:N-1
            f = f + my((j-1)*(N-1)+l)^2/(c((j-1)*N+l)+Gy(j,l)*(c((j-1)*N+l+1)-c((j-1)*N+l)))*dx*dx + dtbeta^2*((c((j-1)*N+l+1)-c((j-1)*N+l))*G1y(j,l)*(log(c((j-1)*N+l+1))-log(c((j-1)*N+l))))*dx;
         end
      end
      for j =1:N
         for l = 1:N
             f = f + 2*dt*D*(c((j-1)*N+l)*log(c((j-1)*N+l))+z*V(j,l)*c((j-1)*N+l))*dx*dx;
          end
       end
   

function C = constraints(u, auxdata)
   [N dx dx_mid dt D z dtbeta V Gx Gy G1x G1y] = deal(auxdata{:});
   [c0 mx my] = deal(u{:});
   temp = diag(-1 * ones(1, N-1),-1) + diag(ones(1,N),0);
   My_temp = temp(:,1:N-1);
   My = My_temp;
   for j = 1:N-1
      My = blkdiag(My, My_temp);
   end
    My = 1/dx_mid .* My;
    Mx_temp = eye(N*N) + diag(-1.*ones(1,N*(N-1)),-N);
    Mx = Mx_temp(:,1:(N-1)*N);
    Mx = 1/dx_mid .* Mx;
    Mx = sparse(Mx);
    My = sparse(My);
   C = [c0 + Mx * mx + My *my; c0];


function g = gradient(u, auxdata)
   [N dx dx_mid dt D z dtbeta V Gx Gy G1x G1y] = deal(auxdata{:});
   [c0 mx my] = deal(u{:});
   c = zeros(N,1);
   for j = 1:N*N
      if(c0(j)<1e-6)
         c(j) = 1e-6;
      else
         c(j) = c0(j);
      end
   end
   g = zeros(N*N+2*N*(N-1),1);
   g(1) = mx(1)^2*(Gx(1,1)-1)*dx*dx/((c(1)+Gx(1,1)*(c(N+1)-c(1)))^2)+my(1)^2*(Gy(1,1)-1)*dx*dx/((c(1)+Gy(1,1)*(c(2)-c(1)))^2)+dtbeta^2*(-G1x(1,1)*(log(c(N+1))-log(c(1)))*dx-G1x(1,1)*(c(N+1)-c(1))/c(1)*dx-G1y(1,1)*(log(c(2))-log(c(1)))*dx-G1y(1,1)*(c(2)-c(1))/c(1)*dx)+2*dt*D*(log(c(1))+1+z*V(1,1))*dx*dx; 

for l = 2:N-1
   g(l) = mx(l)^2*(Gx(1,l)-1)*dx*dx/((c(l)+Gx(1,l)*(c(N+l)-c(l)))^2)+my(l)^2*(Gy(1,l)-1)*dx*dx/((c(l)+Gy(1,l)*(c(l+1)-c(l)))^2)-my(l-1)^2*Gy(1,l-1)*dx*dx/((c(l-1)+Gy(1,l-1)*(c(l)-c(l-1)))^2)+dtbeta^2*(-G1x(1,l)*(log(c(N+l))-log(c(l)))*dx-G1x(1,l)*(c(N+l)-c(l))/c(l)*dx-G1y(1,l)*(log(c(l+1))-log(c(l)))*dx-G1y(1,l)*(c(l+1)-c(l))/c(l)*dx+G1y(1,l-1)*(log(c(l))-log(c(l-1)))*dx+G1y(1,l-1)*(c(l)-c(l-1))/c(l)*dx)+2*dt*D*(log(c(l))+1+z*V(1,l))*dx*dx;
end

   g(N) = mx(N)^2*(Gx(1,N)-1)*dx*dx/((c(N)+Gx(1,N)*(c(N+N)-c(N)))^2)-my(N-1)^2*Gy(1,N-1)*dx*dx/((c(N-1)+Gy(1,N-1)*(c(N)-c(N-1)))^2)+dtbeta^2*(-G1x(1,N)*(log(c(N+N))-log(c(N)))*dx-G1x(1,N)*(c(N+N)-c(N))/c(N)*dx+G1y(1,N-1)*(log(c(N))-log(c(N-1)))*dx+G1y(1,N-1)*(c(N)-c(N-1))/c(N)*dx)+2*dt*D*(log(c(N))+1+z*V(1,N))*dx*dx;

for j = 2:N-1
   g((j-1)*N+1) = mx((j-1)*N+1)^2*(Gx(j,1)-1)*dx*dx/((c((j-1)*N+1)+Gx(j,1)*(c(j*N+1)-c((j-1)*N+1)))^2)+my((j-1)*(N-1)+1)^2*(Gy(j,1)-1)*dx*dx/((c((j-1)*N+1)+Gy(j,1)*(c((j-1)*N+2)-c((j-1)*N+1)))^2)-mx((j-2)*N+1)^2*dx*dx*Gx(j-1,1)/((c((j-2)*N+1)+Gx(j-1,1)*(c((j-1)*N+1)-c((j-2)*N+1)))^2)+dtbeta^2*(-G1x(j,1)*(log(c(j*N+1))-log(c((j-1)*N+1)))*dx-G1x(j,1)*(c(j*N+1)-c((j-1)*N+1))/c((j-1)*N+1)*dx-G1y(j,1)*(log(c((j-1)*N+2))-log(c((j-1)*N+1)))*dx-G1y(j,1)*(c((j-1)*N+2)-c((j-1)*N+1))/c((j-1)*N+1)*dx+G1x(j-1,1)*(log(c((j-1)*N+1))-log(c((j-2)*N+1)))*dx+G1x(j-1,1)*(c((j-1)*N+1)-c((j-2)*N+1))/c((j-1)*N+1)*dx)+2*dt*D*(log(c((j-1)*N+1))+1+z*V(j,1))*dx*dx;
end

for j = 2:N-1
for l = 2:N-1
   g((j-1)*N+l) = mx((j-1)*N+l)^2*(Gx(j,l)-1)*dx*dx/((c((j-1)*N+l)+Gx(j,l)*(c(j*N+l)-c((j-1)*N+l)))^2)+my((j-1)*(N-1)+l)^2*(Gy(j,l)-1)*dx*dx/((c((j-1)*N+l)+Gy(j,l)*(c((j-1)*N+l+1)-c((j-1)*N+l)))^2)-mx((j-2)*N+l)^2*dx*dx*Gx(j-1,l)/((c((j-2)*N+l)+Gx(j-1,l)*(c((j-1)*N+l)-c((j-2)*N+l)))^2)-my((j-1)*(N-1)+l-1)^2*dx*dx*Gy(j,l-1)/((c((j-1)*N+l-1)+Gy(j,l-1)*(c((j-1)*N+l)-c((j-1)*N+l-1)))^2)+dtbeta^2*(-G1x(j,l)*(log(c(j*N+l))-log(c((j-1)*N+l)))*dx-G1x(j,l)*(c(j*N+l)-c((j-1)*N+l))/c((j-1)*N+l)*dx-G1y(j,l)*(log(c((j-1)*N+l+1))-log(c((j-1)*N+l)))*dx-G1y(j,l)*(c((j-1)*N+l+1)-c((j-1)*N+l))/c((j-1)*N+l)*dx+G1x(j-1,l)*(log(c((j-1)*N+l))-log(c((j-2)*N+l)))*dx+G1x(j-1,l)*(c((j-1)*N+l)-c((j-2)*N+l))/c((j-1)*N+l)*dx+G1y(j,l-1)*(log(c((j-1)*N+l))-log(c((j-1)*N+l-1)))*dx+G1y(j,l-1)*(c((j-1)*N+l)-c((j-1)*N+l-1))/c((j-1)*N+l)*dx)+2*dt*D*(log(c((j-1)*N+l))+1+z*V(j,l))*dx*dx;
end
end

for j = 2:N-1
   g((j-1)*N+N) = mx((j-1)*N+N)^2*(Gx(j,N)-1)*dx*dx/((c((j-1)*N+N)+Gx(j,N)*(c(j*N+N)-c((j-1)*N+N)))^2)-mx((j-2)*N+N)^2*dx*dx*Gx(j-1,N)/((c((j-2)*N+N)+Gx(j-1,N)*(c((j-1)*N+N)-c((j-2)*N+N)))^2)-my((j-1)*(N-1)+N-1)^2*dx*dx*Gy(j,N-1)/(c((j-1)*N+N-1)+Gy(j,N-1)*(c((j-1)*N+N)-c((j-1)*N+N-1)))^2+dtbeta^2*(-G1x(j,N)*(log(c(j*N+N))-log(c((j-1)*N+N)))*dx-G1x(j,N)*(c(j*N+N)-c((j-1)*N+N))/c((j-1)*N+N)*dx+G1x(j-1,N)*(log(c((j-1)*N+N))-log(c((j-2)*N+N)))*dx+G1x(j-1,N)*(c((j-1)*N+N)-c((j-2)*N+N))/c((j-1)*N+N)*dx+G1y(j,N-1)*(log(c((j-1)*N+N))-log(c((j-1)*N+N-1)))*dx+G1y(j,N-1)*(c((j-1)*N+N)-c((j-1)*N+N-1))/c((j-1)*N+N)*dx)+2*dt*D*(log(c((j-1)*N+N))+1+z*V(j,N))*dx*dx;
end

  g((N-1)*N+1) = my((N-1)*(N-1)+1)^2*dx*dx*(Gy(N,1)-1)/((c((N-1)*N+1)+Gy(N,1)*(c((N-1)*N+2)-c((N-1)*N+1)))^2)-mx((N-2)*N+1)^2*dx*dx*Gx(N-1,1)/((c((N-2)*N+1)+Gx(N-1,1)*(c((N-1)*N+1)-c((N-2)*N+1)))^2)+dtbeta^2*(-G1y(N,1)*(log(c((N-1)*N+2))-log(c((N-1)*N+1)))*dx-G1y(N,1)*(c((N-1)*N+2)-c((N-1)*N+1))/c((N-1)*N+1)*dx+G1x(N-1,1)*(log(c((N-1)*N+1))-log(c((N-2)*N+1)))*dx+G1x(N-1,1)*(c((N-1)*N+1)-c((N-2)*N+1))/c((N-1)*N+1)*dx)+2*dt*D*(log(c((N-1)*N+1))+1+z*V(N,1))*dx*dx;

  for l = 2:N-1
  g((N-1)*N+l) = my((N-1)*(N-1)+l)^2*dx*dx*(Gy(N,l)-1)/((c((N-1)*N+l)+Gy(N,l)*(c((N-1)*N+l+1)-c((N-1)*N+l)))^2)-mx((N-2)*N+l)^2*dx*dx*Gx(N-1,l)/((c((N-2)*N+l)+Gx(N-1,l)*(c((N-1)*N+l)-c((N-2)*N+l)))^2)-my((N-1)*(N-1)+l-1)^2*Gy(N,l-1)*dx*dx/((c((N-1)*N+l-1)+Gy(N,l-1)*(c((N-1)*N+l)-c((N-1)*N+l-1)))^2)+dtbeta^2*(-G1y(N,l)*(log(c((N-1)*N+l+1))-log(c((N-1)*N+l)))*dx-G1y(N,l)*(c((N-1)*N+l+1)-c((N-1)*N+l))/c((N-1)*N+l)*dx+G1x(N-1,l)*(log(c((N-1)*N+l))-log(c((N-2)*N+l)))*dx+G1x(N-1,l)*(c((N-1)*N+l)-c((N-2)*N+l))/c((N-1)*N+l)*dx+G1y(N,l-1)*(log(c((N-1)*N+l))-log(c((N-1)*N+l-1)))*dx+G1y(N,l-1)*(c((N-1)*N+l)-c((N-1)*N+l-1))/c((N-1)*N+l)*dx)+2*dt*D*(log(c((N-1)*N+l))+1+z*V(N,l))*dx*dx;
  end

  g(N*N) = -mx((N-2)*N+N)^2*dx*dx*Gx(N-1,N)/((c((N-2)*N+N)+Gx(N-1,N)*(c((N-1)*N+N)-c((N-2)*N+N)))^2)-my((N-1)*(N-1)+N-1)^2*Gy(N,N-1)*dx*dx/((c((N-1)*N+N-1)+Gy(N,N-1)*(c((N-1)*N+N)-c((N-1)*N+N-1)))^2)+dtbeta^2*(G1x(N-1,N)*(log(c((N-1)*N+N))-log(c((N-2)*N+N)))*dx+G1x(N-1,N)*(c((N-1)*N+N)-c((N-2)*N+N))/c((N-1)*N+N)*dx+G1y(N,N-1)*(log(c((N-1)*N+N))-log(c((N-1)*N+N-1)))*dx+G1y(N,N-1)*(c((N-1)*N+N)-c((N-1)*N+N-1))/c((N-1)*N+N)*dx)+2*dt*D*(log(c((N-1)*N+N))+1+z*V(N,N))*dx*dx;
   num = N*N;
   for j = 1:N-1
      for l = 1:N
         g(num+(j-1)*N+l) = 2*mx((j-1)*N+l)*dx*dx/(c((j-1)*N+l)+Gx(j,l)*(c(j*N+l)-c((j-1)*N+l)));
      end
   end
   num = N*N+N*(N-1);
   for j = 1:N
      for l = 1:N-1
         g(num+(j-1)*(N-1)+l) = 2*my((j-1)*(N-1)+l)*dx*dx/(c((j-1)*N+l)+Gy(j,l)*(c((j-1)*N+l+1)-c((j-1)*N+l)));
      end
   end

function J = jacobianstructure (auxdata)
   [N dx dx_mid dt D z dtbeta V Gx Gy G1x G1y] = deal(auxdata{:});
   temp = diag(-1 * ones(1, N-1),-1) + diag(ones(1,N),0);
   My_temp = temp(:,1:N-1);
   My = My_temp;
   for j = 1:N-1
      My = blkdiag(My, My_temp);
   end
    Mx_temp = eye(N*N) + diag(-1.*ones(1,N*(N-1)),-N);
    Mx = Mx_temp(:,1:(N-1)*N);
    J = [eye(N*N) Mx My
         eye(N*N) zeros(N*N,(N-1)*N) zeros(N*N,N*(N-1))];
   J = sparse(J);

   
 function J = jacobian (u, auxdata)
   [N dx dx_mid dt D z dtbeta V Gx Gy G1x G1y] = deal(auxdata{:});
   temp = diag(-1 * ones(1, N-1),-1) + diag(ones(1,N),0);
   My_temp = temp(:,1:N-1);
   My = My_temp;
   for j = 1:N-1
      My = blkdiag(My, My_temp);
   end
    My = 1/dx_mid .* My;
    Mx_temp = eye(N*N) + diag(-1.*ones(1,N*(N-1)),-N);
    Mx = Mx_temp(:,1:(N-1)*N);
    Mx = 1/dx_mid .* Mx;
    J = [eye(N*N) Mx My
         eye(N*N) zeros(N*N,(N-1)*N) zeros(N*N,N*(N-1))];
   J = sparse(J);


