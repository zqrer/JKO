function [u, info] = opt_func(cc, m, N, dx, dt, D, z, dtbeta, V, flag)

   b = cc;
   %c0 = zeros(N,1);
   %for j = 1:N
   %   if(2*cc(j) - cc_p(j) < 1e-6)
   %      c0(j) = cc(j);
   %   else
   %      c0(j) = 2*cc(j)-cc_p(j);
   %   end
   %end
   %c0 = cc;
   Dx_mid = dx*ones(N,1);
   Dx_mid(1) = Dx_mid(1);
   Dx_mid(N) = Dx_mid(N);
   if flag ~= 0
      temp = diag(-1 * ones(1, N-1),-1) + diag(ones(1,N),0);
      %M = 1/dx * temp(:,1:N-1);
      M = 1./Dx_mid .* temp(:,1:N-1);
      %M = sparse(M);
      c0 = b - M * m;
   else
      c0 = cc;
   end
   u0 = { c0
         m };

   options.cl = [b; zeros(N,1)];
   options.cu = [b; Inf * ones(N,1)];

   options.auxdata = {N dx dt D z dtbeta V Dx_mid};


    options.ipopt.jac_d_constant   = 'yes';
    options.ipopt.jac_c_constant   = 'yes';
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.derivative_test       = 'first-order';
    options.ipopt.mu_strategy      = 'adaptive';
    options.ipopt.max_iter         = 10000;
    options.ipopt.tol              = 1e-8;


    funcs.objective         = @objective;
    funcs.constraints       = @constraints;
    funcs.gradient          = @gradient;
    funcs.jacobian          = @jacobian;
    funcs.jacobianstructure = @jacobianstructure;

    [u info] = ipopt_auxdata(u0,funcs,options);



% ------------------------------------------------------------------
function f = objective (u, auxdata)
   [N dx dt D z dtbeta V Dx_mid] = deal(auxdata{:});
   f = 0.0;
   [c0 m] = deal(u{:});
   c = zeros(N,1);
   for j = 1:N
      if(c0(j)<1e-6)
         c(j) = 1e-6;
      else
         c(j) = c0(j);
      end
   end
   %dtbeta = 0.0 * (dt)^2;
   for j = 1:N-1
      %f = f + dx *(log(c(j+1))-log(c(j))) * m(j)^2/(c(j+1)-c(j)) + (dt/beta)^2 * (c(j+1)-c(j)) * (log(c(j+1))-log(c(j)))/dx + D * dx * dt * (c(j+1)^2*(log(c(j+1))-0.5)-c(j)^2*(log(c(j))-0.5)) / (c(j+1)-c(j)) + D * z * V(j) * dt * dx * (c(j+1) + c(j));
      f = f + dx *(log(c(j+1))-log(c(j))) * m(j)^2/(c(j+1)-c(j)) + (dtbeta)^2 * (c(j+1)-c(j)) * (log(c(j+1))-log(c(j)))/dx + D * dx * dt * (c(j+1)^2*(log(c(j+1))-0.5)-c(j)^2*(log(c(j))-0.5)) / (c(j+1)-c(j)) + D * z * V(j) * dt * dx * (c(j+1) + c(j));
   end
   


function C = constraints(u, auxdata)
   [N dx dt D z dtbeta V Dx_mid] = deal(auxdata{:});
   %N = auxdata{1};
   %dx = auxdata{2};
   %dt = auxdata{3};
   [c0 m] = deal(u{:});
   temp = diag(-1 * ones(1, N-1),-1) + diag(ones(1,N),0);
   %A = [eye(N) 1/dx*temp(:,1:N-1)];
   %M = 1/dx * temp(:,1:N-1);
   M = 1./Dx_mid .*temp(:,1:N-1);
   %A = sparse(A);
   M = sparse(M);
   %S = [eye(N) zeros(N,N-1)];
   %S = sparse(S);
   %C = [A * u; S * u];
   C = [c0 + M * m; c0];


function g = gradient(u, auxdata)
   [N dx dt D z dtbeta V Dx_mid] = deal(auxdata{:});
   [c0 m] = deal(u{:});
   c = zeros(N,1);
   for j = 1:N
      if(c0(j)<1e-6)
         c(j) = 1e-6;
      else
         c(j) = c0(j);
      end
   end
   g = zeros(2*N-1,1);
   %dtbeta = 0.0 * (dt)^2;
   %g(1) = dx*(log(c(2))-log(c(1)))*m(1)^2/(c(2)-c(1))^2-dx*m(1)^2/((c(2)-c(1))*c(1))-(dt/beta)^2*(c(2)-c(1))/(dx*c(1))-(dt/beta)^2*(log(c(2))-log(c(1)))/dx+dx*dt*D*(c(2)^2*(log(c(2))-0.5)-c(1)^2*(log(c(1))-0.5))/(c(2)-c(1))^2-2*dx*dt*D*c(1)*log(c(1))/(c(2)-c(1))+D*z*V(1)*dx*dt;
   g(1) = dx*(log(c(2))-log(c(1)))*m(1)^2/(c(2)-c(1))^2-dx*m(1)^2/((c(2)-c(1))*c(1))-(dtbeta)^2*(c(2)-c(1))/(dx*c(1))-(dtbeta)^2*(log(c(2))-log(c(1)))/dx+dx*dt*D*(c(2)^2*(log(c(2))-0.5)-c(1)^2*(log(c(1))-0.5))/(c(2)-c(1))^2-2*dx*dt*D*c(1)*log(c(1))/(c(2)-c(1))+D*z*V(1)*dx*dt;
   for k = 2:N-1
      %g(k) = dx*(log(c(k+1))-log(c(k)))*m(k)^2/(c(k+1)-c(k))^2-dx*(log(c(k))-log(c(k-1)))*m(k-1)^2/(c(k)-c(k-1))^2-dx*m(k)^2/((c(k+1)-c(k))*c(k))+dx*m(k-1)^2/((c(k)-c(k-1))*c(k))-(dt/beta)^2*(c(k+1)-2*c(k)+c(k-1))/(dx*c(k))-(dt/beta)^2*(log(c(k+1))-2*log(c(k))+log(c(k-1)))/dx+dx*dt*D*(c(k+1)^2*(log(c(k+1))-0.5)-c(k)^2*(log(c(k))-0.5))/(c(k+1)-c(k))^2-dx*dt*D*(c(k)^2*(log(c(k))-0.5)-c(k-1)^2*(log(c(k-1))-0.5))/(c(k)-c(k-1))^2-2*dx*dt*D*c(k)*log(c(k))*(1/(c(k+1)-c(k))-1/(c(k)-c(k-1)))+dx*dt*D*z*(V(k)+V(k-1));
      g(k) = dx*(log(c(k+1))-log(c(k)))*m(k)^2/(c(k+1)-c(k))^2-dx*(log(c(k))-log(c(k-1)))*m(k-1)^2/(c(k)-c(k-1))^2-dx*m(k)^2/((c(k+1)-c(k))*c(k))+dx*m(k-1)^2/((c(k)-c(k-1))*c(k))-(dtbeta)^2*(c(k+1)-2*c(k)+c(k-1))/(dx*c(k))-(dtbeta)^2*(log(c(k+1))-2*log(c(k))+log(c(k-1)))/dx+dx*dt*D*(c(k+1)^2*(log(c(k+1))-0.5)-c(k)^2*(log(c(k))-0.5))/(c(k+1)-c(k))^2-dx*dt*D*(c(k)^2*(log(c(k))-0.5)-c(k-1)^2*(log(c(k-1))-0.5))/(c(k)-c(k-1))^2-2*dx*dt*D*c(k)*log(c(k))*(1/(c(k+1)-c(k))-1/(c(k)-c(k-1)))+dx*dt*D*z*(V(k)+V(k-1));
   end
   %g(N) = -dx*(log(c(N))-log(c(N-1)))*m(N-1)^2/(c(N)-c(N-1))^2+dx*m(N-1)^2/(c(N)*(c(N)-c(N-1)))+(dt/beta)^2*(log(c(N))-log(c(N-1)))/dx+(dt/beta)^2*(c(N)-c(N-1))/(dx*c(N))-D*dt*dx*(c(N)^2*(log(c(N))-0.5)-c(N-1)^2*(log(c(N-1))-0.5))/(c(N)-c(N-1))^2+D*dx*dt*2*c(N)*log(c(N))/(c(N)-c(N-1))+dt*D*z*V(N-1)*dx;
   g(N) = -dx*(log(c(N))-log(c(N-1)))*m(N-1)^2/(c(N)-c(N-1))^2+dx*m(N-1)^2/(c(N)*(c(N)-c(N-1)))+(dtbeta)^2*(log(c(N))-log(c(N-1)))/dx+(dtbeta)^2*(c(N)-c(N-1))/(dx*c(N))-D*dt*dx*(c(N)^2*(log(c(N))-0.5)-c(N-1)^2*(log(c(N-1))-0.5))/(c(N)-c(N-1))^2+D*dx*dt*2*c(N)*log(c(N))/(c(N)-c(N-1))+dt*D*z*V(N-1)*dx;
   for k = 1:N-1
      g(k+N)=2*dx*(log(c(k+1))-log(c(k)))*m(k)/(c(k+1)-c(k));
   end


function J = jacobianstructure (auxdata)
   N = auxdata{1};
   temp = diag(-1 * ones(1, N-1),-1) + diag(ones(1,N),0);
   %A = [eye(N) temp(:,1:N-1)];
   %A = sparse(A);
   %S = [eye(N) zeros(N,N-1)];
   %S = sparse(S);
   %J = [A; S]; 
   J = [eye(N) temp(:,1:N-1)
        eye(N) zeros(N,N-1)];
   J = sparse(J);

   
 function J = jacobian (u, auxdata)
   [N dx dt D z dtbeta V Dx_mid] = deal(auxdata{:});
   %N = auxdata{1};
   %dx = auxdata{2};
   %dt = auxdata{3};
   temp = diag(-1 * ones(1, N-1),-1) + diag(ones(1,N),0);
   %M = 1/dx * temp(:,1:N-1);
   M = 1./Dx_mid .* temp(:,1:N-1);
   %A = [eye(N) 1/dx*temp(:,1:N-1)];
   %A = sparse(A);
   %S = [eye(N) zeros(N,N-1)];
   %S = sparse(S);
   J = [eye(N) M
        eye(N) zeros(N,N-1)];
   J = sparse(J);
   %J = [A; S];


