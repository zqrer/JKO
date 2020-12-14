function [u, info] = opt_func1(cc, m, N, Dx, Dx_mid, dt, D, z, dtbeta, V, flag, G, G1)

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
   if flag ~= 0
      temp = diag(-1 * ones(1, N-1),-1) + diag(ones(1,N),0);
      %M = 1./Dx * temp(:,1:N-1);
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

   options.auxdata = {N dt D z V Dx Dx_mid G G1};



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
   [N dt D z V Dx Dx_mid G G1] = deal(auxdata{:});
  % dtbeta = dt*ones(N-1,1);
  dtbeta = sqrt(1)*dt*ones(N-1,1);
%   dtbeta = sqrt(50)*dt*ones(N-1,1);
%   dtbeta = zeros(N-1,1);
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
   for j = 1:N-1
      f = f + Dx(j)*m(j)^2/((c(j+1)-c(j))*G(j)+c(j))+(dtbeta(j))^2*(c(j+1)-c(j))*(log(c(j+1))-log(c(j)))*G1(j)+2*dt*D*(c(j)*log(c(j))+z*V(j)*c(j))*Dx_mid(j);
   end
    f = f+2*dt*D*(c(N)*log(c(N))+z*V(N)*c(N))*Dx_mid(N);
   


function C = constraints(u, auxdata)
   [N dt D z V Dx Dx_mid G G1] = deal(auxdata{:});
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
   [N dt D z V Dx Dx_mid G G1] = deal(auxdata{:});
   %dtbeta = zeros(N-1,1);
   dtbeta = sqrt(1)*dt*ones(N-1,1);
   %dtbeta = sqrt(50)*dt*ones(N-1,1);
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
      g(1) = Dx(1)*m(1)^2*(G(1)-1)/(((1-G(1))*c(1)+G(1)*c(2))^2)+(dtbeta(1))^2*(-G1(1)*(log(c(2))-log(c(1)))-G1(1)*(c(2)-c(1))/c(1))+2*dt*D*(log(c(1))+z*V(1)+1)*Dx_mid(1);
   for k = 2:N-1
      g(k) = Dx(k)*m(k)^2*(G(k)-1)/(((1-G(k))*c(k)+G(k)*c(k+1))^2)-G(k-1)*m(k-1)^2*Dx(k-1)/(((1-G(k-1))*c(k-1)+G(k-1)*c(k))^2)+(dtbeta(k))^2*(-G1(k)*(log(c(k+1))-log(c(k)))-G1(k)*(c(k+1)-c(k))/c(k))+(dtbeta(k-1))^2*(G1(k-1)*(log(c(k))-log(c(k-1)))+G1(k-1)*(c(k)-c(k-1))/c(k))+2*dt*D*(log(c(k))+z*V(k)+1)*Dx_mid(k);
   end
      g(N) = -G(N-1)*m(N-1)^2*Dx(N-1)/(((1-G(N-1))*c(N-1)+G(N-1)*c(N))^2)+(dtbeta(N-1))^2*(G1(N-1)*(log(c(N))-log(c(N-1)))+G1(N-1)*(c(N)-c(N-1))/c(N))+2*dt*D*(log(c(N))+1+z*V(N))*Dx_mid(N);
   for k = 1:N-1
      g(k+N)=2*Dx(k)*m(k)/((1-G(k))*c(k)+G(k)*c(k+1));
   end


function J = jacobianstructure (auxdata)
   [N dt D z V Dx Dx_mid G G1] = deal(auxdata{:});
   temp = diag(-1 * ones(1, N-1),-1) + diag(ones(1,N),0);
   J = [eye(N) temp(:,1:N-1)
        eye(N) zeros(N,N-1)];
   J = sparse(J);

   
 function J = jacobian (u, auxdata)
   [N dt D z V Dx Dx_mid G G1] = deal(auxdata{:});
   temp = diag(-1 * ones(1, N-1),-1) + diag(ones(1,N),0);
   M = 1./Dx_mid .* temp(:,1:N-1);
   J = [eye(N) M
        eye(N) zeros(N,N-1)];
   J = sparse(J);


