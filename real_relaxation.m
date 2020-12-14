function [c m status] = real_relaxation(c, m, N, Dx, Dx_mid, dt, D, z, dtbeta, j, t0, g, x_mid, G, G1)
dt1 = dt/10;
flag1 = 0;
j1 = 1;
while j1 <= 10
     T = dt*(j-1)+t0+dt1*j1;
     V = -1.0.*g.*sin(pi.*x)+h.*cos(w.*T);
   [u, info] = opt_func1(c, m, N, Dx, Dx_mid, dt1, D, z, dtbeta, V, flag1, G, G1);
  if info.status < 0 
     if info.status == -2 && flag1 == 0
        fprintf('\n Relaxation!! Relaxation!! Result %d : %d\n', j, j1);
        flag1 = 1;
     elseif info.status == -2 && flag1 == 1
        fprintf('\n Relaxation 3!! Result %d : %d\n', j, j1);
        dt2 = dt1 /10;
        j2 = 1;
        flag2 = 0;
        while j2 <= 10
          T = dt*(j-1)+t0+dt1*(j1-1)+dt2*j2;
          V = -1.0.*g.*sin(pi.*x)+h.*cos(w.*T);
           [u info] = opt_func1(c, m, N, Dx, Dx_mid, dt2, D, z, dtbeta, V, flag2, G, G1);
           if info.status < 0
              if info.status == -2 && flag2 == 0
                  fprintf('\n Relaxation!! Relaxation!! Relaxation!! Result %d : %d : %d\n', j, j1, j2);
                  flag2 = 1;
              else
                 fprintf('\nRelaxation!! Relaxation!! Small step!! Failed %d : %d : %d\n', j, j1, j2);
                  break;% while j2
              end
           else
                  fprintf('\nRelaxation!! Relaxation!! Small step!! Result %d : %d : %d\n', j, j1, j2);
                  [c m] = deal(u{:});
                  flag2 = 0;
                  j2 = j2 + 1;
           end% if under while j2
        end%while
        if j2 > 10
            flag1 = 0;
            j1 = j1 + 1;
        else
            fprintf('\nRelaxation!! Relaxation!! Small step!! Failed %d : %d\n', j, j1);
           break; % while j1
        end
     else
        break; % while j1
     end % if under info.status == -2
  else
     fprintf('\n Relaxation!! Small step!! Result %d : %d\n', j , j1);
     [c m] = deal(u{:});
     flag1 = 0;
     j1 = j1 + 1;
  end % if info.status < 0
end % while j1
if j1 > 10
   status = 1;
else
   status = 0;
end
