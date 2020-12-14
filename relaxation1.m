function [c mx my status] = relaxation1(c, mx, my, N, dx, dx_mid, dt, D, z, dtbeta, V, j, Gx, Gy, G1x, G1y)
dt1 = dt/10;
flag1 = 0;
j1 = 1;
while j1 <= 10
   [u, info] = opt_func1(c, mx, my, N, dx, dx_mid, dt1, D, z, dtbeta, V, flag1, Gx, Gy, G1x, G1y);
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
           [u info] = opt_func1(c, mx, my, N, dx, dx_mid, dt2, D, z, dtbeta, V, flag2, Gx, Gy, G1x, G1y);
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
                  [c mx, my] = deal(u{:});
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
     [c mx, my] = deal(u{:});
     flag1 = 0;
     j1 = j1 + 1;
  end % if info.status < 0
end % while j1
if j1 > 10
   status = 1;
else
   status = 0;
end
