function FunU = FUp2(yU,yB,N,h,D,theta,delta,zeta,dt,yUi)
%FUp2 returns vector function FunU for N nodes in yU for part 2

FunU = zeros(N-1,1);

for i=1:1:N-1
    % If i == 1, boundary condition for yU at x = 0
    if i == 1
        FunU(i) = yUi(i)-yU(i)+(D*(yU(i+1)-2*yU(i))/h^2+delta*yB(i)-zeta*yU(i))*dt;
    elseif i == (N-1)
        FunU(i) = yUi(i)-yU(i)+(D*(yU(i-1)-2*h*theta*yU(i)-2*yU(i)+yU(i-1))/h^2+delta*yB(i)-zeta*yU(i))*dt;
    else
        FunU(i) = yUi(i)-yU(i)+(D*(yU(i+1)-2*yU(i)+yU(i-1))/h^2+delta*yB(i)-zeta*yU(i))*dt;

    end % End if statement

end % End for loop

end