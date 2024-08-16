function FunF = FFp2(yF,yU,N,h,D,zeta,dt,yFi)
%FFp2 returns vector function FunF for N nodes in yF for part 2

FunF = zeros(N-1,1);

for i=1:1:N-1
    % If i == 1, boundary condition for yF at x = 0
    if i == 1
        FunF(i) = yFi(i)-yF(i)+(D*(yF(i+1)-2*yF(i))/h^2+zeta*yU(i))*dt;

    elseif i == (N-1)
        FunF(i) = yFi(i)-yF(i)+(D*(yF(i-1)-2*yF(i)+yF(i-1))/h^2+zeta*yU(i))*dt;

    else
        FunF(i) = yFi(i)-yF(i)+(D*(yF(i+1)-2*yF(i)+yF(i-1))/h^2+zeta*yU(i))*dt;

    end % End if statement
    
end % End for loop

end