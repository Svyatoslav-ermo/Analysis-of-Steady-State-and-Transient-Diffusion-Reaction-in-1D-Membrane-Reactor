function FunU = FU(yU,yB,N,h,D,theta,delta,zeta)
%FU returns vector function FunU for N nodes in yU

FunU = zeros(N-1,1);

for i=1:1:N-1
    % If i == 1, boundary condition for yU at x = 0
    if i == 1
        FunU(i) = D*(yU(i+1)-2*yU(i))/h^2+delta*yB(i)-zeta*yU(i);
    elseif i == (N-1)
        FunU(i) = D*(yU(i-1)-2*h*theta*yU(i)-2*yU(i)+yU(i-1))/h^2+delta*yB(i)-zeta*yU(i);
    else
        FunU(i) = D*(yU(i+1)-2*yU(i)+yU(i-1))/h^2+delta*yB(i)-zeta*yU(i);

    end % End if statement

end % End for loop

end