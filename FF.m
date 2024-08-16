function FunF = FF(yF,yU,N,h,D,zeta)
%FF returns vector function FunF for N nodes in yF

FunF = zeros(N-1,1);

for i=1:1:N-1
    % If i == 1, boundary condition for yF at x = 0
    if i == 1
        FunF(i) = D*(yF(i+1)-2*yF(i))/h^2+zeta*yU(i);

    elseif i == (N-1)
        FunF(i) = D*(yF(i-1)-2*yF(i)+yF(i-1))/h^2+zeta*yU(i);

    else
        FunF(i) = D*(yF(i+1)-2*yF(i)+yF(i-1))/h^2+zeta*yU(i);

    end % End if statement
    
end % End for loop

end