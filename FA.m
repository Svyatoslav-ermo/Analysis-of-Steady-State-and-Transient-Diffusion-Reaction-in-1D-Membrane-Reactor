function FunA = FA(yA,yB,N,h,D,gamma,epsilon)
%FA returns vector function FunA for N nodes in yA

FunA = zeros(N,1);

for i = 1:1:N
    if i == 1
        FunA(i) = D*(yA(i+1)-2*yA(i)+yA(i+1)+2*h*(1-yA(i)))/h^2-yA(i)*yB(i)^2-gamma*yA(i);
    elseif i == N
        FunA(i) = D*(yA(i-1)-2*h*epsilon*yA(i)-2*yA(i)+yA(i-1))/h^2-yA(i)*yB(i)^2-gamma*yA(i);
    else
        FunA(i) = D*(yA(i+1)-2*yA(i)+yA(i-1))/h^2-yA(i)*yB(i)^2-gamma*yA(i);
    end % End if statement
end % End solution for all f

end