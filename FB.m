function FunB = FB(yA,yB,N,h,D,beta,eta,delta)
%FB returns vector function FunB for N nodes in yB

FunB = zeros(N,1);

for i=1:1:N
    if i == 1
        FunB(i) = D*(yB(i+1)-2*yB(i)+yB(i+1)+2*h*(beta-yB(i)))/h^2-2*yA(i)*yB(i)^2-delta*yB(i);
    elseif i == N
        FunB(i) = D*(yB(i-1)-2*h*eta*yB(i)^2-2*yB(i)+yB(i-1))/h^2-2*yA(i)*yB(i)^2-delta*yB(i);
    else
        FunB(i) = D*(yB(i+1)-2*yB(i)+yB(i-1))/h^2-2*yA(i)*yB(i)^2-delta*yB(i);
    end % End if statement
end % End for loop

end