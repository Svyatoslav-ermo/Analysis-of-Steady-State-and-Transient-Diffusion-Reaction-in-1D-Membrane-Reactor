function FunA = FAp2(yA,yB,N,h,D,gamma,epsilon,dt,yAi)
%FAp2 returns vector function FunA for N nodes in yA for part 2

FunA = zeros(N,1);

for i = 1:1:N
    if i == 1
        FunA(i) = yAi(i)-yA(i)+(D*(yA(i+1)-2*yA(i)+yA(i+1)+2*h*(1-yA(i)))/h^2-yA(i)*yB(i)^2-gamma*yA(i))*dt;
    elseif i == N
        FunA(i) = yAi(i)-yA(i)+(D*(yA(i-1)-2*h*epsilon*yA(i)-2*yA(i)+yA(i-1))/h^2-yA(i)*yB(i)^2-gamma*yA(i))*dt;
    else
        FunA(i) = yAi(i)-yA(i)+(D*(yA(i+1)-2*yA(i)+yA(i-1))/h^2-yA(i)*yB(i)^2-gamma*yA(i))*dt;
    end % End if statement
end % End solution for all f

end