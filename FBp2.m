function FunB = FBp2(yA,yB,N,h,D,beta,eta,delta,dt,yBi)
%FBp2 returns vector function FunB for N nodes in yB for part 2

FunB = zeros(N,1);

for i=1:1:N
    if i == 1
        FunB(i) = yBi(i)-yB(i)+(D*(yB(i+1)-2*yB(i)+yB(i+1)+2*h*(beta-yB(i)))/h^2-2*yA(i)*yB(i)^2-delta*yB(i))*dt;
    elseif i == N
        FunB(i) = yBi(i)-yB(i)+(D*(yB(i-1)-2*h*eta*yB(i)^2-2*yB(i)+yB(i-1))/h^2-2*yA(i)*yB(i)^2-delta*yB(i))*dt;
    else
        FunB(i) = yBi(i)-yB(i)+(D*(yB(i+1)-2*yB(i)+yB(i-1))/h^2-2*yA(i)*yB(i)^2-delta*yB(i))*dt;
    end % End if statement
end % End for loop

end