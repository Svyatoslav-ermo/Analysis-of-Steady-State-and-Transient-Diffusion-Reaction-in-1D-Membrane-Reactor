function J = JacobianP2(yA,yB,yU,yF,N,h,D,gamma,epsilon,delta,eta,zeta,theta,dt)
%JAOBIAN Returns Jacobian for part 2

%  Allocate Jacobian segments
JAdyA = zeros(N,N);
JAdyB = zeros(N,N);
JAdyU = zeros(N,N-1);
JAdyF = zeros(N,N-1);

JBdyA = zeros(N,N);
JBdyB = zeros(N,N);
JBdyU = zeros(N,N-1);
JBdyF = zeros(N,N-1);

JUdyA = zeros(N-1,N);
JUdyB = zeros(N-1,N);
JUdyU = zeros(N-1,N-1);
JUdyF = zeros(N-1,N-1);

JFdyA = zeros(N-1,N);
JFdyB = zeros(N-1,N);
JFdyU = zeros(N-1,N-1);
JFdyF = zeros(N-1,N-1);

%  Taking derivative of JA with respect to yA
for row = 1:1:N
    for col = 1:1:N
        if row == 1 && col == 1
            JAdyA(row,col) = (D*(-2-2*h)/h^2-yB(row)^2-gamma)*dt-1;

        elseif row == 1 && col == 2
            JAdyA(row,col) = (D*2/h^2)*dt;

        elseif row == N && col == N-1
            JAdyA(row,col) = (D*2/h^2)*dt;

        elseif row == N && col == N
            JAdyA(row,col) = (D*(-2*h*epsilon-2)/h^2-yB(row)^2-gamma)*dt-1;

        elseif abs(row-col) == 1
            JAdyA(row,col) = (D/h^2)*dt;

        elseif row == col
            JAdyA(row,col) = (-D*2/h^2-yB(row)^2-gamma)*dt-1;

        end

    end
end

%  Taking derivative of JA with respect to yB
for row = 1:1:N
    for col = 1:1:N
        if row == 1 && col == 1
            JAdyB(row,col) = (-2*yA(row)*yB(row))*dt;

        elseif row == N && col == N
            JAdyB(row,col) = (-2*yA(row)*yB(row))*dt;

        elseif row == col
            JAdyB(row,col) = (-2*yA(row)*yB(row))*dt;

        end

    end
end

%  Taking derivative of JB with respect to yA
for row = 1:1:N
    for col = 1:1:N
        if row == 1 && col == 1
            JBdyA(row,col) = (-2*yB(row)^2)*dt;

        elseif row == N && col == N
            JBdyA(row,col) = (-2*yB(row)^2)*dt;

        elseif row == col
            JBdyA(row,col) = (-2*yB(row)^2)*dt;

        end

    end
end

%  Taking derivative of JB with respect to yB
for row = 1:1:N
    for col = 1:1:N
        if row == 1 && col == 1
            JBdyB(row,col) = (D*(-2-2*h)/h^2-4*yA(row)*yB(row)-delta)*dt-1;

        elseif row == 1 && col == 2
            JBdyB(row,col) = (D*2/h^2)*dt;

        elseif row == N && col == N-1
            JBdyB(row,col) = (D*2/h^2)*dt;

        elseif row == N && col == N
            JBdyB(row,col) = (D*(-4*h*eta*yB(row)-2)/h^2-4*yA(row)*yB(row)-delta)*dt-1;

        elseif abs(row-col) == 1
            JBdyB(row,col) = (D/h^2)*dt;

        elseif row == col
            JBdyB(row,col) = (-D*2/h^2-4*yA(row)*yB(row)-delta)*dt-1;

        end

    end
end

%  Taking derivative of JU with respect to yB
for row = 1:1:N-1
    for col = 1:1:N

        if (col-row) == 1
            JUdyB(row,col) = delta*dt;

        end

    end
end

%  Taking derivative of JU with respect to yU
for row = 1:1:N-1
    for col = 1:1:N-1

        if row == (N-1) && col == (N-2)
            JUdyU(row,col) = (2*D/h^2)*dt;

        elseif row == (N-1) && col == (N-1)
            JUdyU(row,col) = (D*(-2*h*theta-2)/h^2-zeta)*dt-1;

        elseif abs(row-col) == 1
            JUdyU(row,col) = (D/h^2)*dt;
        
        elseif row == col
            JUdyU(row,col) = (-2*D/h^2-zeta)*dt-1;
        end

    end
end

%  Taking derivative of JF with respect to yU
for row = 1:1:N-1
    for col = 1:1:N-1

        if row == col
            JFdyU(row,col) = zeta*dt;
        end

    end
end

%  Taking derivative of JF with respect to yF
for row = 1:1:N-1
    for col = 1:1:N-1

        if row == (N-1) && col == (N-2)
            JFdyF(row,col) = (2*D/h^2)*dt;

        elseif abs(row-col) == 1
            JFdyF(row,col) = (D/h^2)*dt;

        elseif row == col
            JFdyF(row,col) = (-2*D/h^2+zeta)*dt-1;
        end

    end
end

J = [JAdyA, JAdyB, JAdyU, JAdyF;
     JBdyA, JBdyB, JBdyU, JBdyF;
     JUdyA, JUdyB, JUdyU, JUdyF;
     JFdyA, JFdyB, JFdyU, JFdyF];

end