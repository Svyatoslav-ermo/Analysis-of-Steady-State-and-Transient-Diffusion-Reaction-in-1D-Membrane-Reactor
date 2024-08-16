function Yf = ImplEuler(yAi,yBi,yUi,yFi,N,dt,h,D,gamma,epsilon,delta,eta,zeta,theta,beta)
%IMPLEULER Returns Yf after one iteration over dt
%   Detailed explanation goes here

%  Make guess vectors
yA = zeros(N,1);
yB = zeros(N,1);
yU = zeros(N,1);
yF = zeros(N,1);

Y0 = zeros(4*N-2,1);

%  Combine vectors without initial conditions
Y0 = [yA; yB; yU(2:N); yF(2:N)];

%  initiate error
error = 100;

while error > 10^-2

    %  Assign vector F(Y0)
    FunA = FAp2(yA,yB,N,h,D,gamma,epsilon,dt,yAi);
    FunB = FBp2(yA,yB,N,h,D,beta,eta,delta,dt,yBi);
    FunU = FUp2(yU,yB,N,h,D,theta,delta,zeta,dt,yUi);
    FunF = FFp2(yF,yU,N,h,D,zeta,dt,yFi);

    %  Combine to construct F solution vector
    F = [FunA; FunB; FunU; FunF];

    %  Get Jacobian
    J = JacobianP2(yA,yB,yU,yF,N,h,D,gamma,epsilon,delta,eta,zeta,theta,dt);

    %  Use LUSolver
    dY = LUSolver(J,-F);

    %  Calculate Y at new iteration
    Y = Y0 + dY;

    %  Calculate 2-norm of vector
    error = norm(dY)^2;

    %  Display the error
    %fprintf('Implicit Euler error: %.2f\n',error)

    %  Assign Y as previous iteration before new iteration
    Y0 = Y;
    %  Obtain vectors yA-yF from Y0
    yA = Y0(1:N);
    yB = Y0(N+1:2*N);
    yU(2:N) = Y0(2*N+1:3*N-1);
    yF(2:N) = Y0(3*N:4*N-2);

end

%  Construct Yf
Yf = [yA; yB; yU; yF];


end