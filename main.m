%% Clear Cache
clear all %#ok<*CLALL>
close all
clc

%% Main Script

%  Set initial percentDiff
percentDiff = 100;

%  Assign number of nodes
    N = 11;

while percentDiff > 0.01

    %  To test different conditions, uncomment one of the parts of the code for
    %  conditions a, b, and c

%     %  Assign conditions for case (a)
%     xi = 0;
%     xf = 1;
%     delta = 0;
%     zeta = 0;
%     D = 0.1;
%     beta = 1.5;
%     gamma = 0.05;
%     epsilon = 0.0;
%     eta = 0.0;
%     theta = 0.0;

%     %  Assign conditions for case (b)
%     xi = 0;
%     xf = 1;
%     delta = 0.05;
%     zeta = 0;
%     D = 0.1;
%     beta = 1.5;
%     gamma = 0.02;
%     epsilon = 0.1;
%     eta = 0.05;
%     theta = 0.1;

    %  Assign conditions for case (c)
    xi = 0;
    xf = 1;
    delta = 0.05;
    zeta = 0.03;
    D = 0.1;
    beta = 1.5;
    gamma = 0.02;
    epsilon = 0.1;
    eta = 0.05;
    theta = 0.1;


    %  Assign vector x
    x = linspace(xi,xf,N);
    x = transpose(x);

    %  Assign distance
    h = (x(N) - x(1))/(N-1);

    %% Part 1

    %  Uncomment proper initial conditions for corresponding part before
    %  running the code
%     % Initial conditions for part (a)
%     yA = zeros(N,1);
%     yB = zeros(N,1);
%     yU = zeros(N,1);
%     yF = zeros(N,1);

    %  Initial conditions for part (b) & (c)
    yA = linspace(0.5,0.3,N);
    yA = transpose(yA);
    yB = linspace(0.5,0.2,N);
    yB = transpose(yB);
    yU = zeros(N,1);
    yF = zeros(N,1);

    %  Combine vectors without initial conditions
    Y0 = [yA; yB; yU(2:N); yF(2:N)];

    %  Initiate error
    error = 100;

    while error > 10^-2

        %  Assign vector F(Y0)
        FunA = FA(yA,yB,N,h,D,gamma,epsilon);
        FunB = FB(yA,yB,N,h,D,beta,eta,delta);
        FunU = FU(yU,yB,N,h,D,theta,delta,zeta);
        FunF = FF(yF,yU,N,h,D,zeta);

        %  Combine to construct F solution vector
        F = [FunA; FunB; FunU; FunF];

        %  Get Jacobian
        J = Jacobian(yA,yB,yU,yF,N,h,D,gamma,epsilon,delta,eta,zeta,theta);

        %  Use LUSolver
        dY = LUSolver(J,-F);

        %  Calculate Y at new iteration
        Y = dY + Y0;

        %  Calculate square of 2-norm of vector
        error = norm(dY)^2;

        %  Display the error
        %fprintf('The error is: %.2f\n', error)

        %  Assign Y as previous iteration before new iteration
        Y0 = Y;
        %  Obtain vectors yA-yF from Y0
        yA = Y0(1:N);
        yB = Y0(N+1:2*N);
        yU(2:N) = Y0(2*N+1:3*N-1);
        yF(2:N) = Y0(3*N:4*N-2);

    end % End while loop

    %  Calculate the test node index
    testNodeIndex = round(0.3 / (1 / (N - 1)) + 1);

    %  If N is greater than 11
    if N > 11
        
        %  Calculate percent difference
        percentDiff = (abs(testValue - yA(testNodeIndex)) / ((testValue + yA(testNodeIndex)) / 2)) * 100;

    end

    fprintf('Percent difference for N = %i is %.4f\n', N, percentDiff)

    %  Save value to be tested for next iteration
    testValue = yA(testNodeIndex);

    %  Increment the number of nodes by 10
    N = N + 10;


end

%  Get least number of nodes for percent difference accuracy to be less
%  than 0.1%
numberOfNodesRequired = N - 10;

%  Decrement number of nodes by 10
N = N - 10;

fprintf('Number of nodes needed for percent difference between values yA for N nodes\nand yA for N-10 nodes at x = 0.3 to be less than 0.01%% is %i\n', numberOfNodesRequired)

%  Copy vector to Y1Result for comparison with next part
Y1result = [yA; yB; yU; yF];

%   Plot solutions
figure(1) % Create figure 1 for the plot
hold on % Draw all on the same figure
grid on % Turn grid on
%  Plot yA
plot(x,yA,'b','LineWidth',3)
%  Plot yB
plot(x,yB,'r','LineWidth',3)
%  Plot yU
plot(x,yU,'g','LineWidth',3)
%  Plot yF
plot(x,yF,'m','LineWidth',3)

title('Chemical Species Concentrations for Part (c)','FontSize',24)
xlabel('X')
ylabel('Y')
%   Set limits on axes
xlim([min(x) max(x)])
ylim([min(Y0) max(Y0)])

%   Set screen position and figure size
set(gcf,'Position',[75 75 1275 600])
%   Format axes lengths and label fonts
set(gca,'LineWidth',3,'FontSize',20)

%   Set figure legend
legend('yA','yB','yU','yF')

%% Part 2

%  Set initial values for all species
yA = zeros(N,1);
yB = zeros(N,1);
yU = zeros(N,1);
yF = zeros(N,1);

%  Set time step
dt = 1;

%  Set current time
t = 0;

%  Set time coordinate vector of N nodes for 3-D scatter plot
tVec = ones(N,1)*t;

%% Implicit Euler

%  Construct vector Yi at time step t
Yi = [yA; yB; yU; yF];

%  Set norm
n = 1;

%  Create figure 2
figure(2)

hold on
grid on

%  Plot first layer of scatter 3-D plot at t = 0
scatter3(x,yA,tVec,'b')
scatter3(x,yB,tVec,'r')
scatter3(x,yU,tVec,'g')
scatter3(x,yF,tVec,'m')

while n > 0.01
    %  Call implicit Euler to advance over 1 iteration
    Yf = ImplEuler(yA,yB,yU,yF,N,dt,h,D,gamma,epsilon,delta,eta,zeta,theta,beta);

    %  Advance current time
    t = t + dt;

    %  Create new time coordinate vector for 3-D scatter plot
    tVec = ones(N,1)*t;

    %  Calculate norm
    n = norm(Y1result-Yf)^2;

    fprintf('Norm at t = %.2f is: %.2f\n', t, n)

    %  Set vectors for next iteration
    yA = Yf(1:N);
    yB = Yf(N+1:2*N);
    yU = Yf(2*N+1:3*N);
    yF = Yf(3*N+1:4*N);

    %  Plot new layer of scatter 3-D plot
    scatter3(x,yA,tVec,'b')
    scatter3(x,yB,tVec,'r')
    scatter3(x,yU,tVec,'g')
    scatter3(x,yF,tVec,'m')

end

xlabel('X')
ylabel('Y')
zlabel('Time t')

title('Chemical Species Evolution Overtime for Part (c)','FontSize',24)

%   Set screen position and figure size
set(gcf,'Position',[75 75 1275 600])

%   Format axes lengths and label fonts
set(gca,'LineWidth',3,'FontSize',20)

%   Set figure legend
legend('yA','yB','yU','yF')

xlim([min(x) max(x)])
ylim([min(Yf) max(Yf)])
zlim([0 t])


% 
% set(gcf,'Position',[75 75 1275 600])
% 
% %% Explicit Euler
% 
% % n = 1;
% % 
% % while n > 0.01
% % 
% %     FunA = FA(yA,yB,N,h,D,gamma,epsilon);
% %     FunB = FB(yA,yB,N,h,D,beta,eta,delta);
% %     FunU = FU(yU,yB,N,h,D,theta,delta,zeta);
% %     FunF = FF(yF,yU,N,h,D,zeta);
% % 
% %     %  Combine to construct F solution vector
% %     F = [FunA; FunB; FunU; FunF];
% % 
% %     %  Calculate Y at new time step
% %     Y = Y0 + F*dt;
% % 
% %     %  Calculate dY
% %     dY = Y - Y1result;
% % 
% %     %  Calculate 2-norm between steady-state and transient
% %     n = norm(dY);
% % 
% %     %  Print norm
% %     fprintf('Norm is: %.2f\n', n)
% %     
% %     %  Square 2-norm
% %     n = n^2;
% % 
% %     %  Set Y to Y0 for next itration
% %     Y0 = Y;
% % 
% %     %  Obtain yA to yF from Y0
% %     yA = Y0(1:N);
% %     yB = Y0(N+1:2*N);
% %     yU(2:N) = Y0(2*N+1:3*N-1);
% %     yF(2:N) = Y0(3*N:4*N-2);
% % 
% % end % End while loop
% % 
