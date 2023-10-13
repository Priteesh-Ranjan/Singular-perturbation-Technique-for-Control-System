
%clear all
close all
format long

% Define the parameters
R = 7.56; L = 0.055; K = 4.23; J = 0.136;
Tl = 14.3; i = 3.38; u = 246.8;
tm = J*R/K^2; te = R/L;
ep = te/tm;

% Define the function representing the system of differential equations
fx = @(t, x) [(K/J)*x(2) - Tl/J; -(K/L)*x(1) - (R/L)*x(2) + u/L];

% Define the time span with the specified step length
tspan = linspace(0,0.5,1000);

% Define the initial conditions
x1_0 =0.1;  % Replace with your initial value for x1
x2_0 = 0.1;  % Replace with your initial value for x2
x0 = [x1_0, x2_0];

% Solve the system of differential equations using ode45
[t, X] = ode45(fx, tspan, x0);

% Extract the solutions for x1 and x2
x1 = X(:, 1);
x2 = X(:, 2);

% Plot the solutions of x1 and x2
figure(1);
plot(t, x1);
xlabel('Time');
ylabel('x1');

figure(2);
plot(t, x2);
xlabel('Time');
ylabel('x2');



fx_spt = @(t,x1)(-x1(1)/tm + u/(tm*K) - Tl/J);
x10 = 0.1;
[t,X1] = ode45(fx_spt,tspan,x10);

x_bar = X1(:,1);


for i = 1:length(tspan)
    z_bar(i) = (u - x_bar(i) *K)/R;
end


figure(3);
subplot(2,1,1);
plot(t, x1,LineWidth=2);
xlabel('Time(in Seconds)','FontSize', 16, 'FontWeight', 'bold');
ylabel('Speed(rad/sec)','FontSize', 16, 'FontWeight', 'bold');
hold on
plot(t, x_bar,LineStyle="--",LineWidth=2);
legend('Response of Full model {\omega}','Response of SPT model of {\omega}(rad/sec)','FontSize', 16, 'FontWeight', 'bold');
title('DC motor Step Response','FontSize', 16, 'FontWeight', 'bold')

subplot(2,1,2);
plot(t, x2,LineWidth=2);
xlabel('Time(in seconds)','FontSize', 16, 'FontWeight', 'bold');
ylabel('Current(A)','FontSize', 16, 'FontWeight', 'bold');
hold on 
plot(t, z_bar(1,:),LineStyle="--",LineWidth=2);
legend('Response of Full model i','Response of SPT model of i(A)','FontSize', 16, 'FontWeight', 'bold');

