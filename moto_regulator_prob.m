%% 
clc
clear all
close all

G = 1.32e-6; L = 0.0136; C = 0.0246; Ra = 7.9; D = 0.0246;
Q = [1 0; 0 1];R = 1;

lmda = 0.5;
A = [0 D/G; -C/(L*lmda) -Ra/(L*lmda)];
B1 = 0; B2 = 1/(L*lmda); B = [B1;B2];

S1 = B1 * R^(-1) * B1;
S2 = B2 * R^(-1) * B2;
S = B1 * R^(-1) * B2;

[P,K,Ei] = icare(A, B, Q, R);
Kr =  R^-1 * B' * P;

x0 = [0.1 0.1];
tspan = linspace(0.001, 0.01,1000);
[t, x] = ode45(@(t, x) (A - B*K) * x, tspan, x0);

figure(1);
subplot(2,1,1);
plot(t,x(:,1));
subplot(2,1,2);
plot(t,x(:,2));


u = -K*x';

figure(2);
plot(t,u);


%%
K3 = 2.77;
E1 = (S-K3*A(1,2))/(A(2,2)-S2*K3);
E2 =  (A(2,1)*K3 + Q(1,2))/(A(2,2)-S2*K3);

Aa = A(1,1) + E1*A(2,1)+S*E2+E1*S1*E2;
Ba = B1+ E1*B2;
Qa = -E2*A(2,1) - A(2,1)*E2 -  E2*S2*E2+ Q(1,1);

k0 = 0.01;
tspan1 = linspace(0.01,0.001,1000);
fx = @(t,k1)(33 * k1 + 53683 * k1.^2 - 1.013);
[t, K1] = ode45(fx, tspan1, k0);
K2(:,1) = K1(:,1)*E1- E2;
figure(4);
plot(tspan,K1);

%%
x_values(1)= 0;
z_values(1) =0;
u_values(1,1) = -0.2;
del = 0.1;

for i = 1:length(tspan)-1
    x_values(i+1) = x_values(i) + del*((Aa*x_values(i) + Ba*u(i)));
    z_values(i+1) = -(1/(A(2, 2) - S2 *2.77)) * (A(2, 1) - S *K1(i) - S2 * K2(i)) * x_values(i);
    u_values(i+1) = -R^(-1) * ((B1 * K1(i) + B2 * K1(i, 1)) * x(i) + (0.5* B1 * K1(i, 1) + B2 * K1(i, 1)) * z_values(i));
end

figure(5);
plot(tspan,(u_values),LineStyle="-", Color='b',LineWidth=1.5);
xlabel('Time(in ms)','FontSize',16);
ylabel('Control input(u)','FontSize',16);
hold on;
plot(tspan,(u),LineStyle="--", Color='r',LineWidth=1.5);
title("u(t) vs t",'FontSize',16)
legend('Regulator Problem using Singular Perturbation', 'As Normal regulator Problem','FontSize',16)

figure(6);
subplot(2,1,1);
plot(tspan/2,z_values,Color = 'b',LineWidth=1.5,LineStyle='-');
xlabel('Time(in ms)','FontSize',16);
ylabel( "Current 'i'",'FontSize',16);
title("States",'FontSize',16)
hold on;
plot(tspan/2,x(:,2),Color = 'r',LineWidth=1.5,LineStyle='--');
legend('Regulator Problem using Singular Perturbation', 'As Normal regulator Problem','FontSize',16)

subplot(2,1,2);
plot(tspan/2,x_values,Color = 'b',LineWidth=1.5,LineStyle='-');
xlabel('Time(in ms)','FontSize',16);
ylabel("\omega ",'FontSize',16);
hold on;
plot(tspan/2,x(:,1),Color = 'r',LineWidth=1.5,LineStyle='--' );
legend('Regulator Problem using Singular Perturbation', 'As Normal regulator Problem','FontSize',16)


