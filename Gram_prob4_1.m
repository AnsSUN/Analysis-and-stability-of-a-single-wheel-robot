%% Calculating input come to initial condition

clear all, clc, close all

%% Initial values of the parameters

J0=4.5*(10^-8);
J1=6.7*(10^-3);
J2=0.9375;
R0=0.025;
R1=0.124;
M1=0.65;
M2=30;
L=0.5;
be=1.85*(10^-3);
g=9.81;

%% Values of the coefficients

h1=J1+J0*((R1^2)/(R0^2))+ (M1+M2)*(R1^2);
h2=M2*L*R1;
h3=be*((R1^2)/(R0^2));
h4=(R1/R0);
h5=J2+(M2*L*L);
h6=-M2*g*L;

%% GENERATING STATE SPACE EQUATION

a22=(h5*h3)/((h1*h5)-(h2*h2));
a23=(h2*h6)/((h1*h5)-(h2*h2));
a42=(h2*h3)/((h1*h5)-(h2*h2));
a43=(h1*h6)/((h1*h5)-(h2*h2));
b21=((h5*h4)-h2)/((h1*h5)-(h2*h2));
b41=(h1-(h2*h4))/((h1*h5)-(h2*h2));

A=[0, 1, 0, 0;
    0, -a22, a23, 0;
    0, 0, 0, 1;
    0, a42, -a43, 0];

B=[0; b21; 0; b41];

C=[1, 0, 0, 0;
    0, 0, 1, 0];
D=0;
t = 0:0.02:2;
%% Creating Create state-space model, convert to state-space model
sys = ss(A,B,C,D);
%{
clf
u = 0.01;
[y t]=step(sys);
plot(t,y);
figure
u = 0.001*t;
[y t]=step(sys);
plot(t,y);
%}
%% Solving the problems for project:
% For transfer function:

syms s
Phi=inv(s*eye(4)-A);
G=C*Phi*B+D;

%% Eigenvalues of the state equation:

lambda= eig(A)
% This system is not stable because one eigenvalue is positive.
% Neither Asymptotically stable nor marginally.

%% Discretization
Ts=0.03;
sysd = c2d(sys,Ts,'zoh');
%{
[Y T]=step(sysd);
plot(T,Y);
%}

%The discrete system which has less sample time represent better output for
%continuous time system.
% It is neither asymptotically stable nor marginally stble.

%% Checking the controllability of the system.

Co = ctrb(sys.A,sys.B);
unco = length(A) - rank(Co) % Determine number of uncontrollable states

% Therefore the system is  controllable.

%% Checking the observability of the system.

Ob = obsv(sys);
unob = length(A) - rank(Ob)
% Therefore the system is also observable.


syms t;
e=exp(1);
X0=[1;0;1;0];
X1=[0;0;0;0];
opt = gramOptions('TimeIntervals',[0 1]);
Wc1=int(e^(A*t)*B*B'*e^(A'*t),0,1);
t=0:0.2:1;
u_gram=-B'.*expm(A'.*(1-t)).*inv(Wc1).*(expm(1*A)*X0);
%[y,t,x]=lsim(sys,u_gram,t,X0);
figure(1)
plot(t,y(:,1),'linewidth',2)
hold on
plot(t,y(:,2),'linewidth',2)
xlabel('t')
ylabel('amplitude')
legend('theta1','theta2')
title('Open loop behavior of initial condition and input')