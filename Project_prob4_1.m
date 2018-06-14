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
x0=[2; 0; 3;0]; % assuming initial position of the ball bot.
%% Creating Create state-space model, convert to state-space model
sys = ss(A,B,C,D);
opt = gramOptions('TimeIntervals',[0 5]);
Wc=gram(sys,'C', opt);
last=inv(Wc)*[2; 0; 3;0]; 
u=-B'*exp((A')*(5-t))*inv(Wc)*(exp(A*5)*x0);
t=[0:.05:5];