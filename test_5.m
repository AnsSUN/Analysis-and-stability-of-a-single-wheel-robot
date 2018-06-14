clc, clear all, close all

tspan = [0:0.2:1];
[t,theta]=ode45(@ nonlin_clloop,tspan, [0,0,0,0]);
figure
plot(tspan,theta(:,1))
hold on
plot(tspan,theta(:,3))
xlabel('Time (sec)')
ylabel('Robot Position (m)')
legend('output theta1','output theta2')

function ncl_new= nonlin_clloop(t, theta)
    

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
K =[-4.1989   -2.4420  -63.8592  -13.7542];
L =[ 21.9805,    5.0916;
    112.0003 , -429.5426;
    8.3776 ,  17.5452;
    85.5629 , 198.9117];
    %% Nonlinearfunction without noise
    %{
    tao=0.01;
    dy=zeros(4,1);
    dy(1)= theta(2);
    dy(2) =(h4/h1)*tao -(h3/h1)*theta(2)- (h2/h1)*((dy(4)*cos(theta(3)))-(theta(4)^2)*sin(theta(3)));
    dy(3)=theta(4);
    dy(4)=-(h6/h5)*sin(theta(3))-(h2/h5)*dy(2)*cos(theta(3)) + tao/h5;
    %}
    
    %% Nonlinearfunction with noise
    
    tao=0.01 ;
    ncl=zeros(4,1);
    nclold=zeros(4,1);
    ncl_new=zeros(4,1);
    ncl(1)= theta(2);
    ncl(2)=(h4/h1)*tao -(h3/h1)*theta(2)- (h2/h1)*((ncl(4)*cos(theta(3)))-(theta(4)^2)*sin(theta(3)));
    ncl(3)= theta(4);
    ncl(4)= -(h6/h5)*sin(theta(3))-(h2/h5)*ncl(2)*cos(theta(3)) + tao/h5 ;
    
    nclold(1)=theta(2);
    nclold(2)=(h4/h1)*tao -(h3/h1)*theta(2)-(h2/h1)*nclold(4);
    nclold(3)=theta(4);
    nclold(4)=tao/h5 - (h2/h5)*nclold(2) - (h6/h5)*theta(3);
    x_dot=ncl-B*K*nclold;
    %x_cap=ncl-L*C*nclold-B*K*nclold+L*C*ncl;
    ncl_new=x_dot;
end