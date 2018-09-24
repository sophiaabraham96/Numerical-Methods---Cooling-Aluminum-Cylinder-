clc
close all 
clear all 
%% Introduction 
% Heat Transfer Project 
disp('Sophia Abraham') 
disp('Project 2')
disp(' ') 
%% Inputs 
disp('INPUTS')
% Diameter of cylinder 
D=45; %mm 
fprintf(' The diameter of the cylinder is %g mm\n',D)
% Length of cylinder 
L=103; %mm 
fprintf(' The length of the cylinder is %g mm\n', L)
% Density of aluminum 
D_A=2700; %kg/m^3
fprintf(' The density of aluminum is %g kg/m^3\n', D_A)
% Specific hear of aluminum 
S=901; %J/(kg-C) 
fprintf(' The specific heat of aluminum is %g J/kg-C\n', S)
% Thermal conductivity of aluminum 
T_c=240; %W/(m=C)
fprintf(' The thermal conductivity of aluminum is %g W/m-C\n\n', T_c)
% Convert units 
D=D/1000; %mm to m 
L=L/1000; %mm to m 
% Experimental Data 
disp('Experimental Data:')
t_s=[0 5 10 15 20 25 30]; 
T_C=[22.9 19.9 17.3 15.5 13.5 11.3 10.8]; 
disp('Data for coefficient of thermal expansion vs. temp for aluminum:')
t_c=[-10 77 127 177 227 277 327 377 427];
therm_ex=[5.8 9.2 13.9 25.5 32.6 34.1 36.1 35.9 39.8];
therm_ex=therm_ex * 0.000001 ;
T_C=T_C + 273.15; %convert to Kelvin 
theta_a=3.8+273.15;
disp('Time in seconds:')
disp(t_s)
disp('Recorded temperatures in degrees Celsius')
disp(T_C)
disp('Coefficient of Thermal Expansion Data:')
disp('Temperature in Celsius')
disp(t_c)
disp('Coefficient of Thermal Expansion in m/m-C')
disp(therm_ex)
%% Part 3 
disp('Part 3')
% Find convective cooling coefficient using point 4 
disp('Find the convective cooling coefficient h:') 
fprintf(' The 4th data point is %g degrees C at %f s\n',T_C(4),t_s(4))
theta_t=T_C(4);
theta_0=T_C(1);
R=D/2; %radius of cylinder 
A= 2*pi*R*L + 2*pi*R*R;
t=t_s(4);
V=pi*R*R*L; %Volume of cylinder 
m=V*D_A;
C=S;
c='(theta_t - theta_a)/(theta_0-theta_a)==exp((-h*A*t)/(m*C))';
syms h 
disp(' ')
disp(' The equation to find cooling coefficient h is:')
disp(c)
%simplifying to solve for h
h=solve(c,h);
disp('h = ') 
disp(h) 
%substitute values to solve for h 
ans=subs(h);
ans=double(ans);
disp(' The convective cooling coefficient h is: ') 
disp(ans) 
%% Part 4 
disp('Part 4')
% Regress the temperature vs. time data to the model 
syms h 
y = T_C - ((theta_a) + (T_C(1)-theta_a)*exp((-h*A*t_s)/(m*C)));
E= y.^2;
P_sum=sum(E);
par_der=diff(P_sum,h);
par_der=vpa(par_der,2);
soln=solve(par_der,h);
fin_soln=subs(soln);
disp(' The convective cooling coefficient h is: ') 
disp(fin_soln)
%% Part 5 Regress with data transformation 
disp('Part 5')
eqn_y=(T_C-theta_a)/(T_C(1)-theta_a);
eqn_z=log(eqn_y);
eqn_num=eqn_z.*t_s;
sum_1=sum(eqn_num);
eqn_den=t_s.^2;
sum_2=sum(eqn_den);
a_1=sum_1/sum_2;
eqn='(-h*A)/(m*C)==a_1';
solv_h=solve(eqn,h);
h_5=subs(solv_h);
h_5=double(h_5);
disp(' The convective cooling coefficient h is: ') 
disp(fin_soln)
%% Part 7 Plot of equations 
figure(1)
% graphing function part 3 
t_s=[0 5 10 15 20 25 30]; 
theta_t=((theta_a) + (T_C(1)-theta_a)*exp((-ans*A*t_s)/(m*C)));
% plotting the graph 
plot(t_s,theta_t,'-');
hold on 
theta_t2= ((theta_a) + (T_C(1)-theta_a)*exp((-fin_soln*A*t_s)/(m*C)));
hold on 
plot(t_s,theta_t2,'r');
theta_t3= ((theta_a) + (T_C(1)-theta_a)*exp((-h_5*A*t_s)/(m*C)));
plot(t_s,theta_t3,'c-.');
hold off 
title('Temperature vs. Time')
xlabel('Time in seconds')
ylabel('Temperature in Kelvin')
legend('Part3','Part4','Part5')
%% Part 8 2nd Order Polynomial Regression 
disp('Part 8')
% solve for coefficient matrix 
n=length(t_c);
Ti=sum(t_c);
Ti2=t_c.^2;
Ti2=sum(Ti2);
Ti3=t_c.^3;
Ti3=sum(Ti3);
Ti4=t_c.^4;
Ti4=sum(Ti4);
coef = [n Ti Ti2;Ti Ti2 Ti3;Ti2 Ti3 Ti4];
% solve for right hand side vector 
alpha=sum(therm_ex);
T_alpha=t_c.*therm_ex;
T_alpha=sum(T_alpha);
T2_alpha=t_c.^2;
T2_alpha=T2_alpha.*therm_ex;
T2_alpha=sum(T2_alpha);
rhs = [alpha;T_alpha;T2_alpha];
% solving for unknowns of regression 
ukn=inv(coef)*rhs;
a_01=ukn(1);
a_11=ukn(2);
a_21=ukn(3);
disp(' The 2nd order polynomial regression model is: ')
fprintf('  alpha=%g + %g*T + %g*T^2\n', a_01,a_11,a_21)
%% Part 9 Use Polyfit Command
disp(' ')
disp('Part 9')
tt=polyfit(t_c,therm_ex,2);
a_0=tt(1);
a_1=tt(2);
a_2=tt(3);
disp(' The polynomial of temperature is: ')
fprintf('  alpha=%g + %g*T + %g*T^2\n', a_0,a_1,a_2)
%% Part 10 Estimate Change of Diameter of Aluminum Cylinder 
disp(' ')
disp('Part 10')
% Placed in ice water for several hours would hit theta_a 
% T= theta_a 
T=theta_a-273.15; 
T_i=T_C(1)-273.15;
T1=T-T_i;
% Find coefficient of thermal expansion 
alpha_ex = (ukn(1)) + (ukn(2)*T1) + (ukn(3)*(T1^2));
alpha_ex=subs(alpha_ex);
alpha_ex=double(alpha_ex);
% Calculate Diameter Change Using Volume Expansion 
delta_V= V*3*alpha_ex*(T-T_i);
syms R_1
eqn_a='-delta_V==(pi)*(R_1^2)*L';
eqn_b=subs(eqn_a);
eqn_c=solve(eqn_b,R_1);
soln=double(eqn_c)*2;
diam=soln*2;
fprintf('The change in diameter would be %g m\n',diam(1))
%% Part 11 Estimate dT/dt 
disp(' ')
disp('Part 11')
% Using central divided difference at 4th point 
a3=t_s(3);
a5=t_s(5);
delta_t=(((t_s(5))-(t_s(3)))/2);
theta_t1=((theta_a) + (T_C(1)-theta_a)*exp((-fin_soln*A*a3)/(m*C)));
theta_t2=((theta_a) + (T_C(1)-theta_a)*exp((-fin_soln*A*a5)/(m*C)));
val1=theta_t1;
val2=theta_t2;
dT_dt=(theta_t2-theta_t1)/(2*delta_t);
disp('The rate of change of temperature with respect to time at 4th point')
disp(dT_dt)
disp('degrees Kelvin per second')
%% Part 12 Estimate Rate of Change of Heat Stored 
disp('Part 11')
sym theta_dt;
func='dT_dt==(m*C*theta_dt)';
func1=solve(func);
rate=subs(func1);
disp(' The rate of change of heat stored in the cylinder:')
disp(rate)
%% Part 13 Estimate Rate of Change of Heat Lost due to Convection 
disp('Part 13')
syms rate_c;
rate_change=exp((-fin_soln*A*t_s(4))/(m*C));
disp(' The rate of change of heat lost due to convection')
disp(rate_change)




























 





















