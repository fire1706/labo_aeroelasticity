clear all
close all
clc
format long
% %Q1
 load('DATAlab3g2.mat');

 Freq_sampling = 201.03; %Hz
 T_sampling = 1/Freq_sampling; %s

 St=0.125; % square section

 s_c = 35*10^(-3); %[m] external side cylinder
 t_cylinder = 2*10^(-3); %[m] thickness cylinder
 l_cylinder = 1.3; %[m] length cylinder

 l_alu=1.06; % [m] length alu
 h_alu=50*10^(-3); %[m] heigth alu
 t_alu = 15*10^(-3); % [m] thickness alu


 
 U = zeros(length(DATA),1); % m.s^-1
 yddot = zeros(length(DATA),1); % g.s

 

 freq_motion = zeros(length(DATA)-1,1);
 
 
% Question 1 
% determintation of the modal properties 
for k=1:length(DATA)
    
     U(k) = DATA(k).U ; % [m/s] tested airspeed
     y = DATA(k).yddot;  % [g.s] Accelerations are expressed in g's.
      %plot(y)
    
    % Creation of time vector for each data
    time_y = ((0:length(y)-1)*T_sampling)';
    
    % Frequency domain
    fft_y = fft(y);  % Fourrier Transform of the signal  
    
    L = length(time_y); % Length of the signal 
    P2 = abs(fft_y/L); % Compute the two-sided spectrum P2
    P1 = P2(1:floor(L/2)+1);  % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                      % floor = round toward negative infinite
    P1(2:end-1) = 2*P1(2:end-1);
    freq = (Freq_sampling * (0:floor(L/2))/L)'; % Frequency domain
    
    [~,P_max] = max(P1(10:end)); % P_max is the maximum value of the FRF [dB].
                       % ind_max gives its correponding index.
    
    freq_motion(k) = freq(P_max);
    
    
    Amplitude(k) = max(abs(-((y)*9.81)/(4*pi^2*freq_motion(k)^2)));
    
end

f_s0=freq_motion(1)

clear y time_y fft_y L P2 P1 P_max
% determintation of the modal properties with elastomer 1
freq_motion1 = zeros(length(DATAadd1)-1,1);
for k=1:length(DATAadd1)
    
     U1(k) = DATAadd1(k).U ; % [m/s] tested airspeed

        y = DATAadd1(k).yddot;  % [g.s] Accelerations are expressed in g's.


    
    % Creation of time vector for each data
    time_y = ((0:length(y)-1)*T_sampling)';
    
    % Frequency domain
    fft_y = fft(y);  % Fourrier Transform of the signal  
    
    L = length(time_y); % Length of the signal 
    P2 = abs(fft_y/L); % Compute the two-sided spectrum P2
    P1 = P2(1:floor(L/2)+1);  % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                      % floor = round toward negative infinite
    P1(2:end-1) = 2*P1(2:end-1);
    freq = (Freq_sampling * (0:floor(L/2))/L)'; % Frequency domain
    
    [~,P_max] = max(P1(10:end)); % P_max is the maximum value of the FRF [dB].
                       % ind_max gives its correponding index.
    
    freq_motion1(k) = freq(P_max);
    
     Amplitude1(k) = max(abs(-(y*9.81)/(4*pi^2*freq_motion1(k)^2)));
end

f_s1=freq_motion1(1)

Amplitude2 = zeros(length(DATAadd2)-1,1);
clear y time_y fft_y L P2 P1 P_max
% determintation of the modal properties with elastomer 2
freq_motion2 = zeros(length(DATAadd2)-1,1);
for k=1:length(DATAadd2)
    
     U2(k) = DATAadd2(k).U ; % [m/s] tested airspeed
     y = DATAadd2(k).yddot;  % [g.s] Accelerations are expressed in g's.
      %plot(y)
    
    % Creation of time vector for each data
    time_y = ((0:length(y)-1)*T_sampling)';
    
    % Frequency domain
    fft_y = fft(y);  % Fourrier Transform of the signal  
    
    L = length(time_y); % Length of the signal 
    P2 = abs(fft_y/L); % Compute the two-sided spectrum P2
    P1 = P2(1:floor(L/2)+1);  % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                      % floor = round toward negative infinite
    P1(2:end-1) = 2*P1(2:end-1);
    freq = (Freq_sampling * (0:floor(L/2))/L)'; % Frequency domain
    %plot(freq,P1)
    [~,P_max] = max(P1(20:end-end/3*2)); % P_max is the maximum value of the FRF [dB].
                       % ind_max gives its correponding index.
    
    freq_motion2(k) = freq(P_max);
    
     Amplitude2(k) = max(abs(-(y*9.81)/(4*pi^2*freq_motion2(k)^2)));
end

f_s2=freq_motion2(1)


figure('name','Maximum amplitude as function of airspeed')
hold on
plot(U(2:end)/(f_s0*s_c),Amplitude(2:end)/s_c,'-o','linewidth',1.5)
plot(U1(2:end)/(f_s0*s_c),Amplitude1(2:end)/s_c,'-o','linewidth',1.5)
plot(U2(2:end)/(f_s0*s_c),Amplitude2(2:end)/s_c,'-o','linewidth',1.5)
% plot(U(2:end),Amplitude(2:end)/s_c,'-o','linewidth',1.5)
% plot(U1(2:end),Amplitude1(2:end)/s_c,'-o','linewidth',1.5)
% plot(U2(2:end),Amplitude2(2:end)/s_c,'-o','linewidth',1.5)
xlabel('$U_r$','FontSize',12,'Interpreter','latex');
ylabel('$A_{max}/D$','FontSize',12,'Interpreter','latex');
lgd = legend('$\xi$=0.0011','$\xi$=0.002','$\xi$=0.0018', 'location', 'northwest');
set(lgd, 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'TickLabelInterpreter','latex','Fontsize',16)
grid on
grid minor



%%damping calculation 
clear yg2 Y_Ug2  yadd2 Y_Uadd2 yadd1 Y_Uadd1 yddot
%% G2
%plot(DATA(1).yddot)
yg2 = DATA(1).yddot;
Y_Ug2 = decrepeak(yg2);
i = 1;
while i<length(Y_Ug2)
    dampg2(i) = 1/(2*pi)*log(Y_Ug2(i)/Y_Ug2(i+1));
    i = i+1;
end
Dampg2 = mean(dampg2)

%%damping calculation with elastomers 1
clear yadd1 Y_Uadd1 yadd2 Y_Uadd2 yddot yg2 Y_Ug2 

U_add = zeros(length(DATAadd1),1); % m.s^-1
yddot_add1 = zeros(length(DATAadd1),1); % g.s

%plot(DATAadd(1).yddot)
yadd1 = DATAadd1(1).yddot;
Y_Uadd1 = decrepeak(yadd1);
i = 1;
while i<length(Y_Uadd1)
    damp_add1(i) = 1/(2*pi)*log(Y_Uadd1(i)/Y_Uadd1(i+1));
    i = i+1;
end
Damp_add1 = mean(damp_add1)

%%damping calculation with elastomers 2
clear yadd2 Y_Uadd2 yadd Y_Uadd yddot yg2 Y_Ug2 

U_add2 = zeros(length(DATAadd2),1); % m.s^-1
yddot_add2 = zeros(length(DATAadd2),1); % g.s

%plot(DATAadd(1).yddot)
yadd2 = DATAadd2(1).yddot;
Y_Uadd2 = decrepeak(yadd2);
i = 1;
while i<length(Y_Uadd2)
    damp_add2(i) = 1/(2*pi)*log(Y_Uadd2(i)/Y_Uadd2(i+1));
    i = i+1;
end
Damp_add2 = mean(damp_add2)

%Question 2

%flexural stiffness
E_alu=69*10^3; %[MPa]

I=h_alu*t_alu^3/12 %  [m^4] inertia of the support

k=3*E_alu*10^6*I/l_alu^3; %[N/m] stiffness of one bar

k_eq= 2*k %[N/m] stiffness equivalent of 2 bars

% equivalent mass

w=f_s0*2*pi; % pulsation  
m=k_eq/w^2 % [kg] mass of the strcuture
m_l=m/l_cylinder %  [kg/m] mass per unit length of the structure

%question 3

%critical velocities VIV

U_VIV=f_s0*s_c/St % strouhal law for our data
U_reduced =1/St % reduced velocity for a square section

%critical velocities galloping

rho=1.225; % [kg/m^3]

%without elastomers
m_r=m/(rho*s_c^2)%((50e-3)^-(50e-3 - 1.5e-3)^2)*(1440e-3)*2710%m_l/(rho*s_c^2) %mass ratio
n=1/(2*m_r)
A=3; %aerodynamic characteristics of a smooth flow


U_glp=2*Dampg2/(n*A) % [m/s] threshold velocity
U_glp_r=U_glp/f_s0/s_c 

% with elastomer 1

U_glp1=2*Damp_add1/(n*A) % [m/s] threshold velocity
U_glp1_r=U_glp1/f_s0/s_c

% with elastomer 2

U_glp2=2*Damp_add2/(n*A) % [m/s] threshold velocity
U_glp2_r=U_glp2/f_s0/s_c


%% Q4
B = Dampg2%/(2*pi)*10;
B1 = Damp_add1%/(2*pi)*10;
B2 = Damp_add2%/(2*pi)*10;
factor = n*A/(2*B);
factor1 = n*A/(2*B1);
factor2 = n*A/(2*B2);
[Uu,Yy] = importfile('UNI_curve.csv', 5, 138);

figure('name','Non dimensionalization')
hold on
plot(U(2:end)/(2*pi*f_s0*s_c)*factor,Amplitude(2:end)/s_c*factor*10,'-o','linewidth',1.5)
plot(U1(2:end)/(2*pi*f_s0*s_c)*factor1,Amplitude1(2:end)/s_c*factor1*10,'-o','linewidth',1.5)
plot(U2(2:end)/(2*pi*f_s0*s_c)*factor2,Amplitude2(2:end)/s_c*factor2*10,'-o','linewidth',1.5)
plot(Uu,Yy,'-o','linewidth',1.5)
% plot(U(2:end),Amplitude(2:end)/s_c,'-o','linewidth',1.5)
% plot(U1(2:end),Amplitude1(2:end)/s_c,'-o','linewidth',1.5)
% plot(U2(2:end),Amplitude2(2:end)/s_c,'-o','linewidth',1.5)
xlabel('$U_r$ * $\frac{n A}{2 \beta}$','FontSize',12,'Interpreter','latex');
ylabel('$A_{max}/D$ * $\frac{n A}{2 \beta}$','FontSize',12,'Interpreter','latex');
lgd = legend('$\xi$=0.0011','$\xi$=0.002','$\xi$=0.0018', 'location', 'northwest');
set(lgd, 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'TickLabelInterpreter','latex','Fontsize',16)
grid on
grid minor