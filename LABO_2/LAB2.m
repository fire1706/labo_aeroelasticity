clear all
close all
clc

% %Q1
 load('DATAg2.mat');
 load('DATAg1.mat');
 Freq_sampling = 1000; %Hz
 T_sampling = 1/Freq_sampling; %s
 D = 0.125; %[m] Cylinder diameter
 Str = 0.2; %[-] Strouhal numbe
 
 
 
 U = zeros(length(DATAg2),1);
 y = zeros(length(DATAg2),1);
 w=zeros(length(DATAg2),1);
 puls_motion = zeros(length(DATAg2)-1,1);
 puls_wake = zeros(length(DATAg2)-1,1);
 freq_motion = zeros(length(DATAg2)-1,1);
 freq_wake = zeros(length(DATAg2)-1,1);

for k=1:length(DATAg2)

%% Question 1    
    if k==1 % determination of the modal properties

        %U(k) = DATAg2(k).U ; 
        y = DATAg2(k).y(14540:18500);
        figure('name','y')
        plot(y)
        %w= DATAg2(k).w ; 

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

        [~,P_max] = max(P1); % P_max is the maximum value of the FRF [dB].
                              % ind_max gives its correponding index.
        fmax = freq(P_max)
        freq_motion(k) = fmax;

        % Plot the signal in the freq. domain
        figure('name','FRF of the acceleration in frequency domain')
        plot(freq,P1)
        xlabel('f [Hz]','FontSize', 18, 'Interpreter', 'latex')
        ylabel('FRF','FontSize', 18, 'Interpreter', 'latex')
        %xlim([4.1,5.2]);
        grid on
        grid minor
      
    else
        U(k) = DATAg2(k).U ; % [m/s] tested airspeed
        y = DATAg2(k).y(5000:15000);  % [m/s^2] horizontal acceleration matrix for each airspeed
        w= DATAg2(k).w(5000:15000) ; %[g] time variation of horizontal component of velocity, in the wake of the cylinder

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

        [~,P_max] = max(P1); % P_max is the maximum value of the FRF [dB].
                           % ind_max gives its correponding index.
        
        freq_motion(k) = freq(P_max);
        puls_motion(k) = freq_motion(k)*2*pi;
        
        

        clear f P1 P2 L 

        %U(k) = DATAg2(k).U ; % [m/s] tested airspeed
        %y = DATAg2(k).y;  % [m/s^2] horizontal acceleration matrix for each airspeed
        %w= DATAg2(k).w ; %[g] time variation of horizontal component of velocity, in the wake of the cylinder

        % Creation of time vector for each data
        time_w = ((0:length(w)-1)*T_sampling)';
    
        % Frequency domain
        fft_w = fft(w);  % Fourrier Transform of the signal  

        L = length(time_w); % Length of the signal 
        P2 = abs(fft_w/L); % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1);  % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                          % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        freq = (Freq_sampling * (0:floor(L/2))/L)'; % Frequency domain

        [~,P_max] = max(P1); % P_max is the maximum value of the FRF [dB].
                              % ind_max gives its correponding index.
        
        freq_wake(k) = freq(P_max);
        puls_wake(k) = freq_wake(k)*2*pi;

    end

 % Q4

    acceleration = y * 9.81; 
    amplitude = - acceleration/(4*pi^2*freq_motion(k)^2);   
    A_max(k) = max(abs(amplitude))/D; % Maximum amplitude /D


end 

%% Q2 mass computation
w = 2*pi*fmax;
N = 40800;
w = w*w;
Mass =  N/w;

%% Q4

[U, sortU] = sort(U); 
A_max = A_max(sortU);
Ampl_rms = A_max/sqrt(2);

figure('name','Maximum amplitude as function of airspeed')
plot(U(2:end),A_max(2:end), '-o', 'linewidth',1.5)
hold on
plot(U(2:end),Ampl_rms(2:end), '-o', 'linewidth',1.5)
xlabel('$U_r$ [m/s]','FontSize', 12, 'Interpreter', 'latex');
ylabel('$\frac{A_{max}}{D}$ ','FontSize', 12, 'Interpreter', 'latex');
lgd = legend('Max amplitude','RMS amplitude', 'location', 'northwest');
set(lgd, 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'TickLabelInterpreter','latex','Fontsize',16)
grid on
grid minor


%% Q3
%  load('DATAe_3.mat');
%  
%  Freq_sampling = 1000; %Hz
%  T_sampling = 1/Freq_sampling; %s
%  D = 0.125; %[m] Cylinder diameter
%  Str = 0.2; %[-] Strouhal number
%  U = zeros(length(DATAe_3),1);
%  puls = zeros(length(DATAe_3)-1,1);
%  y = zeros(length(DATAe_3),1);
% 
% for k=1:length(DATAe_3)
% 
%     
%     if k==1 % determination of the modal properties
% 
%         %U(k) =DATAe_3(k).U ; % [m/s] tested airspeed
%         y = DATAe_3(k).y;  % [m/s^2] horizontal acceleration matrix for each airspeed
%         %puls= DATAe_3(k).w ; %[g] time variation of horizontal component of velocity, in the wake of the cylinder
% 
%         % Creation of time vector for each data
%         time_y = ((0:length(y)-1)*T_sampling)';
%     
%         % Frequency domain
%         fft_y = fft(y);  % Fourrier Transform of the signal  
% 
%         L = length(time_y); % Length of the signal 
%         P2 = abs(fft_y/L); % Compute the two-sided spectrum P2
%         P1 = P2(1:floor(L/2)+1);  % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
%                           % floor = round toward negative infinite
%         P1(2:end-1) = 2*P1(2:end-1);
%         freq = (Freq_sampling * (0:floor(L/2))/L)'; % Frequency domain
% 
%         [~,P_max] = max(P1); % P_max is the maximum value of the FRF [dB].
%                               % ind_max gives its correponding index.
%         fmax = freq(P_max)
% 
%         % Plot the signal in the freq. domain
%         figure('name','FRF of the acceleration in frequency domain')
%         plot(freq,P1)
%         xlabel('f [Hz]','FontSize', 18, 'Interpreter', 'latex')
%         ylabel('FRF','FontSize', 18, 'Interpreter', 'latex')
%         grid on
%         grid minor
%        
%  
%     
%         
%     end

%end

%% Q5

Str_law = (Str/D) * U;

figure('name','Q4: Variation of the motion and wake frequencies')
plot(U(2:end),freq_wake(2:end),'o', 'linewidth',1.5)
hold on
plot(U(2:end),freq_motion(2:end),'+', 'linewidth',1.5)
plot(U,Str_law,'-o', 'linewidth',1.5)
plot(U,linspace(fmax,fmax,length(U)),'-','color','k','linewidth',1.5)

xlabel('Flow velocity [m/s]','FontSize', 12, 'Interpreter', 'latex');
ylabel('Frequency [Hz]','FontSize', 12, 'Interpreter', 'latex');
lgd = legend('Wake frequency $f_{VS}$','Motion frequency $f_{S}$', 'Strouhal law','Location','northwest');
set(lgd, 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'TickLabelInterpreter','latex','Fontsize',16)
grid on
grid minor


%% Q6
% the critical velocity is the one when VIV appear, it will be computed by
% the fs*D/Str 
freq_viv = mean(freq_motion(6:9));
U_cr = freq_viv*D/Str
U_r_exp = U_cr/D/fmax
U_r_lockin = 1/Str

% on remarque que lees deux Ur sont fort proche , it is good

%% Q7 Effect of the damping on VIV

 dampingf
 
 %% Q8
 
 





