%% script damping f

 Freq_sampling = 1000; %Hz
 T_sampling = 1/Freq_sampling; %s
 D = 0.125; %[m] Cylinder diameter
 Str = 0.2; %[-] Strouhal number
 
 %% DaTA 1

 load('DATAe_1.mat');
 

 U_1 = zeros(length(DATAe_1),1);

 y_1 = zeros(length(DATAe_1),1);

for k=1:length(DATAe_1)

    
    if k==1 % determination of the modal properties
        clear P1 P2 L freq
        
        y_1 = DATAe_1(k).y(25550:30000);  % [m/s^2] horizontal acceleration matrix for each airspeed
%         figure('name','ye1 [mm]')
%         plot(y_1)
        % Creation of time vector for each data
        time_y_1 = ((0:length(y_1)-1)*T_sampling)';
    
        % Frequency domain
        fft_y = fft(y_1);  % Fourrier Transform of the signal  

        L = length(time_y_1); % Length of the signal 
        P2 = abs(fft_y/L); % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1);  % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                          % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        freq = (Freq_sampling * (0:floor(L/2))/L)'; % Frequency domain

        [~,P_max] = max(P1); % P_max is the maximum value of the FRF [dB].
                              % ind_max gives its correponding index.
        freq_motion_1(k) = freq(P_max)
        
%                         figure('name','FRF e1 of the acceleration in frequency domain')
%         plot(freq,P1)
%         xlabel('f [Hz]','FontSize', 18, 'Interpreter', 'latex')
%         ylabel('FRF','FontSize', 18, 'Interpreter', 'latex')
%         %xlim([4.1,5.2]);
%         grid on
%         grid minor
  
    
    else
        clear P1 P2 L freq
        U_1(k) = DATAe_1(k).U ; % [m/s] tested airspeed
        y_1 = DATAe_1(k).y(5000:15000);  % [m/s^2] horizontal acceleration matrix for each airspeed
        w_1= DATAe_1(k).w(5000:15000) ; %[g] time variation of horizontal component of velocity, in the wake of the cylinder

        % Creation of time vector for each data
        time_y_1 = ((0:length(y_1)-1)*T_sampling)';
    
        % Frequency domain
        fft_y = fft(y_1);  % Fourrier Transform of the signal  

        L = length(time_y_1); % Length of the signal 
        P2 = abs(fft_y/L); % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1);  % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                          % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        freq = (Freq_sampling * (0:floor(L/2))/L)'; % Frequency domain

        [~,P_max] = max(P1); % P_max is the maximum value of the FRF [dB].
                           % ind_max gives its correponding index.
        
        freq_motion_1(k) = freq(P_max);
        



        clear P1 P2 L freq
        % Creation of time vector for each data
        time_w = ((0:length(w_1)-1)*T_sampling)';
    
        % Frequency domain
        fft_w = fft(w_1);  % Fourrier Transform of the signal  

        L = length(time_w); % Length of the signal 
        P2 = abs(fft_w/L); % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1);  % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                          % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        freq = (Freq_sampling * (0:floor(L/2))/L)'; % Frequency domain

        [~,P_max] = max(P1); % P_max is the maximum value of the FRF [dB].
                              % ind_max gives its correponding index.
        
        freq_wake_1(k) = freq(P_max);

    end
%     acceleration = y_1 * 9.81; 
%     amplitude = - acceleration/(4*pi^2*freq_motion_1(k)^2);  
    amplitude = y_1(2:end);%./(sin(2*pi*freq_motion_1(k)*time_y_1(2:end)))/1000;
    A_max_1(k) = max(abs(amplitude)); % Maximum amplitude /D

end

[U_1, sortU] = sort(U_1); 
A_max_1 = A_max_1(sortU);

%% data e2
load('DATAe_2.mat');
 

 U_2 = zeros(length(DATAe_2),1);

 y_2 = zeros(length(DATAe_2),1);

for k=1:length(DATAe_2)

    
    if k==1 % determination of the modal properties
        clear P1 P2 L freq
        
        y_2 = DATAe_2(k).y(21840:26000);  % [m/s^2] horizontal acceleration matrix for each airspeed
%         figure('name','ye2 [mm]')
%         plot(y_2)
        % Creation of time vector for each data
        time_y_2 = ((0:length(y_2)-1)*T_sampling)';
    
        % Frequency domain
        fft_y = fft(y_2);  % Fourrier Transform of the signal  

        L = length(time_y_2); % Length of the signal 
        P2 = abs(fft_y/L); % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1);  % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                          % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        freq = (Freq_sampling * (0:floor(L/2))/L)'; % Frequency domain

        [~,P_max] = max(P1); % P_max is the maximum value of the FRF [dB].
                              % ind_max gives its correponding index.
        freq_motion_2(k) = freq(P_max)
%         
%                 figure('name','FRF e2 of the acceleration in frequency domain')
%         plot(freq,P1)
%         xlabel('f [Hz]','FontSize', 18, 'Interpreter', 'latex')
%         ylabel('FRF','FontSize', 18, 'Interpreter', 'latex')
%         %xlim([4.1,5.2]);
%         grid on
%         grid minor
%   
    
    else
        clear P1 P2 L freq
        U_2(k) = DATAe_2(k).U ; % [m/s] tested airspeed
        y_2 = DATAe_2(k).y(5000:15000);  % [m/s^2] horizontal acceleration matrix for each airspeed
        w_2= DATAe_2(k).w(5000:15000) ; %[g] time variation of horizontal component of velocity, in the wake of the cylinder

        % Creation of time vector for each data
        time_y_2 = ((0:length(y_2)-1)*T_sampling)';
    
        % Frequency domain
        fft_y = fft(y_2);  % Fourrier Transform of the signal  

        L = length(time_y_2); % Length of the signal 
        P2 = abs(fft_y/L); % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1);  % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                          % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        freq = (Freq_sampling * (0:floor(L/2))/L)'; % Frequency domain

        [~,P_max] = max(P1); % P_max is the maximum value of the FRF [dB].
                           % ind_max gives its correponding index.
        
        freq_motion_2(k) = freq(P_max);

        clear P1 P2 L freq
        % Creation of time vector for each data
        time_w = ((0:length(w_2)-1)*T_sampling)';
    
        % Frequency domain
        fft_w = fft(w_2);  % Fourrier Transform of the signal  

        L = length(time_w); % Length of the signal 
        P2 = abs(fft_w/L); % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1);  % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                          % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        freq = (Freq_sampling * (0:floor(L/2))/L)'; % Frequency domain

        [~,P_max] = max(P1); % P_max is the maximum value of the FRF [dB].
                              % ind_max gives its correponding index.
        
        freq_wake_2(k) = freq(P_max);

    end
%     acceleration = y_2 * 9.81; 
%     amplitude = - acceleration/(4*pi^2*freq_motion_2(k)^2);   
    amplitude = y_2(2:end);%./(sin(2*pi*freq_motion_2(k)*time_y_2(2:end)))/1000;
    A_max_2(k) = max(abs(amplitude)); % Maximum amplitude /D

end

[U_2, sortU] = sort(U_2); 
A_max_2 = A_max_2(sortU);

%% data e3
load('DATAe_3.mat');
 

 U_3 = zeros(length(DATAe_3),1);

 y_3 = zeros(length(DATAe_3),1);

for k=1:length(DATAe_3)

    
    if k==1 % determination of the modal properties
        clear P1 P2 L freq
        
        y_3 = DATAe_3(k).y(18260:21000);  % [m/s^2] horizontal acceleration matrix for each airspeed
%         figure('name','ye3 [mm]')
%         plot(y_3)
        % Creation of time vector for each data
        time_y_3 = ((0:length(y_3)-1)*T_sampling)';
    
        % Frequency domain
        fft_y = fft(y_3);  % Fourrier Transform of the signal  

        L = length(time_y_3); % Length of the signal 
        P2 = abs(fft_y/L); % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1);  % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                          % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        freq = (Freq_sampling * (0:floor(L/2))/L)'; % Frequency domain

        [~,P_max] = max(P1); % P_max is the maximum value of the FRF [dB].
                              % ind_max gives its correponding index.
        freq_motion_3(k) = freq(P_max)
  
%         figure('name','FRF e3 of the acceleration in frequency domain')
%         plot(freq,P1)
%         xlabel('f [Hz]','FontSize', 18, 'Interpreter', 'latex')
%         ylabel('FRF','FontSize', 18, 'Interpreter', 'latex')
%         %xlim([4.1,5.2]);
%         grid on
%         grid minor
    
    else
        clear P1 P2 L freq
        U_3(k) = DATAe_3(k).U ; % [m/s] tested airspeed
        y_3 = DATAe_3(k).y(5000:15000);  % [m/s^2] horizontal acceleration matrix for each airspeed
        w_3= DATAe_3(k).w(5000:15000) ; %[g] time variation of horizontal component of velocity, in the wake of the cylinder

        % Creation of time vector for each data
        time_y_3 = ((0:length(y_3)-1)*T_sampling)';
    
        % Frequency domain
        fft_y = fft(y_3);  % Fourrier Transform of the signal  

        L = length(time_y_3); % Length of the signal 
        P2 = abs(fft_y/L); % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1);  % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                          % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        freq = (Freq_sampling * (0:floor(L/2))/L)'; % Frequency domain

        [~,P_max] = max(P1); % P_max is the maximum value of the FRF [dB].
                           % ind_max gives its correponding index.
        
        freq_motion_3(k) = freq(P_max);
  

        clear P1 P2 L freq
        % Creation of time vector for each data
        time_w = ((0:length(w_3)-1)*T_sampling)';
    
        % Frequency domain
        fft_w = fft(w_3);  % Fourrier Transform of the signal  

        L = length(time_w); % Length of the signal 
        P2 = abs(fft_w/L); % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1);  % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                          % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        freq = (Freq_sampling * (0:floor(L/2))/L)'; % Frequency domain

        [~,P_max] = max(P1); % P_max is the maximum value of the FRF [dB].
                              % ind_max gives its correponding index.
        
        freq_wake_3(k) = freq(P_max);

    end
    %acceleration = y_3 * 9.81; 
    %amplitude = - acceleration/(4*pi^2*freq_motion_3(k)^2);
    amplitude = y_3(2:end);%./(sin(2*pi*freq_motion_3(k)*time_y_3(2:end)))/1000;
    A_max_3(k) = max(abs(amplitude)); % Maximum amplitude /D

end

[U_3, sortU] = sort(U_3); 
A_max_3 = A_max_3(sortU);

% Compute the damping for each Datae_123

figure('name','Maximum amplitude as function of airspeed and damping')
hold on
plot(U_1(2:end)/D/freq_motion_1(1),A_max_1(2:end)/D, '-o', 'linewidth',1.5)
plot(U_2(2:end)/D/freq_motion_2(1),A_max_2(2:end)/D, '-o', 'linewidth',1.5)
plot(U_3(2:end)/D/freq_motion_3(1),A_max_3(2:end)/D, '-o', 'linewidth',1.5)
plot(U(2:end)/D/freq_motion(1),A_max(2:end)/D, '-o', 'linewidth',1.5)
xlabel('$U$/$f_S$D','FontSize', 12, 'Interpreter', 'latex');
ylabel('$A_{max}/D$ ','FontSize', 12, 'Interpreter', 'latex');
lgd = legend('$\xi$=0.010127','$\xi$=0.011374','$\xi$=0.017113','$\xi$=0.008474');
set(lgd, 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'TickLabelInterpreter','latex','Fontsize',16)
grid on
grid minor


Str_law = (Str/D) * U_3;

figure('name',' Variation of the motion and wake frequencies depending on damping ratio')
hold on
plot(U_1(2:end),freq_wake_1(2:end),'o', 'linewidth',1.5)
plot(U_1(4:end),freq_motion_1(4:end),'+', 'linewidth',1.5)

plot(U_2(2:end),freq_wake_2(2:end),'o', 'linewidth',1.5)
plot(U_2(3:end),freq_motion_2(3:end),'+', 'linewidth',1.5)

plot(U_3(2:end),freq_wake_3(2:end),'o', 'linewidth',1.5)
plot(U_3(3:end),freq_motion_3(3:end),'+', 'linewidth',1.5)

plot(U_3,Str_law,'-', 'linewidth',1.5)
plot(U_3,linspace(fmax,fmax,length(U_3)),'-','color','k','linewidth',1.5)
xlabel('Flow velocity [m/s]','FontSize', 12, 'Interpreter', 'latex');
ylabel('Frequency [Hz]','FontSize', 12, 'Interpreter', 'latex');
lgd = legend('$\xi$=0.010127 $f_{VS}$','$\xi$=0.010127 $f_S$','$\xi$=0.011374 $f_{VS}$','$\xi$=0.011374 $f_{S}$','$\xi$=0.017113 $f_{VS}$','$\xi$=0.017113 $f_{S}$');
set(lgd, 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'TickLabelInterpreter','latex','Fontsize',16)
grid on
grid minor
