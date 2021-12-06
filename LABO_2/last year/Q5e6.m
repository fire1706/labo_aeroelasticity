%%%% MAIN - Aeroelasticity Lab 2 - Vortex Induced Vibration %%%%
clearvars
close all


%%%%%%%%%%%%%%%%%%%%%%%ù
% The code used are made with the same method than Elena for the natural
% frequencies and damping and Maf for the frequency matching
% The unique interesting method is the Skop-Griffin graph.

% Other things are a huge hardcode, sorry


% DATA 

graph_displ = 0;

load('DATAg1.mat');
datag1 = data;          % Data of our experiments
load('DATAe_1.mat');
datae1 = data;          % Data with elastomer 1
load('DATAe_2.mat');
datae2 = data;          % Data with elastomer 2
fs_acc = 201.03; % [Hz]. Constant sampling frequency of horizontal acceleration
dt_acc = 1/fs_acc; %[s] Sampling time 
fs_acc = 201.03; % [Hz]. Constant sampling frequency 
fs_vel = 250; % [Hz]. Constant sampling frequency of flow velocity
dt_vel = 1/fs_vel; %[s] Sampling time of flow velocity

D = 0.1;            % [m] external diameter of the cylinder

Stiffg1 = 6155;     % [N/m] The stiffness in flexion in the horizontal direction

Ug1 = zeros(length(datag1),1);
Ue1 = zeros(length(datae1),1);
Ue2 = zeros(length(datae2),1);
% ydd = zeros(length(data),1);
% w_velocity = zeros(length(data),1);

omega = 7.1237; % Natural frequency found by Elena for OUR Data

% Airspeed you want to analyze

%%%%%%%%%%   Question 6   %%%%%%%%%%%%
% Without elastomer
for i = 1: length(datag1)
% data = 12

    Ug1(i) = datag1(i).U; % [m/s] tested airspeed
    ydd= datag1(i).yddot; % [m/s^2] horizontal acceleration matrix for each airspeed
    wg1= datag1(i).w; %[g] time variation of horizontal component of velocity, in the wake of the cylinder
    A(i) = (max(abs(ydd))*9.81)/(4*pi^2*omega^2);
    if i == 1
        ydd_0 = ydd;
    end
    
    if i > 1
        k = i;
        % Motion
        L = length(ydd);
        fft_ydd = fft(ydd);
        P2 = abs(fft_ydd/L);     % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1); % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                                 % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        f = (fs_acc * (0:floor(L/2))/L)';
        [frequency_Peak,frequency_Peak_index] = findpeaks(P1(2:end)); % Return local maxima and corresponding index 
        [max_Peak,frequency_max_index] = max(frequency_Peak);
        frequency_max_index = frequency_Peak_index(frequency_max_index);
        frequency_motiong1(k) = f(frequency_max_index);
        pulsation_motiong1(k) = frequency_motiong1(k)*2*pi;
        
        clearvars f P1 P2 L frequency_Peak frequency_Peak_index max_Peak frequency_max_index
        
        % Wake
        L = length(wg1);
        fft_w = fft(wg1);
        P2 = abs(fft_w/L);       % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1); % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                                 % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        f = (fs_vel * (0:floor(L/2))/L)';
        [frequency_Peak,frequency_Peak_index] = findpeaks(P1(2:end)); % Return local maxima and corresponding index 
        [max_Peak,frequency_max_index] = max(frequency_Peak);
        frequency_max_index = frequency_Peak_index(frequency_max_index);
        frequency_wakeg1(k) = f(frequency_max_index);
        pulsation_wakeg1(k) = frequency_wakeg1(k)*2*pi;
    end

end


% With elastomer 1
for i = 1 : 16 %length(datae1)

    Ue1(i) = datae1(i).U;
    ydde1= datae1(i).yddot;
    w_vel_e1= datae1(i).w;
    Time_ydde1 = ((0 : length(ydde1)-1) * dt_acc)';
    if i == 1
        ydd_0e1 = ydde1;
        
        % Frequency domain
        fft_ydd = fft(ydde1);  % Fourrier Transform of the signal  

        L = length(Time_ydde1); % Length of the signal 
        P2 = abs(fft_ydd/L); % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1);  % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                                  % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        freq = (fs_acc * (0:floor(L/2))/L)'; % Frequency domain


        % -- Computation of the natural frequency -- %
        % 1st method
        [P_max,P_max] = max(P1); % P_max is the maximum value of the FRF [dB].
                                      % ind_max gives its correponding index.
        w_e1 = freq(P_max); 
        
    end

    if i > 1
        k = i;
        % Motion
        L = length(ydde1);
        fft_ydd = fft(ydde1);
        P2 = abs(fft_ydd/L);     % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1); % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                                 % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        f = (fs_acc * (0:floor(L/2))/L)';
        [frequency_Peak,frequency_Peak_index] = findpeaks(P1(2:end)); % Return local maxima and corresponding index 
        [max_Peak,frequency_max_index] = max(frequency_Peak);
        frequency_max_index = frequency_Peak_index(frequency_max_index);
        frequency_motione1(k) = f(frequency_max_index);
        pulsation_motione1(k) = frequency_motione1(k)*2*pi;
        
        clearvars f P1 P2 L frequency_Peak frequency_Peak_index max_Peak frequency_max_index
        
        % Wake
        L = length(w_vel_e1);
        fft_w = fft(w_vel_e1);
        P2 = abs(fft_w/L);       % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1); % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                                 % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        f = (fs_vel * (0:floor(L/2))/L)';
        [frequency_Peak,frequency_Peak_index] = findpeaks(P1(2:end)); % Return local maxima and corresponding index 
        [max_Peak,frequency_max_index] = max(frequency_Peak);
        frequency_max_index = frequency_Peak_index(frequency_max_index);
        frequency_wakee1(k) = f(frequency_max_index);
        pulsation_wakee1(k) = frequency_wakee1(k)*2*pi;
    end
    
    Ae1(i) = (max(abs(ydde1))*9.81)/(4*pi^2*w_e1^2);

end


% With elastomer 2
for i = 1 : 11 %length(datae2)
    
    Ue2(i) = datae2(i).U;
    ydde2= datae2(i).yddot;
    w_vel_e2= datae2(i).w;
    Time_ydde2 = ((0 : length(ydde2)-1) * dt_acc)';
    
    if i == 1
        ydd_0e2 = ydde2;
        
        % Elastomer 2
        % Frequency domain
        fft_ydd = fft(ydde2);  % Fourrier Transform of the signal  

        L = length(Time_ydde2); % Length of the signal 
        P2 = abs(fft_ydd/L); % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1);  % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                                  % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        freq = (fs_acc * (0:floor(L/2))/L)'; % Frequency domain


        % -- Computation of the natural frequency -- %
        % 1st method
        [P_max,P_max] = max(P1); % P_max is the maximum value of the FRF [dB].
                                      % ind_max gives its correponding index.
        w_e2 = freq(P_max); 
    end
    
    if i > 1
        k = i;
        % Motion
        L = length(ydde2);
        fft_ydd = fft(ydde2);
        P2 = abs(fft_ydd/L);     % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1); % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                                 % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        f = (fs_acc * (0:floor(L/2))/L)';
        [frequency_Peak,frequency_Peak_index] = findpeaks(P1(2:end)); % Return local maxima and corresponding index 
        [max_Peak,frequency_max_index] = max(frequency_Peak);
        frequency_max_index = frequency_Peak_index(frequency_max_index);
        frequency_motione2(k) = f(frequency_max_index);
        pulsation_motione2(k) = frequency_motione2(k)*2*pi;
        
        clearvars f P1 P2 L frequency_Peak frequency_Peak_index max_Peak frequency_max_index
        
        % Wake
        L = length(w_vel_e2);
        fft_w = fft(w_vel_e2);
        P2 = abs(fft_w/L);       % Compute the two-sided spectrum P2
        P1 = P2(1:floor(L/2)+1); % Compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
                                 % floor = round toward negative infinite
        P1(2:end-1) = 2*P1(2:end-1);
        f = (fs_vel * (0:floor(L/2))/L)';
        [frequency_Peak,frequency_Peak_index] = findpeaks(P1(2:end)); % Return local maxima and corresponding index 
        [max_Peak,frequency_max_index] = max(frequency_Peak);
        frequency_max_index = frequency_Peak_index(frequency_max_index);
        frequency_wakee2(k) = f(frequency_max_index);
        pulsation_wakee2(k) = frequency_wakee2(k)*2*pi;
    end
    
    Ae2(i) = (max(abs(ydde2))*9.81)/(4*pi^2*w_e2^2);

end


%% Find the damping of the new elastomers

% Elastomer initial
[ydd_peaks, ydd_peaks_index] = findpeaks(ydd_0);  % Return local maxima and corresponding index 
peak_wrong = find(ydd_peaks<0);
ydd_peaks(peak_wrong) = [];

n_peaks = length(ydd_peaks); % No. of peaks of interest

damp = zeros(n_peaks,1);
for i = 1:n_peaks-1
    damp(i) = (log(ydd_peaks(i)/ydd_peaks(i+1)))/(2*pi);   
end
damp_mean = mean(damp);
zeta_logdec = 1/sqrt(1+((2*pi/(damp_mean))^2)); %Assesses Damping Constant

% Elastomer 1
[ydd_peaks, ydd_peaks_index] = findpeaks(ydd_0e1);  % Return local maxima and corresponding index 
peak_wrong = find(ydd_peaks<=0);
ydd_peaks(peak_wrong) = [];

n_peaks = length(ydd_peaks); % No. of peaks of interest

damp = zeros(n_peaks,1);
for i = 1:n_peaks-1
    damp(i) = (log(ydd_peaks(i)/ydd_peaks(i+1)))/(2*pi);  
end
damp(isnan(damp)) = [];
damp_mean1 = mean(damp);
zeta_logdec1 = 1/sqrt(1+((2*pi/(damp_mean1))^2)); %Assesses Damping Constant

% Elastomer 2
[ydd_peaks, ydd_peaks_index] = findpeaks(ydd_0e2);  % Return local maxima and corresponding index 
peak_wrong = find(ydd_peaks<=0);
ydd_peaks(peak_wrong) = [];

n_peaks = length(ydd_peaks); % No. of peaks of interest

damp = zeros(n_peaks,1);
for i = 1:n_peaks-1
    damp(i) = (log(ydd_peaks(i)/ydd_peaks(i+1)))/(2*pi); 
end
damp(isnan(damp)) = [];
damp_mean2 = mean(damp);
zeta_logdec2 = 1/sqrt(1+((2*pi/(damp_mean2))^2)); %Assesses Damping Constant

%% Plots

Ue1 = Ue1(1:16);
Ue2 = Ue2(1:11);
Ue1 = Ue1 - 0.4;        % Blockage effect of 0.4 m/s
Ue2 = Ue2 - 0.4;
A(1) = 0; Ae1(1) = 0; Ae2(1) = 0;
[Ug1,sortUg1] = sort(Ug1);      % To have an order in the velocities and amplitude
A = A(sortUg1);
[Ue1,sortUe1] = sort(Ue1);
Ae1 = Ae1(sortUe1);
[Ue2,sortUe2] = sort(Ue2);
Ae2 = Ae2(sortUe2);

damp = figure('name','Effect of the damping on VIV');
hold on;
plot(Ug1/(omega*D),A/D,'-o','linewidth',2);
plot(Ue1/(w_e1*D),Ae1/D,'-^','linewidth',2);
plot(Ue2/(w_e2*D),Ae2/D,'-s','linewidth',2);
grid on;
grid minor;
xlabel('Reduced velocity $U_r$ [-]','interpreter','latex');
ylabel('$A/D$ [-]','interpreter','latex');
set(gca,'fontsize',16);
leg = legend('show',{'$8.34\cdot 10^{-4}$','$2.7\cdot 10^{-3}$','$3.4 \cdot 10^{-3}$'},'fontsize',16,'Interpreter','latex','location','northwest');
title(leg,'$\zeta$ [-]','interpreter','latex')

%% Skop-Griffin
% SG = 4*pi^2*St^2*Sc
% Sc = pi/2 * (1+m_r)*eta_s
% St = 0.2 (cylinder)
% m_r = m_s/m_f = m_s/(rho*pi*D^2/4) = 3.1/(1.225*pi*0.1^2/4)
% eta = damping factor = E_dissipated_per_cycle/(4*pi*E_struct) = damp_mean

m_s = Stiffg1/((omega*2*pi)^2*1.3);
m_r = m_s*4/(1.225*pi*D^2);
St = 0.2;

Scg1 = pi/2 * (1+m_r)*damp_mean;
Sce1 = pi/2 * (1+m_r)*damp_mean1;
Sce2 = pi/2 * (1+m_r)*damp_mean2;

SGg1 = 4*pi^2*St^2*Scg1;
SGe1 = 4*pi^2*St^2*Sce1;
SGe2 = 4*pi^2*St^2*Sce2;

SG = [SGg1 SGe1 SGe2];

A_SG = [max(A) max(Ae1) max(Ae2)]/D;

SG_figure = figure('name','Effect of SG on VIV');
semilogx(SG,A_SG,'-o','linewidth',2);
xlabel('SG','interpreter','latex');
ylabel('$A/D$ [-]','interpreter','latex');
set(gca,'fontsize',16);
grid on;
grid minor;
xlim([0.1 10]);


%% Frequency vs. airspeed

figure
h = plot(1:10,1:10,1:10,2:11,1:5,1:5);
c = get(h,'Color');


Str_law = (St/D) * Ue1(2:16);

figure('name','Variation of the motion and wake frequencies')
hold on
p1 = plot(Ug1(2:11)/(omega*D),frequency_wakeg1(2:11),'x','color',c{1},'MarkerSize',12,'linewidth',2);
p5 = plot(Ug1(2:11)/(omega*D),frequency_wakeg1(2:11),'x','color',c{1},'MarkerSize',12,'linewidth',2);
p2 = plot(Ug1(2:11)/(omega*D),frequency_motiong1(2:11),'o','color',c{1},'MarkerSize',12,'linewidth',2);
p3 = plot(Ue1(2:16)/(w_e1*D),frequency_wakee1(2:16),'x','color',c{2},'MarkerSize',12,'linewidth',2);
p4 = plot(Ue2(2:11)/(w_e2*D),frequency_wakee2(2:11),'x','color',c{3},'MarkerSize',12,'linewidth',2);


plot(Ue1(2:16)/(w_e1*D),frequency_motione1(2:16),'o','color',c{2},'MarkerSize',12,'linewidth',2)
plot(Ue2(2:11)/(w_e2*D),frequency_motione2(2:11),'o','color',c{3},'MarkerSize',12,'linewidth',2)
plot(Ue1(2:16)/(w_e1*D),Str_law(:),'--','color','k','linewidth',1.5)
plot(Ug1(2:11)/(omega*D),linspace(omega,omega,length(Ug1)-1),'--','color','k','linewidth',1.5)
% Arrow with two head at both end and text between
%an = annotation('doublearrow',[0.4825 0.775],[0.35 0.35],'linewidth',1.5);
plot([Ug1(4)/(omega*D) Ug1(4)/(omega*D)],[0 max(frequency_wakeg1)],':','color','k','linewidth',1.5)
plot([Ug1(10)/(omega*D) Ug1(10)/(omega*D)],[0 max(frequency_wakeg1)],':','color','k','linewidth',1.5)
xlabel('Reduced velocity [-]','FontSize', 12, 'Interpreter', 'latex');
ylabel('Frequency [Hz]','FontSize', 12, 'Interpreter', 'latex');
grid on;
grid minor; 
a=axes('position',get(gca,'position'),'visible','off');
lgd = legend([p1,p2],'Wake frequency $f_{VS}$','Motion frequency $f_{S}$','Location','northwest');
set(lgd, 'Interpreter', 'latex', 'FontSize', 14)
leg = legend(a,[p5,p3,p4],'$8.34\cdot 10^{-4}$','$2.7\cdot 10^{-3}$','$3.4 \cdot 10^{-3}$','interpreter','latex','location','southeast');
title(leg,'$\zeta$ [-]','interpreter','latex')
set(leg, 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'TickLabelInterpreter','latex','Fontsize',16)








