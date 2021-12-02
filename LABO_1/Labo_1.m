%% This code will initialize the data from the labo 1.
%clc;
clear all;

Freq_sampling = 201.03; %Hz
T_sampling = 1/Freq_sampling; %s

%% Plot of all the data
pitch = zeros(7,11000);
plunge = zeros(7,11000);
x = 1:11000;

DATA = load('DATAG2.mat');
for i = 1:7
    airspeed(i) = DATA.exp_data_G2(i).airspeed;
    for j=1:length(DATA.exp_data_G2(i).pitch)
        pitch(i,j) = DATA.exp_data_G2(i).pitch(j);
        plunge(i,j) = DATA.exp_data_G2(i).plunge(j);
    end
    
%     figure('name',['Figure for airspeed equal to ', num2str(airspeed(i)),' ']);
%     hold on
%     xlabel('Time of sampling')
%     ylabel('Acceleration in m/s^2')
%     plot(x*T_sampling,pitch(i,:),x*T_sampling,plunge(i,:));
%     %plot(x,pitch(i,:),x,plunge(i,:));
%     legend('pitch','plunge');
%     grid on
end

%% Initialization of the data we want from all the peak

% The index corresponding to all the peak 
% there is a general first index variable and a bonus one if the first one
% is not enough

% % % % number on the first collumn are the start value whereas the one on second
% collumn are the end value

index = [5820,2950,6000,3780,2450,4650,750;...
         6350,3500,6400,4300,2900,5500,5280];

index_bonus = [6980,4250,7150,5550,3150,5900,750;...
               7500,4800,7600,6000,3600,6800,5280];

%% Computation of the frequency and damping

% Using the half point method

% compute of the frequency range : 


for i = 1:7
    windowidth = (index(2,i)-index(1,i))*T_sampling;
    t = 0:T_sampling:windowidth;
    
    smooth_plunge =plunge(i,index_bonus(1,i):index_bonus(2,i));%sgolayfilt(plunge(i,index_bonus(1,i):index_bonus(2,i)),9,27);
    smooth_pitch =pitch(i,index_bonus(1,i):index_bonus(2,i));%sgolayfilt(pitch(i,index_bonus(1,i):index_bonus(2,i)),9,27);
    
    fft_plunge = fft(smooth_plunge);
    fft_pitch = fft(smooth_pitch);
    
    L = length(smooth_plunge);
    P2_plunge = abs(fft_plunge/L);
    P1_plunge = P2_plunge(1:L/2+1);
    P1_plunge(2:end-1) = 2*P1_plunge(2:end-1);
    f = Freq_sampling*(0:(L/2))/L;
    
    P2_pitch = abs(fft_pitch/L);
    P1_pitch = P2_pitch(1:L/2+1);
    P1_pitch(2:end-1) = 2*P1_pitch(2:end-1);
     
    %freq_range = linspace(1*2*pi,201.03*2*pi,length(fft_plunge));
    
    [a,dplunge(:,:)] = rfp(P1_plunge,f*2*pi,4);
    [b,dpitch(:,:)] = rfp(P1_pitch,f*2*pi,4);
%      [Max_plunge,ii] = max(P1_plunge);
%      val_int_plunge = Max_plunge/sqrt(2)
%      wn = f(ii)
%      plot(f,P1_plunge);
%      grid on
%      [Max_pitch,ii] = max(P1_pitch);
%      val_int_pitch = Max_pitch/sqrt(2)
%      wn = f(ii)
%      plot(f,P1_pitch);
%      grid on
%     j =1;
%     for t=1:length(P1_plunge)
%         if abs(val_int_plunge-P1_p lunge(t))<=0.00001
%             w(j) = f(t);
%             j = j+1;
%         end
%         if j==3
%             break;
%         end
%     end
%     damp_plunge(i) = (w(2)-w(1))/(2*wn);
    
    
    
    freq_plunge(i,:) = dplunge(1,:)/2/pi% in Hz
    damp_plunge(i,:) = dplunge(2,:)% the two other quantity are seems to be useless

    freq_pitch(i,:) = dpitch(:,1)/2/pi% in Hz
    damp_pitch(i,:) = dpitch(:,2)
end
Damp_HP_plunge = [0.05604,0.0871,0.1172,0.05961,0.0561,0.05104,0.00557];
Damp_HP_pitch = [0.06132,0.0779,0.08229,0.0933,0.0598,0.0629,0.0106];
Freq_HP_plunge = [3.7859,3.6485,4.0106,3.8585,4.0117,4.0159,4.2149];
Freq_HP_pitch = [5.6788,5.8375,4.0106,3.86    ,4.017,4.016,4.1706];

%plot the damp and freq
figure('name','frequency')
hold on
xlabel('Velocity')
ylabel('Frequency')
% plot(airspeed,Freq_HP_plunge(:))
% plot(airspeed,Freq_HP_pitch(:))
plot(airspeed,freq_plunge(:,1))
plot(airspeed,freq_pitch(:,1))
grid on
legend('plunge','pitch')


figure('name','damping')
hold on
xlabel('Velocity')
ylabel('Damping')
% plot(airspeed,Damp_HP_plunge(:))
% plot(airspeed,Damp_HP_pitch(:))
plot(airspeed,damp_plunge(:,1))
plot(airspeed,damp_pitch(:,1))
grid on
legend('plunge','pitch')


%% ---------------------------
%Partie Raph
%%-------------

%% Stability criterion part

% Stability criterion of lab data -----------------------------------------
S_pitch = zeros(7,1);
S_plunge = zeros(7,1);
for i = 1:7
    y_pitch = pitch(i,index_bonus(1,i):index_bonus(2,i)); % [rad/s?] Pitch signal
    y_plunge = plunge(i,index_bonus(1,i):index_bonus(2,i)); % [m/s?] Plunge signal

    Y_pitch = fft(y_pitch); % Y(omega) FFT for pitch
    Y_plunge = fft(y_plunge); % Y(omega) FFT for plunge
    yh_pitch = ifft(imag(Y_pitch) - 1i*real(Y_pitch)); % yh(t) Hilbert Transform for pitch
    yh_plunge = ifft(imag(Y_plunge) - 1i*real(Y_plunge)); % yh(t) Hilbert Transform for plunge

    E_pitch = sqrt(y_pitch.^2 - yh_pitch.^2); % E(t) Envelope function for pitch
    E_plunge = sqrt(y_plunge.^2 - yh_plunge.^2); % E(t) Envelope function for plunge
    
    t = (1:length(E_pitch))*T_sampling; % [s] Time steps
    t1 = max(t); % [s] End time 

    t_bar_pitch = sum(E_pitch.*t*T_sampling)/sum(E_pitch*T_sampling); % Time centroid for pitch
    t_bar_plunge = sum(E_plunge.*t*T_sampling)/sum(E_plunge*T_sampling); % Time centroid for plunge

    S_pitch(i) = real(1/t_bar_pitch - 2/t1);
    S_plunge(i) = real(1/t_bar_plunge - 2/t1);
end

% 2nd order interpolations ------------------------------------------------
coeffs_pitch = polyfit(airspeed', S_pitch, 2);
coeffs_plunge = polyfit(airspeed', S_plunge, 2);
xPos = (0:0.1:21);
interpolVal_pitch = polyval(coeffs_pitch, (0:0.1:21));
interpolVal_plunge = polyval(coeffs_plunge, (0:0.1:21));

% Plots -------------------------------------------------------------------
figure;
plot(airspeed,S_pitch,'bo');
hold on;
plot(airspeed,S_plunge,'ro');

plot((0:0.1:21),interpolVal_pitch,'b');
plot((0:0.1:21),interpolVal_plunge,'r');
zeroLine = yline(0,':');

xlabel('Airspeed $U$ [m/s]','Interpreter','latex','Fontsize',24);
ylabel('Envelope function stability criterion $S$','Interpreter','latex','Fontsize',24);
legend('Pitch data', 'Plunge data', 'Pitch $2^{nd}$ order fit','Plunge $2^{nd}$ order fit','Interpreter','latex')
set(gca,'fontsize',24,'TickLabelInterpreter','latex')
set(zeroLine,'linewidth',1.5);
grid on;
