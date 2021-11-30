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

% number on the first collumn are the start value whereas the one on second
% collumn are the end value

index = [5820,2950,6000,3780,2450,4650,750;...
         6350,3500,6400,4300,2900,5500,5280];

Index_bonus = [6980,4250,7150,5550,3150,5900,750;...
               7500,4800,7600,6000,3600,6800,5280];

%% Computation of the frequency and damping

% Using the half point method

for i = 1:1
    windowidth = (index(2,i)-index(1,i))*T_sampling;
    t = 0:T_sampling:windowidth;
    
    smooth_plunge = sgolayfilt(plunge(i,index(1,i):index(2,i)),9,27);
    smooth_pitch = sgolayfilt(pitch(i,index(1,i):index(2,i)),9,27);
    
    fftplunge = fft(smooth_plunge);
    fftpitch = fft(smooth_pitch);
     
    [a,dplunge(i,:,:)] = rfp(fftplunge,smooth_plunge,4);
    [b,dpitch(i,:,:)] = rfp(fftpitch,smooth_pitch,4);

end

%plot the damp and freq
%figure('name','frequency')



