%% This code will initialize the data from the labo 1.
clc;
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
    
    figure('name',['Figure for airspeed equal to ', num2str(airspeed(i)),' ']);
    hold on
    xlabel('Time of sampling')
    ylabel('Acceleration in m/s^2')
    plot(x*T_sampling,pitch(i,:),x*T_sampling,plunge(i,:));
    legend('pitch','plunge');
    grid on
end

