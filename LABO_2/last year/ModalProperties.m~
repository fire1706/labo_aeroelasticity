%% Code for labo 2
clear all
close all


Freq_sampling = 201.03; %Hz
T_sampling = 1/Freq_sampling; %s

%% Initialization of the data
 load('DATAg2.mat');
   x = 1:20480;
   puls=zeros(11,20480);
   L = length(puls);


for i = 2:11
    airspeed(i-1) = DATAg2(i).U;
      for j=1:length(DATAg2(i).w)
         puls(i-1,j)=DATAg2(i).w(j);
      end
    
    figure(i-1);
    
    xlabel('Time of sampling')
    ylabel('Acceleration in m/s^2')
    plot(x*T_sampling,puls(i-1,:));
    %plot(x,puls(i-1,:));
    
    legend('puls');
    grid on
end

%%

