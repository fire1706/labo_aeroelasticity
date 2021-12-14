%%damping calculation 

clear y y_1 y_2 y_3 Y_Ug2 Y_Ue1 Y_Ue2 Y_Ue3
%% G2
yg2 = DATAg2(1).y(14540:18500);
% Y_Ug2 = [0.003548,0.003262,0.003158,0.002738,0.002725,0.002402,0.002395,0.002198,...
%     0.001996,0.001981,0.001729,0.001784,0.001489,0.001576,0.001227,0.001382,...
%     0.001227,0.001252,0.001169,0.0009699,0.001036];
Y_Ug2 = decrepeak(yg2);%(1:3000));
i = 1;
while i<length(Y_Ug2)
    dampg2(i) = 1/(2*pi)*log(Y_Ug2(i)/Y_Ug2(i+1));
    i = i+1;
end
Dampg2 = mean(dampg2)


%% E1
y_1 = DATAe_1(1).y(25550:30000);
% Y_Ue1 = [0.004474,0.004173,0.003934,0.003582,0.003411,0.003196,0.002963,0.002867,...
%     0.002627,0.002485,0.002296,0.002199,0.002077,0.001922,0.001808,0.001689,...
%     0.001656,0.001618,0.001471,0.001415,0.001261,0.001179,0.001089];
Y_Ue1 = decrepeak(y_1);%(1:2001));
i = 1;
while i<length(Y_Ue1)
    dampe1(i) = 1/(2*pi)*log(Y_Ue1(i)/Y_Ue1(i+1));
    i = i+1;
end
Dampe1 = mean(dampe1)


%% E2
y_2 = DATAe_2(1).y(21840:26000);
% Y_Ue2 = [0.004015,0.004076,0.003412,0.003415,0.002954,0.002795,0.002583,0.002354,...
%     0.002324,0.001958,0.001998,0.001699,0.001695,0.001458,0.001358,0.001302,...
%     0.001173,0.001182,0.0009603];
% Y_Ue2 = [0.004076,0.003415,0.002954,0.002795,0.002583,0.002354,...
%     0.001958,0.001699,0.001458,0.001358,0.001302,...
%     0.001173,0.0009603];
Y_Ue2 = decrepeak(y_2);%(1:2000));
i = 1;
while i<length(Y_Ue2)
    dampe2(i) = 1/(2*pi)*log(Y_Ue2(i)/Y_Ue2(i+1));
    i = i+1;
end
Dampe2 = mean(dampe2)


%% E3
% 
y_3 = DATAe_3(1).y(18260:21000);
% Y_Ue3 = [0.004125,0.004072,0.003414,0.003387,0.002919,0.002696,0.002528,0.002143,0.002105,0.001804,0.001636,0.001444,0.001402,0.001102,0.001148,0.000883,...
%     0.0008209,0.0006681,0.0006236,0.0005554,];

%  Y_Ue3 = [0.004125,0.004072,0.003414,0.003387,0.002919,0.002696,0.002528,...
%      0.002143,0.002105,0.001804,0.001636,0.001444,0.001402,0.001102,0.000883,...
%     0.0008209,0.0006681,0.0006236,0.0005554,];
Y_Ue3 = decrepeak(y_3);%(1:1000));
i = 1;
while i<length(Y_Ue3)
    dampe3(i) = 1/(2*pi)*log(Y_Ue3(i)/Y_Ue3(i+1));
    i = i+1;
end
Dampe3 = mean(dampe3)
