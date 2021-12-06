%code for FIV 1977


%% data for damp = 0.00145
[U_f_1,f_1] = importfile('UF_1.csv', 5, 25);
[U_A_1,A_1] = importfile('UA_1.csv', 5, 34);



%% data for damp = 0.00181
[U_f_2,f_2] = importfile('UF_2.csv', 5, 28);
[U_A_2,A_2] = importfile('UA_2.csv', 5, 28);


%% Graph

figure('name','Maximum amplitude as for different value , undimensionalized')
hold on
plot(U(2:end)/(fmax*D),A_max(2:end)/D, '-o', 'linewidth',1.5)
plot(U_A_1,A_1, 'o', 'linewidth',1.5)
plot(U_A_2,A_2, 'o', 'linewidth',1.5)
xlabel('$U/fD$ [m/s]','FontSize', 12, 'Interpreter', 'latex');
ylabel('$A_{max}/D$ ','FontSize', 12, 'Interpreter', 'latex');
lgd = legend('our case','0.00145','0.00181');
set(lgd, 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'TickLabelInterpreter','latex','Fontsize',16)
grid on
grid minor


Str_law = (Str/D) * U;
figure('name',' Variation of the motion and wake frequencies depending on damping ratio')
hold on
plot(U(2:end)/(fmax*D),freq_wake(2:end)/fmax,'o', 'linewidth',1.5)
plot(U_f_1,f_1,'o', 'linewidth',1.5)
plot(U_f_2,f_2,'o', 'linewidth',1.5)
plot(U/(fmax*D),Str_law/fmax,'-', 'linewidth',1.5)
plot(U/(fmax*D),linspace(fmax,fmax,length(U))/fmax,'-','color','k','linewidth',1.5)
xlabel('Flow velocity [m/s]','FontSize', 12, 'Interpreter', 'latex');
ylabel('Frequency [Hz]','FontSize', 12, 'Interpreter', 'latex');
lgd = legend('our case','0.00145','0.00181');
set(lgd, 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'TickLabelInterpreter','latex','Fontsize',16)
grid on
grid minor


