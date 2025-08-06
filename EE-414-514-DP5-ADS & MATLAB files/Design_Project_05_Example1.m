%% Devon Overgaard 
%----------------------------------------------------------------------
% EE 414
%  Design Project 05
%----------------------------------------------------------------------
clearvars
clc
close all
%----------------------------------------------------------------------
% Units 
%----------------------------------------------------------------------
G = 10^9;
Meg = 10^6;
k = 10^+3;
c = 10^-2;
m = 10^-3;
u = 10^-6;
n = 10^-9;
%----------------------------------------------------------------------
% 
%----------------------------------------------------------------------
IFigure = 0;
NF = 32;
j = 1*j;
dfreq = 1;
df = 1*Meg;
f_min = 1*G;
f_max = 5*G;
freq = f_min : df : f_max;
freq = sort(freq);
freq = freq';
%----------------------------------------------------------------------
% port 4 is isolated
%----------------------------------------------------------------------
Z0 = 50;
f0 = 3*G;
theta_P = 75;

C_dB = 20;
C = 10^-(C_dB/20);
er = 2.2;

Z0e = Z0*sqrt((1+abs(C))/(1-abs(C)));
Z0o = Z0*sqrt((1-abs(C))/(1+abs(C)));
Print_Real_Unit('Z0e',Z0e,'Ohms')
Print_Real_Unit('Z0o',Z0o,'Ohms')
theta_e = 90;
theta_o = 90;

phi_e = (1/2)*theta_e;
phi_o = (1/2)*theta_o;

Zin_ee = -j*Z0e*cotd(phi_e);
Zin_eo = +j*Z0e*tand(phi_e);

Zin_oe = -j*Z0o*cotd(phi_o);
Zin_oo = +j*Z0o*tand(phi_o);


Print_Polar_Unit('Zin_ee',Zin_ee,'ohms')
Print_Polar_Unit('Zin_oe',Zin_oe,'ohms')
Print_Polar_Unit('Zin_eo',Zin_eo,'ohms')
Print_Polar_Unit('Zin_oo',Zin_oo,'ohms')


S11_ee=(Zin_ee - Z0)/(Zin_ee + Z0);
S11_eo=(Zin_eo - Z0)/(Zin_eo + Z0);
S11_oe=(Zin_oe - Z0)/(Zin_oe + Z0);
S11_oo=(Zin_oo - Z0)/(Zin_oo + Z0);


S11 = (1/4)*(S11_ee + S11_eo + S11_oe + S11_oo);
S21 = (1/4)*(S11_ee - S11_eo + S11_oe - S11_oo);
S31 = (1/4)*(S11_ee + S11_eo - S11_oe - S11_oo);
S41 = (1/4)*(S11_ee - S11_eo - S11_oe + S11_oo);
Ep = Polar_2_Rect(1, -theta_P);
u = eye(4,4);

s = [S11 S21 S31 S41;
      S21 S11 S41 S31;
      S31 S41 S11 S21;
      S41 S31 S21 S11];


Ep = Ep * u;
sp = Ep * s * Ep;

S11 = sp(1,1);
S21 = sp(2,1);
S31 = sp(3,1);
S41 = sp(4,1);
Print_Polar('S11',S11)
Print_Polar('S21',S21)
Print_Polar('S31',S31)
Print_Polar('S41',S41)
Print_Break
%%
RL = -20*log10((S11));
IL = -20*log10((S21));
C = -20*log10((S31));
I = -20*log10((S41));
D = I - C;
%
N_Freq = length(freq);
S_Filter_TRL0 = zeros(N_Freq,4,4);
for kk = 1 : N_Freq
fk = freq(kk);
u = eye(4,4);
theta_e = (fk / f0) * 90;
theta_o = (fk / f0) * 90;
TRL0 = EE414_Zparam_CPL_4_Port(Z0e,Z0o,theta_e,theta_o);
Zoc0 = TRL0;
S_Filter_TRL0(kk, :, :) = Z_2_S(Zoc0,[Z0] );
end
%
S11_Filter_TRL0 = S_Filter_TRL0(:,1,1);
S11_Filter_Mag_TRL0 = abs(S11_Filter_TRL0);
S11_Filter_dB_TRL0 = 20*log10(S11_Filter_Mag_TRL0);

S21_Filter_TRL0 = S_Filter_TRL0(:,2,1);
S21_Filter_Mag_TRL0 = abs(S21_Filter_TRL0);
S21_Filter_dB_TRL0 = 20*log10(S21_Filter_Mag_TRL0);

S31_Filter_TRL0 = S_Filter_TRL0(:,3,1);
S31_Filter_Mag_TRL0 = abs(S31_Filter_TRL0);
S31_Filter_dB_TRL0 = 20*log10(S31_Filter_Mag_TRL0);

S41_Filter_TRL0 = S_Filter_TRL0(:,4,1);
S41_Filter_Mag_TRL0 = abs(S41_Filter_TRL0);
S41_Filter_dB_TRL0 = 20*log10(S41_Filter_Mag_TRL0);

Print_Real_Unit('RL',RL,'dB')
Print_Real_Unit('IL',IL,'dB')
Print_Real_Unit('C',C,'dB')
Print_Real_Unit('I',I,'dB')
Print_Real_Unit('D',D,'dB')

Print_Polar('S11',S11)
Print_Polar('S21',S21)
Print_Polar('S31',S31)
Print_Polar('S41',S41)
%
% Case 2
% Given
f0 = 3*G;
theta_e_2 = 94.0415;
theta_o_2 = 85.9434; 
lambda_fs = physconst('LightSpeed')/f0;
er_eff_e = 5.0240;
er_eff_o = 4.1960;
%  Calc LC
lambda_ge = lambda_fs/sqrt(er_eff_e);
LC = (theta_e_2 * lambda_ge) / 360;
Print_Real_Unit('LC',LC,'m')
Print_Real_Unit('lambda_fs',lambda_fs,'m')



%  Calc LC
%lambda_ge = lambda_fs/sqrt(er_eff_e);
lambda_go = lambda_fs/sqrt(er_eff_o);

% Calc S-Parameters
phi_e = (1/2)*theta_e_2;
phi_o = (1/2)*theta_o_2;

Zin_ee = -j*Z0e*cotd(phi_e);
Zin_eo = +j*Z0e*tand(phi_e);

Zin_oe = -j*Z0o*cotd(phi_o);
Zin_oo = +j*Z0o*tand(phi_o);


Print_Polar_Unit('Zin_ee',Zin_ee,'ohms')
Print_Polar_Unit('Zin_oe',Zin_oe,'ohms')
Print_Polar_Unit('Zin_eo',Zin_eo,'ohms')
Print_Polar_Unit('Zin_oo',Zin_oo,'ohms')


S11_ee=(Zin_ee - Z0)/(Zin_ee + Z0);
S11_eo=(Zin_eo - Z0)/(Zin_eo + Z0);
S11_oe=(Zin_oe - Z0)/(Zin_oe + Z0);
S11_oo=(Zin_oo - Z0)/(Zin_oo + Z0);


S11 = (1/4)*(S11_ee + S11_eo + S11_oe + S11_oo);
S21 = (1/4)*(S11_ee - S11_eo + S11_oe - S11_oo);
S31 = (1/4)*(S11_ee + S11_eo - S11_oe - S11_oo);
S41 = (1/4)*(S11_ee - S11_eo - S11_oe + S11_oo);
Ep = Polar_2_Rect(1, -theta_P);
u = eye(4,4);

s = [S11 S21 S31 S41;
      S21 S11 S41 S31;
      S31 S41 S11 S21;
      S41 S31 S21 S11];


Ep = Ep * u;
sp = Ep * s * Ep;

S11 = sp(1,1);
S21 = sp(2,1);
S31 = sp(3,1);
S41 = sp(4,1);
Print_Polar('S11',S11)
Print_Polar('S21',S21)
Print_Polar('S31',S31)
Print_Polar('S41',S41)
%%
RL = -20*log10((S11));
IL = -20*log10((S21));
C = -20*log10((S31));
I = -20*log10((S41));
D = I - C;

Print_Real_Unit('RL',RL,'dB')
Print_Real_Unit('IL',IL,'dB')
Print_Real_Unit('C',C,'dB')
Print_Real_Unit('I',I,'dB')
Print_Real_Unit('D',D,'dB')

Print_Polar('S11',S11)
Print_Polar('S21',S21)
Print_Polar('S31',S31)
Print_Polar('S41',S41)
% Plot Attempt
% function [Z] =EE414_Zparam_CPL_4_Port(ZOe,ZOo,theta_e,theta_o)


N_Freq = length(freq);
S_Filter_TRL = zeros(N_Freq,4,4);
for kk = 1 : N_Freq
fk = freq(kk);
u = eye(4,4);
theta_e = (fk / f0) * theta_e_2;
theta_o = (fk / f0) * theta_o_2;
TRL1 = EE414_Zparam_CPL_4_Port(Z0e,Z0o,theta_e,theta_o);
Zoc = TRL1;
S_Filter_TRL(kk, :, :) = Z_2_S(Zoc,[Z0] );
end

S11_Filter_TRL = S_Filter_TRL(:,1,1);
S11_Filter_Mag_TRL = abs(S11_Filter_TRL);
S11_Filter_dB_TRL = 20*log10(S11_Filter_Mag_TRL);

S21_Filter_TRL = S_Filter_TRL(:,2,1);
S21_Filter_Mag_TRL = abs(S21_Filter_TRL);
S21_Filter_dB_TRL = 20*log10(S21_Filter_Mag_TRL);

S31_Filter_TRL = S_Filter_TRL(:,3,1);
S31_Filter_Mag_TRL = abs(S31_Filter_TRL);
S31_Filter_dB_TRL = 20*log10(S31_Filter_Mag_TRL);

S41_Filter_TRL = S_Filter_TRL(:,4,1);
S41_Filter_Mag_TRL = abs(S41_Filter_TRL);
S41_Filter_dB_TRL = 20*log10(S41_Filter_Mag_TRL);

% %----------------------------------------------------------------------
% % ADS V3
% %----------------------------------------------------------------------
FileName_ADS = 'CPL_Example_TRL_V2.s4p';
Folder_ADS = 'G:\EE 414 Microwave Engineering\Design Project 5\';
File_Loc_ADS_Amp = [Folder_ADS,FileName_ADS];
[f_ADS, S_ADS, Mult_ADS] = ...
Read_SParam_sNp(File_Loc_ADS_Amp);
f_ADS = f_ADS * Mult_ADS;
I_f0_ADS = f_ADS == f0;

% %----------------------------------------------------------------------
% %S ADS Parameters
% %----------------------------------------------------------------------
S11_ADS = S_ADS(:, 1, 1);
S11_ADS_Mag = abs(S11_ADS);
S11_ADS_dB = 20*log10(S11_ADS_Mag);

S21_ADS = S_ADS(:, 2, 1);
S21_ADS_Mag = abs(S21_ADS);
S21_ADS_dB = 20*log10(S21_ADS_Mag);

S31_ADS = S_ADS(:, 3, 1);
S31_ADS_Mag = abs(S31_ADS);
S31_ADS_dB = 20*log10(S31_ADS_Mag);

S41_ADS = S_ADS(:, 4, 1);
S41_ADS_Mag = abs(S41_ADS);
S41_ADS_dB = 20*log10(S41_ADS_Mag);
Print_Polar('S_11',S11_ADS(201))
RL_2 = abs(S11_ADS_dB(201));
IL_2 = abs(S21_ADS_dB(201));
C_2 = abs(S31_ADS_dB(201));
I_2 = abs(S41_ADS_dB(201));
D_2 = I_2 - C_2;
Print_Real_Unit('RL_2',RL_2,'dB')
Print_Real_Unit('IL_2',IL_2,'dB')
Print_Real_Unit('C_2',C_2,'dB')
Print_Real_Unit('I_2',I_2,'dB')
Print_Real_Unit('D_2',D_2,'dB')
%%
% %----------------------------------------------------------------------
% % ADS V2
% %----------------------------------------------------------------------
FileName_ADS3 = 'CPL_Example_TRL_V3.s4p';
Folder_ADS3 = 'G:\EE 414 Microwave Engineering\Design Project 5\';
File_Loc_ADS_Amp3 = [Folder_ADS,FileName_ADS3];
[f_ADS3, S_ADS3, Mult_ADS3] = ...
Read_SParam_sNp(File_Loc_ADS_Amp3);
f_ADS3 = f_ADS3 * Mult_ADS3;


% %----------------------------------------------------------------------
% %S ADS Parameters
% %----------------------------------------------------------------------
S11_ADS3 = S_ADS3(:, 1, 1);
S11_ADS_Mag3 = abs(S11_ADS3);
S11_ADS_dB3 = 20*log10(S11_ADS_Mag3);

S21_ADS3 = S_ADS3(:, 2, 1);
S21_ADS_Mag3 = abs(S21_ADS3);
S21_ADS_dB3 = 20*log10(S21_ADS_Mag3);

S31_ADS3 = S_ADS3(:, 3, 1);
S31_ADS_Mag3 = abs(S31_ADS3);
S31_ADS_dB3 = 20*log10(S31_ADS_Mag3);

S41_ADS3 = S_ADS3(:, 4, 1);
S41_ADS_Mag3 = abs(S41_ADS3);
S41_ADS_dB3 = 20*log10(S41_ADS_Mag3);

RL_3 = abs(S11_ADS_dB3(201));
IL_3 = abs(S21_ADS_dB3(201));
C_3 = abs(S31_ADS_dB3(201));
I_3 = abs(S41_ADS_dB3(201));
D_3 = I_3 - C_3;
Print_Real_Unit('RL_3',RL_3,'dB')
Print_Real_Unit('IL_3',IL_3,'dB')
Print_Real_Unit('C_3',C_3,'dB')
Print_Real_Unit('I_3',I_3,'dB')
Print_Real_Unit('D_3',D_3,'dB')
%
% %----------------------------------------------------------------------
% % S PARAMETER f0 
% %----------------------------------------------------------------------

IFigure = IFigure + 1;
figure_max(IFigure)
plot(freq/G, S11_Filter_dB_TRL0, 'r', 'linewidth', 6)
hold on
plot(freq/G, S11_Filter_dB_TRL, 'b', 'linewidth', 6)
plot(f_ADS/G, S11_ADS_dB, 'm', 'linewidth', 6);
plot(f_ADS/G, S11_ADS_dB3, 'c', 'linewidth', 6);
hold off
grid on
grid minor 
xlabel('{\itf} (GHz) ')
ylabel('| {\itS}_{11} | ( dB ) ', ...
'VerticalAlignment', 'bottom')
set(gca, 'FontName', 'times new roman', 'FontSize', NF)
title('S11')
legend('MATLAB (V1)','MATLAB (V2)','ADS (MStrip V1)','ADS (MStrip V2)')
axis([1*G/G, 5*G/G, -60 0])
set(gca, 'xtick', f_min/G : 0.5 : f_max/G);
set(gca, 'linewidth', 2.0)

IFigure = IFigure + 1;
figure_max(IFigure)
plot(freq/G, S21_Filter_dB_TRL0, 'r', 'linewidth', 6)
hold on
plot(freq/G, S21_Filter_dB_TRL, 'b', 'linewidth', 6);
plot(f_ADS/G, S21_ADS_dB, 'm', 'linewidth', 6);
plot(f_ADS/G, S21_ADS_dB3, 'c', 'linewidth', 8);
hold off
grid on
grid minor 
xlabel('{\itf} (GHz) ')
ylabel('| {\itS}_{21} | ( dB ) ', ...
'VerticalAlignment', 'bottom')
set(gca, 'FontName', 'times new roman', 'FontSize', NF)
title('S21')
legend('MATLAB (V1)','MATLAB (V2)','ADS (MStrip V1)','ADS (MStrip V2)')
axis([1*G/G, 5*G/G, -0.3, 0])
set(gca, 'xtick', f_min/G : 0.5 : f_max/G);
set(gca, 'linewidth', 2.0)
%
IFigure = IFigure + 1;
figure_max(IFigure)
plot(freq/G, S31_Filter_dB_TRL0, 'r', 'linewidth', 6)
hold on
plot(freq/G, S31_Filter_dB_TRL, 'b', 'linewidth', 6);
plot(f_ADS/G, S31_ADS_dB, 'm', 'linewidth', 6);
plot(f_ADS/G, S31_ADS_dB3, 'c', 'linewidth', 8);
hold off
grid on
grid minor 
xlabel('{\itf} (GHz) ')
ylabel('| {\itS}_{31} | ( dB ) ', ...
'VerticalAlignment', 'bottom')
set(gca, 'FontName', 'times new roman', 'FontSize', NF)
title('S31')
legend('MATLAB (V1)','MATLAB (V2)','ADS (MStrip V1)','ADS (MStrip V2)')
axis([1*G/G, 5*G/G, -26, -19])
set(gca, 'xtick', f_min/G : 0.5 : f_max/G);
set(gca, 'linewidth', 2.0)


IFigure = IFigure + 1;
figure_max(IFigure)
plot(freq/G, S41_Filter_dB_TRL0, 'r', 'linewidth', 6)
hold on
plot(freq/G, S41_Filter_dB_TRL, 'b', 'linewidth', 6);
plot(f_ADS/G, S41_ADS_dB, 'm', 'linewidth', 6);
plot(f_ADS/G, S41_ADS_dB3, 'c', 'linewidth', 8);
hold off
grid on
grid minor 
xlabel('{\itf} (GHz) ')
ylabel('| {\itS}_{41} | ( dB ) ', ...
'VerticalAlignment', 'bottom')
set(gca, 'FontName', 'times new roman', 'FontSize', NF)
title('S41')
legend('MATLAB (V1)','MATLAB (V2)','ADS (MStrip V1)','ADS (MStrip V2)')
axis([1*G/G, 5*G/G, -40, 0])
set(gca, 'xtick', f_min/G : 0.5 : f_max/G);
set(gca, 'linewidth', 2.0)