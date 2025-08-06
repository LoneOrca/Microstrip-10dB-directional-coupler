%--------------------------------------------------------------------
%% EE-414-514-sn620094 - Design Project 5 -  Coupled Line Ideal TRL
%----------------------------------------------------------------------
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
f_min = 3*G;
f_max = 7*G;
freq = f_min : df : f_max;
freq = sort(freq);
freq = freq';
%----------------------------------------------------------------------
% port 4 is isolated
%----------------------------------------------------------------------
Z0 = 50;
f0 = 5*G;
theta_P = 60;

C_dB = 10;
C = 10^-(C_dB/20);
er = 10.7;

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
% %----------------------------------------------------------------------
% % ADS V0
% %----------------------------------------------------------------------
FileName_ADS = 'Ideal_TRL.s4p';
Folder_ADS = 'E:\Design Project-5\EE-414-514-DP5-ADS & MATLAB files\';
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
RL_0 = abs(S11_ADS_dB(201));
IL_0 = abs(S21_ADS_dB(201));
C_0 = abs(S31_ADS_dB(201));
I_0 = abs(S41_ADS_dB(201));
D_0 = I_0 - C_0;
Print_Real_Unit('RL_0',RL_0,'dB')
Print_Real_Unit('IL_0',IL_0,'dB')
Print_Real_Unit('C_0',C_0,'dB')
Print_Real_Unit('I_0',I_0,'dB')
Print_Real_Unit('D_0',D_0,'dB')
%%
% %----------------------------------------------------------------------
% % S PARAMETER f0 
% %----------------------------------------------------------------------
IFigure = IFigure + 1;
figure_max(IFigure)
plot(f_ADS/G, S11_ADS_dB, 'r', 'linewidth', 6)
hold on
plot(f_ADS/G, S21_ADS_dB, 'g', 'linewidth', 6);
plot(f_ADS/G, S31_ADS_dB, 'b', 'linewidth', 6);
plot(f_ADS/G, S41_ADS_dB, 'm', 'linewidth', 6);
hold off
grid on
grid minor 
xlabel('{\itf} (GHz) ')
ylabel(' ( dB ) ', ...
'VerticalAlignment', 'bottom')
set(gca, 'FontName', 'times new roman', 'FontSize', NF)
title('Ideal TRL')
legend('S11','S21','S31','S41')
axis([3*G/G, 7*G/G, -200, 0])
set(gca, 'xtick', f_min/G : 0.5 : f_max/G);
set(gca, 'linewidth', 2.0)

