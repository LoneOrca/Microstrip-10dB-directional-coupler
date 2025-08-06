%--------------------------------------------------------------------

function [freq, S, Mult, S_Type] = Read_SParam_sNp(FileName, ~)

%--------------------------------------------------------------------
% Read_SParam_snp
%--------------------------------------------------------------------
%   Written by John J. Burke 08/01/2017
%   Modified by John J. Burke 08/01/2017
%   Modified by John J. Burke 03/23/2018
%   Modified by John J. Burke 05/05/2018
%   Modified by John J. Burke 03/30/2020
%   Modified by John J. Burke 05/10/2022
%--------------------------------------------------------------------

FileName1 = 'temp1.txt';
fid = fopen(FileName);
fid_1 = fopen(FileName1, 'wt');
IFlag = 0;
ILine = 0;

%--------------------------------------------------------------------

while ~feof(fid) % loop over the following until the end of the file is reached.
    line = fgets(fid); % read in one line
    IFlag_1 = 0;
    IFlag_2 = 0;
    if contains(line, '#')
        [Mult] = freq_mult(line);
        [S_Type] = Sparam_Type(line);
        IFlag_1 = 1;
    end
    if contains(line, '!')
        IFlag_2 = 1;
    end
    if (IFlag_1 == 0 && IFlag_2 == 0)
        Data = str2num(line);
        N_Line = length(Data);
        if (N_Line ~= 0)
            ILine = ILine + 1;
            fprintf(fid_1, '%s\n', line(1:end-1));
            if (IFlag == 0)
                NPorts = ( N_Line - 1 ) / 2;
                IFlag = 1;
            end
        end
    end
end
fclose(fid);
fclose(fid_1);

%--------------------------------------------------------------------

Nfreq = ILine / NPorts;
fid_1 = fopen(FileName1);
freq = zeros(NPorts, 1);
C1 = zeros(NPorts, NPorts);
C2 = zeros(NPorts, NPorts);
S = zeros(Nfreq, NPorts, NPorts);

%--------------------------------------------------------------------

for kk = 1:Nfreq
    for jj = 1:NPorts
        line = fgets(fid_1); % read in one line
        Data = str2num(line);
        if (jj == 1)
            freq(kk) = Data(1);
            Data = Data(2:end);
        end
        C1(jj, 1:NPorts) = Data(1:2:2*NPorts);
        C2(jj, 1:NPorts) = Data(2:2:2*NPorts);
        fgets(fid_1); % read in one line (Note this read doesn't make sense)
    end
    [S(kk, :, :)] = Sparam_Extract(C1, C2, S_Type);
end

%--------------------------------------------------------------------

fclose(fid_1);
dos('del temp1.txt');

%--------------------------------------------------------------------











%--------------------------------------------------------------------

function [Mult] = freq_mult(line)

%--------------------------------------------------------------------

T = 10^+12;
G = 10^+9;
M = 10^+6;
k = 10^+3;

%--------------------------------------------------------------------

line = lower(line);
if contains(line, 'thz')
    Mult = T;
elseif contains(line, 'ghz')
    Mult = G;
elseif contains(line, 'mhz')
    Mult = M;
elseif contains(line, 'khz')
    Mult = k;
else
    Mult = 1;
end

%--------------------------------------------------------------------










%--------------------------------------------------------------------

function [S_Type] = Sparam_Type(line)

%--------------------------------------------------------------------

line = lower(line);
if contains(line, 'ma')
    S_Type = 'MA';
elseif contains(line, 'db')
    S_Type = 'dB';
elseif contains(line, 'ri')
    S_Type = 'RI';
else
    S_Type = 'MA';
end

%--------------------------------------------------------------------









%--------------------------------------------------------------------

function [S] = Sparam_Extract(C1, C2, S_Type)

%--------------------------------------------------------------------

switch S_Type
    case 'MA'
        S = Polar_2_Rect(C1, C2);
    case 'RI'
        S = C1 + 1j*C2;
    case 'dB'
        C1 = 10.^(C1/20);
        S = Polar_2_Rect(C1, C2);
end

%--------------------------------------------------------------------
