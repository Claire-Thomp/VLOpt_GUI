%Neutral Zones - choose values from literature 
Fx_Nz_pos = 0.01; %mm
Fx_Nz_neg = -0.01;  %mm
Fy_Nz_pos = 0.01; %mm
Fy_Nz_neg = -0.01; %mm
Fz_Nz_pos = 0.01; %mm
Fz_Nz_neg = -0.01; %mm
Mx_Nz_pos = 0.01; %degrees
Mx_Nz_neg = -0.01;  %degrees
My_Nz_pos = 0.01; %degrees
My_Nz_neg = -0.01; %degrees
Mz_Nz_pos = 0.01; %degrees
Mz_Nz_neg = -0.01; %degrees

file = "N:\VLOptCURRENT\ConversionMatrix.txt";
Unit_conversion = readtable(file);
Unit_conversion = table2array(Unit_conversion);

file = "N:\VLOptCURRENT\TransformMatrix.txt";
Transform = readtable(file);
Transform = table2array(Transform);

%Stiffness matrix 
file = "N:\VLOptCURRENT\StiffnessMatrix250pre.txt";
stiffMat = readtable(file);
stiffMat.Properties.VariableNames = {'TX', 'TY', 'TZ', 'RX', 'RY', 'RZ'};
stiffMat = table2array(stiffMat);
stiffMat = mtimes(stiffMat, Unit_conversion); %Changes from radians to degrees for all rotations
stiffMat = mtimes(stiffMat, Transform); %changes stiffness matrix coordinates to VIVO coordinates

%VIVO File
file2 = "N:\VLOptCURRENT\EmptyVivo.txt";
vivo_file = readtable(file2, VariableNamingRule="preserve");
%ensure vivo file has even height
if ~mod(height(vivo_file),2)
    last_row = vivo_file(end, :);
    repeated_rows = repmat(last_row, 1, 1);
    vivo_file = [vivo_file; repeated_rows];
end

%% Disp Targets

Xdeg_pos = 10; %deg
Xdeg_neg = 10;
%Ydeg_pos = 0;
%Ydeg_neg = 2;
%Zdeg_pos = 3;
%Zdeg_neg = 3;

%XMm_pos = 0; %mm
%XMm_neg = 20;
%YMm_pos = 0;
%YMm_neg = 10;
ZMm_pos = 2;
ZMm_neg = 2;

%% Create Waves
DispWave = zeros(height(vivo_file),6);

%triangle wave
% Define a time array
fs = (height(vivo_file)-1);  % Sampling frequency
duration = 1; % Duration in seconds
t = 0:1/fs:duration;
% Generate a triangle wave(s) in desired column(s) of DispWave
DispWave(:,4) = (Xdeg_pos)*sawtooth(2*pi*1*t+(pi/2), 0.5);  % 1 Hz triangle wave
%DispWave(:,3) = (ZMm_pos)*sawtooth(2*pi*1*t+(pi/2), 0.5);
%linear wave
% t2 = 0:1/fs:0;
% DispWave(:,3) = (t2+ZMm_pos);

% Parameters: number of points in each section
% n_up = round(height(vivo_file)/4);    % Number of points from 0 to 0.75
% n_hold = round(height(vivo_file)/2);  % Number of points holding at 0.75
% n_down = round(height(vivo_file)/4);  % Number of points from 0.75 back to 0

% % Create each section
% up = linspace(0, 0.75, n_up);
% hold_val = -0.75 * ones(1, n_hold);
% down = linspace(0.75, 0, n_down);
% 
% % Combine all sections
% DispWave(:,3) = [down, hold_val, up];

% sine wave            
% x = linspace(0, 2*pi, height(vivo_file))';       % Create input values (e.g., 0 to 2π over N points)
% amplitude = Xdeg_pos;         % Scale the wave
% frequency = 1;         % Number of full cycles over the range
% phase = 2*pi;          % Phase shift
% sineWave = amplitude * sin(frequency * x + phase);                              % Calculate sine values
% DispWave(:,1) = sineWave;                       % Add the sine wave as a new column

%Create 6DOF Displacement values for VIVO file 
for i=1:(height(vivo_file))
    vivo_file.Xmm(i) = DispWave(i,1);
    vivo_file.Ymm(i) = DispWave(i,2);
    vivo_file.Zmm(i) = DispWave(i,3);
    vivo_file.Xdeg(i) = DispWave(i,4);
    vivo_file.Ydeg(i) = DispWave(i,5);
    vivo_file.Zdeg(i) = DispWave(i,6);
end

%% Creat Cut Vivo File
outputPath = "N:\VLOptCURRENT\1DOF_NZ_Vivo_Cut.txt";  % Change this to your desired path/filename
writetable(vivo_file, outputPath, 'Delimiter', '\t', 'FileType', 'text');

%% Create Intact Vivo File

DispMat = [vivo_file.Xmm, vivo_file.Ymm,vivo_file.Zmm,vivo_file.Xdeg,vivo_file.Ydeg,vivo_file.Zdeg];
ForceMat = zeros(length(DispMat),6);
for i = 1:length(DispMat)
    
    DispMatRow = DispMat(i,:);
    ForceMat(i,:) = stiffMat*DispMatRow';
    %Neutral Zones: if the displacement is within the NZ the force is 0
    % if DispMat(i,4) > Mx_Nz_neg && DispMat(i,4) < Mx_Nz_pos
    %     ForceMat(i,4) = 0;
    % elseif DispMat(i,6) > Mz_Nz_neg && DispMat(i,6) < Mz_Nz_pos
    %     ForceMat(i,6) = 0;
    % elseif DispMat(i,3) > Fz_Nz_neg && DispMat(i,3) < Fz_Nz_pos
    %    ForceMat(i,3) = 0;
    % else
    %     ForceMat(i,:) = stiffMat*DispMatRow';
    % end
end

for k = 1:height(vivo_file)
    vivo_file.Fx(k) = ForceMat(k,1);
    vivo_file.Fy(k) = ForceMat(k,2);
    vivo_file.Fz(k) = -ForceMat(k,3);
    vivo_file.Mx(k) = ForceMat(k,4)/1000;
    vivo_file.My(k) = -ForceMat(k,5)/1000;
    vivo_file.Mz(k) = -ForceMat(k,6)/1000;
end

outputPath = "N:\VLOptCURRENT\1DOF_noNZ_Vivo_Intact.txt";  % Change this to your desired path/filename
writetable(vivo_file, outputPath, 'Delimiter', '\t', 'FileType', 'text');

