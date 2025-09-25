cutoff_freq = 2;

% Set the path where the files are located
% pathname = "C:\MATLAB Process Results App\MATLAB Process Results App\Sample results\Intact\";
pathname = "";
%filename = "LigAB_downandML_0_00002.txt";
filename = "ForceControlIntact_0_00002.txt";

pname = char(pathname);
dname = 'Displacement_waveforms';
mkdir(fullfile(pname, dname)); % Create the directory if it doesn't exist

fullpathname = fullfile(pathname, filename);
iter = 1;

for i = 1:1
    [filepath, filename1, ext] = fileparts(fullpathname);
    
    % Construct the full path for the displacement waveform file
    fullpathnameDisp = fullfile(filepath, 'Displacement_waveforms', strcat(filename1, '_Disp'));
    
    % Call the Disp2Wave_v2 function with the constructed paths
    % Disp2Wave_v2(fullpathname, fullpathnameDisp, cutoff_freq);
    Disp2Wave_v2(fullpathname, fullpathnameDisp); % Uncomment the above line if cutoff_freq is needed
    iter = iter + 1;
end
