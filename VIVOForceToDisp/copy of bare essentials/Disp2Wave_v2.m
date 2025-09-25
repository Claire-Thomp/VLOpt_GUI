function []=Disp2Wave_v2(input,output)
    
    if ~ischar(input) && ~isstring(input)
        error('Input file name must be a character vector or string');
    end

    if ~ischar(output) && ~isstring(output)
        error('Output file name must be a character vector or string');
    end

    % Convert to character vector if needed
    if isstring(input)
        input = char(input);
    end
    if isstring(output)
        output = char(output);
    end
fname=input;
new_sample_rate=1024;
[data smooth_data]=cyclicsmooth_forces_last_cycle(fname,new_sample_rate);
disp('confirm this is actually kinematics!')
smooth_kinematics=smooth_data(:,7:12);
CreateWaveFile(output,smooth_kinematics');
disp(output)
