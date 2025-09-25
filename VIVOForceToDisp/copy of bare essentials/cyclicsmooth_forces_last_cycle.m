function [data new_data]=cyclicsmooth_forces_last_cycle(fname,new_sample_rate)
% Smooth data


% Get sample and cycle frequencies
fid=fopen(fname,'r');
for i=1:5
    tline=fgetl(fid);
end
fclose(fid);

[temp]=sscanf(tline,'%*s %*s %f %*s %*s %f');
samp_freq=temp(1);
cycle_freq=temp(2);

% Get overall file length
fid=fopen(fname,'r');
B=0;
while 1
    tline=fgetl(fid);
    if tline==-1
        break
    end
    B=B+1;
end

fclose(fid);


% Pull out the actual data
data=dlmread(fname,'\t',[9 0 B-1 46]);
n_seconds=(B-9)/samp_freq;
n_cycles=n_seconds*cycle_freq;
sample_per_cycle=(B-9)/n_cycles;

n_extracted_cycles=1;
last_cycle=[(size(data,1)-round(sample_per_cycle))+1:1:size(data,1)];


assert(all(mod(last_cycle, 1) == 0), 'last_cycle contains non-integer values');

data=[data(last_cycle,1:12) data(last_cycle,42:47) data(last_cycle,13:18) data(last_cycle,25:30)];



for j=1:size(data,2)
    size(data(:,j));
    repeated_data=[data(:,j); data(:,j); data(:,j)];
    size(repeated_data);
    
    filtered_data(:,1)=butterworth(samp_freq,samp_freq/100,2,repeated_data(:,1));
    interp_data(:,1)=interp1(1/sample_per_cycle:1/sample_per_cycle:n_extracted_cycles*3,filtered_data(:,1),1/new_sample_rate:1/new_sample_rate:n_extracted_cycles*3);
    
    
    % extract sample from middle
    
    for i=1:n_extracted_cycles*3
        start_range=(i-1)*new_sample_rate+1;
        end_range=start_range+new_sample_rate-1;
        
        interp_sample{i}=interp_data(start_range:end_range,:);
    end
    % Only use the middle sample; the others have "end effects".... 
    new_data(:,j)=interp_sample{2};

end

%csvwrite(fout,new_data);

% sample 1 = 1 --> 600+1
% sample 2 = 602 --> 600+602
% sample 3