function []=CreateWaveFile(filename,wavedata)

fname=[filename '.wvp'];
copyfile('backup/Template.wvp',fname)

binfile = memmapfile(fname, 'offset',94808, 'writable', true, 'format', {'double' [8 1024] 'Amplitude'});

if size(wavedata)==[6 1024]
    binfile.Data.Amplitude(2:7,:)=wavedata;
else
    disp('Provide wavedata is not appropriate size (must contain 6 rows, 1024 numbers long)')
end

clear binfile

% function CreateWaveFile(output, wavedata)
%     % Ensure 'output' is a character vector or string
%     if ~ischar(output) && ~isstring(output)
%         error('Output file name must be a character vector or string');
%     end
% 
%     % Convert to character vector if needed
%     if isstring(output)
%         output = char(output);
%     end
% 
%     fname=[output '.wvp'];
%     % Create the wave file from the template
%     copyfile('backup\Template.wvp', output);
% 
%     binfile = memmapfile(fname, 'offset',94808, 'writable', true, 'format', {'double' [8 1024] 'Amplitude'});
% 
%     if size(wavedata)==[6 1024]
%         binfile.Data.Amplitude(2:7,:)=wavedata;
%     else
%         disp('Provide wavedata is not appropriate size (must contain 6 rows, 1024 numbers long)')
%     end
% 
%     clear binfile
% 
% end
