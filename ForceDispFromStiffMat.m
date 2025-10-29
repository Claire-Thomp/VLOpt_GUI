function [vivo_file] = ForceDispFromStiffMat(MomentWave,stiffMat,tag)
%Takes a stiffness matrix and displacement waveform and creates a force
%displacement curve to optimize IVD virtual ligaments to 

%% Setup
%Stiffness matrix
stiffMat = readtable(stiffMat);
stiffMat.Properties.VariableNames = {'TX', 'TY', 'TZ', 'RX', 'RY', 'RZ'};
stiffMat = table2array(stiffMat);

convMat = "/ConversionMatrix.txt";
Unit_conversion = readtable(convMat);
Unit_conversion = table2array(Unit_conversion);

%Put stiffness matrix into correct units (based on Stokes paper) 
stiffMat_base = stiffMat .* Unit_conversion;

%Transform stiffness matrix from paper coordinates to VIVO coordinates
%Rotate about z
phi = -90; %degrees
Zrot = [cosd(phi) -sind(phi) 0;sind(phi) cosd(phi) 0;0 0 1];
Zrot_6by6 = blkdiag(Zrot, Zrot);
stiffMat_base_z = Zrot_6by6' * stiffMat_base * Zrot_6by6;

%Reflect about y (mirror X axis)
theta = 180; %degrees

Refx = [cosd(2*theta) 0 sind(2*theta); 0 1 0; sind(2*theta) 0 -cosd(2*theta)];
Refx_6by6 = blkdiag(Refx, Refx);
stiffMat_base_z_y = Refx_6by6' * stiffMat_base_z * Refx_6by6;

%Reflect about x (mirror Y axis)
zeta = 180; %degrees

Refy = [1 0 0; 0 cosd(2*zeta) sind(2*zeta); 0 sind(2*zeta) -cosd(2*zeta)];
Refy_6by6 = blkdiag(Refy, Refy);
stiffMat_base_z_y_x = Refy_6by6' * stiffMat_base_z_y * Refy_6by6;

stiffMat_inv = inv(stiffMat_base_z_y_x);

%VIVO File
emptyVivo = "/EmptyVivo.txt";
vivo_file = readtable(emptyVivo, VariableNamingRule="preserve");

%% 

% %if input is displacement
% DispWave = table2array(DispWave);
% 
% % Matrix Math (create optimization input)
% if tag == "FE"
%     for i=1:(height(vivo_file))
%         vivo_file.Xdeg(i) = DispWave(1,i);
%     end
% elseif tag == "LB"
%     for i=1:(height(vivo_file))
%         vivo_file.Ydeg(i) = DispWave(1,i);
%     end
% else
%     for i=1:(height(vivo_file))
%         vivo_file.Zdeg(i) = DispWave(1,i);
%     end
% end 
% % Create intact file 
% DispMat = [vivo_file.Xmm, vivo_file.Ymm,vivo_file.Zmm,vivo_file.Xdeg,vivo_file.Ydeg,vivo_file.Zdeg];
% ForceMat = zeros(length(DispMat),6);
% for i = 1:length(DispMat)
%     DispMatRow = DispMat(i,:);
%     ForceMat(i,:) = stiffMat*DispMatRow';
% end
% for k = 1:height(vivo_file)
%     vivo_file.Fx(k) = ForceMat(k,1);
%     vivo_file.Fy(k) = ForceMat(k,2);
%     vivo_file.Fz(k) = ForceMat(k,3);
%     vivo_file.Mx(k) = ForceMat(k,4);
%     vivo_file.My(k) = ForceMat(k,5);
%     vivo_file.Mz(k) = ForceMat(k,6);
% end
%% 

%if input is force
MomentWave = table2array(MomentWave);
MomentWave = MomentWave';
% Matrix Math (create optimization input)
if tag == "FE"
    for i=1:(height(vivo_file))
        vivo_file.Mx(i) = MomentWave(1,i);
    end
elseif tag == "LB"
    for i=1:(height(vivo_file))
        vivo_file.My(i) = MomentWave(1,i);
    end
else
    for i=1:(height(vivo_file))
        vivo_file.Mz(i) = MomentWave(1,i);
    end
end  

% Create intact file 
ForceMat = [vivo_file.Fx, vivo_file.Fy,vivo_file.Fz,vivo_file.Mx,vivo_file.My,vivo_file.Mz];
DispMat = zeros(length(ForceMat),6);
for i = 1:length(ForceMat)
    ForceMatRow = ForceMat(i,:);
    DispMat(i,:) = stiffMat_inv*ForceMatRow';
end

for k = 1:height(vivo_file)
    vivo_file.Xmm(k) = DispMat(k,1);
    vivo_file.Ymm(k) = DispMat(k,2);
    vivo_file.Zmm(k) = DispMat(k,3);
    vivo_file.Xdeg(k) = DispMat(k,4);
    vivo_file.Ydeg(k) = DispMat(k,5);
    vivo_file.Zdeg(k) = DispMat(k,6);
end
% % Create cut file **I'm not sure we need this**
% outputPath = "/VLOptCURRENT/6DOF_Vivo_Cut.txt";  % Change this to your desired path/filename
% writetable(vivo_file, outputPath, 'Delimiter', '\t', 'FileType', 'text');


% outputPath = "N:\VLOptCURRENT\1DOF_noNZ_Vivo_Intact.txt";  % Change this to your desired path/filename
% writetable(vivo_file, outputPath, 'Delimiter', '\t', 'FileType', 'text');

end