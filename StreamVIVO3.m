function [trans_results]=StreamVIVO3(vdat)

[fem_transformation, tib_transformation]=vdat_to_trans(vdat);
[flexion_transformation, gimbal_transformation, platen_transformation]=ActuatorOffsets(vdat,fem_transformation,tib_transformation);

trans_results{1}=fem_transformation;
trans_results{2}=tib_transformation;
trans_results{3}=flexion_transformation;
trans_results{4}=gimbal_transformation;
trans_results{5}=platen_transformation;

% f_transform.matrix = fem_transformation;
% f_transform.timestamp = igtlTimestampNow()-startTime;
% igtlSendTransform(igtlConnection, f_transform);
% 
% t_transform.matrix = tib_transformation;
% t_transform.timestamp = igtlTimestampNow()-startTime;
% igtlSendTransform(igtlConnection, t_transform);
% 
% flexion_transform.matrix = flexion_transformation;
% flexion_transform.timestamp = igtlTimestampNow()-startTime;
% igtlSendTransform(igtlConnection, flexion_transform);
% 
% gimbal_transform.matrix = gimbal_transformation;
% gimbal_transform.timestamp = igtlTimestampNow()-startTime;
% igtlSendTransform(igtlConnection, gimbal_transform);
% 
% platen_transform.matrix = platen_transformation;
% platen_transform.timestamp = igtlTimestampNow()-startTime;
% igtlSendTransform(igtlConnection, platen_transform);
