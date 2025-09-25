function [flexion_transformation, gimbal_transformation, platen_transformation]=ActuatorOffsets(vdat,fem_transformation,tib_transformation)

% This assumes that reference pose is all zeros (normal practice)

points_T=[[0 0 0]' [1 0 0]' [0 1 0]' [0 0 1]' [0 0 0]' [0 0 1]'];
points_F_arm=[[0 0 0]' [1 0 0]' [0 1 0]' [0 0 1]'];
points_G_arm=[[0 0 0]' [1 0 0]' [0 1 0]' [0 0 1]'];

offsets=vdat(1:6,2);

%first offset G arm
gimble_offset_rad=offsets(5)*pi()/180;
R0=[cos(gimble_offset_rad) 0 sin(gimble_offset_rad);
    0 1 0;
    -1*sin(gimble_offset_rad) 0 cos(gimble_offset_rad)];
points_G_arm=R0*points_G_arm;

%Now offset flexion arm and G arm (for flexion)
flexion_offset_rad=offsets(4)*pi()/180;
R0=[1 0 0;
    0 cos(flexion_offset_rad) -1*sin(flexion_offset_rad);
    0 sin(flexion_offset_rad) cos(flexion_offset_rad)];

points_F_arm=R0*points_F_arm;
points_G_arm=R0*points_G_arm;

%Now, we don't want the flexion arm to move the same as the gimbal (there
%is a joint), so instead we calculate the projected flexion of the gimbal
%to tease out just that part:
temp=fem_transformation(1:3,1:3)*points_F_arm(:,2:4);

y_vector=temp(:,2);
y_vector(1)=0;
flexion_angle_new=-1*(acosd(dot(y_vector,[0;0;1]))-90);

flexion_offset_rad_new=flexion_angle_new*pi()/180-flexion_offset_rad;
R0=[1 0 0;
    0 cos(flexion_offset_rad_new) -1*sin(flexion_offset_rad_new);
    0 sin(flexion_offset_rad_new) cos(flexion_offset_rad_new)];


points_F_arm=R0*points_F_arm;

% the gimbal, however, receives the normal femoral transformation (note we
% can ignore translations, cause always zero). 
points_G_arm=fem_transformation(1:3,1:3)*points_G_arm;


% for the tibia, first we rotate as per offset:
ie_offset_rad=offsets(6)*pi()/180;
R0=[cos(ie_offset_rad) -1*sin(ie_offset_rad) 0;
    sin(ie_offset_rad) cos(ie_offset_rad) 0;
    0 0 1];
points_T=R0*points_T;

% now apply translation offsets
offsets(3)=-1*offsets(3);
points_T=points_T-offsets(1:3);

% now apply actual motion transformation:
points_T=tib_transformation(1:3,1:3)*points_T;
points_T=points_T+tib_transformation(:,4);

flexion_trans=points_F_arm(1:3,1);
flexion_rot=points_F_arm(1:3,2:4)-flexion_trans;
flexion_transformation=[flexion_rot flexion_trans];

gimbal_trans=points_G_arm(1:3,1);
gimbal_rot=points_G_arm(1:3,2:4)-gimbal_trans;
gimbal_transformation=[gimbal_rot gimbal_trans];

platen_trans=points_T(1:3,1);
platen_rot=points_T(1:3,2:4)-platen_trans;
platen_transformation=[platen_rot platen_trans];
