function [points]=rot_2_vivo3(points,FE,AD,rot_offset)

pivot=[0;1;0];
offset=-20+rot_offset;
offset_rad=offset*pi()/180;
R0=[1 0 0;
    0 cos(offset_rad) -1*sin(offset_rad);
    0 sin(offset_rad) cos(offset_rad)];
pivot=R0*pivot;

%FE_rad=FE*pi()/180;
% AD_rad=AD*pi()/180;

R1=[1 0 0;
    0 cos(FE) -1*sin(FE);
    0 sin(FE) cos(FE)];

R2=[cos(AD) 0 sin(AD);
    0 1 0;
    -1*sin(AD) 0 cos(AD)];

points=R1*points;
points=R2*points;
pivot=R1*pivot;
pivot=R2*pivot;

% find how far pivot has rotated from [0;1;0] if projected onto the table;
pivot_projected=(pivot.*[1;1;0])/norm(pivot.*[1;1;0]);
extra_ie_rot=acos(dot(pivot_projected,[0;1;0]));

if pivot_projected(1)<0
    extra_ie_rot=-1*extra_ie_rot;
end

R3=[cos(extra_ie_rot) -1*sin(extra_ie_rot) 0;
    sin(extra_ie_rot) cos(extra_ie_rot) 0;
    0 0 1];

points=R3*points;
