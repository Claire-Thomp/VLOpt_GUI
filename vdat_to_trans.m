function [fem_transformation, tib_transformation]=vdat_to_trans(vdat)


pose=vdat(1:6,1);
Roffset=vdat(4:6,2)+[20 0 0]';
T=[0 0 0];
[T_fem,T_tib]=VIVO_Positions_rotmatrix(pose(4),pose(5),pose(6),pose(2),pose(1),pose(3),Roffset,T,'r');

fem_trans=T_fem(1:3,1);
fem_rot=T_fem(1:3,2:4)-fem_trans;
fem_transformation=[fem_rot fem_trans];

tib_trans=T_tib(1:3,1);
tib_rot=T_tib(1:3,2:4)-tib_trans;
tib_transformation=[tib_rot tib_trans];


