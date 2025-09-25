function [T_fem,T_tib]=VIVO_Positions_rotmatrix(FE,AD,IE,AP,ML,VL,Roffset,T,chirality)

points_T=[[0 0 0]' [1 0 0]' [0 1 0]' [0 0 1]' [0 0 0]' [0 0 1]'];
points_F=[[0 0 0]' [1 0 0]' [0 1 0]' [0 0 1]'];

T_lower=[T(1) T(2) -T(3)];


% Calculate the alpha beta gammas. 
if chirality=='r'
    gamma=(IE);
else
    gamma=-(IE);
end

beta=(AD+(pi()/2));
alpha=(FE);

% Calculate translations in terms of S (grood-suntay directions) and then H
% (tibia wrt femur)
q1=ML;
q2=AP;
q3=VL;

S3=(q1*cos(beta)+q3)/((cos(beta)*cos(beta))-1);
S1=((-1*S3)-q3)/(cos(beta));
S2=q2;

S=[S1;S2;S3];

U=[1 0 cos(beta);
    0 cos(alpha) sin(alpha)*sin(beta);
    0 -1*sin(alpha) cos(alpha)*sin(beta)];

H=U*S;

% Calculate the transposed rotation, R^T, givin in appendix of grood-suntay
R_trans_1_1=cos(gamma)*sin(beta);
R_trans_1_2=sin(gamma)*sin(beta);
R_trans_1_3=cos(beta);

R_trans_2_1=-cos(alpha)*sin(gamma)-cos(gamma)*sin(alpha)*cos(beta);
R_trans_2_2=cos(alpha)*cos(gamma)-sin(gamma)*sin(alpha)*cos(beta);
R_trans_2_3=sin(beta)*sin(alpha);

R_trans_3_1=sin(alpha)*sin(gamma)-cos(gamma)*cos(alpha)*cos(beta);
R_trans_3_2=-cos(gamma)*sin(alpha)-cos(alpha)*cos(beta)*sin(gamma);
R_trans_3_3=cos(alpha)*sin(beta);

R_trans=[R_trans_1_1 R_trans_2_1 R_trans_3_1;
    R_trans_1_2 R_trans_2_2 R_trans_3_2;
    R_trans_1_3 R_trans_2_3 R_trans_3_3];

% Calculate the rotation / translation of the tibia wrt to femur. The
% equation is given as R=[R]r+H where r is coordinates in tibial frame and
% R is coordinates in femoral frame, R is the matrix we just calculated the
% transpose of above, and H is the transformation vector from femur to
% tibia
R=transpose(R_trans);

points_T2=R*points_T;
points_T2=points_T2+H;

%disp(points_T2)
% this option should be "as on VIVO in real life"...
[points_F]=rot_2_vivo3(points_F,FE,AD,Roffset(1));
[points_T3]=rot_2_vivo3(points_T2,FE,AD,Roffset(1));
    
T_fem=points_F;
T_tib=points_T3;

    % end
