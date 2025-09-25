%VIVO Coordinates (aka vivo pose) 
testc = {[10,0,0,20,0,0]}; %ML, AP, VL, FE, ABD, IE
%Ligament parameters as inserted into the VIVO 
ligmdl = {[0,0,0,0,-10,10,50,50]}; % Px, PY, PZ, DX, DY, DZ, K, RefStrain
%reference pose 
refpos = [0,0,0,0,0,0];

%%number of ligs 
nligs = length(ligmdl);
%number of poses
ntc = length(testc);
ligmdlN = cell(size(nligs));
%scale the ligament parameters to meters and decimals 
ligscale = [.001 .001 .001 .001 .001 .001 1 .01];
%scale each lig and find its slack length 
for ilg = 1:nligs
	ligmdlN{ilg} = ligmdl{ilg} .* ligscale ;  % normalize units
	ligmdlN{ilg} = [ligmdlN{ilg}, cptL0(ligmdlN{ilg})]; % add L0 as 9th element
end

function L0 = cptL0(lig)
    % return L0 given ligament insertion site positions at reference and reference strain
    P = lig(1:3);
    D = lig(4:6);
    dlig = P - D;
    lrf = sqrt(dlig * dlig');
    L0 = lrf/(lig(8)+1);  %lig(8) = reference strain
end

testcN = cell(size(testc));
%scale the coordinates to Metres and Radians 
coordScale = [.001 .001 .001  pi/180 pi/180 pi/180];
% add in 90 degree offset 
coordOffs = [0 0 0 0 pi/2 0]; 
for itc = 1:ntc
	testcN{itc} = (testc{itc} .* coordScale ) + coordOffs; % normalize units and add G&S beta offset
end
%scale the ref pose and add the offset 
refposN = (refpos .* coordScale) + coordOffs;

% Remap ligament distal sites to T-referenced coordinates using reference
% pose --> Only needed if the refpose is not all zeros
% refcell = num2cell(refposN);
% tZfRef = tZfgs(refcell{:});
% for ilg = 1:nligs
% 	dlig = [ 1; (ligmdlN{ilg}(4:6))']; % convert to column
% 	dlig = tZfRef * dlig;
% 	ligmdlN{ilg}(4:6) = (dlig(2:4))'; % back to row
% end

% Preallocate output arrays
LigF = cell(ntc, nligs);
NetF = cell(ntc, 1);

%this assumes chirality is R 
rotsign = [-1 -1 1];
%% 

%Uses Ryan's VIVO transform streaming and lig code 
for itc = 1:ntc
    pose = num2cell(testc{itc} .* coordScale);
    Roffset = [0; 0; 0] + [20 0 0]';
    T = [0 0 0];
    [T_fem, T_tib] = VIVO_Positions_rotmatrix(pose{4}, pose{5}, pose{6}, pose{2}, pose{1}, pose{3},Roffset, T, 'r');

    fem_trans = T_fem(1:3, 1);
    fem_rot = T_fem(1:3, 2:4) - fem_trans;
    fem_transformation = [fem_rot fem_trans];

    tib_trans = T_tib(1:3, 1);
    tib_rot = T_tib(1:3, 2:4) - tib_trans;
    tib_transformation = [tib_rot tib_trans];

    %[flexion_transformation, gimbal_transformation, platen_transformation] = ActuatorOffsets(vdat, fem_transformation, tib_transformation);

    % transformation matrices 
    trans_results{1} = fem_transformation;
    trans_results{2} = tib_transformation;
    %trans_results{3} = flexion_transformation;
    %trans_results{4} = gimbal_transformation;
    %trans_results{5} = platen_transformation;

    % Compute current ligament insertion positions
    F = ((trans_results{1}(1:3, 1:3)) * (ligmdlN{1}(1:3))' + (trans_results{1}(1:3, 4)))';
    T = ((trans_results{2}(1:3, 1:3)) * (ligmdlN{1}(4:6))' + (trans_results{2}(1:3, 4)))';

    %Inverted transforms 
    F_tib = inv([trans_results{2}; 0 0 0 1])*[F';1];
    T_tib = inv([trans_results{2}; 0 0 0 1])*[T';1];
    Ft = (F_tib(1:3, 1))';
    Tt = (T_tib(1:3, 1))';
    L0 = ligmdlN{1}(9);
    K = ligmdlN{1}(7);
    D = ligmdlN{1}(4:6);
    % Compute ligament vectors and lengths
    dlig = Ft - Tt;
    dlig =dlig.*rotsign;
    llig = sqrt( dlig * dlig'); % Compute lengths along rows
    if llig < L0   % detect ligament slack
	    fT = zeros(1,6);
	    return
    end
    nlig = dlig / llig; % unit vector aligned with distal to proximal line
		    %; gives direction of ligament force wrt T

    %compute force using Blankevoort model
    e = (llig-L0) / L0;
    et = .03;
    if e < 2*et
	    flig = ( K * .25 * e * e / et) * nlig;
    else
	    flig = ( K * (e-et) ) * nlig;
    end

    mlig = cross((rotsign.*Tt), flig); % moment wrt T origin
    fT = [flig, (rotsign.*mlig)];
   % tforce = LForce(ligmdlN{ilig}, trans_results{2}.*pose{:});
end
%% 

%Uses VIVO internal code in addition to Ryans VIVO streaming code (use
%either this or the above lines) 
for itc = 1:ntc %loop over all supplied test poses
    pose = num2cell(testc{itc} .* coordScale);
    Roffset = [0 0 0] + [20 0 0]';
    T = [0 0 0];
    [T_fem, T_tib] = VIVO_Positions_rotmatrix(pose{4}, pose{5}, pose{6}, pose{2}, pose{1}, pose{3},Roffset, T, 'r');

    fem_trans = T_fem(1:3, 1);
    fem_rot = T_fem(1:3, 2:4) - fem_trans;
    fem_transformation = [fem_rot fem_trans];

    tib_trans = T_tib(1:3, 1);
    tib_rot = T_tib(1:3, 2:4) - tib_trans;
    tib_transformation = [tib_rot tib_trans];

    %[flexion_transformation, gimbal_transformation, platen_transformation] = ActuatorOffsets(vdat, fem_transformation, tib_transformation);

    trans_results{1} = fem_transformation;
    trans_results{2} = tib_transformation;
    %trans_results{3} = flexion_transformation;
    %trans_results{4} = gimbal_transformation;
    %trans_results{5} = platen_transformation;

    % Compute current ligament insertion positions
    F = ((trans_results{1}(1:3, 1:3)) * (ligmdlN{1}(1:3))' + (trans_results{1}(1:3, 4)))';
    T = ((trans_results{2}(1:3, 1:3)) * (ligmdlN{1}(4:6))' + (trans_results{2}(1:3, 4)))';
    F_tib = inv([trans_results{2}; 0 0 0 1])*[F';1];
    T_tib = inv([trans_results{2}; 0 0 0 1])*[T';1];
	Ft = (F_tib(1:3, 1))';
    Tt = (T_tib(1:3, 1))';

%pull the pose from the input list.

	% T-to-Q coordinate conversion (det(∑) ~= 1; not a rotation)
	[beta, gamma] = pose{5:6};
	qXt = [	cos(gamma)*sin(beta)		sin(beta)*sin(gamma)	cos(beta);
			-sin(gamma)					cos(gamma)				0;
			0							0						-1];
	NetF{itc} = zeros(1,6); % init resultant

    % compute forces referred to G-S basis vectors for each ligament fiber
    for ilig = 1:nligs
		tforce = LForce(ligmdlN{ilig}, Ft, Tt); %force and moment referred to T
			% Constraint forces are multiplied by -1 for reporting in the VivoControl UI.
			% Reproduce that here so that results look consistent.
		LigF{itc, ilig} = -1*[(qXt * tforce(1:3)')', (qXt * tforce(4:6)')' .* rotsign]; % compute projections along G-S bases
		NetF{itc} = NetF{itc} + LigF{itc, ilig}; % accumulate in resultant
    end
end


function fT = LForce(lig, F, T)
    % lig is a [1 x 9] specifying the ligament: [Px Py Pz Dx Dy Dz K eRf L0]
    % The distal insertion site has already been statically transformed to the T frame
    % The proximal point is in F and must be converted to T using the dynamic xform passed as tZf.

    P = lig(1:3); % proximal (in F)
    D = lig(4:6); % distal (in T)
    K = lig(7);
    %eRf = lig(8) % unused
    L0 = lig(9);


    dlig = F-T; % vector from distal to proximal
    llig = sqrt( dlig * dlig'); % current length of ligament
    if llig < L0   % detect ligament slack
	    fT = zeros(1,6);
	    return
    end
    nlig = dlig / llig; % unit vector aligned with distal to proximal line
		    %; gives direction of ligament force wrt T

    %compute force using Blankevoort model
    e = (llig-L0) / L0;
    et = .03;
    if e < 2*et
	    flig = ( K * .25 * e * e / et) * nlig;
    else
	    flig = ( K * (e-et) ) * nlig;
    end
    mlig = cross(D, flig); % moment wrt T origin
    fT = [flig, mlig];
end