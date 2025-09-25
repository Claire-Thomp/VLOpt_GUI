function [WorkingLig] = LigamentResultantForces(intact, cutLig)
 %Find the delta forces in 3 DoF

    DeltaFx = intact.Fx-cutLig.Fx; % leave directions? Ask Stewart
    DeltaFy = intact.Fy-cutLig.Fy;
    DeltaFz = intact.Fz-cutLig.Fz; 
    Resultant = sqrt(DeltaFx.^2 + DeltaFy.^2 + DeltaFz.^2);

    DeltaMx = intact.Mx-cutLig.Mx;
    DeltaMy = intact.My-cutLig.My;
    DeltaMz = intact.Mz-cutLig.Mz;
 
    %Add resultant forces to working table
    WorkingLig = table(DeltaFx,DeltaFy,DeltaFz,Resultant,DeltaMx,DeltaMy,DeltaMz,'VariableNames', {'DeltaFx', 'DeltaFy', 'DeltaFz','Resultant','DeltaMx','DeltaMy','DeltaMz'});
end
