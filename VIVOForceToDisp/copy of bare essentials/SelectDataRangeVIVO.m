FileName = 'N:\VLOptTest3\Tension_Disp_Q3_0_00003.txt';


TempFile = readtable(FileName);

FinalTable = TempFile(2000:5000, :);


txtFileName = "ShortenedcutTension.txt";

writetable(FinalTable, txtFileName, 'WriteVariableNames', true, 'Delimiter', '\t');





