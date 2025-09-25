function [intact, cut] = FileLengthDiscrep(intact, cut)
    % This repeats the last row of the shorter file until both are the same length
    len_intact = height(intact);
    len_cut = height(cut);
    
    if len_intact > len_cut
        rows_to_add = len_intact - len_cut;
        last_row = cut(end, :);
        repeated_rows = repmat(last_row, rows_to_add, 1);
        cut = [cut; repeated_rows];
    elseif len_cut > len_intact
        rows_to_add = len_cut - len_intact;
        last_row = intact(end, :);
        repeated_rows = repmat(last_row, rows_to_add, 1);
        intact = [intact; repeated_rows];
    end
end
