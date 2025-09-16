function lbl = pretty_label(resp)
    lbl = strrep(char(resp), '_','-');
    lbl = strrep(lbl, '^',' ');
end
