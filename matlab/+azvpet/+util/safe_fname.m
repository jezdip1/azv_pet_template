function s = safe_fname(txt)
% Nahraď ne-alfa-numerické znaky podtržítkem (pro název složky/souboru)
    if ~ischar(txt) && ~isstring(txt), txt = char(string(txt)); end
    s = regexprep(txt,'[^A-Za-z0-9\-]+','_');
    if isempty(s), s = 'case'; end
end
