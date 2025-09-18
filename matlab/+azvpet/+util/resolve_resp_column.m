function col = resolve_resp_column(resp, T, opts)
% Vrátí skutečný název sloupce v tabulce T pro "hezký" název resp.
% Preferuje opts.name_map z ensure_valid_varnames.

    resp = string(resp);
    col  = resp;  % fallback

    % Pokud máme name_map
    if nargin >= 3 && isfield(opts,'name_map') && ~isempty(opts.name_map)
        nm = opts.name_map;
        if isa(nm,'containers.Map')
            if isKey(nm, char(resp))
                col = string(nm(char(resp)));
            end
        elseif isstruct(nm)
            if isfield(nm, char(resp))
                col = string(nm.(char(resp)));
            end
        end
    end

    % Pokud col pořád není ve sloupcích T, zkus safe variantu
    vnames = string(T.Properties.VariableNames);
    if ~ismember(col, vnames)
        safe = string(azvpet.util.safe_resp_name(resp));
        if ismember(safe, vnames)
            col = safe;
        end
    end
end
