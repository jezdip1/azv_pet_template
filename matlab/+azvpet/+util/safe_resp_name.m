function s = safe_resp_name(resp)
% Convert table variable name to a valid identifier for formulas.
% Keep suffix; replace '-' with '_' etc.
    s = matlab.lang.makeValidName(char(resp), 'ReplacementStyle', 'delete');
    % makeValidName z 'I-V' udělá 'IV', což je OK – hlavně ať je to validní
end
