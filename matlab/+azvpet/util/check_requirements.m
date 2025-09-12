function check_requirements
  v = ver('MATLAB'); assert(~isempty(v), 'MATLAB not found');
  assert(~verLessThan('MATLAB','9.14'), 'Use MATLAB R2023a+ for best compatibility');
end
