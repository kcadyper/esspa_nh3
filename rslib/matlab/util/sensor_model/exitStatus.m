function exitStatus(code)
% function to return exit code. It is designed to work with
% upper-level script (e.g., python).
%
% History: 
%   May 20, 2005, First Version.................J. Galantowicz, AER
%
% Copyright, 2005, AER, Inc. All rights reserved.
%-
if (isempty(code)) code = 0; end
if (code == 0)
  disp('<PYIEXE:OK>');
else
  disp(['<PYIEXE ERROR:' num2str(code) '>']);
end
return