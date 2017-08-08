function [t,f,phase] = instf(y,t)

xt = hilbert(y);
phase = unwrap(angle(xt));
f = gradient(phase,t)/(2*pi);
phase = (angle(xt));

% if nargout == 2
%     varargout(1) = {f};
% elseif nargout ==3
%     varargout(1) = {f};
%     varargout(2) = {phase};
% else
%     disp('no output arguments assigned');
% end



