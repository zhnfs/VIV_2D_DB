% Make contour plot of various lift and drag coefficients
% Jason Dahl, 6/29/07

% phase order {'vr6ph0' 'vr6ph135' 'vr6ph135n' 'vr6ph180n' 'vr6ph45'
% 'vr6ph45n' 'vr6ph90' 'vr6ph90n'};

clear all

phase = [0 135 -135 -180 45 -45 90 -90];

Yad = [0.5 0.5 0.5 0.5 0.5 0.75 0.75 0.75 0.75 0.75 1.0 1.0 1.0 1.0 1.0 ...
    1.25 1.25 1.25 1.25 1.25 1.5 1.5 1.5 1.5 1.5];

Xad = [0 0.15 0.3 0.45 0.6 0 0.15 0.3 0.45 0.6 0 0.15 0.3 0.45 0.6 ...
    0 0.15 0.3 0.45 0.6 0 0.15 0.3 0.45 0.6];

