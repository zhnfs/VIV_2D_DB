

theta = [0:0.01:2*pi];

y = sin(theta);
x = cos(theta);

phase_JMD = [-180 -135 -90 -45 0 45 90 135 180]./180.*pi;

% phase_HZ = [-180:22.5:-45, -30:15:30, 45:22.5:180]./180.*pi;
phase_HZ =  [-180:22.5:-135, -120:15:-45, -22.5:22.5:180]./180.*pi;


figure
plot(x,y,'-')
axis equal
hold on;
plot(cos(phase_JMD), sin(phase_JMD),'og','MarkerSize',24)
hold on
plot(cos(phase_HZ), sin(phase_HZ),'+r','MarkerSize',24)