function plot_testMatrix_com()

X_vector = [0:0.15:0.75 [0:0.15:0.75] [0:0.15:0.75] [0:0.15:0.75] [0:0.15:0.75] [0:0.15:0.75]];
Y_vector = [0.25*ones(1,6) 0.5*ones(1,6) 0.75*ones(1,6) 1*ones(1,6) 1.25*ones(1,6) 1.5*ones(1,6)];
outputPath = '/Users/haining/VIV/src/2012SpringSmalltank/outputJason/';

figure
plot(X_vector, Y_vector, 'o')
xlim([min(0,0.9*min(X_vector)-0.1) 1.1*max(X_vector)])
ylim([min(0,0.9*min(Y_vector)-0.1) 1.1*max(Y_vector)])

hold on

Y_vector = [0.15 0.15 0.25 0.25 0.25 0.5 0.5 0.5 0.5 0.75 0.75 0.75 0.75 0.75 0.75 1 1 1 1 1 1];
X_vector = [0.05 0.1 0.05 0.1 0.15 0.05 0.1 0.15 0.2 0.05 0.1 0.15 0.2 0.25 0.3 0.05 0.1 0.15 0.2 0.25 0.3 ];
plot(X_vector, Y_vector, '*r')

ylabel('Y/D')
xlabel('X/D')


set(gca, 'XTick',[0:0.15:0.75])
set(gca, 'YTick',[0.25:0.25:1.5])

legend('Original Test Points', 'New Test Points')
title('Test Matrix')
saveas(gcf, [outputPath 'testMatrixComp.jpg'])


