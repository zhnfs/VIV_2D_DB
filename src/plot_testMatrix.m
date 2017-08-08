function plot_testMatrix(X_vector, Y_vector, outputPath)

figure
plot(X_vector, Y_vector, 'o')
ylabel('Y/D')
xlabel('X/D')
xlim([0 1.1*max(X_vector)])
ylim([0 1.1*max(Y_vector)])
legend('Test Points')
title('Test Matrix')
saveas(gcf, [outputPath 'testMatrix.jpg'])


