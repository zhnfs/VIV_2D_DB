clear all;

outputname = 'output_';
% numruns = 36;
numruns = 17;
numiters = 1;
k = 1;

for i = 1:numruns
    for j = 1:numiters
        files{k} = [outputname num2str(i+10) '_' num2str(j) '_' num2str(k) '.mat'];
        k = k+1;
    end
end

% allzerodiff = zeros(36,6);
allzerodiff = zeros(17,6);

for i = 1:numruns
    eval([ 'load ' char(files(i)) ';']);
    allzerodiff(i,:) = zerodiff;
end

figure(4);clf;
plot(allzerodiff);
xlabel('Run Number');
ylabel('Voltage Diff.');
legend('Fx','Fy','Fz','Mx','My','Mz');
hold on;
plot([1 36],[.004 .004], 'k');
set(gcf,'Position',[0 0 1024 695]);

% clear A Yad files names t_long B forcehz numiters vel C calmat 
% clear gatherperiod numruns vel_long CD_long channame gathertime 
% clear outputname wn1 CL_long cofm h p wn2 D colnames hexfilename
% clear pathname xpos_long F_CD cylD i probruns ypos_long
% clear F_CL cylaccx_long j psd_CD zerodiff F_x cylaccy_long k psd_CL
% clear zeroind1 F_y cylvelx_long lengthtosensor psd_x zeroind2 N1
% clear cylvely_long m psd_y N2 endzero1 mass servocycle Xad endzero2
% clear momarm span
% 
% save allzerodiff.mat
