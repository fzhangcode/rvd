clc;clear;close all
%h5disp('control.hdf5')
d=[0.1,0.3,1,10,100];
muControl_s= h5read('control.hdf5','/mu');
J=size(muControl_s,2);
t=zeros(5,J);
for i=1:5
    dilution=d(i);
    if dilution==0.1
        caseFile='case0_1.hdf5';
    elseif dilution==0.3
        caseFile='case0_3.hdf5';
    elseif dilution==1
        caseFile='case1_0.hdf5';
    elseif dilution==10
        caseFile='case10_0.hdf5';
    else
        caseFile='case100_0.hdf5';
    end

    muCase_s= h5read(caseFile,'/mu');


    %% plot histogram
    postion1=50;
    postion2=245;
    figure,
    subplot(2,1,1),    
    [nb,xb]=hist(muCase_s(:,postion1));
    bh=bar(xb,nb,'hist');
    set(bh,'facecolor',[1 0 0]);
    hold on
    [nb,xb]=hist(muControl_s(:,postion1));
    bh=bar(xb,nb,'hist');
    set(bh,'facecolor',[0 0 1]);
    legend('Case','Control')
    title(['Histogram of mu when ', 'dilution=',num2str(dilution),' position=',num2str(postion1)])
    ylabel('t')
    xlabel('mu')
    hold off
    
    subplot(2,1,2),    
    [nb,xb]=hist(muCase_s(:,postion2));
    bh=bar(xb,nb,'hist');
    set(bh,'facecolor',[1 0 0]);
    hold on
    [nb,xb]=hist(muControl_s(:,postion2));
    bh=bar(xb,nb,'hist');
    set(bh,'facecolor',[0 0 1]);
    legend('Case','Control')     
    title(['Histogram of mu when ', 'dilution=',num2str(dilution),' position=',num2str(postion2)])
    xlabel('mu')
    ylabel('t')
    hold off

%% h = ttest2(x,y) returns a test decision for the null hypothesis that
%%the data in vectors xand y comes from independent random samples from 
%%normal distributions with equal means and equal but unknown 
%%variances, using the two-sample t-test.
    [h,p,~,stats] = ttest2(muCase_s,muControl_s,0.01,'right','equal',1);
    
    t(i,:)=stats.tstat;
end

inx_pos=1:J;
inx_mpos=85:20:345;
figure,

for i=1:5
    subplot(3,2,i)
    stem(inx_pos,t(i,:));
    hold on
    plot(inx_mpos,t(i,inx_mpos),'linestyle','none','marker','*','color','r')
    title(['dilution=',num2str(d(i))])
    xlabel('position')
    ylabel('t')
end

