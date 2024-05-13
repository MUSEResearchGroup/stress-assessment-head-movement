%% t-test for NASA-TLX results

%This script performes a t-test to evaluate whether the NASA-TLX scores
%obtained from the Stroop Color Word Test are significantly
%different from thos obtained from the Mental Arithmetic Test.

stroop = [36,56.6,75,30.3,45.6,58.3,69.3,59.9,57.15,54.6, 80.3,54.3,21.6,74.66,16.6,78.6,74.6,76.6,33.43,41.6];
MA = [61.3,48.1,97.3,78.2,58.3,32.6,92.3,40.6,87.3,61,59.6,67.3,47.3,81.6,35.5,76,77,89.6,60.3,65];

[h1,p1,ci1,stats1] = ttest(stroop,MA,'Tail','both','Alpha',0.05);

[h2,p2,ci2,stats2] = ttest(stroop,MA,'Tail','right','Alpha',0.05);

f=figure;
h=boxplot([stroop MA], [zeros(1,20) ones(1,20)]);

ylabel('Scores', 'FontSize', 18, 'FontName', 'Times New Roman')

set(gca, 'FontSize', 80, 'FontName', 'Times New Roman')
xticklabels(["Stroop test" "MA test"])

set(h,{'linew'},{3})