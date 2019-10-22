%PLOT 2_2

x2_2 = [0 1 2 3 4 5 6 7 8 9 10];

y2_2_CS_1 = [1991.865 1978.263 1975.783 1948.077 1977.055 1982.428 ...
    1949.393 1984.967 1985.964 1944.201 1988.437]; %artdata0.5
y2_2_NMI_1 = [.322 .3 .3 .302 .308 .324 .311 .304 .3 .3 .316]; %artdata0.5

figure
hold on
title('CS and NMI for Artdata0.5 Dataset')
xlabel('Run Index') % x-axis label

yyaxis left
plot(x2_2,y2_2_CS_1);
ylim([0 2000])
ylabel('Cluster Scatter') % y-axis label

yyaxis right
plot(x2_2,y2_2_NMI_1,'--');
ylim([0 1])
ylabel('NMI') % y-axis label
legend('CS','NMI','Location','southeast')

hold off

y2_2_CS_2 = [2164.341 2149.963 2148.585 2166.38 2185.879 2205.834 ...
    2150.606 2147.777 2177.559 2148.205 2196.802]; %artdata1
y2_2_NMI_2 = [.612 .653 .638 .582 .562 .571 .602 .605 .62 .646 .568]; ...
    %artdata1

figure
hold on
title('CS and NMI for Artdata1 Dataset')
xlabel('Run Index') % x-axis label

yyaxis left
plot(x2_2,y2_2_CS_2);
ylim([0 2500])
ylabel('Cluster Scatter') % y-axis label

yyaxis right
plot(x2_2,y2_2_NMI_2,'--');
ylim([0 1])
ylabel('NMI') % y-axis label
legend('CS','NMI','Location','southeast')

hold off

y2_2_CS_3 = [2457.405 3232.451 2457.405 3226.773 2457.405 2457.405 ... 
    2457.405 3237.746 3237.746 2457.423 2457.423]; %artdata2
y2_2_NMI_3 = [.948 .751 .948 .757 .948 .948 .948 .754 .754 .941 ...
    .941]; %artdata2

figure
hold on
title('CS and NMI for Artdata2 Dataset')
xlabel('Run Index') % x-axis label

yyaxis left
plot(x2_2,y2_2_CS_3);
ylim([0 3500])
ylabel('Cluster Scatter') % y-axis label

yyaxis right
plot(x2_2,y2_2_NMI_3,'--');
ylim([0 1])
ylabel('NMI') % y-axis label
legend('CS','NMI','Location','southeast')

hold off

y2_2_CS_4 = [4532.358 4573.328 4558.953 2530.738 4557.628 2530.738 ...
    2530.738 4554.471 2530.738 2530.738 2530.738]; %artdata3
y2_2_NMI_4 = [.847 .806 .809 1.0 .809 1.0 1.0 .823 1.0 1.0 1.0]; %artdata3

figure
hold on
title('CS and NMI for Artdata3 Dataset')
xlabel('Run Index') % x-axis label

yyaxis left
plot(x2_2,y2_2_CS_4);
ylim([0 5000])
ylabel('Cluster Scatter') % y-axis label

yyaxis right
plot(x2_2,y2_2_NMI_4,'--');
ylim([0 1])
ylabel('NMI') % y-axis label
legend('CS','NMI','Location','southeast')

hold off

y2_2_CS_5 = [6360.72 2502.347 6360.7 6204.969 6399.995 6375.145 ...
    6760.896 2502.347 6206.851 6206.789 2502.347]; %artdata4
y2_2_NMI_5 = [.86 1.0 .86 .852 .859 .852 .79 1.0 .849 .849 1.0]; %artdata4

figure
hold on
title('CS and NMI for Artdata4 Dataset')
xlabel('Run Index') % x-axis label

yyaxis left
plot(x2_2,y2_2_CS_5);
ylim([0 7000])
ylabel('Cluster Scatter') % y-axis label

yyaxis right
plot(x2_2,y2_2_NMI_5,'--');
ylim([0 1])
ylabel('NMI') % y-axis label
legend('CS','NMI','Location','southeast')

hold off

y2_2_CS_6 = [9060.096 9060.096 9060.096 9060.096 9060.133 9060.096 ...
    9060.096 9060.133 9060.133 9060.096 9060.096]; %ionosphere
y2_2_NMI_6 = [.138 .138 .138 .138 .142 .138 .138 .142 .142 .138 .138]; ...
    %ionosphere

figure
hold on
title('CS and NMI for Ionosphere Dataset')
xlabel('Run Index') % x-axis label

yyaxis left
plot(x2_2,y2_2_CS_6);
ylim([0 10000])
ylabel('Cluster Scatter') % y-axis label

yyaxis right
plot(x2_2,y2_2_NMI_6,'--');
ylim([0 1])
ylabel('NMI') % y-axis label
legend('CS','NMI','Location','southeast')

hold off

y2_2_CS_7 = [140.279 140.213 140.279 140.279 140.213 140.213 140.026 ...
    190.757 140.213 140.213 140.029 ]; %iris
y2_2_NMI_7 = [0.76 0.733 0.76 0.76 0.733 0.733 0.748 0.689 0.733 0.733 ...
    0.741]; %iris

figure
hold on
title('CS and NMI for Iris Dataset')
xlabel('Run Index') % x-axis label

yyaxis left
plot(x2_2,y2_2_CS_7);
ylim([0 200])
ylabel('Cluster Scatter') % y-axis label

yyaxis right
plot(x2_2,y2_2_NMI_7,'--');
ylim([0 1])
ylabel('NMI') % y-axis label
legend('CS','NMI','Location','southeast')

hold off

y2_2_CS_8 = [3810.572 3746.932 3797.343 3864.68 3874.768 3806.134 ...
    3792.325 3780.988 3814.0 3847.293 3763.695 ]; %soybean-processed
y2_2_NMI_8 = [0.921 0.944 0.958 0.961 0.945 0.987 0.938 0.965 0.966 ...
    0.985 1.0]; %soybean-processed

figure
hold on
title('CS and NMI for Soybean-Processed Dataset')
xlabel('Run Index') % x-axis label

yyaxis left
plot(x2_2,y2_2_CS_8);
ylim([0 4000])
ylabel('Cluster Scatter') % y-axis label

yyaxis right
plot(x2_2,y2_2_NMI_8,'--');
ylim([0 1])
ylabel('NMI') % y-axis label
legend('CS','NMI','Location','southeast')

hold off
 
%PLOT 2_3

fileID1 = fopen('artdata0.5output23.txt','r');
fileID2 = fopen('artdata1output23.txt','r');
fileID3 = fopen('artdata2output23.txt','r');
fileID4 = fopen('artdata3output23.txt','r');
fileID5 = fopen('artdata4output23.txt','r');
fileID6 = fopen('ionosphereoutput23.txt','r');
fileID7 = fopen('irisoutput23.txt','r');
fileID8 = fopen('soybean-processedoutput23.txt','r');

formatSpec = '%d %f';
sizeA = [2 Inf];

A1 = fscanf(fileID1,formatSpec,sizeA);
A2 = fscanf(fileID2,formatSpec,sizeA);
A3 = fscanf(fileID3,formatSpec,sizeA);
A4 = fscanf(fileID4,formatSpec,sizeA);
A5 = fscanf(fileID5,formatSpec,sizeA);
A6 = fscanf(fileID6,formatSpec,sizeA);
A7 = fscanf(fileID7,formatSpec,sizeA);
A8 = fscanf(fileID8,formatSpec,sizeA);

fclose(fileID1);
fclose(fileID2);
fclose(fileID3);
fclose(fileID4);
fclose(fileID5);
fclose(fileID6);
fclose(fileID7);
fclose(fileID8);

A1 = A1';
A2 = A2';
A3 = A3';
A4 = A4';
A5 = A5';
A6 = A6';
A7 = A7';
A8 = A8';

y2_3_CS_1 = A1(:,2); %artdata0.5
y2_3_CS_2 = A2(:,2); %artdata1
y2_3_CS_3 = A3(:,2); %artdata2
y2_3_CS_4 = A4(:,2); %artdata3
y2_3_CS_5 = A5(:,2); %artdata4
y2_3_CS_6 = A6(:,2); %ionosphere
y2_3_CS_7 = A7(:,2); %iris
y2_3_CS_8 = A8(:,2); %soybean-processed

x2_3 = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22];

figure
hold on
plot(x2_3,y2_3_CS_1);
title('CS for Different K Values for Artdata0.5 Dataset')
xlabel('k') % x-axis label
ylabel('Cluster Scatter') % y-axis label
hold off

figure
hold on
plot(x2_3,y2_3_CS_2);
title('CS for Different K Values for Artdata1 Dataset')
xlabel('k') % x-axis label
ylabel('Cluster Scatter') % y-axis label
hold off

figure
hold on
plot(x2_3,y2_3_CS_3);
title('CS for Different K Values for Artdata2 Dataset')
xlabel('k') % x-axis label
ylabel('Cluster Scatter') % y-axis label
hold off

figure
hold on
plot(x2_3,y2_3_CS_4);
title('CS for Different K Values for Artdata3 Dataset')
xlabel('k') % x-axis label
ylabel('Cluster Scatter') % y-axis label
hold off

figure
hold on
plot(x2_3,y2_3_CS_5);
title('CS for Different K Values for Artdata4 Dataset')
xlabel('k') % x-axis label
ylabel('Cluster Scatter') % y-axis label
hold off

figure
hold on
plot(x2_3,y2_3_CS_6);
title('CS for Different K Values for Ionosphere Dataset')
xlabel('k') % x-axis label
ylabel('Cluster Scatter') % y-axis label
hold off

figure
hold on
plot(x2_3,y2_3_CS_7);
title('CS for Different K Values for Iris Dataset')
xlabel('k') % x-axis label
ylabel('Cluster Scatter') % y-axis label
hold off

figure
hold on
plot(x2_3,y2_3_CS_8);
title('CS for Different K Values for Soybean-Processed Dataset')
xlabel('k') % x-axis label
ylabel('Cluster Scatter') % y-axis label
hold off

