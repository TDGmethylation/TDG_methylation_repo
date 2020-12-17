%% This code is to finalize the TDG methylation on liver right lobe & HepG2 cell lines.

MutS_MM = zeros(1, 1)
MutS_WT = zeros(1, 1)
hist(MutS_MM, 30, 'auto')
hold on
hist(MutS_WT, 30,  'auto')
hold off


Data = zeros(1, 4)

CpN_bins_mean = zeros(32, 5);
CpN_bins_median = zeros(32, 5);

j = 1

for i = 1:32:1024
    CpN_bins_mean(j, 2) = mean(Data(i:(i+31), 1));
    CpN_bins_mean(j, 3) = mean(Data(i:(i+31), 2));
    CpN_bins_mean(j, 4) = mean(Data(i:(i+31), 3));
    CpN_bins_mean(j, 5) = mean(Data(i:(i+31), 4));
    CpN_bins_mean(j, 1) = j
    
    CpN_bins_median(j, 2) = median(Data(i:(i+31), 1));
    CpN_bins_median(j, 3) = median(Data(i:(i+31), 2));
    CpN_bins_median(j, 4) = median(Data(i:(i+31), 3));
    CpN_bins_median(j, 5) = median(Data(i:(i+31), 4));
    CpN_bins_median(j, 1) = j

    j = j + 1
end

coef_CpN_mean = zeros(4, 1)
coef_CpN_median = zeros(4, 1)

%%Regress TDG_FI on #bin

for k = 2:5
    R_CpN_mean = corrcoef(CpN_bins_mean(:, k), CpN_bins_mean(:, 1));
    coef_CpN_mean((k - 1), 1) = R_CpN_mean(1, 2);

    R_CpN_median = corrcoef(CpN_bins_median(:, k), CpN_bins_median(:, 1));
    coef_CpN_median((k - 1), 1) = R_CpN_median(1, 2);
end