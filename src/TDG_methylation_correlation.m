CpA_table = readtable("/Volumes/Genome_data/Qbic/All_Data_Result/Matlab_input/CpA.csv");
CpC_table = readtable("/Volumes/Genome_data/Qbic/All_Data_Result/Matlab_input/CpC.csv");
CpG_table = readtable("/Volumes/Genome_data/Qbic/All_Data_Result/Matlab_input/CpG.csv");
CpT_table = readtable("/Volumes/Genome_data/Qbic/All_Data_Result/Matlab_input/CpT.csv");

%Here, all the CpN files shoule be sorted by methylation level in
%ASCENDING order. Reason is that we hope to get the R to be positive when
%methylation level increase with protein binding FI increase.

CpA_TFs = table2array(CpA_table(:, 6: 671));
CpC_TFs = table2array(CpC_table(:, 6: 671));
CpG_TFs = table2array(CpG_table(:, 6: 671));
CpT_TFs = table2array(CpT_table(:, 6: 671));

CpA_bins = zeros(32, 666);
CpG_bins = zeros(32, 666);
CpC_bins = zeros(32, 666);
CpT_bins = zeros(32, 666);

j = 1
for i = 1:32:1024
    CpA_bins(j, :) = mean(CpA_TFs(i:(i+31), :));
    CpC_bins(j, :) = mean(CpC_TFs(i:(i+31), :));
    CpG_bins(j, :) = mean(CpG_TFs(i:(i+31), :));
    CpT_bins(j, :) = mean(CpT_TFs(i:(i+31), :));
    j  = j + 1
end

%CpN_x = zeros(32, 1)
%Here, for matlab regression, you need to put :1 col filled with 1, and :2
%col filled with your X. In this case, x=range(32) in reverse order

coef_CpA = zeros(3, 666);
coef_CpC = zeros(3, 666);
coef_CpG = zeros(3, 666);
coef_CpT = zeros(3, 666);

for k = 1:666
    CpA_y = CpA_bins(:, k);
    coef_CpA(1:2, k) = CpN_x\CpA_y;
    R_CpA = corrcoef(CpN_x(:, 2), CpA_y);
    coef_CpA(3, k) = R_CpA(1, 2);
    
    CpC_y = CpC_bins(:, k);
    coef_CpC(1:2, k) = CpN_x\CpC_y;
    R_CpC = corrcoef(CpN_x(:, 2), CpC_y);
    coef_CpC(3, k) = R_CpC(1, 2);
    
    CpG_y = CpG_bins(:, k);
    coef_CpG(1:2, k) = CpN_x\CpG_y;
    R_CpG = corrcoef(CpN_x(:, 2), CpG_y);
    coef_CpG(3, k) = R_CpG(1, 2);
    
    CpT_y = CpT_bins(:, k);
    coef_CpT(1:2, k) = CpN_x\CpT_y;
    R_CpT = corrcoef(CpN_x(:, 2), CpT_y);
    coef_CpT(3, k) = R_CpT(1, 2);
end

% %TDG_FI = zeros(1, 1) 
% %Here, just copy from excel.
% coef_TDG = zeros(3, 4)
% 
% for m = 1:4
%     TDG_y = log(TDG_FI(:, m))
%     coef_TDG(1:2, m) = CpN_x\TDG_y
%     R_TDG = corrcoef(CpN_x(:, 2), TDG_y)
%     coef_TDG(3, m) = R_TDG(1, 2)^2
% end

%R_boxplot = zeros(1, 1)
result  = boxplot(R_boxplot.')


count = zeros(4, 1)
A = zeros(1, 1)
C = zeros(1, 1)
G = zeros(1, 1)
T = zeros(1, 1)

for t = 1:666
    if coef_CpA(3, t) >= 0.71
        count(1, 1) = count(1, 1) + 1;
        fprintf('A%i\n', t)
    end
        
    if coef_CpC(3, t) >= 0.57
        count(2, 1) = count(2, 1) + 1;
        fprintf("C%i\n", t)
    end    
    
    if coef_CpG(3, t) >= 0.77
        count(3, 1) = count(3, 1) + 1;
        fprintf("G%i\n", t)
    end    
    
    if coef_CpT(3, t) >= 0.74
        count(4, 1) = count(4, 1) + 1;
        fprintf("T%i\n", t)
    end    
end

