library(boot)
library(caTools)
library(glmnet)
library(ggplot2)
library(ggpubr)

#The codes do bootstrapping in TDG positional LASSO. The plot is done by comparing
#p-value between groups.
#Author: Yuze Hou

set.seed(101)

rsq <- function (k, n_mer){
  xn_factor <- paste("x", as.character(n_mer), "_factor", sep = "")
  ratio = 1 - 1/k
  splitted_index <- sample.split(rownames(get(xn_factor)), SplitRatio=ratio)
  train_set <- subset(get(xn_factor), splitted_index == TRUE)
  test_set <- subset(get(xn_factor), splitted_index == FALSE)  
  TDG_test_set <- subset(Intensity, splitted_index == FALSE)
  TDG_train_set <- subset(Intensity, splitted_index == TRUE)
  #Here, I haven't write alpha = 1 in previous LASSO but the default is 1.
  #Here is just for clarification.
  cvfitn <- cv.glmnet(train_set, TDG_train_set, alpha = 1)

  LASSO_predictn <- predict(cvfitn, newx = test_set, s = "lambda.1se")
  #LASSO_predictn <- predict(cvfitn, newx = test_set, s = "lambda.min")
  
  R_calc_datan = data.frame(x = LASSO_predictn, y = TDG_test_set)
  predict_lmn = lm(X1 ~ y, data = R_calc_datan)
  return(summary(predict_lmn)$r.squared)
}
 
#Change the # of iteration for bootstrapping.
iteration = 20

result = matrix(0L, nrow = iteration, ncol = 6)

for (n in c(2, 3, 4, 52, 5, 6)){
  for (r in 1 : iteration){
    set.seed(101 + r + 100 * n)
    if(n == 52){
      result[r, 4] <- rsq(4, n)
    }
    else if(n < 5){
      result[r, n - 1] <- rsq(4, n)
    }
    else{
      result[r, n] <- rsq(4, n)
    }
  }
}

colnames(result) <- c("2mer", "3mer", "4mer", "5mer1", "5mer2", "6mer")

#r = boxplot(result)
#hist(result[, 5])

alternative_param = "less"

p1 = wilcox.test(result[, 1], result[, 2], alternative = alternative_param)$p.value
p2 = wilcox.test(result[, 2], result[, 3], alternative = alternative_param)$p.value
p3 = wilcox.test(result[, 3], result[, 4], alternative = alternative_param)$p.value
p4 = wilcox.test(result[, 3], result[, 5], alternative = alternative_param)$p.value
p5 = wilcox.test(result[, 4], result[, 6], alternative = alternative_param)$p.value
p6 = wilcox.test(result[, 5], result[, 6], alternative = alternative_param)$p.value

#Here, to use ggplot2 package, we have to reshape the result.

for (i in 1:5){
  if (i == 1){
    gresult = rbind(as.matrix(result[, i]), as.matrix(result[, i + 1]))
  }
  else{
    gresult = rbind(gresult, as.matrix(result[, i + 1]))
  }
}

total = iteration * 6
Nmer = matrix("NA", nrow = iteration * 6, ncol = 1)
for (i in 1:total){
  if(i <= iteration * 3) {name = paste(as.character(ceiling(i/iteration) + 1), "mer", sep="")}
  else if(i >= iteration * 3 + 1 && i <= iteration * 4) {name = "5mer1"}
  else if(i >= iteration * 4 + 1 && i <= iteration * 5) {name = "5mer2"}
  else {name = "6mer"}
  Nmer[i, 1] = name
}

boxplot_result_forggplot = as.data.frame(cbind(Nmer, gresult))
boxplot_result_forggplot$V1 = as.factor(boxplot_result_forggplot$V1)
boxplot_result_forggplot$V2 = as.numeric(paste(boxplot_result_forggplot$V2))

#ggplot(boxplot_result_forggplot, aes(x=V1, y=V2)) + geom_boxplot()

stat.test <- tibble::tribble(
  ~group1, ~group2,   ~p.adj,
  "2mer",     "3mer", signif(p1, 3),
  "3mer",     "4mer", signif(p2, 3),
  "4mer",     "5mer1", signif(p3, 3),
  "4mer",     "5mer2", signif(p4, 3),
  "5mer1",     "6mer", signif(p5, 3),
  "5mer2",     "6mer", signif(p6, 3)
)

colnames(boxplot_result_forggplot) <- c("N_mer", "R_square")

ggboxplot(boxplot_result_forggplot, x = "N_mer", y = "R_square") +
  stat_pvalue_manual(
    stat.test, 
    y.position = 0.95, step.increase = 0.1,
    label = "p.adj"
  )