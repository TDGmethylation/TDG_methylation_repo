library(glmnet)

#Input using TDG binding data
TDG_context <- read.csv("/Users/gerogehou/Desktop/R/PBM/LoopContext_v1/chamber1_1st_use_TDG/Multivariate_done.csv")
TDG_input <- read.csv("/Users/gerogehou/Desktop/R/PBM/LoopContext_v1/chamber1_1st_use_TDG/Input.csv")

#Change categorical variables to dummies for LASSO
p0_nt1 = as.factor(c(TDG_context$X0))
p1_nt1 = as.factor(c(TDG_context$X1))
p3_nt1 = as.factor(c(TDG_context$X3))
p4_nt1 = as.factor(c(TDG_context$X4))
p5_nt1 = as.factor(c(TDG_context$X5))
p6_nt1 = as.factor(c(TDG_context$X6))

p0_nt2 = as.factor(c(TDG_context$X7))
p1_nt2 = as.factor(c(TDG_context$X8))
p2_nt2 = as.factor(c(TDG_context$X9))
p3_nt2 = as.factor(c(TDG_context$X10))
p4_nt2 = as.factor(c(TDG_context$X11))
p5_nt2 = as.factor(c(TDG_context$X12))

p0_nt3 = as.factor(c(TDG_context$X13))
p1_nt3 = as.factor(c(TDG_context$X14))
p2_nt3 = as.factor(c(TDG_context$X15))
p3_nt3 = as.factor(c(TDG_context$X16))
p4_nt3 = as.factor(c(TDG_context$X17))

Intensity <- TDG_context$Intensity

#######################N-mer_test_TDGcontext_boundary##########################
#Here, we run 2-6 mer LASSO and predict the binding. 2mer is also for Kmax to
#substantiate 2mer is not enough for modelling TDG binding. 

x1_factor <- model.matrix(Intensity ~ p3_nt1 + p2_nt2)[, -1]
LASSO_RE1 <- glmnet(x1_factor, Intensity, alpha = 1)
cvfit1 <- cv.glmnet(x1_factor, Intensity)
LASSO_predict1 <- predict(cvfit1, newx = x1_factor, s = "lambda.min")
plot(LASSO_predict1, TDG_context$Intensity)

R_calc_data1 = data.frame(x = LASSO_predict1, y = TDG_context$Intensity)
predict_lm1 = lm(X1 ~ y, data = R_calc_data1)
summary(predict_lm1)

#######################N-mer_test_TDGcontext_boundary_end######################

x2_factor <- model.matrix(Intensity ~ p3_nt1 + p4_nt1 + p2_nt2 + p3_nt2 + p2_nt3)[, -1]
LASSO_RE2 <- glmnet(x2_factor, Intensity, alpha = 1)
cvfit2 <- cv.glmnet(x2_factor, Intensity)
LASSO_predict2 <- predict(cvfit2, newx = x2_factor, s = "lambda.min")
LASSO_predict2_test <- predict(cvfit2, newx = x2_factor, s = "lambda.1se")
plot(LASSO_predict2, LASSO_predict2_test)

#Here, discuss with Raluca, which lambda to use, default is 1se.We could shrink variables
#using 1se but not lambda.min.
plot(LASSO_predict2, TDG_context$Intensity)
plot(LASSO_predict2_test, TDG_context$Intensity)
#

R_calc_data2 = data.frame(x = LASSO_predict2, y = TDG_context$Intensity)
predict_lm2 = lm(X1 ~ y, data = R_calc_data2)
summary(predict_lm2)
#######################N-mer_test_TDGcontext_boundary_end######################

x3_factor <- model.matrix(Intensity ~ p1_nt1 + p3_nt1 + p4_nt1 + 
                            p1_nt2 + p2_nt2 + p3_nt2 + 
                            p1_nt3 + p2_nt3)[, -1]
LASSO_RE3 <- glmnet(x3_factor, Intensity, alpha = 1)
cvfit3 <- cv.glmnet(x3_factor, Intensity)
LASSO_predict3 <- predict(cvfit3, newx = x3_factor, s = "lambda.min")
plot(LASSO_predict3, TDG_context$Intensity)

R_calc_data3 = data.frame(x = LASSO_predict3, y = TDG_context$Intensity)
predict_lm3 = lm(X1 ~ y, data = R_calc_data3)
summary(predict_lm3)

#######################N-mer_test_TDGcontext_boundary_end######################

x4_factor <- model.matrix(Intensity ~ p1_nt1 + p3_nt1 + p4_nt1 + p5_nt1 + 
                            p1_nt2 + p2_nt2 + p3_nt2 + p4_nt2 +
                            p1_nt3 + p2_nt3 + p3_nt3)[, -1]
LASSO_RE4 <- glmnet(x4_factor, Intensity, alpha = 1)
cvfit4 <- cv.glmnet(x4_factor, Intensity)
LASSO_predict4 <- predict(cvfit4, newx = x4_factor, s = "lambda.min")
plot(LASSO_predict4, TDG_context$Intensity)

R_calc_data4 = data.frame(x = LASSO_predict4, y = TDG_context$Intensity)
predict_lm4 = lm(X1 ~ y, data = R_calc_data4)
summary(predict_lm4)

#######################N-mer_test_TDGcontext_boundary_end######################

x5_factor <- model.matrix(Intensity ~ p0_nt1 + p1_nt1 + p3_nt1 + p4_nt1 + p5_nt1 + 
                            p0_nt2 + p1_nt2 + p2_nt2 + p3_nt2 + p4_nt2 +
                            p0_nt3 + p1_nt3 + p2_nt3 + p3_nt3)[, -1]
LASSO_RE5 <- glmnet(x5_factor, Intensity, alpha = 1)
cvfit5 <- cv.glmnet(x5_factor, Intensity)
LASSO_predict5 <- predict(cvfit5, newx = x5_factor, s = "lambda.min")
plot(LASSO_predict5, TDG_context$Intensity)

R_calc_data5 = data.frame(x = LASSO_predict5, y = TDG_context$Intensity)
predict_lm5 = lm(X1 ~ y, data = R_calc_data5)
summary(predict_lm5)

#######################N-mer_test_TDGcontext_boundary_end######################

x6_factor <- model.matrix(Intensity ~ p0_nt1 + p1_nt1 + p3_nt1 + p4_nt1 + p5_nt1 + p6_nt1 +
                            p0_nt2 + p1_nt2 + p2_nt2 + p3_nt2 + p4_nt2 + p5_nt2 +
                            p0_nt3 + p1_nt3 + p2_nt3 + p3_nt3 + p4_nt3)[, -1]
LASSO_RE6 <- glmnet(x6_factor, Intensity, alpha = 1)
cvfit6 <- cv.glmnet(x6_factor, Intensity) #Here, default alpha = 1, we will use the cv-ed one.
LASSO_predict6 <- predict(cvfit6, newx = x6_factor, s = "lambda.min")
LASSO_predict6_test <- predict(cvfit6, newx = x6_factor, s = "lambda.1se")
plot(LASSO_predict6, TDG_context$Intensity)

R_calc_data6 = data.frame(x = LASSO_predict6, y = TDG_context$Intensity)
predict_lm6 = lm(X1 ~ y, data = R_calc_data6)
summary(predict_lm6)

###6mer is the longest model, here we output the result###
write.table(as.matrix(cvfit6[["glmnet.fit"]][["beta"]]), 
            file = "/Users/gerogehou/Desktop/R/TDG_context_analysis/LASSO_beta.tsv", sep = "\t")

#######################N-mer_test_TDGcontext_boundary_end######################

x52_factor <- model.matrix(Intensity ~ p1_nt1 + p3_nt1 + p4_nt1 + p5_nt1 + p6_nt1 +
                            p1_nt2 + p2_nt2 + p3_nt2 + p4_nt2 + p5_nt2 +
                            p1_nt3 + p2_nt3 + p3_nt3 + p4_nt3)[, -1]
LASSO_RE52 <- glmnet(x52_factor, Intensity, alpha = 1)
cvfit52 <- cv.glmnet(x52_factor, Intensity)
LASSO_predict52 <- predict(cvfit52, newx = x52_factor, s = "lambda.min")
plot(LASSO_predict52, TDG_context$Intensity)

R_calc_data52 = data.frame(x = LASSO_predict52, y = TDG_context$Intensity)
predict_lm52 = lm(X1 ~ y, data = R_calc_data52)
summary(predict_lm52)