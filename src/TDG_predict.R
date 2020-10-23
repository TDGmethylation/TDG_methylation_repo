library(glmnet)

#Input using TDG binding data
TDG_context <- read.csv("/Users/gerogehou/Desktop/R/PBM/LoopContext_v1/Multivariate_done.csv")
TDG_input <- read.csv("/Users/gerogehou/Desktop/R/PBM/LoopContext_v1/Input.csv")

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

x_factor <- model.matrix(Intensity ~ p0_nt1 + p1_nt1 + p3_nt1 + p4_nt1 + p5_nt1 + p6_nt1
                         + p0_nt2 + p1_nt2 + p2_nt2 + p3_nt2 + p4_nt2 + p5_nt2 +
                           p0_nt3 + p1_nt3 + p2_nt3 + p3_nt3 + p4_nt3)[, -1]

Intensity <- TDG_context$Intensity

#Run LASSO
LASSO_RE <- glmnet(x_factor, Intensity, alpha = 1)
plot(LASSO_RE)

#Cross validation for lambda choice
cvfit <- cv.glmnet(x_factor, Intensity)
plot(cvfit)
coef_minMSE_sparse = coef(cvfit, s = "lambda.min")
coef_minMSE <- summary(coef_minMSE_sparse)
write.csv(coef_minMSE, file = "/Users/gerogehou/Desktop/R/PBM/LoopContext_v1/LASSO_Coef_minMSE.csv", 
         fileEncoding = "UTF-8")

LASSO_predict = predict(cvfit, newx = x_factor, s = "lambda.min")
LASSO_predict2 = predict(LASSO_RE, newx = x_factor)

write.csv(LASSO_predict2, file = "/Users/gerogehou/Desktop/R/PBM/LoopContext_v1/LASSO_Coed_predy2.csv")

write.csv(LASSO_predict, file = "/Users/gerogehou/Desktop/R/PBM/LoopContext_v1/LASSO_Coed_predy.csv")
write.csv(x_factor, file = "/Users/gerogehou/Desktop/R/PBM/LoopContext_v1/LASSO_Coed_predx.csv")

#plot
plot(LASSO_predict, TDG_context$Intensity)

#Data summary for interpretation and save results
beta <- as.matrix(summary(LASSO_RE$beta))

write.csv(beta, file = "/Users/gerogehou/Desktop/R/PBM/LoopContext_v1/LASSO_Coef.csv",
          fileEncoding = "UTF-8")

Dimnames_x = LASSO_RE$beta@Dimnames[[1]]
write.csv(Dimnames_x, file = "/Users/gerogehou/Desktop/R/PBM/LoopContext_v1/dimname1.csv",
     fileEncoding = "UTF-8")

Dimnames_y = LASSO_RE$beta@Dimnames[[2]]
write.csv(Dimnames_y, file = "/Users/gerogehou/Desktop/R/PBM/LoopContext_v1/dimname2.csv",
          fileEncoding = "UTF-8")



