TDG_context <- read.csv("/Users/gerogehou/Desktop/R/PBM/LoopContext_v1/Multivariate_done.csv")
TDG_input <- read.csv("/Users/gerogehou/Desktop/R/PBM/LoopContext_v1/Input.csv")

TDG_reg = lm(Intensity ~ X0 + X1 + X3 + X4 + X5 + X6
             + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14
             + X15 + X16 + X17, data = TDG_context)

coef1 <- summary(TDG_reg)$coefficients

write.csv(coef1, file = "/Users/gerogehou/Desktop/R/PBM/LoopContext_v1/Coef.csv")

rSquared <- summary(TDG_reg)$r.squared


pred_TDG = predict.lm(TDG_reg, TDG_input)

write.csv(pred_TDG, file = "/Users/gerogehou/Desktop/R/PBM/LoopContext_v1/Pred.csv")

typeof(pred_TDG)

TDG_align <- read.csv("/Users/gerogehou/Desktop/R/PBM/LoopContext_v1/Pred.csv")
plot(TDG_align$X0, TDG_align$X1)
plot(TDG_align$X0, TDG_align$X2)
plot(TDG_align$X2, TDG_align$X1)

#TDG_multi = lm(Intensity ~ X, data = TDG_context)

#summary(TDG_multi)
