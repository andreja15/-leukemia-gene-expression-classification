library(glmnet)
library(randomForest)
library(pROC)
library(readxl)
library(ggplot2)

train_raw <- read_excel("train dataset.xlsx")
indep_raw <- read_excel("independent dataset.xlsx")
actual    <- read_excel("actual.csv.xlsx")

# Expression matrix
get_expr <- function(df) {
  df           <- as.data.frame(df)
  gene_ids     <- df[["Gene Accession Number"]]
  drop         <- c("Gene Description", "Gene Accession Number")
  df           <- df[, !colnames(df) %in% drop]
  call_idx     <- grep("call", colnames(df), ignore.case = TRUE)
  df           <- df[, -call_idx]
  df           <- apply(df, 2, as.numeric)
  rownames(df) <- gene_ids
  return(df)
}

train_expr <- get_expr(train_raw)
indep_expr <- get_expr(indep_raw)

# Preprocessing
X_train_raw <- t(train_expr)
X_test_raw  <- t(indep_expr)

control_idx <- grep("^AFFX", colnames(X_train_raw))
X_train_raw <- X_train_raw[, -control_idx]
X_test_raw  <- X_test_raw[,  -control_idx]

sd_train   <- apply(X_train_raw, 2, sd)
keep_genes <- sd_train > 100
X_train_f  <- X_train_raw[, keep_genes]
X_test_f   <- X_test_raw[,  keep_genes]

train_means <- colMeans(X_train_f)
train_sds   <- apply(X_train_f, 2, sd)
X_train     <- scale(X_train_f)
X_test      <- scale(X_test_f, center = train_means, scale = train_sds)

train_ids <- as.integer(rownames(X_train))
test_ids  <- as.integer(rownames(X_test))
Y_train   <- ifelse(actual$cancer[match(train_ids, actual$patient)] == "AML", 1, 0)
Y_test    <- ifelse(actual$cancer[match(test_ids,  actual$patient)] == "AML", 1, 0)

# Metrics
get_metrics <- function(pred_prob, true_labels) {
  pred_class  <- ifelse(pred_prob >= 0.5, 1, 0)
  TP          <- sum(pred_class == 1 & true_labels == 1)
  TN          <- sum(pred_class == 0 & true_labels == 0)
  FP          <- sum(pred_class == 1 & true_labels == 0)
  FN          <- sum(pred_class == 0 & true_labels == 1)
  sensitivity <- round(TP / (TP + FN), 4)
  specificity <- round(TN / (TN + FP), 4)
  accuracy    <- round((TP + TN) / length(true_labels), 4)
  return(c(sensitivity = sensitivity,
           specificity = specificity,
           accuracy    = accuracy))
}

# Ridge regression
set.seed(1)
ridge_cv          <- cv.glmnet(X_train, Y_train, alpha = 0, family = "binomial")
best_lambda_ridge <- ridge_cv$lambda.min
ridge_final       <- glmnet(X_train, Y_train, alpha = 0, family = "binomial",
                             lambda = best_lambda_ridge)
ridge_pred        <- predict(ridge_final, X_test, type = "response")
ridge_roc         <- roc(Y_test, as.vector(ridge_pred))
ridge_auc         <- auc(ridge_roc)
ridge_m           <- get_metrics(as.vector(ridge_pred), Y_test)

ridge_coef <- coef(ridge_final)
ridge_df   <- data.frame(gene = rownames(ridge_coef), coef = as.vector(ridge_coef))
ridge_df   <- ridge_df[ridge_df$gene != "(Intercept)", ]
ridge_top  <- head(ridge_df[order(abs(ridge_df$coef), decreasing = TRUE), ], 10)

cat("RIDGE REGRESSION\n")
cat("Lambda:      ", round(best_lambda_ridge, 4), "\n")
cat("AUC:         ", round(ridge_auc, 4), "\n")
cat("Sensitivity: ", ridge_m["sensitivity"], "\n")
cat("Specificity: ", ridge_m["specificity"], "\n")
cat("Accuracy:    ", ridge_m["accuracy"],    "\n")
print(ridge_top)

# Lasso regression
set.seed(1)
lasso_cv          <- cv.glmnet(X_train, Y_train, alpha = 1, family = "binomial")
best_lambda_lasso <- lasso_cv$lambda.min
lasso_final       <- glmnet(X_train, Y_train, alpha = 1, family = "binomial",
                             lambda = best_lambda_lasso)
lasso_pred        <- predict(lasso_final, X_test, type = "response")
lasso_roc         <- roc(Y_test, as.vector(lasso_pred))
lasso_auc         <- auc(lasso_roc)
lasso_m           <- get_metrics(as.vector(lasso_pred), Y_test)

lasso_coef <- coef(lasso_final)
lasso_df   <- data.frame(gene = rownames(lasso_coef), coef = as.vector(lasso_coef))
lasso_df   <- lasso_df[lasso_df$coef != 0 & lasso_df$gene != "(Intercept)", ]
lasso_df   <- lasso_df[order(abs(lasso_df$coef), decreasing = TRUE), ]

cat("LASSO REGRESSION\n")
cat("Lambda:         ", round(best_lambda_lasso, 4), "\n")
cat("Genes Selected: ", nrow(lasso_df), "\n")
cat("AUC:            ", round(lasso_auc, 4), "\n")
cat("Sensitivity:    ", lasso_m["sensitivity"], "\n")
cat("Specificity:    ", lasso_m["specificity"], "\n")
cat("Accuracy:       ", lasso_m["accuracy"],    "\n")
print(lasso_df)

# Elastic net
set.seed(1)
enet_cv          <- cv.glmnet(X_train, Y_train, alpha = 0.9, family = "binomial")
best_lambda_enet <- enet_cv$lambda.min
enet_final       <- glmnet(X_train, Y_train, alpha = 0.9, family = "binomial",
                            lambda = best_lambda_enet)
enet_pred        <- predict(enet_final, X_test, type = "response")
enet_roc         <- roc(Y_test, as.vector(enet_pred))
enet_auc         <- auc(enet_roc)
enet_m           <- get_metrics(as.vector(enet_pred), Y_test)

enet_coef <- coef(enet_final)
enet_df   <- data.frame(gene = rownames(enet_coef), coef = as.vector(enet_coef))
enet_df   <- enet_df[enet_df$coef != 0 & enet_df$gene != "(Intercept)", ]
enet_df   <- enet_df[order(abs(enet_df$coef), decreasing = TRUE), ]

cat("ELASTIC NET\n")
cat("Lambda:         ", round(best_lambda_enet, 4), "\n")
cat("Alpha:           0.9\n")
cat("Genes Selected: ", nrow(enet_df), "\n")
cat("AUC:            ", round(enet_auc, 4), "\n")
cat("Sensitivity:    ", enet_m["sensitivity"], "\n")
cat("Specificity:    ", enet_m["specificity"], "\n")
cat("Accuracy:       ", enet_m["accuracy"],    "\n")
print(enet_df)

# Random forest
set.seed(1)
train_df <- data.frame(leukemia = as.factor(Y_train), X_train)
test_df  <- data.frame(leukemia = as.factor(Y_test),  X_test)
rf_final <- randomForest(leukemia ~ ., data = train_df)
rf_pred  <- predict(rf_final, test_df, type = "prob")[, 2]
rf_roc   <- roc(Y_test, rf_pred)
rf_auc   <- auc(rf_roc)
rf_m     <- get_metrics(rf_pred, Y_test)

importance_df      <- as.data.frame(importance(rf_final))
importance_df$gene <- rownames(importance_df)
importance_df      <- importance_df[order(importance_df$MeanDecreaseGini,
                                          decreasing = TRUE), ]
rf_top <- head(importance_df[, c("gene", "MeanDecreaseGini")], 10)

cat("RANDOM FOREST\n")
cat("AUC:         ", round(rf_auc, 4), "\n")
cat("Sensitivity: ", rf_m["sensitivity"], "\n")
cat("Specificity: ", rf_m["specificity"], "\n")
cat("Accuracy:    ", rf_m["accuracy"],    "\n")
print(rf_top)

# Variable importance plot
ggplot(rf_top, aes(x = reorder(gene, MeanDecreaseGini),
                   y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  xlab("Gene") +
  ylab("Mean Decrease Gini") +
  ggtitle("Random Forest Variable Importance") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9))

# SIS + Lasso
sis_screen <- function(X, Y, dn) {
  cors              <- abs(cor(X, Y))
  cors[is.na(cors)] <- 0
  return(order(cors, decreasing = TRUE)[1:dn])
}

set.seed(1)
folds               <- numeric(length(Y_train))
folds[Y_train == 0] <- sample(rep(1:5, length.out = sum(Y_train == 0)))
folds[Y_train == 1] <- sample(rep(1:5, length.out = sum(Y_train == 1)))

dn_values <- c(20, 100, 200, 500, 1000)
cv_errors <- rep(0, length(dn_values))

for (i in seq_along(dn_values)) {
  dn          <- dn_values[i]
  fold_errors <- rep(0, 5)
  for (k in 1:5) {
    Xft <- X_train[folds != k, ]; Xfv <- X_train[folds == k, ]
    Yft <- Y_train[folds != k];   Yfv <- Y_train[folds == k]
    cc  <- apply(Xft, 2, function(x) all(is.finite(x)) & sd(x) > 0)
    Xft <- Xft[, cc]; Xfv <- Xfv[, cc]
    tg  <- sis_screen(Xft, Yft, dn)
    Xtt <- Xft[, tg]; Xtv <- Xfv[, tg]
    fc  <- apply(Xtt, 2, function(x) all(is.finite(x)) & sd(x) > 0)
    Xtt <- Xtt[, fc]; Xtv <- Xtv[, fc]
    lasso_cv_fold  <- cv.glmnet(Xtt, Yft, alpha = 1, family = "binomial")
    pred           <- predict(lasso_cv_fold, Xtv, s = "lambda.min", type = "response")
    fold_errors[k] <- auc(roc(Yfv, as.vector(pred), quiet = TRUE))
  }
  cv_errors[i] <- mean(fold_errors)
  cat("dn =", dn, "| CV AUC =", round(cv_errors[i], 4), "\n")
}

best_dn    <- dn_values[which.max(cv_errors)]
best_genes <- sis_screen(X_train, Y_train, best_dn)
best_lambda <- cv.glmnet(X_train[, best_genes], Y_train,
                          alpha = 1, family = "binomial")$lambda.min
sis_model  <- glmnet(X_train[, best_genes], Y_train, alpha = 1,
                      family = "binomial", lambda = best_lambda)
sis_pred   <- predict(sis_model, X_test[, best_genes], type = "response")
sis_roc    <- roc(Y_test, as.vector(sis_pred))
sis_auc    <- auc(sis_roc)
sis_m      <- get_metrics(as.vector(sis_pred), Y_test)

sis_coef <- coef(sis_model)
sis_df   <- data.frame(gene = rownames(sis_coef), coef = as.vector(sis_coef))
sis_df   <- sis_df[sis_df$coef != 0 & sis_df$gene != "(Intercept)", ]
sis_df   <- sis_df[order(abs(sis_df$coef), decreasing = TRUE), ]

cat("SIS + LASSO\n")
cat("Best dn:        ", best_dn, "\n")
cat("Lambda:         ", round(best_lambda, 4), "\n")
cat("Genes Selected: ", nrow(sis_df), "\n")
cat("AUC:            ", round(sis_auc, 4), "\n")
cat("Sensitivity:    ", sis_m["sensitivity"], "\n")
cat("Specificity:    ", sis_m["specificity"], "\n")
cat("Accuracy:       ", sis_m["accuracy"],    "\n")
print(sis_df)

# Bivariate analysis
top_genes <- c("Y12670_at", "U50136_rna1_at", "X95735_at", "U22376_cds2_s_at",
               "D49950_at", "M19507_at", "U82759_at", "M23197_at",
               "M37435_at", "L08246_at")

all_idx <- which(Y_train == 0)
aml_idx <- which(Y_train == 1)

p_vals <- sapply(top_genes, function(g) {
  t.test(X_train[aml_idx, g], X_train[all_idx, g])$p.value
})

p_labels <- sapply(p_vals, function(p) {
  if (p < 0.001) "p < 0.001"
  else paste0("p = ", round(p, 3))
})

cat("BIVARIATE ANALYSIS P-VALUES\n")
for (i in seq_along(top_genes)) {
  cat(top_genes[i], ":", p_labels[i], "\n")
}

plot_df <- data.frame(
  expression = as.vector(X_train[, top_genes]),
  gene       = rep(top_genes, each = nrow(X_train)),
  type       = rep(ifelse(Y_train == 0, "ALL", "AML"), length(top_genes))
)

label_df <- data.frame(
  gene  = top_genes,
  label = p_labels
)

ggplot(plot_df, aes(x = type, y = expression, fill = type)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  geom_text(data = label_df,
            aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE,
            vjust = 1.5, size = 3) +
  xlab("Leukemia Type") +
  ylab("Standardized Expression") +
  ggtitle("Differential Expression of Top Consensus Genes Between ALL and AML Patients") +
  scale_fill_manual(values = c("ALL" = "blue", "AML" = "red")) +
  theme_minimal() +
  theme(strip.text   = element_text(size = 8),
        legend.title = element_blank())

# ROC curves
par(mfrow = c(2, 3))
plot(ridge_roc, col = "blue",      lwd = 2, legacy.axes = TRUE,
     main = paste("Ridge\nAUC =", round(ridge_auc, 4)))
plot(lasso_roc, col = "red",       lwd = 2, legacy.axes = TRUE,
     main = paste("Lasso\nAUC =", round(lasso_auc, 4)))
plot(enet_roc,  col = "green", lwd = 2, legacy.axes = TRUE,
     main = paste("Elastic Net\nAUC =", round(enet_auc, 4)))
plot(rf_roc,    col = "purple",    lwd = 2, legacy.axes = TRUE,
     main = paste("Random Forest\nAUC =", round(rf_auc, 4)))
plot(sis_roc,   col = "orange",    lwd = 2, legacy.axes = TRUE,
     main = paste("SIS + Lasso\nAUC =", round(sis_auc, 4)))
par(mfrow = c(1, 1))
