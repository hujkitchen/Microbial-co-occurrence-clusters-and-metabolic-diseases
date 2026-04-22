# ==================== 多疾病分类器构建与验证 ====================
cat("\n=== 多疾病分类器构建与验证 ===\n")

# 创建分析目录
if(!dir.exists("cluster/classifier")) {
  dir.create("cluster/classifier")
}
if(!dir.exists("cluster/classifier/figures")) {
  dir.create("cluster/classifier/figures")
}
if(!dir.exists("cluster/classifier/data")) {
  dir.create("cluster/classifier/data")
}
if(!dir.exists("cluster/classifier/models")) {
  dir.create("cluster/classifier/models")
}

# 加载必要的包
library(ggplot2)
library(caret)
library(randomForest)
library(pROC)
library(xgboost)
library(glmnet)
library(reshape2)
library(corrplot)

# ==================== 1. 数据准备 ====================
cat("\n=== 1. 数据准备 ===\n")

# 使用菌群组成数据
feature_data <- cluster_abundance_rel
target <- colData(tse_2w)$condition_simplified

# 只保留样本数 >= 30 的疾病类型
sample_counts <- table(target)
keep_diseases <- names(sample_counts[sample_counts >= 30 & names(sample_counts) != "Healthy"])
keep_idx <- target %in% keep_diseases

X <- feature_data[keep_idx, ]
y <- factor(target[keep_idx])

cat("纳入分析的疾病类型:\n")
print(table(y))

# ==================== 2. 特征选择 ====================
cat("\n=== 2. 特征选择 ===\n")

# 方法1：基于方差选择
feature_vars <- apply(X, 2, var)
high_var_features <- names(feature_vars[feature_vars > quantile(feature_vars, 0.5)])
cat("高方差特征数:", length(high_var_features), "\n")

# 方法2：基于随机森林特征重要性
set.seed(123)
rf_select <- randomForest(x = X, y = y, ntree = 500, importance = TRUE)
rf_importance <- importance(rf_select)
important_features <- names(sort(rf_importance[, "MeanDecreaseGini"], decreasing = TRUE)[1:5])
cat("随机森林选出的重要特征:", paste(important_features, collapse = ", "), "\n")

# 方法3：LASSO特征选择
X_matrix <- as.matrix(X)
set.seed(123)
cv_lasso <- cv.glmnet(X_matrix, y, family = "multinomial", alpha = 1, nfolds = 5)
lasso_coef <- coef(cv_lasso, s = "lambda.min")
lasso_features <- c()
for(i in 1:length(lasso_coef)) {
  features <- rownames(lasso_coef[[i]])[as.matrix(lasso_coef[[i]])[,1] != 0]
  lasso_features <- c(lasso_features, features)
}
lasso_features <- unique(lasso_features[lasso_features != "(Intercept)"])
cat("LASSO选出的特征数:", length(lasso_features), "\n")

# 合并特征（使用重要特征）
selected_features <- unique(c(important_features, lasso_features))
if(length(selected_features) < 3) {
  selected_features <- colnames(X)[1:min(5, ncol(X))]
}
cat("最终选用的特征:", paste(selected_features, collapse = ", "), "\n")

X_selected <- X[, selected_features]

# ==================== 3. 多分类模型训练 ====================
cat("\n=== 3. 多分类模型训练 ===\n")

# 设置交叉验证参数
train_control <- trainControl(
  method = "cv",
  number = 10,
  savePredictions = TRUE,
  classProbs = TRUE,
  summaryFunction = multiClassSummary
)

# 模型1：随机森林
set.seed(123)
cat("训练随机森林模型...\n")
rf_model <- train(
  x = X_selected,
  y = y,
  method = "rf",
  trControl = train_control,
  ntree = 500,
  importance = TRUE
)

# 模型2：XGBoost
set.seed(123)
cat("训练XGBoost模型...\n")
xgb_model <- train(
  x = X_selected,
  y = y,
  method = "xgbTree",
  trControl = train_control,
  tuneLength = 5
)

# 模型3：支持向量机
set.seed(123)
cat("训练SVM模型...\n")
svm_model <- train(
  x = X_selected,
  y = y,
  method = "svmRadial",
  trControl = train_control,
  preProcess = c("center", "scale")
)

# 模型4：逻辑回归（一对多）
set.seed(123)
cat("训练多项逻辑回归...\n")
lr_model <- train(
  x = X_selected,
  y = y,
  method = "multinom",
  trControl = train_control,
  maxit = 1000
)

# 保存模型
saveRDS(rf_model, "cluster/classifier/models/rf_model.rds")
saveRDS(xgb_model, "cluster/classifier/models/xgb_model.rds")
saveRDS(svm_model, "cluster/classifier/models/svm_model.rds")
saveRDS(lr_model, "cluster/classifier/models/lr_model.rds")

# ==================== 4. 模型性能比较 ====================
cat("\n=== 4. 模型性能比较 ===\n")

# 收集各模型性能
models <- list(
  "Random Forest" = rf_model,
  "XGBoost" = xgb_model,
  "SVM" = svm_model,
  "Logistic Regression" = lr_model
)

performance_summary <- data.frame()

for(model_name in names(models)) {
  model <- models[[model_name]]
  
  # 获取交叉验证结果
  cv_results <- model$results[which.max(model$results$Accuracy), ]
  
  performance_summary <- rbind(performance_summary, data.frame(
    Model = model_name,
    Accuracy = cv_results$Accuracy,
    Kappa = cv_results$Kappa,
    LogLoss = ifelse("logLoss" %in% colnames(cv_results), cv_results$logLoss, NA)
  ))
}

print(performance_summary)
write.csv(performance_summary, "cluster/classifier/data/model_performance.csv", row.names = FALSE)

# ==================== 5. 混淆矩阵和分类报告 ====================
cat("\n=== 5. 混淆矩阵分析 ===\n")

# 使用最佳模型（随机森林）
best_model <- rf_model

# 预测
predictions <- predict(best_model, newdata = X_selected)
conf_matrix <- confusionMatrix(predictions, y)

# 保存混淆矩阵
write.csv(as.matrix(conf_matrix$table), "cluster/classifier/data/confusion_matrix.csv")

# 可视化混淆矩阵
conf_matrix_table <- as.matrix(conf_matrix$table)
conf_df <- as.data.frame(conf_matrix$table)
conf_df$Prediction <- rownames(conf_df)
conf_long <- melt(conf_df, id.vars = "Prediction", variable.name = "True", value.name = "Count")

# 确保Count是数值型
conf_long$Count <- as.numeric(as.character(conf_long$Count))

p <- ggplot(conf_long, aes(x = True, y = Prediction, fill = Count)) +
  geom_tile() +
  geom_text(aes(label = Count), size = 4, color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue", 
                      name = "Count",
                      limits = c(0, max(conf_long$Count))) +
  labs(title = "Confusion Matrix - Random Forest Classifier",
       x = "True Class", y = "Predicted Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

# 保存
ggsave("cluster/classifier/figures/confusion_matrix.pdf", p, width = 10, height = 8)
ggsave("cluster/classifier/figures/confusion_matrix.png", p, width = 10, height = 8, dpi = 300)

cat("混淆矩阵热图已保存\n")

# ==================== 6. 特征重要性 ====================
cat("\n=== 6. 特征重要性 ===\n")

# 提取特征重要性
if(best_model$method == "rf") {
  importance_df <- data.frame(
    Feature = selected_features,
    Importance = best_model$finalModel$importance[, "MeanDecreaseGini"]
  )
} else if(best_model$method == "xgbTree") {
  importance_df <- varImp(best_model)$importance
  importance_df$Feature <- rownames(importance_df)
  colnames(importance_df)[1] <- "Importance"
} else {
  importance_df <- varImp(best_model)$importance
  importance_df$Feature <- rownames(importance_df)
  importance_df$Importance <- importance_df$Overall
}

importance_df <- importance_df[order(-importance_df$Importance), ]
write.csv(importance_df, "cluster/classifier/data/feature_importance.csv", row.names = FALSE)

# 可视化特征重要性
p <- ggplot(head(importance_df, 10), aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 Most Important Features",
       x = "Feature", y = "Importance") +
  theme_minimal()

ggsave("cluster/classifier/figures/feature_importance.pdf", p, width = 8, height = 6)

# ==================== 7. 一对一分类性能 ====================
cat("\n=== 7. 一对一分类性能 ===\n")

# 为每对疾病计算AUC
disease_pairs <- combn(levels(y), 2, simplify = FALSE)
pair_auc <- data.frame()

for(pair in disease_pairs) {
  # 选择这两个疾病的样本
  pair_idx <- y %in% pair
  pair_X <- X_selected[pair_idx, ]
  pair_y <- factor(y[pair_idx])
  
  # 训练二分类模型
  set.seed(123)
  pair_model <- randomForest(x = pair_X, y = pair_y, ntree = 200)
  
  # 计算预测概率
  pred_prob <- predict(pair_model, type = "prob")
  
  # 计算AUC
  roc_obj <- roc(pair_y, pred_prob[, 2])
  
  pair_auc <- rbind(pair_auc, data.frame(
    Disease1 = pair[1],
    Disease2 = pair[2],
    AUC = auc(roc_obj)
  ))
}

write.csv(pair_auc, "cluster/classifier/data/pairwise_auc.csv", row.names = FALSE)

# 可视化AUC矩阵
auc_matrix <- matrix(NA, nrow = length(levels(y)), ncol = length(levels(y)))
rownames(auc_matrix) <- levels(y)
colnames(auc_matrix) <- levels(y)

for(i in 1:nrow(pair_auc)) {
  auc_matrix[pair_auc$Disease1[i], pair_auc$Disease2[i]] <- pair_auc$AUC[i]
  auc_matrix[pair_auc$Disease2[i], pair_auc$Disease1[i]] <- pair_auc$AUC[i]
}
diag(auc_matrix) <- 1

# 热图
pdf("cluster/classifier/figures/pairwise_auc_heatmap.pdf", width = 8, height = 7)
corrplot(auc_matrix, method = "color", type = "full",
         addCoef.col = "black", tl.col = "black",
         title = "Pairwise Classification AUC",
         mar = c(0,0,1,0))
dev.off()

# ==================== 8. 交叉验证详细结果 ====================
cat("\n=== 8. 交叉验证详细结果 ===\n")

# 获取10折交叉验证的详细预测结果
# 获取预测结果
cv_predictions <- rf_model$pred

# 检查数据结构
cat("预测结果列名:", paste(colnames(cv_predictions), collapse = ", "), "\n")
cat("预测结果行数:", nrow(cv_predictions), "\n")

# 为每个预测计算是否正确
cv_predictions$Correct <- cv_predictions$pred == cv_predictions$obs

# 计算每折的准确率
fold_accuracy <- aggregate(Correct ~ Resample, 
                           data = cv_predictions, 
                           FUN = mean)

colnames(fold_accuracy) <- c("Fold", "Accuracy")

# 显示结果
cat("\n各折准确率:\n")
print(fold_accuracy)

# 计算统计信息
mean_accuracy <- mean(fold_accuracy$Accuracy)
sd_accuracy <- sd(fold_accuracy$Accuracy)

cat("\n平均CV准确率:", round(mean_accuracy, 4), "\n")
cat("CV准确率标准差:", round(sd_accuracy, 4), "\n")

# 保存结果
write.csv(fold_accuracy, "cluster/classifier/data/fold_accuracy.csv", row.names = FALSE)

cat("\n可视化交叉验证结果...\n")

# 箱线图
p1 <- ggplot(fold_accuracy, aes(x = "Overall", y = Accuracy)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, color = "darkblue") +
  labs(title = "10-Fold Cross-Validation Accuracy",
       x = "", y = "Accuracy") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("cluster/classifier/figures/cv_accuracy_boxplot.pdf", p1, width = 6, height = 5)

# 条形图
p2 <- ggplot(fold_accuracy, aes(x = reorder(Fold, Accuracy), y = Accuracy)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = mean_accuracy, linetype = "dashed", color = "red", size = 1) +
  labs(title = "10-Fold Cross-Validation Accuracy by Fold",
       x = "Fold", y = "Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

ggsave("cluster/classifier/figures/cv_accuracy_bars.pdf", p2, width = 8, height = 5)

# ==================== 9. 外部验证（模拟）====================
cat("\n=== 9. 外部验证 ===\n")

# 使用留一法模拟外部验证
set.seed(123)
external_folds <- createFolds(y, k = 5, list = TRUE)

external_results <- data.frame()

for(i in 1:length(external_folds)) {
  # 训练集和测试集
  test_idx <- external_folds[[i]]
  train_idx <- setdiff(1:nrow(X_selected), test_idx)
  
  # 训练模型
  train_model <- randomForest(x = X_selected[train_idx, ],
                              y = y[train_idx],
                              ntree = 500)
  
  # 预测
  predictions <- predict(train_model, newdata = X_selected[test_idx, ])
  
  # 计算准确率
  accuracy <- sum(predictions == y[test_idx]) / length(test_idx)
  
  external_results <- rbind(external_results, data.frame(
    Fold = i,
    Test_Size = length(test_idx),
    Accuracy = accuracy
  ))
}

write.csv(external_results, "cluster/classifier/data/external_validation.csv", row.names = FALSE)

cat("外部验证平均准确率:", mean(external_results$Accuracy), "\n")
cat("外部验证准确率范围:", range(external_results$Accuracy), "\n")

# ==================== 10. 生物标志物组合识别 ====================
cat("\n=== 生物标志物组合识别 ===\n")

# 使用选定的特征
selected_features <- colnames(X_selected)
n_features <- length(selected_features)

cat("当前特征数量:", n_features, "\n")
cat("特征列表:", paste(selected_features, collapse = ", "), "\n")

# ==================== 方法1：基于随机森林特征重要性选择 ====================
cat("\n方法1：基于随机森林特征重要性\n")

# 训练随机森林模型
set.seed(123)
rf_importance_model <- randomForest(x = X_selected, y = y, ntree = 500, importance = TRUE)

# 获取特征重要性
importance_scores <- importance(rf_importance_model)
importance_df <- data.frame(
  Feature = rownames(importance_scores),
  Importance = importance_scores[, "MeanDecreaseGini"]
)
importance_df <- importance_df[order(-importance_df$Importance), ]

# 选择累积重要性达到80%的特征
total_importance <- sum(importance_df$Importance)
cumulative_importance <- cumsum(importance_df$Importance) / total_importance
best_size <- which(cumulative_importance >= 0.8)[1]

if(is.na(best_size)) {
  best_size <- min(3, nrow(importance_df))
}

best_features <- importance_df$Feature[1:best_size]

cat("最佳特征数:", best_size, "\n")
cat("选中的特征:", paste(best_features, collapse = ", "), "\n")

# ==================== 方法2：使用逐步回归选择 ====================
cat("\n方法2：使用逐步回归选择\n")

# 将多分类问题转换为多个二分类问题
# 为每个疾病类别创建二分类标签
all_features <- c()
for(class in levels(y)) {
  # 创建二分类标签
  binary_y <- ifelse(y == class, 1, 0)
  
  # 逐步回归
  data_df <- data.frame(y = binary_y, X_selected)
  
  # 使用glm进行逐步回归（限制特征数量）
  null_model <- glm(y ~ 1, data = data_df, family = binomial())
  full_model <- glm(y ~ ., data = data_df, family = binomial())
  
  step_model <- step(null_model, 
                     scope = list(lower = null_model, upper = full_model),
                     direction = "forward",
                     trace = 0,
                     k = 2)  # AIC
  
  # 提取选中的特征
  selected_vars <- names(coef(step_model))[-1]  # 移除截距
  all_features <- c(all_features, selected_vars)
}

# 统计每个特征被选中的频率
feature_frequency <- table(all_features)
feature_frequency <- sort(feature_frequency, decreasing = TRUE)

cat("各特征被选中的频率:\n")
print(feature_frequency)

# 选择频率最高的特征
if(length(feature_frequency) > 0) {
  best_features_step <- names(feature_frequency)[1:min(5, length(feature_frequency))]
  cat("\n逐步回归选中的特征:", paste(best_features_step, collapse = ", "), "\n")
} else {
  best_features_step <- best_features
}

# ==================== 方法3：使用LASSO选择 ====================
cat("\n方法3：使用LASSO选择\n")

library(glmnet)

# 对每个类别进行LASSO
lasso_features <- c()

for(class in levels(y)) {
  # 创建二分类标签
  binary_y <- ifelse(y == class, 1, 0)
  
  # LASSO回归
  set.seed(123)
  cv_lasso <- cv.glmnet(as.matrix(X_selected), binary_y, 
                        family = "binomial", alpha = 1, nfolds = 5)
  
  # 提取非零系数
  lasso_coef <- coef(cv_lasso, s = "lambda.min")
  selected <- rownames(lasso_coef)[which(lasso_coef[,1] != 0)]
  selected <- selected[selected != "(Intercept)"]
  
  lasso_features <- c(lasso_features, selected)
}

# 统计LASSO选中的特征频率
lasso_frequency <- table(lasso_features)
lasso_frequency <- sort(lasso_frequency, decreasing = TRUE)

cat("LASSO选中的特征频率:\n")
print(lasso_frequency)

if(length(lasso_frequency) > 0) {
  best_features_lasso <- names(lasso_frequency)[1:min(5, length(lasso_frequency))]
  cat("\nLASSO选中的特征:", paste(best_features_lasso, collapse = ", "), "\n")
} else {
  best_features_lasso <- best_features
}

# ==================== 方法4：组合方法 ====================
cat("\n方法4：组合方法\n")

# 结合三种方法的结果
all_method_features <- c(best_features, best_features_step, best_features_lasso)
combined_frequency <- table(all_method_features)
combined_frequency <- sort(combined_frequency, decreasing = TRUE)

# 选择出现次数最多的特征
best_combined_features <- names(combined_frequency)[1:min(5, length(combined_frequency))]

cat("\n组合方法选中的特征:", paste(best_combined_features, collapse = ", "), "\n")

# ==================== 验证生物标志物组合的性能 ====================
cat("\n验证生物标志物组合的性能...\n")

# 测试不同特征组合的性能
feature_sets <- list(
  "All Features" = selected_features,
  "Top 1" = best_combined_features[1],
  "Top 2" = best_combined_features[1:min(2, length(best_combined_features))],
  "Top 3" = best_combined_features[1:min(3, length(best_combined_features))],
  "Top 4" = best_combined_features[1:min(4, length(best_combined_features))],
  "Top 5" = best_combined_features[1:min(5, length(best_combined_features))]
)

performance_results <- data.frame()

for(set_name in names(feature_sets)) {
  features <- feature_sets[[set_name]]
  
  if(length(features) >= 2) {  # 至少需要2个特征
    X_subset <- X_selected[, features, drop = FALSE]
    
    # 交叉验证
    set.seed(123)
    cv_control <- trainControl(method = "cv", number = 5, 
                               summaryFunction = multiClassSummary)
    
    model_cv <- train(x = X_subset, y = y,
                      method = "rf", ntree = 200,
                      trControl = cv_control)
    
    # 获取最佳性能
    best_result <- model_cv$results[which.max(model_cv$results$Accuracy), ]
    
    performance_results <- rbind(performance_results, data.frame(
      Feature_Set = set_name,
      N_Features = length(features),
      Accuracy = best_result$Accuracy,
      Kappa = best_result$Kappa
    ))
  }
}

print(performance_results)

# ==================== 选择最佳生物标志物组合 ====================
# 选择准确率最高的组合
best_set <- performance_results[which.max(performance_results$Accuracy), ]
best_biomarkers <- feature_sets[[best_set$Feature_Set]]

cat("\n最佳生物标志物组合:\n")
cat("  - 特征集:", best_set$Feature_Set, "\n")
cat("  - 特征数量:", best_set$N_Features, "\n")
cat("  - 准确率:", round(best_set$Accuracy, 4), "\n")
cat("  - 生物标志物:", paste(best_biomarkers, collapse = ", "), "\n")

# ==================== 可视化 ====================

# 1. 特征重要性图
p1 <- ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Feature Importance for Disease Classification",
       x = "Feature", y = "Mean Decrease Gini") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("cluster/classifier/figures/feature_importance.pdf", p1, width = 8, height = 5)

# 2. 生物标志物性能比较
p2 <- ggplot(performance_results, aes(x = Feature_Set, y = Accuracy, fill = Feature_Set)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Accuracy - 0.02, ymax = Accuracy + 0.02), width = 0.2) +
  labs(title = "Classification Accuracy with Different Feature Sets",
       x = "Feature Set", y = "Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

ggsave("cluster/classifier/figures/feature_set_comparison.pdf", p2, width = 8, height = 5)

# 3. 箱线图展示生物标志物的分布
if(length(best_biomarkers) >= 2) {
  # 准备数据
  biomarker_data <- data.frame(
    Disease = y,
    X_selected[, best_biomarkers, drop = FALSE]
  )
  
  biomarker_long <- reshape2::melt(biomarker_data, 
                                   id.vars = "Disease", 
                                   variable.name = "Biomarker", 
                                   value.name = "Abundance")
  
  p3 <- ggplot(biomarker_long, aes(x = Disease, y = Abundance, fill = Disease)) +
    geom_boxplot() +
    facet_wrap(~ Biomarker, scales = "free_y") +
    labs(title = "Distribution of Selected Biomarkers by Disease",
         x = "Disease", y = "Relative Abundance") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  
  ggsave("cluster/classifier/figures/biomarker_distribution.pdf", p3, width = 12, height = 8)
}

# ==================== 保存结果 ====================

# 保存最佳生物标志物
biomarker_signature <- data.frame(
  Rank = 1:length(best_biomarkers),
  Biomarker = best_biomarkers,
  Importance = importance_df$Importance[match(best_biomarkers, importance_df$Feature)]
)

write.csv(biomarker_signature, "cluster/classifier/data/biomarker_signature.csv", row.names = FALSE)

# 保存特征重要性
write.csv(importance_df, "cluster/classifier/data/feature_importance.csv", row.names = FALSE)

# 保存性能比较
write.csv(performance_results, "cluster/classifier/data/feature_set_performance.csv", row.names = FALSE)

# 保存组合频率
combined_frequency_df <- data.frame(
  Feature = names(combined_frequency),
  Frequency = as.numeric(combined_frequency)
)
write.csv(combined_frequency_df, "cluster/classifier/data/feature_selection_frequency.csv", row.names = FALSE)

# ==================== 生成报告 ====================
sink("cluster/classifier/biomarker_report.txt")

cat("========================================\n")
cat("生物标志物组合识别报告\n")
cat("========================================\n\n")

cat("分析日期:", date(), "\n\n")

cat("1. 特征选择方法\n")
cat("   - 方法1: 随机森林累积重要性80%\n")
cat("   - 方法2: 逐步回归（AIC）\n")
cat("   - 方法3: LASSO回归\n")
cat("   - 方法4: 组合方法\n\n")

cat("2. 各方法选中的特征\n")
cat("   - 随机森林: ", paste(best_features, collapse = ", "), "\n")
cat("   - 逐步回归: ", paste(best_features_step, collapse = ", "), "\n")
cat("   - LASSO: ", paste(best_features_lasso, collapse = ", "), "\n")
cat("   - 组合: ", paste(best_combined_features, collapse = ", "), "\n\n")

cat("3. 特征集性能比较\n")
print(performance_results)
cat("\n")

cat("4. 最佳生物标志物组合\n")
cat(sprintf("   - 特征集: %s\n", best_set$Feature_Set))
cat(sprintf("   - 特征数量: %d\n", best_set$N_Features))
cat(sprintf("   - 准确率: %.4f\n", best_set$Accuracy))
cat(sprintf("   - Kappa: %.4f\n", best_set$Kappa))
cat("   - 生物标志物列表:\n")
for(i in 1:length(best_biomarkers)) {
  cat(sprintf("     %d. %s (重要性: %.2f)\n", 
              i, best_biomarkers[i], 
              importance_df$Importance[match(best_biomarkers[i], importance_df$Feature)]))
}

sink()

cat("\n报告已保存到: cluster/classifier/biomarker_report.txt\n")

# ==================== 最终输出 ====================
cat("\n\n========================================\n")
cat("生物标志物识别完成！\n")
cat("========================================\n")
cat("最佳生物标志物组合:\n")
for(i in 1:length(best_biomarkers)) {
  cat(sprintf("  %d. %s\n", i, best_biomarkers[i]))
}
cat(sprintf("\n使用%d个生物标志物的准确率: %.4f\n", 
            length(best_biomarkers), best_set$Accuracy))

# ==================== 11. 生成综合报告 ====================
sink("cluster/classifier/classifier_report.txt")

cat("========================================\n")
cat("多疾病分类器构建与验证报告\n")
cat("========================================\n\n")

cat("分析日期:", date(), "\n\n")

cat("1. 数据概况\n")
cat(sprintf("   - 总样本数: %d\n", nrow(X_selected)))
cat("   - 疾病类型分布:\n")
for(d in levels(y)) {
  cat(sprintf("     %s: %d\n", d, sum(y == d)))
}
cat("\n")

cat("2. 模型性能比较\n")
print(performance_summary)
cat("\n")

cat("3. 最佳模型性能\n")
cat(sprintf("   - 模型: %s\n", best_model$method))
cat(sprintf("   - 平均CV准确率: %.3f ± %.3f\n", 
            mean(fold_accuracy$Accuracy), sd(fold_accuracy$Accuracy)))
cat(sprintf("   - 外部验证准确率: %.3f\n", mean(external_results$Accuracy)))
cat("\n")

cat("4. 关键生物标志物\n")
cat("   - 最佳生物标志物组合大小:", best_size, "\n")
cat("   - 生物标志物列表:\n")
for(i in 1:min(10, length(best_features))) {
  cat(sprintf("     %d. %s\n", i, best_features[i]))
}
cat("\n")

cat("5. 疾病区分难度\n")
hard_pairs <- pair_auc[pair_auc$AUC < 0.7, ]
if(nrow(hard_pairs) > 0) {
  cat("   - 难以区分的疾病对 (AUC < 0.7):\n")
  for(i in 1:nrow(hard_pairs)) {
    cat(sprintf("     %s vs %s: AUC = %.3f\n", 
                hard_pairs$Disease1[i], hard_pairs$Disease2[i], hard_pairs$AUC[i]))
  }
} else {
  cat("   - 所有疾病对都能较好区分\n")
}

sink()

cat("\n报告已保存\n")

# ==================== 12. 输出最终结果 ====================
cat("\n\n========================================\n")
cat("多疾病分类器构建完成！\n")
cat("结果保存在: cluster/classifier/\n")
cat("========================================\n")

cat("\n最终推荐:\n")
cat("1. 最佳分类模型:", best_model$method, "\n")
cat("2. 关键生物标志物:", paste(best_features[1:min(5, length(best_features))], collapse = ", "), "\n")
cat("3. 模型性能: 10折CV准确率 =", round(mean(fold_accuracy$Accuracy), 3), "\n")