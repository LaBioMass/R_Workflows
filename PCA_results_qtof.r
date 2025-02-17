# Carregar os pacotes necessários
library(MSnbase)
library(xcms)
library(ggplot2)

# Ler e processar os dados com xcms
data_dir <- "C:/Users/vinic/Desktop/INFUSÕES/AQUOSO/NEGATIVO"
files <- sort(list.files(data_dir, pattern = "\\.mzXML$", full.names = TRUE))
print(files)
raw_data <- readMSData(files, mode = "onDisk")

#Ver a quantidade de espectros
table(msLevel(raw_data))

#ver o perfil de alguns espectros
first_spectrum <- spectra(raw_data)[[1]]
some_spectrum <- spectra(raw_data)[[555]]

#informações dos espectros
View(first_spectrum)
View(some_spectrum)

plot(mz(first_spectrum), intensity(first_spectrum), type = "h",
     main = "Primeiro Espectro Original", xlab = "m/z", ylab = "Intensidade")
plot(mz(some_spectrum), intensity(first_spectrum), type = "h",
     main = "Algum Espectro Original", xlab = "m/z", ylab = "Intensidade")

# Converter para modo centroid (caso necessário)
raw_data <- pickPeaks(raw_data, refineMz = "kNeighbors", k = 2)

first_spectrum_centroid <- spectra(raw_data)[[1]]

# Plotar novamente para conferir
plot(mz(first_spectrum_centroid), intensity(first_spectrum_centroid), type = "h",
     main = "Espectro Após pickPeaks()", xlab = "m/z", ylab = "Intensidade")

View(first_spectrum_centroid)

# Criar o Base Peak Chromatogram (BPC)
bpis <- chromatogram(raw_data, aggregationFun = "max")

#Plote um cromatograma total (TIC)
tic <- chromatogram(raw_data, aggregationFun = "sum")
plot(tic)

# Criar uma lista com os picos base de todas as amostras
base_peak_list <- lapply(seq_along(bpis), function(i) {
  idx <- which.max(intensity(bpis[[i]]))  
  data.frame(
    sample = i,
    base_peak_rt = rtime(bpis[[i]])[idx], 
    base_peak_intensity = intensity(bpis[[i]])[idx]
  )
})

# Juntar tudo em um único data frame
base_peak_df <- do.call(rbind, base_peak_list)

# Calcular a média dos tempos de retenção e intensidades
mean_base_peak_rt <- mean(base_peak_df$base_peak_rt)
mean_base_peak_intensity <- mean(base_peak_df$base_peak_intensity)
head(base_peak_df)

# Exibir os resultados
print(mean_base_peak_rt)
print(mean_base_peak_intensity)

#Ajustar os Parâmetros Automaticamente
cwp <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 10)
xset <- findChromPeaks(raw_data, param = cwp)

# Filtrar picos com intensidade > 51988.83 (baseado no mean peak base)
chrom_peaks <- chromPeaks(xset)
filtered_peaks <- chrom_peaks[chrom_peaks[, "into"] > 51988.83, ]
head(chrom_peaks)

# Agrupar picos cromatográficos
xset_grouped <- groupChromPeaks(xset, param = PeakDensityParam(sampleGroups = rep(1, length(files))))

# Verificar as definições de features após o agrupamento
feature_defs <- featureDefinitions(xset_grouped)
head(feature_defs)

# Criar a matriz de m/z por amostras

feature_matrix <- featureValues(xset_grouped, value = "into")
feature_matrix[is.na(feature_matrix)] <- 0
feature_matrix <- feature_matrix[!apply(feature_matrix, 1, function(x) all(x == 0)), ]

#Renomeando (de acordo com sua necessidade)
new_sample_names <- rep(LETTERS[1:13], each = 3)  
new_sample_names <- paste0(new_sample_names, rep(1:3, 13))  

# Verifique se o número de amostras é 39
length(new_sample_names)  

# Substituir os nomes das amostras (colunas 10 a 48)
colnames(feature_matrix)[10:48] <- new_sample_names

# Verificar se a alteração foi aplicada corretamente
head(colnames(feature_matrix)[10:48])

# Obter os valores de m/z das features
mz_values <- featureDefinitions(xset_grouped)$mzmed
rownames(feature_matrix) <- as.character(round(mz_values, 4))
rownames(feature_matrix) <- gsub("^X", "", rownames(feature_matrix))

# Verificar o tipo de mz_values
class(mz_values)

# Visualizar a matriz no RStudio
View(feature_matrix)

# Supondo que as 6 primeiras amostras são os brancos
blank_values <- feature_matrix[,1:6]  # Extrair os brancos
View(blank_values)

# Calcular o maior valor dos brancos para cada m/z (linha)
max_blank_values <- apply(blank_values, 1, max)

# Subtrair os valores máximos dos brancos das amostras e do controle de qualidade
feature_matrix_adjusted <- feature_matrix[, 7:ncol(feature_matrix)] - max_blank_values

# Substituir valores negativos por zero
feature_matrix_adjusted[feature_matrix_adjusted < 0] <- 0

# Substituir valores menores que 1000 por zero
feature_matrix_adjusted[feature_matrix_adjusted < 1000] <- 0

# Remover as linhas que possuem apenas zeros
feature_matrix_adjusted <- feature_matrix_adjusted[!apply(feature_matrix_adjusted, 1, function(x) all(x == 0)), ]

# Verificar a matriz ajustada (verifique se os valores de m/z agora têm as intensidades ajustadas)
View(feature_matrix_adjusted)

# Calcular a média de cada coluna, ignorando os zeros
column_means <- apply(feature_matrix_adjusted, 2, function(x) mean(x[x > 0], na.rm = TRUE))

# Substituir os valores zero pela média da coluna correspondente
feature_matrix_normalized <- feature_matrix_adjusted
for (j in seq_along(column_means)) {
  feature_matrix_normalized[feature_matrix_normalized[, j] == 0, j] <- column_means[j]}

# Aplicar o PCA (com controle de qualidade)
pca_matrix_qc <- t(feature_matrix_normalized)
pca_result_qc <- prcomp(pca_matrix_qc, center = TRUE, scale. = TRUE)

# Ver o resumo dos scores para todos os PCs
summary(pca_result_qc)

# Criar DataFrame com os scores do PCA
scores_df_qc <- as.data.frame(pca_result_qc$x)
scores_df_qc$Sample <- rownames(scores_df_qc)

# Plotar o PCA
ggplot(scores_df_qc, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 4, color = "blue") +
  geom_text(vjust = -1) +
  theme_minimal() +
  ggtitle("PCA - Análise com Controle de Qualidade (QC)") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_qc)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_qc)$importance[2, 2], 1), "%)"))

#Criar o Dendograma/k-means 

library(cluster)

set.seed(123)
wcss <- sapply(1:10, function(k) {
  kmeans(scores_df_qc[, c("PC1", "PC2")], centers = k, nstart = 25)$tot.withinss
})

print(wcss)  # Checar se os valores são todos zeros
plot(1:10, wcss, type = "b", pch = 19, col = "blue",
     xlab = "Número de clusters", ylab = "Soma dos Quadrados Dentro dos Clusters (WCSS)",
     main = "Método do Cotovelo")

set.seed(123) # Garantir reprodutibilidade
k <- 4 # escolhemos clusters com base no método do cotovelo
kmeans_result <- kmeans(scores_df_qc[, c("PC1", "PC2")], centers = k, nstart = 25)

# Adicionar os clusters ao DataFrame
scores_df_qc$Cluster <- as.factor(kmeans_result$cluster)

ggplot(scores_df_qc, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1) +
  stat_ellipse(geom = "polygon", alpha = 0.2)+
  theme_minimal() +
  ggtitle("K-means Clustering no PCA (Sem Controle de Qualidade)") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_qc)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_qc)$importance[2, 2], 1), "%)")) +
  scale_color_manual(values = c("red", "blue", "green","yellow"))  # Ajuste de cores conforme necessário


# Remover as amostras de controle de qualidade
feature_matrix_without_qc <- feature_matrix_normalized[, -c(1:3)]
View(feature_matrix_without_qc)

# Aplicar o PCA (sem controle de qualidade)
pca_matrix_no_qc <- t(feature_matrix_without_qc)
pca_result_no_qc <- prcomp(pca_matrix_no_qc, center = TRUE, scale. = TRUE)

# Criar DataFrame com os scores do PCA
scores_df_no_qc <- as.data.frame(pca_result_no_qc$x)
scores_df_no_qc$Sample <- rownames(scores_df_no_qc)

# Plotar o PCA
ggplot(scores_df_no_qc, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 4, color = "blue") +
  geom_text(vjust = -1) +
  theme_minimal() +
  ggtitle("PCA - Análise Sem Controle de Qualidade (QC)") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_no_qc)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_no_qc)$importance[2, 2], 1), "%)"))

# Calcular a soma dos quadrados dentro dos clusters para diferentes k
library(cluster)

set.seed(123)
wcss <- sapply(1:10, function(k) {
  kmeans(scores_df_no_qc[, c("PC1", "PC2")], centers = k, nstart = 25)$tot.withinss
})

print(wcss)  # Checar se os valores são todos zeros
plot(1:10, wcss, type = "b", pch = 19, col = "blue",
     xlab = "Número de clusters", ylab = "Soma dos Quadrados Dentro dos Clusters (WCSS)",
     main = "Método do Cotovelo")

set.seed(123) # Garantir reprodutibilidade
k <- 4 # escolhemos clusters com base no método do cotovelo
kmeans_result <- kmeans(scores_df_no_qc[, c("PC1", "PC2")], centers = k, nstart = 25)

# Adicionar os clusters ao DataFrame
scores_df_no_qc$Cluster <- as.factor(kmeans_result$cluster)

ggplot(scores_df_no_qc, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1) +
  stat_ellipse(geom = "polygon", alpha = 0.2)+
  theme_minimal() +
  ggtitle("K-means Clustering no PCA (Sem Controle de Qualidade)") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_no_qc)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_no_qc)$importance[2, 2], 1), "%)")) +
  scale_color_manual(values = c("red", "blue", "green","yellow"))  # Ajuste de cores conforme necessário

# Calcular a matriz de distâncias (usando a matriz PCA reduzida)
dist_matrix <- dist(scores_df_no_qc[, c("PC1", "PC2", "PC3", "PC4")])

# Aplicar o método de agrupamento hierárquico (usando "ward.D2")
hclust_result <- hclust(dist_matrix, method = "ward.D2")

# Plotar o dendrograma
plot(hclust_result, main = "Dendrograma - Agrupamento Hierárquico", 
     xlab = "Amostras", sub = "", cex = 0.8)

rect.hclust(hclust_result, k = 4, border = "red")

loadings_df <- as.data.frame(pca_result_no_qc$rotation)  # Loadings dos PCs
loadings_df$Feature <- rownames(loadings_df)  # Adiciona nomes das variáveis
head(loadings_df)

loadings_pc1_pc2 <- loadings_df[, c("PC1", "PC2", "Feature")]

biplot(pca_result_no_qc, scale = 0, col = c("blue", "red"))

ggplot(loadings_pc1_pc2, aes(x = PC1, y = PC2, label = Feature)) +
  geom_point(color = "red", size = 3) +
  geom_text(vjust = -1, size = 3) +
  theme_minimal() +
  ggtitle("Loadings - PC1 vs PC2") +
  xlab("PC1 Loadings") +
  ylab("PC2 Loadings")

library(ggplot2)
library(ggrepel)

# Criar um dataframe com os loadings
loadings_df <- as.data.frame(pca_result_no_qc$rotation)
loadings_df$Feature <- rownames(loadings_df)

# Selecionar apenas as variáveis com maior impacto (exemplo: top 10)
top_loadings <- loadings_df[order(abs(loadings_df$PC1) + abs(loadings_df$PC2), decreasing = TRUE), ]
top_loadings <- head(top_loadings, 10)  # Pegar as 10 principais variáveis

# Criar o biplot simplificado
loadings_df <- as.data.frame(pca_result_no_qc$rotation[, 1:2]) # Apenas PC1 e PC2
loadings_df$Variable <- rownames(loadings_df)

# Criar o biplot
biplot(pca_result_no_qc, scale = 0, col = c("blue", "red"))
ggplot(loadings_pc1_pc2, aes(x = PC1, y = PC2, label = Feature)) +
  geom_point(color = "red", size = 3) +
  geom_text(vjust = -1, size = 3) +
  theme_minimal() +
  ggtitle("Loadings - PC1 vs PC2") +
  xlab("PC1 Loadings") +
  ylab("PC2 Loadings")

#PCA pela torra
View(feature_matrix_normalized)
feature_matrix_without_cb <- feature_matrix_normalized[, -c(1:12)]
View(feature_matrix_without_cb)
feature_matrix_org_values_real <- feature_matrix_without_cb[,-c(28:30)]
View(feature_matrix_org_values_real)

# Aplicar o PCA (sem controle de qualidade)
pca_org_values <- t(feature_matrix_org_values_real)
pca_result_org <- prcomp(pca_org_values, center = TRUE)

# Criar DataFrame com os scores do PCA
scores_df_org <- as.data.frame(pca_result_org$x)
scores_df_org$Sample <- rownames(scores_df_org)

# Plotar o PCA
ggplot(scores_df_org, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 4, color = "blue") +
  geom_text(vjust = -1) +
  theme_minimal() +
  ggtitle("PCA - Análise Sem Controle de Qualidade (QC)") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_org)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_org)$importance[2, 2], 1), "%)"))

# Calcular a soma dos quadrados dentro dos clusters para diferentes k

set.seed(123)
wcss <- sapply(1:10, function(k) {
  kmeans(scores_df_org[, c("PC1", "PC2")], centers = k, nstart = 25)$tot.withinss
})

print(wcss)  # Checar se os valores são todos zeros
plot(1:10, wcss, type = "b", pch = 19, col = "blue",
     xlab = "Número de clusters", ylab = "Soma dos Quadrados Dentro dos Clusters (WCSS)",
     main = "Método do Cotovelo")

set.seed(123) # Garantir reprodutibilidade
k <- 3 # Suponha que escolhemos 3 clusters com base no método do cotovelo
kmeans_result <- kmeans(scores_df_org[, c("PC1", "PC2")], centers = k, nstart = 25)

# Adicionar os clusters ao DataFrame
scores_df_org$Cluster <- as.factor(kmeans_result$cluster)

ggplot(scores_df_org, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1) +
  theme_minimal() +
  stat_ellipse(geom = "polygon", alpha = 0.2)+
  ggtitle("K-means Clustering no PCA (Sem Controle de Qualidade)") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_no_qc)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_no_qc)$importance[2, 2], 1), "%)")) +
  scale_color_manual(values = c("red", "blue", "green","yellow"))  # Ajuste de cores conforme necessário
