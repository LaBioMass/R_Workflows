# Carregar os pacotes necessários
library(MSnbase)
library(xcms)
library(ggplot2)

# Ler e processar os dados com o pacote xcms que vamos usar
data_dir <- "C:/Users/vinic/Desktop/INFUSÕES/AQUOSO/NEGATIVO/mzML"
files <- sort(list.files(data_dir, pattern = "\\.mzML", full.names = TRUE))
print(files)
raw_data <- readMSData(files, mode = "onDisk")

#Ver a quantidade de espectros e o perfil de cada (1 = ms1; 2 = ms2 e assim por diante)
table(msLevel(raw_data))

#ver o perfil de alguns espectros 
#se estiver na dúvida se for centroid ou não (pular para o bpis se for centroid)

first_spectrum <- spectra(raw_data)[[1]]
X = #algum numero que queria ver
some_spectrum <- spectra(raw_data)[[X]]

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

# Calcular a média dos tempos de retenção e intensidades (caso queira usar as medias como filtro)
mean_base_peak_rt <- mean(base_peak_df$base_peak_rt)
mean_base_peak_intensity <- mean(base_peak_df$base_peak_intensity)
head(base_peak_df)

# Exibir os resultados
print(mean_base_peak_rt)
print(mean_base_peak_intensity)

# Filtrar apenas as amostras de branco (as 6 primeiras amostras)
#é recomendável que na hora de fazer seu diretorio os brancos sejam as primeiras amostras
base_peak_branco <- base_peak_df[base_peak_df$sample <= 6, ]

# Calcular a média do tempo de retenção e intensidade apenas dos brancos
mean_base_peak_rt_branco <- mean(base_peak_branco$base_peak_rt)
mean_base_peak_intensity_branco <- mean(base_peak_branco$base_peak_intensity)

# Exibir os resultados
print(mean_base_peak_rt_branco)
print(mean_base_peak_intensity_branco)

# Visualizar as primeiras linhas da tabela com os brancos
head(base_peak_branco)

#Ajustar os Parâmetros COM TESTES 

#primeiro testar o snthresh
# ** O parâmetro snthresh define o limiar mínimo da relação sinal-ruído (S/N) para um pico ser detectado 
# Apenas picos com S/N maior que snthresh serão considerados válidos.

cwp_low <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 5, fitgauss = TRUE)
cwp_mid <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 10, fitgauss = TRUE)
cwp_high <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 20, fitgauss = TRUE)

xset_low <- findChromPeaks(raw_data, param = cwp_low)
xset_mid <- findChromPeaks(raw_data, param = cwp_mid)
xset_high <- findChromPeaks(raw_data, param = cwp_high)

nrow(chromPeaks(xset_low))  # Com snthresh = 5
nrow(chromPeaks(xset_mid))  # Com snthresh = 10
nrow(chromPeaks(xset_high)) # Com snthresh = 20

#voce pode testar duas configurações até o fim e ver como isso afeta seu resultado estatistico
#mas no caso achei pertinente 10

#TESTAR O NOISE (USAR O THRESH CORRETO DAI)
#ajustar valores de noise conforme sua necessidade e o que tem no cromatograma

cwp_default <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 10, fitgauss = TRUE)
cwp_noise1 <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 10, noise = 500, fitgauss = TRUE)
cwp_noise2 <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 10, noise = 5000, fitgauss = TRUE)

xset_default <- findChromPeaks(raw_data, param = cwp_default)
xset_noise1 <- findChromPeaks(raw_data, param = cwp_noise1)
xset_noise2 <- findChromPeaks(raw_data, param = cwp_noise2)

nrow(chromPeaks(xset_default)) # Sem filtro de noise
nrow(chromPeaks(xset_noise1))  # Com noise = 500
nrow(chromPeaks(xset_noise2))  # Com noise = 5000

#TESTAR A LARGURA DO PICO
#USAR O THRESH CORRETO E O NOISE CORRETO
# Definir diferentes configurações de peakwidth
cwp_narrow <- CentWaveParam(ppm = 10, peakwidth = c(2, 10), snthresh = 10, noise = 500 , fitgauss = TRUE)
cwp_medium <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 10, noise = 500, fitgauss = TRUE)
cwp_wide   <- CentWaveParam(ppm = 10, peakwidth = c(10, 50), snthresh = 10, noise = 500, fitgauss = TRUE)

# Aplicar nos dados
xset_narrow <- findChromPeaks(raw_data, param = cwp_narrow)
xset_medium <- findChromPeaks(raw_data, param = cwp_medium)
xset_wide   <- findChromPeaks(raw_data, param = cwp_wide)

# Comparar número de picos detectados
cat("Picos detectados (peakwidth 2-10s):", nrow(chromPeaks(xset_narrow)), "\n")
cat("Picos detectados (peakwidth 5-20s):", nrow(chromPeaks(xset_medium)), "\n")
cat("Picos detectados (peakwidth 10-50s):", nrow(chromPeaks(xset_wide)), "\n")

#ACHO QUE UM AJUSTE MAIS FINO DE 4 A 25
#OU ALGO ENTRE 5 E 30 PODE SER ÚTIL

#prefiltragem 

# Diferentes configurações de prefilter
cwp_low_prefilter  <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 10, noise = 500, prefilter = c(3, 100), fitgauss = TRUE)
cwp_medium_prefilter <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 10, noise = 500, prefilter = c(5, 500), fitgauss = TRUE)
cwp_strict_prefilter <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 10, noise = 500, prefilter = c(10, 1000),fitgauss = TRUE)

# Aplicar nos dados
xset_low_prefilter  <- findChromPeaks(raw_data, param = cwp_low_prefilter)
xset_medium_prefilter <- findChromPeaks(raw_data, param = cwp_medium_prefilter)
xset_strict_prefilter <- findChromPeaks(raw_data, param = cwp_strict_prefilter)

# Comparar número de picos detectados
cat("Picos detectados (prefilter 3,100):", nrow(chromPeaks(xset_low_prefilter)), "\n")
cat("Picos detectados (prefilter 5,500):", nrow(chromPeaks(xset_medium_prefilter)), "\n")
cat("Picos detectados (prefilter 10,1000):", nrow(chromPeaks(xset_strict_prefilter)), "\n")

#tempo de retenção (nao esquecer de carregar o setup escolhido)
n <- 3
I <- 100
cwp_setup <- CentWaveParam(ppm = 10, peakwidth = c(5, 20), snthresh = 10, prefilter = c(n, I), noise = 100,  fitgauss = FALSE)
xset_setup <- findChromPeaks(raw_data, param = cwp_setup)
nrow(chromPeaks(xset_setup))
table(chromPeaks(xset_setup)[, "sample"])

str(xset_setup)
head(chromPeaks(xset_setup))

#Agrupar os grupos 

pdp <- PeakDensityParam(
  sampleGroups = rep(1, length(sampleNames(raw_data))), 
  bw = 15,
  minFraction = 0.01,
  minSamples = 1)

# Agrupar picos
xset_grouped <- groupChromPeaks(xset_setup, param = pdp)
head(xset_grouped)

# Verifique se os picos estão sendo agrupados corretamente
table(chromPeaks(xset_grouped)[, "sample"])

hist(chromPeaks(xset_grouped)[, "rt"], breaks = 50, main = "Distribuição do Tempo de Retenção", xlab = "Tempo de Retenção (s)")

pdp_adjust <- PeakGroupsParam(
  minFraction = 0.5,  # Mínimo de 50% das amostras precisam ter o pico para ser usado no ajuste
  smooth = "loess",   # Ajuste suavizado LOESS
  span = 0.5          # Parâmetro de suavização
)

# Aplicar a correção do tempo de retenção
xset_adjusted <- adjustRtime(xset_grouped, param = pdp_adjust)

# Visualizar os ajustes no tempo de retenção
head(adjustedRtime(xset_adjusted))


# Filtrar picos com intensidade > 5571.167 (baseado no mean base peak do branco)
chrom_peaks <- chromPeaks(xset_grouped)
filtered_peaks <- chrom_peaks[chrom_peaks[, "into"] > 5571.167, ]
head(chrom_peaks)

# Verificar as definições de features após o agrupamento
feature_defs <- featureDefinitions(xset_grouped)
head(feature_defs)

#A partir daqui temos a parte de quimiometria
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

#Normalização por função loess
install.packages(xMSanalyzer)
library(xMSanalyzer)
feature_matrix_normalized <- loess.fit(feature_matrix, qc_samples = 7:9)

qc_median <- apply(feature_matrix[, 7:9], 1, median, na.rm = TRUE)
feature_matrix_normalized <- sweep(feature_matrix, 1, qc_median, "/")
feature_matrix_normalized[is.na(feature_matrix_normalized)] <- 0

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
  ylab(paste0("PC2 (", round(100 * summary(pca_result_qc)$importance[2, 2], 1), "%)")) +
  theme(
    plot.title = element_text(hjust = 0.5),  
    panel.grid = element_blank()             
  )

#Criar o Kmeans 

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
k <- 3 # escolhemos clusters com base no método do cotovelo
kmeans_result <- kmeans(scores_df_qc[, c("PC1", "PC2")], centers = k, nstart = 25)

# Adicionar os clusters ao DataFrame
scores_df_qc$Cluster <- as.factor(kmeans_result$cluster)

ggplot(scores_df_qc, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1) +
  stat_ellipse(geom = "polygon", alpha = 0.2)+
  theme_minimal() +
  ggtitle("K-means Clustering no PCA (Com Controle de Qualidade)") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_qc)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_qc)$importance[2, 2], 1), "%)")) +
  scale_color_manual(values = c("red", "blue", "yellow"))  

# Calcular a matriz de distâncias (usando a matriz PCA reduzida)
dist_matrix <- dist(scores_df_qc[, c("PC1", "PC2", "PC3", "PC4")])

# Aplicar o método de agrupamento hierárquico (usando "ward.D2")
hclust_result <- hclust(dist_matrix, method = "ward.D2")

# Plotar o dendrograma
plot(hclust_result, main = "Dendrograma - Agrupamento Hierárquico", 
     xlab = "Amostras", sub = "", cex = 0.8)

rect.hclust(hclust_result, k = 3, border = "red")

#Loadings

loadings_df <- as.data.frame(pca_result_qc$rotation)  # Loadings dos PCs
loadings_df$Feature <- rownames(loadings_df)  # Adiciona nomes das variáveis
head(loadings_df)

loadings_pc1_pc2 <- loadings_df[, c("PC1", "PC2", "Feature")]

biplot(pca_result_qc, scale = 0, col = c("blue", "red"))

ggplot(loadings_pc1_pc2, aes(x = PC1, y = PC2, label = Feature)) +
  geom_point(color = "red", size = 3) +
  geom_text(vjust = -1, size = 3) +
  theme_minimal() +
  ggtitle("Loadings - PC1 vs PC2") +
  xlab("PC1 Loadings") +
  ylab("PC2 Loadings")

# Selecionar os 10 m/z com maior contribuição (pela soma dos valores absolutos de PC1 e PC2)
loadings_pc1_pc2$Importance <- abs(loadings_pc1_pc2$PC1) + abs(loadings_pc1_pc2$PC2)
top10_loadings <- loadings_pc1_pc2[order(-loadings_pc1_pc2$Importance), ][1:10, ]

# Criar o gráfico com apenas os 10 m/z mais importantes
ggplot(top10_loadings, aes(x = PC1, y = PC2, label = Feature)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm")), color = "red") +
  geom_text(vjust = -0.5, size = 4) +
  theme_minimal() +
  ggtitle("Loadings - Top 10 Contribuintes (PC1 vs PC2)") +
  xlab("PC1 Loadings") +
  ylab("PC2 Loadings")


#ANALISE SEM O QC

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
k <- 3 # escolhemos clusters com base no método do cotovelo
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
  scale_color_manual(values = c("red", "blue", "green"))  # Ajuste de cores conforme necessário

# Calcular a matriz de distâncias (usando a matriz PCA reduzida)
dist_matrix <- dist(scores_df_no_qc[, c("PC1", "PC2", "PC3", "PC4")])

# Aplicar o método de agrupamento hierárquico (usando "ward.D2")
hclust_result <- hclust(dist_matrix, method = "ward.D2")

# Plotar o dendrograma
plot(hclust_result, main = "Dendrograma - Agrupamento Hierárquico", 
     xlab = "Amostras", sub = "", cex = 0.8)

rect.hclust(hclust_result, k = 3, border = "red")

#Loadings

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

# Selecionar os 10 m/z com maior contribuição (pela soma dos valores absolutos de PC1 e PC2)
loadings_pc1_pc2$Importance <- abs(loadings_pc1_pc2$PC1) + abs(loadings_pc1_pc2$PC2)
top10_loadings <- loadings_pc1_pc2[order(-loadings_pc1_pc2$Importance), ][1:10, ]

# Criar o gráfico com apenas os 10 m/z mais importantes
ggplot(top10_loadings, aes(x = PC1, y = PC2, label = Feature)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm")), color = "red") +
  geom_text(vjust = -0.5, size = 4) +
  theme_minimal() +
  ggtitle("Loadings - Top 10 Contribuintes (PC1 vs PC2)") +
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
pca_result_org <- prcomp(pca_org_values, center = TRUE, scale = TRUE)

# Criar DataFrame com os scores do PCA
scores_df_org <- as.data.frame(pca_result_org$x)
scores_df_org$Sample <- rownames(scores_df_org)

# Plotar o PCA
ggplot(scores_df_org, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 4, color = "blue") +
  geom_text(vjust = -1) +
  theme_minimal() +
  ggtitle("PCA - Análise Amostras Orgânicas") +
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
  ggtitle("K-means Clustering no PCA (Amostras Orgânicas)") +
  xlab(paste0("PC1 (", round(100 * summary(pca_result_no_qc)$importance[2, 1], 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca_result_no_qc)$importance[2, 2], 1), "%)")) +
  scale_color_manual(values = c("red", "blue", "yellow"))  # Ajuste de cores conforme necessário

# Calcular a matriz de distâncias (usando a matriz PCA reduzida)
dist_matrix <- dist(scores_df_org[, c("PC1", "PC2", "PC3", "PC4")])

# Aplicar o método de agrupamento hierárquico (usando "ward.D2")
hclust_result <- hclust(dist_matrix, method = "ward.D2")

# Plotar o dendrograma
plot(hclust_result, main = "Dendrograma - Agrupamento Hierárquico", 
     xlab = "Amostras", sub = "", cex = 0.8)

rect.hclust(hclust_result, k = 3, border = "red")

#Loadings

loadings_df <- as.data.frame(pca_result_org$rotation)  # Loadings dos PCs
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

# Selecionar os 10 m/z com maior contribuição (pela soma dos valores absolutos de PC1 e PC2)
loadings_pc1_pc2$Importance <- abs(loadings_pc1_pc2$PC1) + abs(loadings_pc1_pc2$PC2)
top10_loadings <- loadings_pc1_pc2[order(-loadings_pc1_pc2$Importance), ][1:10, ]

# Criar o gráfico com apenas os 10 m/z mais importantes
ggplot(top10_loadings, aes(x = PC1, y = PC2, label = Feature)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm")), color = "red") +
  geom_text(vjust = -0.5, size = 4) +
  theme_minimal() +
  ggtitle("Loadings - Top 10 Contribuintes (PC1 vs PC2)") +
  xlab("PC1 Loadings") +
  ylab("PC2 Loadings")

