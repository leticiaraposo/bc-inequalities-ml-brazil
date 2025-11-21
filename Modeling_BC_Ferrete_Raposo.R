# ==========================================================================================
# Título: Spatially Varying Socioeconomic Determinants of Breast Cancer Mortality in Brazil
# Autores: Livia Ferrete, Letícia Raposo
# Instituição: Universidade Federal do Estado do Rio de Janeiro (UNIRIO)
# Etapa 4: modelagens RF, GRF e análise SHAP
# Data: Setembro de 2025
# ==========================================================================================

# --- INSTALAÇÃO DE PACOTES ----
install.packages(c("sp", "sf", "randomForest", "caret", "SpatialML",
                   "tmap", "ggplot2", "dplyr", "GWmodel", "readxl"))



# --- CARREGAMENTO DE PACOTES ----
library(sp)          # Para objetos espaciais (SpatialPointsDataFrame)
library(sf)          # Para Simple Features (modelagem espacial, compatível com sp)
library(randomForest) # Para modelos Random Forest
library(caret)       # Para machine learning
library(SpatialML)   # Para aprendizado de máquina espacial (GRF)
library(tmap)        # Para visualização de dados espaciais (mapas)
library(ggplot2)     # Para visualização de dados em geral
library(dplyr)       # Para manipulação de dados
library(GWmodel)     # Outro pacote para modelos geograficamente ponderados
library(readxl)      # Para leitura de arquivos Excel



# --- IMPORTAÇÃO DA BASE ---
dados <- readRDS("base_cancer_mama_imputados_sem_correlacao_alta_ri.rds")



# --- PADRONIZAÇÃO (Z-SCORE) DAS VARIÁVEIS NUMÉRICAS EXPLICATIVAS ---

# Identificação das colunas numéricas a serem padronizadas (ajustar conforme necessário)
colunas_para_padronizar <- colnames(dados)[5:27]

# Cálculo dos parâmetros de padronização (média e desvio padrão) com base no TREINO
parametros_padronizacao <- lapply(colunas_para_padronizar, function(col) {
  list(media = mean(dados[[col]], na.rm = TRUE),
       dp    = sd(dados[[col]], na.rm = TRUE))
})
names(parametros_padronizacao) <- colunas_para_padronizar

# Padronização 
dados_padronizado <- dados
for (col in colunas_para_padronizar) {
  media <- parametros_padronizacao[[col]]$media
  dp    <- parametros_padronizacao[[col]]$dp
  dados_padronizado[[col]] <- (dados[[col]] - media) / dp
}



# --- PREPARAÇÃO PARA MODELAGEM ESPACIAL  ---

# Conversão para objeto espacial (SpatialPointsDataFrame)
dados_sp <- dados_padronizado
coordinates(dados_sp) <- ~X + Y
proj4string(dados_sp) <- CRS("+proj=longlat +datum=WGS84")

# Coordenadas separadas para métodos que exigem matriz de coordenadas
coords <- cbind(dados_padronizado$Y, dados_padronizado$X)



# --- DEFINIÇÃO DA FÓRMULA DO MODELO ---

# Definir variável resposta
variavel_resposta <- "TAXA_AJUSTADA"

# Excluir variáveis não preditoras da fórmula
variaveis_para_excluir <- c("regiao_geografica_imediata", "X", "Y", variavel_resposta)

# Identificação das variáveis preditoras
variaveis_preditoras <- setdiff(colnames(dados_padronizado), variaveis_para_excluir)

# Construção dinâmica da fórmula para modelagem
formula_string <- paste(variavel_resposta, "~", paste(variaveis_preditoras, collapse = " + "))
form <- as.formula(formula_string)





# =============================================================================
# RANDOM FOREST GLOBAL - OTIMIZAÇÃO DE PARÂMETROS (ntree x mtry)
# Com paralelização, validação cruzada e seleção por RMSE
# =============================================================================

library(ranger)
library(doParallel)
library(caret)

# --- DEFINIÇÃO DOS PARÂMETROS A SEREM AVALIADOS ---

# Número de árvores a ser testado
ntree_vals <- c(100, 500, 1000, 1500, 2000)

# Número de preditores disponíveis (excluindo a variável resposta)
num_preditores <- length(all.vars(form)) - 1

# Valores de mtry (número de variáveis sorteadas por árvore)
mtry_vals_to_test <- seq(1, num_preditores, by = 2)
mtry_vals_to_test <- unique(sort(pmin(mtry_vals_to_test, num_preditores)))

message("Número de preditores:", num_preditores)
message("Variações de mtry a testar: ", paste(mtry_vals_to_test, collapse = ", "))
message("Variações de ntree a testar: ", paste(ntree_vals, collapse = ", "))

# --- CONFIGURAÇÃO DE PARALELIZAÇÃO ---
num_cores <- max(1, floor(0.75 * parallel::detectCores(logical = FALSE)))
cl <- makeCluster(num_cores)
registerDoParallel(cl)
message(paste("Utilizando", num_cores, "núcleos para paralelização."))

# --- OBJETO PARA ARMAZENAR RESULTADOS ---
resultados_rf <- data.frame()

# --- LAÇO DE OTIMIZAÇÃO VIA FOREACH PARALELIZADO ---
resultados_lista <- foreach(nt = ntree_vals, .combine = rbind, .packages = c("caret", "ranger")) %dopar% {
  
  message(paste0("Rodando ntree = ", nt))
  
  grid <- expand.grid(
    mtry = mtry_vals_to_test,
    splitrule = "variance",
    min.node.size = 5
  )
  
  control <- trainControl(
    method = "cv",
    number = 5,
    allowParallel = TRUE,
    verboseIter = FALSE
  )
  
  set.seed(123)
  modelo <- train(
    form,
    data = dados_padronizado,
    method = "ranger",
    trControl = control,
    tuneGrid = grid,
    num.trees = nt,
    num.threads = 1,
    importance = "impurity"
  )
  
  modelo$results$ntree <- nt
  return(modelo$results)
}

# --- ENCERRAMENTO DO CLUSTER ---
stopCluster(cl)
registerDoSEQ()
message("Paralelização encerrada.")

# --- RESULTADOS FINAIS ---
resultados_rf <- resultados_lista

# Seleção do melhor modelo com menor RMSE
melhor_rf <- resultados_rf[which.min(resultados_rf$RMSE), ]
best_ntree <- melhor_rf$ntree
best_mtry <- melhor_rf$mtry

# Exibição do resultado ótimo
message("===== MELHOR MODELO RANDOM FOREST =====")
print(melhor_rf)

# Visualização gráfica (opcional)
library(ggplot2)
ggplot(resultados_rf, aes(x = factor(ntree), y = RMSE, color = factor(mtry))) +
  geom_point(size = 3) +
  geom_line(aes(group = mtry), linewidth = 0.7) +
  labs(
    title = "RMSE por combinação de ntree e mtry (Random Forest)",
    x = "Número de Árvores (ntree)",
    y = "RMSE",
    color = "mtry"
  ) +
  theme_minimal()

# --- LIBERAÇÃO DE MEMÓRIA ---
gc()





# =============================================================================
# GRF – GEOGRAPHICALLY WEIGHTED RANDOM FOREST
# Otimização de Bandwidth via Validação Cruzada (RMSE)
# =============================================================================

# --- PACOTES NECESSÁRIOS ---
library(SpatialML)
library(dplyr)
library(sp)

# --- DEFINIÇÃO DOS BANDWIDTHS A SEREM TESTADOS ---
bws <- c(50, 100, 150, 200, 250, 300, 350, 400)
num_preditores <- length(all.vars(form)) - 1
bws <- bws[bws > num_preditores]

if (length(bws) == 0) stop("A sequência de valores de 'bw' está inadequada. Ajuste para valores maiores.")

# --- VALIDAÇÃO CRUZADA (5-fold) ---
num_folds_cv <- 5
set.seed(123)
n <- nrow(dados_sp)
folds_indices <- sample(rep(1:num_folds_cv, length.out = n))

# --- PREPARAÇÃO DAS VARIÁVEIS E COORDENADAS ---
response_var_name <- all.vars(form)[1]
coords <- coordinates(dados_sp)
dados_sp@data$X <- coords[, 1]
dados_sp@data$Y <- coords[, 2]

variaveis_usadas <- unique(c(all.vars(form), "X", "Y"))
dados_sp@data <- dados_sp@data[, variaveis_usadas, drop = FALSE]
dados_sp <- SpatialPointsDataFrame(coords = coords, data = dados_sp@data)

# --- OTIMIZAÇÃO DA LARGURA DE BANDA (BANDWIDTH) ---
resultados_grf_cv_list <- list()

for (bw_val in bws) {
  message(sprintf("▶ Avaliando GRF com Bandwidth = %d", bw_val))
  tempo_inicio <- Sys.time()
  rmse_folds <- numeric(num_folds_cv)
  
  for (k in 1:num_folds_cv) {
    treino_fold <- dados_sp[folds_indices != k, ]
    validacao_fold <- dados_sp[folds_indices == k, ]
    
    form_local <- as.formula(form)
    dados_modelo <- tryCatch(
      model.frame(form_local, data = treino_fold@data),
      error = function(e) {
        message(sprintf("⚠ Erro em model.frame (BW = %d, Fold = %d): %s", bw_val, k, e$message))
        return(NULL)
      }
    )
    
    if (is.null(dados_modelo)) {
      rmse_folds[k] <- NA
      next
    }
    
    modelo_grf <- tryCatch(
      grf(
        formula = form_local,
        dframe = dados_modelo,
        coords = coordinates(treino_fold),
        bw = bw_val,
        kernel = "adaptive",
        ntree = best_ntree,
        mtry = best_mtry,
        print.results = FALSE
      ),
      error = function(e) {
        message(sprintf("⚠ Erro no ajuste do GRF (BW = %d, Fold = %d): %s", bw_val, k, e$message))
        return(NULL)
      }
    )
    
    if (is.null(modelo_grf)) {
      rmse_folds[k] <- NA
      next
    }
    
    preds <- tryCatch(
      predict.grf(modelo_grf, new.data = validacao_fold@data, x.var.name = "X", y.var.name = "Y"),
      error = function(e) {
        message(sprintf("⚠ Erro na predição (BW = %d, Fold = %d): %s", bw_val, k, e$message))
        return(rep(NA, nrow(validacao_fold)))
      }
    )
    
    observado <- validacao_fold@data[[response_var_name]]
    rmse_folds[k] <- sqrt(mean((observado - preds)^2, na.rm = TRUE))
    
    rm(treino_fold, validacao_fold, modelo_grf, preds, dados_modelo)
    gc()
  }
  
  resultados_grf_cv_list[[as.character(bw_val)]] <- data.frame(
    bw = bw_val,
    RMSE_medio_CV = mean(rmse_folds, na.rm = TRUE)
  )
  
  tempo_fim <- Sys.time()
  duracao <- as.numeric(difftime(tempo_fim, tempo_inicio, units = "mins"))
  message(sprintf("✓ Bandwidth = %d finalizado em %.2f minutos", bw_val, duracao))
}

# --- RESULTADOS FINAIS ---
resultados_grf_cv <- do.call(rbind, resultados_grf_cv_list)
rownames(resultados_grf_cv) <- NULL

message("=== Resultados Finais da Validação Cruzada ===")
print(resultados_grf_cv)

melhor_bw <- resultados_grf_cv$bw[which.min(resultados_grf_cv$RMSE_medio_CV)]
message(sprintf("✔ Melhor Bandwidth (menor RMSE): %d", melhor_bw))

# --- VISUALIZAÇÃO DOS RESULTADOS ---
plot(resultados_grf_cv$bw, resultados_grf_cv$RMSE_medio_CV,
     type = "b", pch = 16, col = "blue",
     xlab = "Bandwidth (Número de Vizinhos)",
     ylab = "RMSE Médio (Validação Cruzada)",
     main = "Otimização de Bandwidth para GRF")
abline(v = melhor_bw, col = "red", lty = 2)

# =============================================================================
# AJUSTE FINAL DO MODELO GRF COM BANDWIDTH OTIMIZADO
# =============================================================================

grf_model <- grf(
  formula = form,
  dframe = dados_padronizado,
  coords = coords,
  bw = melhor_bw,
  kernel = "adaptive",
  ntree = best_ntree,
  mtry = best_mtry
)

message("Modelo GRF final ajustado com largura de banda otimizada.")

# Importância média por variável
imp_media_grf <- apply(grf_model$Local.Variable.Importance, 2, mean) %>% 
  sort(decreasing = TRUE)





# =============================================================================
# MAPAS DE IMPORTÂNCIA LOCAL - IncMSE por variável
# =============================================================================

# Carregar pacotes necessários
library(SpatialML)
library(sf)
library(dplyr)
library(tmap)       # ou use ggplot2 se preferir
library(viridis)    # para escalas de cores perceptualmente uniformes
library(tmaptools)

library(geobr)
regioes_imediatas <- read_immediate_region(year = 2019)  # ou o ano desejado
regioes_imediatas$code_immediate <- as.character(regioes_imediatas$code_immediate)

# 1. Extração da importância local das variáveis do modelo GRF
# ------------------------------------------------------------
importancia_local <- grf_model$Local.Variable.Importance  # matriz de IncMSE (linhas = regiões, colunas = variáveis)

# 2. Construir dataframe com coordenadas
# ------------------------------------------------------------
importancia_local_df <- as.data.frame(importancia_local)
importancia_local_df$X <- dados_padronizado$X
importancia_local_df$Y <- dados_padronizado$Y
importancia_local_df$regiao_imediata <- dados_padronizado$regiao_geografica_imediata

library(dplyr)
library(sf)

# Suponha que sua base tenha a coluna 'code_immediate'
importancia_local_df <- regioes_imediatas %>%
  left_join(importancia_local_df, by = c("code_immediate" = "regiao_imediata"))

regioes_com_dados <- st_join(regioes_imediatas, importancia_local_df, join = st_contains)

# Pegando as 9 variáveis com maiores médias de importância da GRF
variaveis_importancia <- names(imp_media_grf)[1:9]

library(stringr)

# Gerar e salvar os mapas
dir.create("mapas_importancia", showWarnings = FALSE)

for (variavel in variaveis_importancia) {
  
  p <- ggplot(regioes_com_dados) +
    geom_sf(aes_string(fill = variavel), color = NA) +
    scale_fill_viridis_c(option = "D", name = "Importância local") +
    theme_minimal() +
    labs(
      title = paste("Importância local -", variavel),
      subtitle = "Modelo GRF por região imediata"
    )
  
  nome_arquivo <- paste0("mapas_importancia/", str_replace_all(variavel, "[^\\w]", "_"), ".png")
  
  ggsave(nome_arquivo, plot = p, width = 8, height = 6, dpi = 300)
}

# Pacotes necessários
library(geobr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(forcats)
library(stringr)
library(patchwork)

# Transformar os dados para formato longo
dados_long <- regioes_com_dados %>%
  pivot_longer(cols = all_of(variaveis_importancia),
               names_to = "variavel", values_to = "importancia")

# Escalar todos os valores entre 0–100 com base no maior valor global
valor_max_global <- max(dados_long$importancia, na.rm = TRUE)

dados_long <- dados_long %>%
  mutate(
    importancia_esc = (importancia / valor_max_global) * 100
  )

# Classificar os valores em categorias discretas
dados_long <- dados_long %>%
  mutate(
    categoria = cut(importancia_esc,
                    breaks = c(-Inf, 5, 15, 35, 100),
                    labels = c("Low (≤5%)", 
                               "Medium low (>5 a 15%)", 
                               "Medium high (>15 a 35%)", 
                               "High (>35 a 100%)")),
    categoria = fct_rev(categoria)
  )

dados_long <- dados_long[-which(is.na(dados_long$categoria)),]


# Paleta de cores
paleta_cores <- c(
  "Low (≤5%)" = "#b4d0e8",
  "Medium low (>5 a 15%)" = "#f0edf4",
  "Medium high (>15 a 35%)" = "#de66b0", 
  "High (>35 a 100%)" = "#cf1157"
)

# Filtrando apenas as 10 variáveis definidas anteriormente
dados_long2 <- dados_long[dados_long$variavel %in% variaveis_importancia, ]

library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)  # para wrap_plots()
library(stringr)

# Criando a lista de mapas
lista_mapas <- lapply(unique(dados_long2$variavel), function(var) {
  ggplot(dados_long2 %>% filter(variavel == var)) +
    
    geom_sf(aes(fill = categoria, geometry = geom), color = "gray70", size = 0.1) +
    
    geom_sf(data = regiao, fill = NA, color = "#000", size = 0.1) + 
    scale_fill_manual(values = paleta_cores, drop = FALSE) +
    labs(title = str_wrap(var, 30), fill = "Importância relativa") +
    theme_void() +
    theme(
      plot.title      = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.title    = element_text(size = 9),
      legend.text     = element_text(size = 8)
    )
})


# Títulos em português
# novos_titulos_imp <- (c(
#   "Ensino médio completo", "Densidade demográfica", "Domicílios chefiados por mulheres",
#   "PIB per capita", "Médicos especialistas", "Domicílios com entorno pavimentado", 
#   "Domicílios com esgoto ou fossa séptica", "Trabalho excedente", "Domicílios com entorno arborizado"
# ))


# Títulos em inglês
novos_titulos_imp_eng <- (c(
  "Completed secondary education", "Population density", "Female-headed households",
  "Per capita GDP", "Specialist physicians", "Paved surroundings", 
  "Households with sewer or septic tank", "Excessive working hours", "Tree-covered surroundings"
))

for (i in seq_along(lista_mapas)) {
  lista_mapas[[i]] <- lista_mapas[[i]] +
    labs(title = novos_titulos_imp_eng[i])
}


# Extrai **todas** as caixas de legenda do primeiro mapa
legs_all <- cowplot::get_plot_component(
  lista_mapas[[3]], # antes = 1; não extraía legenda
  "guide-box",
  return_all = TRUE
)

# Seleciona **somente** a primeira legenda
legenda_unica <- legs_all[[3]] # antes == 1; não extraía legenda

# Remove legendas de TODOS os mapas
lista_sem_legenda <- lapply(lista_mapas, function(p) {
  p + theme(legend.position = "none")
})

# Combina os mapas (sem legenda)
mapas_comb <- wrap_plots(lista_sem_legenda, ncol = 3)

# Monta a figura final: mapas + legenda única abaixo
final <- cowplot::plot_grid(
  mapas_comb,
  legenda_unica,
  ncol = 1,
  rel_heights = c(1, 0.1)   # ajusta a altura da área da legenda
)

# Exibe e salva
print(final)
ggsave("Local_Importance.png", width = 20, height = 10, dpi = 300)


# IMPORTANCIA DE VARIAVEIS COM RANDOM FOREST
# Vendo as variáveis mais importantes
rf_model <- randomForest::randomForest(TAXA_AJUSTADA ~ .,
                                       mtry = best_mtry, 
                                       ntree = best_ntree,
                                       importance = T, 
                                       data = dados_padronizado[, 4:27])

# Vendo as variáveis por região
dados_padronizado_com_regiao <- regioes_imediatas[,c(1,7)] %>%
  left_join(dados_padronizado, by = c("code_immediate" = "regiao_geografica_imediata"))
dados_padronizado_com_regiao <- dados_padronizado_com_regiao[!is.na(dados_padronizado_com_regiao$TAXA_AJUSTADA), ]
dados_padronizado_com_regiao <- as.data.frame(dados_padronizado_com_regiao)

# NORTE
library(sf)

dados_norte <- dados_padronizado_com_regiao %>%
  filter(name_region == "Norte") %>%
  st_drop_geometry() %>%
  select(-geom)  # remove a coluna geom que permanece como lista

rf_model_norte <- randomForest::randomForest(
  TAXA_AJUSTADA ~ .,
  mtry = best_mtry,
  ntree = best_ntree,
  importance = T,
  data = dados_norte[, 5:28]
)

varImpPlot(rf_model_norte)

# NORDESTE
library(sf)

dados_nordeste <- dados_padronizado_com_regiao %>%
  filter(name_region == "Nordeste") %>%
  st_drop_geometry() %>%
  select(-geom)  # remove a coluna geom que permanece como lista

rf_model_nordeste <- randomForest::randomForest(
  TAXA_AJUSTADA ~ .,
  mtry = best_mtry,
  ntree = best_ntree,
  importance = T,
  data = dados_nordeste[, 5:28]
)

varImpPlot(rf_model_nordeste)

# CENTRO OESTE
library(sf)

dados_centro_oeste <- dados_padronizado_com_regiao %>%
  filter(name_region == "Centro Oeste") %>%
  st_drop_geometry() %>%
  select(-geom)  # remove a coluna geom que permanece como lista

rf_model_centro_oeste <- randomForest::randomForest(
  TAXA_AJUSTADA ~ .,
  mtry = best_mtry,
  ntree = best_ntree,
  importance = T,
  data = dados_centro_oeste[, 5:28]
)

varImpPlot(rf_model_centro_oeste)

# SUDESTE
library(sf)

dados_sudeste <- dados_padronizado_com_regiao %>%
  filter(name_region == "Sudeste") %>%
  st_drop_geometry() %>%
  select(-geom)  # remove a coluna geom que permanece como lista

rf_model_sudeste <- randomForest::randomForest(
  TAXA_AJUSTADA ~ .,
  mtry = best_mtry,
  ntree = best_ntree,
  importance = T,
  data = dados_sudeste[, 5:28]
)

varImpPlot(rf_model_sudeste)

# SUL
library(sf)

dados_sul <- dados_padronizado_com_regiao %>%
  filter(name_region == "Sul") %>%
  st_drop_geometry() %>%
  select(-geom)  # remove a coluna geom que permanece como lista

rf_model_sul <- randomForest::randomForest(
  TAXA_AJUSTADA ~ .,
  mtry = best_mtry,
  ntree = best_ntree,
  importance = T,
  data = dados_sudeste[, 5:28]
)

varImpPlot(rf_model_sul)

library(dplyr)
library(ggplot2)
library(randomForest)
library(tidytext)

extrair_importancia <- function(modelo, regiao) {
  imp <- importance(modelo, type = 1)  # %IncMSE
  df <- data.frame(
    Variavel = rownames(imp),
    Importancia = imp[, 1],
    Regiao = regiao
  )
  return(df)
}

# Gerando data.frame com todas as regiões
importancias <- bind_rows(
  extrair_importancia(rf_model, "Brasil"),
  extrair_importancia(rf_model_norte, "Norte"),
  extrair_importancia(rf_model_nordeste, "Nordeste"),
  extrair_importancia(rf_model_centro_oeste, "Centro-Oeste"),
  extrair_importancia(rf_model_sudeste, "Sudeste"),
  extrair_importancia(rf_model_sul, "Sul")
)

importancias_norm <- importancias %>%
  group_by(Regiao) %>%
  mutate(ImportanciaRelativa = 100 * Importancia / max(Importancia)) %>%
  ungroup()

top_vars <- importancias_norm %>%
  group_by(Regiao) %>%
  slice_max(ImportanciaRelativa, n = 15) %>%
  ungroup()

# Criando dicionário de variáveis para modificação dos gráficos ----
# dicionario_vars <- c("mp_cobertura_esf_2010" = "Cobertura de ESF",
#                      "mp_densidade_demografica" = "Densidade demográfica",
#                      "mp_mais_de_0.5_a_1_morador" = "Densidade domiciliar (até 1 morador por cômodo)",
#                      "mp_desemprego" = "Desemprego",
#                      "mp_domicilio_entorno_arborizado" = "Domicílios com entorno arborizado",
#                      "mp_domicilio_mulher_responsavel" = "Domicílios chefiados por mulheres",
#                      "mp_domicilio_entorno_pavimentado" = "Domicílios com entorno pavimentado",
#                      "mp_enfermeiros_2010" = "Enfermeiros",
#                      "mp_esgoto_ou_fossa" = "Domicílios com esgoto ou fossa séptica",
#                      "mp_evangelicos" = "Evangélicos",
#                      "mp_gini" = "Índice de Gini",
#                      "mp_homicidios_2010" = "Homicídios",               
#                      "mp_leitos_hosp_2010" = "Leitos hospitalares",
#                      "mp_mamografos_2010" = "Mamógrafos",
#                      "mp_medicos_gyn_ob_2010" = "Médicos especialistas",          
#                      "mp_medio_completo" = "Ensino médio completo",
#                      "mp_PCDs" = "Deficiência",
#                      "mp_pib_per_capita" = "PIB per capita",                
#                      "mp_raca_preta" = "Raça preta",
#                      "mp_recursos_destinados_saude_2010" = "Recursos destinados à saúde",
#                      "mp_suicidios_2010" = "Suicídios",
#                      "mp_trabalha_mais_44h" = "Trabalho excedente",
#                      "mp_ultrassons_2010" = "Ultrassons")


# Dicionário em inglês
dicionario_vars_eng <- c("mp_cobertura_esf_2010" = "FHS coverage",
                     "mp_densidade_demografica" = "Population density",
                     "mp_mais_de_0.5_a_1_morador" = "Low household crowding",
                     "mp_desemprego" = "Unemployment",
                     "mp_domicilio_entorno_arborizado" = "Tree-covered surroundings",
                     "mp_domicilio_mulher_responsavel" = "Female-headed households",
                     "mp_domicilio_entorno_pavimentado" = "Paved surroundings",
                     "mp_enfermeiros_2010" = "Nurses",
                     "mp_esgoto_ou_fossa" = "Households with sewer or septic tank",
                     "mp_evangelicos" = "Evangelical population",
                     "mp_gini" = "Gini Index",
                     "mp_homicidios_2010" = "Homicide rate",               
                     "mp_leitos_hosp_2010" = "Hospital beds",
                     "mp_mamografos_2010" = "Mammography units",
                     "mp_medicos_gyn_ob_2010" = "Specialist physicians",          
                     "mp_medio_completo" = "Completed secondary education",
                     "mp_PCDs" = "Disability",
                     "mp_pib_per_capita" = "Per capita GDP",                
                     "mp_raca_preta" = "Black population",
                     "mp_recursos_destinados_saude_2010" = "Own-source health expenditure",
                     "mp_suicidios_2010" = "Suicide rate",
                     "mp_trabalha_mais_44h" = "Excessive working hours",
                     "mp_ultrassons_2010" = "Ultrasound units")


top_vars <- top_vars %>% 
  mutate(Variavel = recode(Variavel, !!!dicionario_vars_eng)) # operador bang-bang-bang desempacota o vetor nomeado

top_vars$Regiao <- factor(top_vars$Regiao,
                          levels = c("Brasil", "Norte", "Nordeste", "Centro-Oeste", "Sudeste", "Sul"))

levels(top_vars$Regiao) <- c("Brazil", "North", "Northeast", "Central-West", "Southeast", "South")



# Facet wrap em 2 colunas
importancia <- ggplot(top_vars, aes(x = ImportanciaRelativa, 
                     y = reorder_within(Variavel, ImportanciaRelativa, Regiao))) +
  geom_point(size = 2) +
  facet_wrap(~ Regiao, scales = "free_y", ncol = 2) +
  scale_y_reordered() +  # função do pacote tidytext
  labs(
    x = "Relative importance (%)",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_blank()
  )


png("Variable_Importance.png", units = "cm", res = 500, width = 20, height = 30)
importancia
dev.off()





# ===============================================
# MAPA DA TAXA AJUSTADA
# ===============================================
library(geobr)
library(sf)
library(ggplot2)
library(dplyr)
library(classInt)

# Carregar as regioes com geometria
regiao <- read_region(year = 2010)

dados_mapa <- regioes_imediatas %>%
  left_join(dados, by = c("code_immediate" = "regiao_geografica_imediata"))
dados_mapa <- dados_mapa[!is.na(dados_mapa$TAXA_AJUSTADA), ]

# Criar categorias de quartis
quarts <- classIntervals(dados_mapa$TAXA_AJUSTADA, n = 4, style = "quantile")
dados_mapa$faixa <- cut(dados_mapa$TAXA_AJUSTADA,
                        breaks = quarts$brks,
                        include.lowest = TRUE)

mapa_taxa <- ggplot() +
  geom_sf(data = dados_mapa, aes(fill = faixa), color = NA) +
  geom_sf(data = regiao, fill = NA, color = "black", size = 0.8) +  # contorno das 5 regiões
  scale_fill_manual(
    name = "Age-adjusted female breast cancer\n mortality rate by 100.000 inhabitants",
    values = c("#f0edf4", "#d4b6d8", "#de66b0", "#e72989"),
    na.value = "grey80"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",             # legenda na parte debaixo
    legend.box.background = element_rect(   # borda preta ao redor
      colour = "black",
      fill = NA,
      size = 0.3
    ),
    legend.title = element_text(size = 8)
  )

png("Adjusted_Mortality_Rate.png", units = "cm", res = 500, width = 15, height = 15)
mapa_taxa
dev.off()





# ---------------------------------------------------
# MAPA DA PREVALÊNCIA + IMPORTÂNCIA DAS VARIÁVEIS
# ---------------------------------------------------

library(sf)
library(ggplot2)
library(dplyr)
library(classInt)
library(patchwork)
library(cowplot)

# Criando lista para armazenar os novos mapas combinados
lista_mapas_combinados <- list()

# Duplicando o banco para criar as faixas de prevalência
dados_mapa_prevalencia <- dados_mapa

# --- INÍCIO DO LOOP COMBINADO ---
for (variavel_atual in variaveis_importancia) {
  
  # --- Parte 1: lógica do mapa base --- #
  # Seguindo a lógica do mapa da taxa ajustada
  
  # Cria o nome da nova coluna de faixas
  novo_nome_faixa <- paste0("faixas_", variavel_atual)
  
  # Calcula os quantis para a variável atual
  quarts <- classIntervals(
    na.omit(dados_mapa_prevalencia[[variavel_atual]]), 
    n = 4, 
    style = "quantile"
  )
  
  # Cria a coluna de faixas no dataframe
  dados_mapa_prevalencia[[novo_nome_faixa]] <- cut(
    dados_mapa_prevalencia[[variavel_atual]],
    breaks = quarts$brks,
    include.lowest = TRUE
  )
  
  # Cria o objeto ggplot com a camada base azul
  mapa_base_prev <- ggplot() +
    
    # Camada de preenchimento azul (prevalência)
    geom_sf(data = dados_mapa_prevalencia, aes(fill = .data[[novo_nome_faixa]]), color = NA) +
    
    # Camada do contorno geral do Brasil + regiões
    geom_sf(data = regiao, fill = NA, color = "#000", size = 0.1) +
    
    # Define a paleta de cores azul
    scale_fill_manual(
      name = "Variable prevalence",
      values = c("#F2F2F2", "#C9DDEB", "#89B8D8", "#2166AC"),
      labels = c("Low", "Medium-Low", "Medium-High", "High"),
      na.value = "grey80"
    ) +
    labs(title = variavel_atual) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text = element_blank(), axis.title = element_blank(),
      axis.ticks = element_blank(), panel.grid = element_blank(),
      legend.position = "bottom",
      legend.box.background = element_rect(   # borda preta ao redor
       colour = "black",
       fill = NA,
       size = 0.3
      )
    )
  
  
  # # --- Parte 2: lógica do contorno vermelho --- #
  # Filtra dados_long2 para a variável da vez e depois para a categoria "Alta"
  regioes_alta <- dados_long2 %>%
    filter(variavel == variavel_atual,
           categoria == "High")

  # Verifica se existem regiões na categoria "Alta" para esta variável
  if (nrow(regioes_alta) > 0) {

    # Dissolve as geometrias para criar um contorno externo único
    contorno_externo <- regioes_alta %>%
      summarise(geometry = st_union(.))

    # Adiciona a camada do contorno vermelho ao mapa base que já criamos
    mapa_base_prev <- mapa_base_prev +
      geom_sf(data = contorno_externo,
              fill = NA,          # Sem preenchimento
              color = "red",      # Cor do contorno
              size = 0.7)         # Espessura do contorno
  }

  # Guarda o mapa final (azul + vermelho) na lista
  lista_mapas_combinados[[variavel_atual]] <- mapa_base_prev
}


# Trocando o título dos mapas, como feito nos mapas de importância local
for (i in seq_along(lista_mapas_combinados)) {
  lista_mapas_combinados[[i]] <- lista_mapas_combinados[[i]] +
    labs(title = novos_titulos_imp_eng[i])
}

# Extrai todas as caixas de legenda do primeiro mapa
legs_prevalencia <- cowplot::get_plot_component(
  lista_mapas_combinados[[3]],
  "guide-box",
  return_all = TRUE
)

# Seleciona **somente** a legenda
legenda_prevalencia <- legs_prevalencia[[3]]

# Remove as legendas individuais de todos os mapas
lista_sem_legenda <- lapply(lista_mapas_combinados, function(p) {
  p + theme(legend.position = "none")
})

# Combina os mapas (sem legenda) em uma grade
mapas_combinados_prev <- wrap_plots(lista_sem_legenda, ncol = 3)

# Monta a figura final: grade de mapas + legenda única abaixo
mapas_prev_final <- cowplot::plot_grid(
  mapas_combinados_prev,
  legenda_prevalencia,
  ncol = 1,
  rel_heights = c(1, 0.1) # Ajusta a altura relativa da legenda
)

# Exibe e salva
print(mapas_prev_final)
ggsave("Prevalence_Importance_Var.png", 
       plot = mapas_prev_final, 
       width = 20, height = 10, dpi = 300)


# Devido a impurezas no shapefile das RGIs brasileiras, foi necessário remover
# os resíduos dos contornos vermelhos manualmente, via Photoshop (v. 22.0)






# ----------------------------------------------
# SHAP DO RF DO BRASIL E CADA GRANDE REGIÃO
# ----------------------------------------------
# Calculando valores SHAP com base nos RFs já feitos para cada grande região
library(shapviz)
library(fastshap)
library(ggplot2)
library(patchwork)

# Wrapper de previsão
pred_wrapper <- function(object, newdata) {
  predict(object, newdata = newdata)
}

# SHAP BRASIL
set.seed(42)
shap_values_rf <- fastshap::explain(
  object       = rf_model,
  X            = dados_padronizado[, 5:27],  # background (distribuição de referência)
  pred_wrapper = pred_wrapper,               # função wrapper
  newdata      = dados_padronizado[, 5:27],   # instâncias a explicar
  nsim         = 1000                        # nº de permutações/amostragens
)

# Informar o "baseline" (valor esperado)
baseline_rf <- mean(pred_wrapper(rf_model, dados_padronizado)) # prob. média no background

# Construir objeto shapviz
sv_ind <- shapviz(
  object = shap_values_rf,
  X = dados_padronizado[, 5:27],
  y = dados_padronizado[, 4], 
  baseline = baseline_rf
)

# Gráfico
shap_ind <- sv_importance(sv_ind, kind = "beeswarm") + 
  ggtitle("Brazil") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) + 
  scale_y_discrete(labels = dicionario_vars_eng)

# Extraindo legenda
legenda_shap <- cowplot::get_legend(shap_ind)

# Removendo legenda e eixo X
shap_brasil <- shap_ind + 
  theme(legend.position = "none",
        axis.title.x = element_blank())



# SHAP NORTE
set.seed(42)
shap_values_rf_norte <- fastshap::explain(
  object       = rf_model_norte,
  X            = dados_norte[, 6:28],  
  pred_wrapper = pred_wrapper,              
  newdata      = dados_norte[, 6:28],  
  nsim         = 1000
)

# Informar o "baseline" (valor esperado)
baseline_rf_norte <- mean(pred_wrapper(rf_model_norte, dados_norte))  # prob. média no background

# Construir objeto 'sv' com SHAP
sv_ind_norte <- shapviz(
    object    = shap_values_rf_norte,
    X         = dados_norte[, 6:28], # apenas preditores
    y         = dados_norte[, 5], # apenas o desfecho
    baseline  = baseline_rf_norte
)

# Gráfico
shap_ind_norte <-  sv_importance(sv_ind_norte, kind = "beeswarm") + 
  ggtitle("North") +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)) + 
  scale_y_discrete(labels = dicionario_vars_eng)



# SHAP NORDESTE
set.seed(42)
shap_values_rf_nordeste <- fastshap::explain(
  object       = rf_model_nordeste,
  X            = dados_nordeste[, 6:28],  
  pred_wrapper = pred_wrapper,              
  newdata      = dados_nordeste[, 6:28],  
  nsim         = 1000
)

# Informar o "baseline" (valor esperado)
baseline_rf_nordeste <- mean(pred_wrapper(rf_model_nordeste, dados_nordeste))  # prob. média no background

# Construir objeto 'sv' com SHAP
sv_ind_nordeste <- shapviz(
  object    = shap_values_rf_nordeste,
  X         = dados_nordeste[, 6:28], # apenas preditores
  y         = dados_nordeste[, 5], # apenas o desfecho
  baseline  = baseline_rf_nordeste
)

# Gráfico
shap_ind_nordeste <-  sv_importance(sv_ind_nordeste, kind = "beeswarm") + 
  ggtitle("Northeast") +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)) + 
  scale_y_discrete(labels = dicionario_vars_eng)



# SHAP CENTRO-OESTE
set.seed(42)
shap_values_rf_centro_oeste <- fastshap::explain(
  object       = rf_model_centro_oeste,
  X            = dados_centro_oeste[, 6:28],  
  pred_wrapper = pred_wrapper,              
  newdata      = dados_centro_oeste[, 6:28],  
  nsim         = 1000
)

# Informar o "baseline" (valor esperado)
baseline_rf_centro_oeste <- mean(pred_wrapper(rf_model_centro_oeste, dados_centro_oeste))  # prob. média no background

# Construir objeto 'sv' com SHAP
sv_ind_centro_oeste <- shapviz(
  object    = shap_values_rf_centro_oeste,
  X         = dados_centro_oeste[, 6:28], # apenas preditores
  y         = dados_centro_oeste[, 5], # apenas o desfecho
  baseline  = baseline_rf_centro_oeste
)

# Gráfico
shap_ind_centro_oeste <-  sv_importance(sv_ind_centro_oeste, kind = "beeswarm") + 
  ggtitle("Central-West") +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)) + 
  scale_y_discrete(labels = dicionario_vars_eng)



# SHAP SUDESTE
set.seed(42)
shap_values_rf_sudeste <- fastshap::explain(
  object       = rf_model_sudeste,
  X            = dados_sudeste[, 6:28],  
  pred_wrapper = pred_wrapper,              
  newdata      = dados_sudeste[, 6:28],  
  nsim         = 1000
)

# Informar o "baseline" (valor esperado)
baseline_rf_sudeste <- mean(pred_wrapper(rf_model_sudeste, dados_sudeste))  # prob. média no background

# Construir objeto 'sv' com SHAP
sv_ind_sudeste <- shapviz(
  object    = shap_values_rf_sudeste,
  X         = dados_sudeste[, 6:28], # apenas preditores
  y         = dados_sudeste[, 5], # apenas o desfecho
  baseline  = baseline_rf_sudeste
)

# Gráfico
shap_ind_sudeste <-  sv_importance(sv_ind_sudeste, kind = "beeswarm") + 
  ggtitle("Southeast") +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)) + 
  scale_y_discrete(labels = dicionario_vars_eng)



# SHAP SUL
set.seed(42)
shap_values_rf_sul <- fastshap::explain(
  object       = rf_model_sul,
  X            = dados_sul[, 6:28],  
  pred_wrapper = pred_wrapper,              
  newdata      = dados_sul[, 6:28],  
  nsim         = 1000
)

# Informar o "baseline" (valor esperado)
baseline_rf_sul <- mean(pred_wrapper(rf_model_sul, dados_sul))  # prob. média no background

# Construir objeto 'sv' com SHAP
sv_ind_sul <- shapviz(
  object    = shap_values_rf_sul,
  X         = dados_sul[, 6:28], # apenas preditores
  y         = dados_sul[, 5], # apenas o desfecho
  baseline  = baseline_rf_sul
)

# Gráfico
shap_ind_sul <-  sv_importance(sv_ind_sul, kind = "beeswarm") + 
  ggtitle("South") +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)) + 
  scale_y_discrete(labels = dicionario_vars_eng)



# Montando grid de SHAP individuais
grid_shap_ind <- (shap_brasil + shap_ind_norte) /
  (shap_ind_nordeste + shap_ind_centro_oeste) / 
  (shap_ind_sudeste + shap_ind_sul)


# Combinando o grid de gráficos à legenda extraída
shap_com_legenda <- grid_shap_ind | legenda_shap

# Adiciona o título do eixo X como uma anotação de rodapé (caption)
# e ajusta o espaço da legenda com plot_layout
grid_shap_ind_final <- shap_com_legenda +
  plot_annotation(
    caption = 'SHAP Value',
    theme = theme(plot.caption = element_text(hjust = 0.5, size = 12)) # Centraliza e formata
  ) +
  plot_layout(
    widths = c(1, 0.1) # Ajusta a largura relativa do grid de gráficos vs. a legenda
  )

# Exibe e salva
print(grid_shap_ind_final)
ggsave(
  "SHAP_Brazil_Major_Regions.PNG",
  plot = grid_shap_ind_final,
  width = 15, height = 20, dpi = 300
)







# ---------------------------------------------------------------------------
# SHAP GRID TOP 3 DEPENDENCE PARA BRASIL E GRANDES REGIÕES
# ---------------------------------------------------------------------------
# Os gráficos da coluna do meio recebem como titulo os nomes das regiões da 
# linha seguinte, para melhor interpretacao do grafico

# Brasil
dep_bra_pib <- sv_dependence(sv_ind_original, "mp_pib_per_capita", 
                             color_var = "mp_pib_per_capita") +
  labs(x = "Per capita GDP") +
  theme(legend.position = 'none')


dep_bra_esgoto <- sv_dependence(sv_ind_original, "mp_esgoto_ou_fossa",
                                color_var = "mp_esgoto_ou_fossa")  +
  ggtitle("Brazil") +
  labs(x = "Households with sewer or septic tank") +
  theme(legend.position = 'none',
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15))


dep_bra_mulher <- sv_dependence(sv_ind_original, "mp_domicilio_mulher_responsavel",
                                color_var = "mp_domicilio_mulher_responsavel")  +
  labs(x = "Female-headed households") +
  theme(legend.position = 'none')



# Norte
dep_norte_enf <- sv_dependence(sv_ind_original_norte, "mp_enfermeiros_2010", 
                               color_var = "mp_enfermeiros_2010") +
  labs(x = "Nurses") +
  theme(legend.position = 'none')


dep_norte_dens_demo <- sv_dependence(sv_ind_original_norte, "mp_densidade_demografica", 
                                     color_var = "mp_densidade_demografica") +
  ggtitle("North") +
  labs(x = "Population density") +
  theme(legend.position = 'none',
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15))


dep_norte_pib <- sv_dependence(sv_ind_original_norte, "mp_pib_per_capita", 
                               color_var = "mp_pib_per_capita") +
  labs(x = "Per capita GDP") +
  theme(legend.position = 'none')




# Nordeste
dep_nordeste_medio <- sv_dependence(sv_ind_original_nordeste, "mp_medio_completo", 
                                    color_var = "mp_medio_completo") +
  labs(x = "Completed secondary education") +
  theme(legend.position = 'none')


dep_nordeste_pib <- sv_dependence(sv_ind_original_nordeste, "mp_pib_per_capita", 
                                  color_var = "mp_pib_per_capita") +
  ggtitle("Northeast") + 
  labs(x = "Per capita GDP") +
  theme(legend.position = 'none',
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15))


dep_nordeste_dens_demo <- sv_dependence(sv_ind_original_nordeste, "mp_densidade_demografica", 
                                        color_var = "mp_densidade_demografica") +
  labs(x = "Population density") +
  theme(legend.position = 'none')




# Centro-Oeste
dep_centro_oeste_medio <- sv_dependence(sv_ind_original_centro_oeste, "mp_medio_completo", 
                                        color_var = "mp_medio_completo") +
  labs(x = "Completed secondary education") +
  theme(legend.position = 'none')


dep_centro_oeste_dens_demo <- sv_dependence(sv_ind_original_centro_oeste, "mp_densidade_demografica", 
                                            color_var = "mp_densidade_demografica") +
  ggtitle("Central-West") +
  labs(x = "Population density") +
  theme(legend.position = 'none',
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15))


dep_centro_oeste_esf <- sv_dependence(sv_ind_original_centro_oeste, "mp_cobertura_esf_2010", 
                                      color_var = "mp_cobertura_esf_2010") +
  labs(x = "FHS coverage") +
  theme(legend.position = 'none')




# Sudeste
dep_sudeste_medio <- sv_dependence(sv_ind_original_sudeste, "mp_medio_completo", 
                                   color_var = "mp_medio_completo") +
  labs(x = "Completed secondary education") +
  theme(legend.position = 'none')

dep_sudeste_mulher <- sv_dependence(sv_ind_original_sudeste, "mp_domicilio_mulher_responsavel",
                                    color_var = "mp_domicilio_mulher_responsavel")  +
  ggtitle("Southeast") +
  labs(x = "Female-headed households") +
  theme(legend.position = 'none',
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15))

dep_sudeste_medicos <- sv_dependence(sv_ind_original_sudeste, "mp_medicos_gyn_ob_2010", 
                                     color_var = "mp_medicos_gyn_ob_2010") +
  labs(x = "Specialist physicians") +
  theme(legend.position = 'none')



# Sul
dep_sul_medio <- sv_dependence(sv_ind_original_sul, "mp_medio_completo", 
                               color_var = "mp_medio_completo") +
  labs(x = "Completed secondary education") +
  theme(legend.position = 'none')

dep_sul_mulher <- sv_dependence(sv_ind_original_sul, "mp_domicilio_mulher_responsavel",
                                color_var = "mp_domicilio_mulher_responsavel")  +
  ggtitle("South") +
  labs(x = "Female-headed households") +
  theme(legend.position = 'none',
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15))

dep_sul_dens_demo <- sv_dependence(sv_ind_original_sul, "mp_densidade_demografica", 
                                   color_var = "mp_densidade_demografica") +
  labs(x = "Population density") +
  theme(legend.position = 'none')



# Montando grid de SHAPs
grid_dep <- (dep_bra_pib + dep_bra_esgoto + dep_bra_mulher) / # Brasil
  (dep_norte_enf + dep_norte_dens_demo + dep_norte_pib) / # Norte
  (dep_nordeste_medio + dep_nordeste_pib + dep_nordeste_dens_demo) / # Nordeste
  (dep_centro_oeste_medio + dep_centro_oeste_dens_demo + dep_centro_oeste_esf) / # Centro-Oeste
  (dep_sudeste_medio + dep_sudeste_mulher + dep_sudeste_medicos) / # Sudeste
  (dep_sul_medio + dep_sul_mulher + dep_sul_dens_demo) # Sul


# Combinando o grid de gráficos à legenda extraída
grid_dep_leg <- grid_dep | legenda_shap

# Adiciona o título do eixo X como uma anotação de rodapé (caption)
# e ajusta o espaço da legenda com plot_layout
grid_dep_final <- grid_dep_leg +
  plot_annotation(
    theme = theme(plot.caption = element_text(hjust = 0.5, size = 12)) # Centraliza e formata
  ) +
  plot_layout(
    widths = c(1, 0.1) # Ajusta a largura relativa do grid de gráficos vs. a legenda
  )

ggsave(
  "SHAP_dep_grid.PNG",
  plot = grid_dep_final,
  width = 10, height = 18, dpi = 300
)






# ---------------------------------------------------------------
# MEDIDAS RESUMO DAS VARIÁVEIS PELAS GRANDES REGIÕES BRASILEIRAS
# ---------------------------------------------------------------
library(gtsummary)
library(flextable)
library(dplyr)

# Agregando regiao ao banco, para divisao da tabela
dados_com_regiao <- regioes_imediatas[,c(1,7)] %>%
  left_join(dados, by = c("code_immediate" = "regiao_geografica_imediata"))
dados_com_regiao <- dados_com_regiao[!is.na(dados_com_regiao$TAXA_AJUSTADA), ]
dados_com_regiao <- as.data.frame(dados_com_regiao)


# Reordenando as colunas para facilitar seleção no tbl_summary
dados_com_regiao <- dados_com_regiao %>% 
  relocate(name_region, .before = TAXA_AJUSTADA)

# Configurando a tabela
theme_gtsummary_language(
  language = "en",     # Define o idioma para Inglês
  decimal.mark = ",",    # Define a vírgula como separador decimal
  big.mark = ".",        # Define o ponto como separador de milhares
  iqr.sep = "-",         # Define o hífen como separador para intervalos interquartis
  ci.sep = "-",          # Define o hífen como separador para intervalos de confiança
  set_theme = TRUE       # Aplica essas configurações como tema padrão para as tabelas
)

# Aplica a função personalizada de formatação de porcentagens como tema padrão
list("tbl_summary-fn:percent_fun" = function(x) sprintf(x * 100, fmt='%#.2f')) %>% 
  set_gtsummary_theme()

# Criando tabela por região
tbl_summary(
  dados_com_regiao[, 5:29], 
  by = name_region,
  
  statistic = list(all_continuous() ~ "{mean} ± {sd} | {median} ({IQR})"),
  
  label = list(
    TAXA_AJUSTADA ~ "Age-adjusted female breast cancer mortality rate",
    mp_cobertura_esf_2010 ~ "FHS coverage",
    mp_densidade_demografica	~ "Population density",
    mp_mais_de_0.5_a_1_morador	~ "Low household crowding",
    mp_desemprego	~ "Unemployment",
    mp_domicilio_entorno_arborizado ~ "Tree-covered surroundings", 
    mp_domicilio_mulher_responsavel	~ "Female-headed households",
    mp_domicilio_entorno_pavimentado	~ "Paved surroundings",
    mp_enfermeiros_2010	~ "Nurses",
    mp_esgoto_ou_fossa	~ "Households with sewer or septic tank",
    mp_evangelicos	~ "Evangelical population",
    mp_gini	~ "Gini Index",
    mp_homicidios_2010	~ "Homicide rate",
    mp_leitos_hosp_2010	~ "Hospital beds",
    mp_mamografos_2010	~ "Mammography units",
    mp_medicos_gyn_ob_2010	~ "Specialist physicians",
    mp_medio_completo	~ "Completed secondary education",
    mp_PCDs	~ "Disability",
    mp_pib_per_capita	~ "Per capita GDP",
    mp_raca_preta	~ "Black population",
    mp_recursos_destinados_saude_2010 ~ "Own-source health expenditure",
    mp_suicidios_2010	~ "Suicide rate",
    mp_trabalha_mais_44h	~ "Excessive working hours",
    mp_ultrassons_2010 ~ "Ultrasound units"
  )
) %>% 
  add_overall(last = T) %>%
  bold_labels() %>% 
  modify_header(label ~ "**Variable**") %>% 
  modify_caption("eTable 1. Summmary statistics of study variables by brazilian major geographic region.") %>% 
  
  # Exportando para .docx
  # (em um tamanho bizarramente pequeno, mas deu pra ajeitar na mão até o momento)
  as_flex_table() %>% 
  fontsize(size = 8, part = "all") %>%
  padding(padding = 2, part = "all") %>%  # Reduzir padding
  set_table_properties(layout = "autofit", width = 1) %>%  # Layout automático
  fit_to_width(max_width = 10) %>% 
  save_as_docx(
    path = "Summary_statistics_study_variables.docx",
    pr_section = officer::prop_section(
      page_size = officer::page_size(orient = "landscape", width = 11, height = 8.5),
      type = "continuous",
      page_margins = officer::page_mar(bottom = 0.5, top = 0.5, right = 0.5, left = 0.5)
    )
  )





#--------------------------------------------
# GRÁFICO DO PSEUDO R² DO GRF
#--------------------------------------------

# Pegando os pseudo R² dos submodelos locais de GRF
dados_R2 <- grf_model$LGofFit

# Os bancos possuem a mesma ordem de obs
dados_R2_mapa <- bind_cols(dados_R2, dados_mapa)

# Pegando valores maximo e minimo de pseudo R² para
# determinacao dos breaks

min_R2 <- min(dados_R2_mapa$LM_Rsq100)
max_R2 <- max(dados_R2_mapa$LM_Rsq100)
mid_R2 <- round(min_R2 + max_R2)

mapa_R2 <- ggplot(data = dados_R2_mapa) +
  geom_sf(aes(fill = LM_Rsq100, geometry = geom), color = "gray50", lwd = 0.1) +
  geom_sf(data = regiao, fill = NA, color = "#000", size = 0.5) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    breaks = c(min_R2, mid_R2, max_R2),
    limits = c(min_R2, max_R2),
    
    # limitando a 2 casas decimais
    labels = function(x) sprintf("%.3f", x)
  ) +
  labs(fill = "Pseudo R squared (%)") +
  theme_void()
  
# Salvando o grafico
ggsave(
  "pseudo_r2.png",
  plot = mapa_R2,
  width = 6, height = 6, dpi = 300
)






# Preparação dos Dados para Modelagem Preditiva e Espacial
# Com a base analítica consolidada, completa (sem valores ausentes) e depurada 
# quanto à multicolinearidade, procedeu-se à preparação estatística para as etapas
# de modelagem preditiva e espacial. Esta fase contemplou: (i) divisão amostral 
# entre conjuntos de treino e teste; (ii) padronização das variáveis numéricas 
# explicativas; (iii) estruturação dos dados em objetos espaciais; e (iv) formulação 
# algébrica do modelo estatístico.
# 
# 1. Divisão dos Dados em Conjuntos de Treino e Teste
# Para garantir a validade externa dos modelos e permitir avaliação de desempenho
# em dados não utilizados no ajuste, a base foi aleatoriamente dividida em dois
# subconjuntos: 80% das observações foram destinadas ao treinamento e 20% ao teste,
# respeitando-se uma semente fixa de aleatoriedade (set.seed = 123) 
# para assegurar reprodutibilidade. A divisão foi realizada sem estratificação, 
# uma vez que a variável dependente (taxa de mortalidade ajustada) é contínua e 
# geograficamente distribuída.
# 
# 2. Padronização das Variáveis Explicativas (Z-score)
# As variáveis explicativas numéricas foram submetidas à padronização por escore-z,
# transformando seus valores de forma que cada variável tenha média zero e desvio 
# padrão igual a um. Esta transformação é essencial para:
#   
# evitar que variáveis com escalas maiores dominem a modelagem;
# garantir a comparabilidade entre coeficientes estimados;
# permitir a convergência estável de algoritmos sensíveis à escala.
# 
# Importante destacar que os parâmetros de padronização (média e desvio padrão) 
# foram calculados exclusivamente com base no conjunto de treino, conforme 
# recomendações de validação cruzada. Esses mesmos parâmetros foram utilizados 
# posteriormente para padronizar o conjunto de teste, assegurando que os dados 
# de validação não influenciassem o treinamento do modelo.
# 
# 3. Estruturação dos Dados Espaciais
# Para aplicação de métodos de regressão espacial, o conjunto de treino foi convertido
# para o formato SpatialPointsDataFrame, com sistema de projeção geográfica WGS84 
# (EPSG:4326). Esse formato é compatível com pacotes especializados em análise espacial
# como spgwr, GWmodel e SpatialML, permitindo a modelagem com base na vizinhança
# geográfica e no cálculo de kernels espaciais adaptativos ou fixos.
# 
# Adicionalmente, foi construída uma matriz com as coordenadas (longitude e latitude) 
# das Regiões Geográficas Imediatas, a ser utilizada como argumento em métodos que 
# exigem matriz de localização explícita (e.g., GWR clássica ou Random Forest espacial).
# 
# 4. Formulação do Modelo Preditivo
# A fórmula estatística que representa a estrutura funcional dos modelos foi gerada de 
# forma programática. Definiu-se como variável dependente a TAXA_AJUSTADA (taxa de 
# mortalidade por câncer de mama ajustada por idade). As variáveis explicativas foram
# selecionadas automaticamente, excluindo-se identificadores espaciais (nome da região,
# coordenadas geográficas) e a própria variável dependente.
# 
# Modelagem Global com Random Forest: Otimização de Hiperparâmetros
# A etapa de modelagem preditiva foi conduzida por meio do algoritmo Random Forest 
# (RF), uma técnica de aprendizado de máquina baseada em ensemble de árvores de 
# decisão (Breiman, 2001). O RF foi selecionado por sua robustez frente à multicolinearidade,
# capacidade de modelar relações não lineares e elevada acurácia preditiva mesmo em 
# conjuntos de dados com estrutura complexa.
# 
# 1. Formulação do Modelo e Conjunto de Treinamento
# Utilizou-se como variável dependente a taxa de mortalidade ajustada por idade 
# por câncer de mama nas Regiões Geográficas Imediatas. As variáveis explicativas
# foram definidas a partir da base padronizada e pré-processada (dados_treino_padronizado),
# excluindo-se identificadores espaciais e coordenadas geográficas. 
# 
# 2. Otimização de Hiperparâmetros via Validação Cruzada
# Foram avaliadas diversas combinações dos dois principais hiperparâmetros do modelo RF:
#   
# ntree: número de árvores da floresta, testado nos valores de 100, 500, 1000, 1500 e 2000;
# mtry: número de variáveis sorteadas aleatoriamente para divisão em cada nó da árvore,
# variando de 1 até o número total de preditores, com incrementos de 2 unidades.
# 
# O processo de otimização foi conduzido com validação cruzada k-fold (k = 5), 
# utilizando-se o pacote caret. Em cada repetição, o modelo foi ajustado a 80% dos 
# dados (treinamento) e validado em 20%, com métricas de desempenho 
# (RMSE – Root Mean Square Error) calculadas e registradas para cada combinação testada.
# 
# 3. Paralelização Computacional
# A otimização dos hiperparâmetros foi realizada de forma paralelizada por meio do
# pacote doParallel, utilizando até 75% dos núcleos físicos disponíveis no sistema,
# com o objetivo de reduzir o tempo computacional e acelerar as simulações. 
# Para garantir reprodutibilidade dos resultados, uma semente fixa foi utilizada (set.seed(123)).
# 
# 4. Critério de Seleção e Avaliação do Modelo
# Ao final do processo, foi selecionado o modelo com o menor valor de RMSE médio
# entre as dobras da validação cruzada. Os hiperparâmetros ótimos (mtry*, ntree*)
# foram identificados, e a performance do modelo foi visualmente inspecionada por 
# meio de gráficos de desempenho.
# 
# Referências Sugeridas
# Breiman, L. (2001). Random forests. Machine Learning, 45(1), 5–32.
# 
# Kuhn, M., & Johnson, K. (2013). Applied Predictive Modeling. Springer.
# 
# James, G., Witten, D., Hastie, T., & Tibshirani, R. (2021). An Introduction to
# Statistical Learning (2nd ed.). Springer.
# 
# Liaw, A., & Wiener, M. (2002). Classification and regression by randomForest.
# R News, 2(3), 18–22.
# 
# Modelagem Espacial por Regressão Geograficamente Ponderada (GWR)
# Com o objetivo de explorar variações espaciais na associação entre determinantes 
# contextuais e a taxa de mortalidade ajustada por idade por câncer de mama, foi
# aplicada a técnica de Regressão Geograficamente Ponderada (GWR). A GWR é um 
# modelo estatístico espacial que permite estimar coeficientes de regressão localmente,
# ou seja, em cada unidade espacial, considerando a vizinhança geográfica como fator 
# de ponderação (Fotheringham, Brunsdon & Charlton, 2002).
# 
# Justificativa para o Uso da GWR
# A abordagem GWR foi escolhida por sua capacidade de capturar heterogeneidade 
# espacial não modelada por regressões tradicionais de efeitos fixos.
# Isso é especialmente relevante em estudos de saúde coletiva, nos quais fatores
# sociais, econômicos e demográficos podem exercer efeitos diferenciados sobre 
# desfechos de interesse conforme a região geográfica analisada.
# 
# Otimização da Largura de Banda (Bandwidth)
# A largura de banda (bandwidth) — parâmetro que define a extensão geográfica da
# vizinhança considerada para o ajuste local — foi otimizada com base no Critério 
# de Informação de Akaike corrigido (AICc). Esse procedimento foi realizado 
# exclusivamente sobre o conjunto de dados de treino, a fim de evitar vazamento 
# de informação (data leakage) no processo de seleção.
# 
# Foi adotada uma largura de banda adaptativa, o que implica que a vizinhança de 
# cada unidade é composta por um número fixo de vizinhos mais próximos, em vez de
# uma distância fixa. Tal abordagem é particularmente recomendada quando há 
# desigualdade na distribuição espacial das unidades de análise.
# 
# O kernel de ponderação utilizado foi do tipo bisquare, comumente empregado na 
# literatura por atribuir pesos decrescentes com o aumento da distância e reduzir 
# a influência de observações distantes.
# 
# Ajuste e Interpretação do Modelo
# O modelo foi ajustado utilizando a função gwr.basic() do pacote GWmodel, com os 
# parâmetros otimizados previamente. O resultado do modelo compreende:
# 
# coeficientes de regressão locais para cada unidade espacial (Região Geográfica 
# Imediata),
# valores locais de R², que expressam a qualidade do ajuste em nível local,
# estatísticas diagnósticas globais, como AICc, log-verossimilhança e resíduos.
# 
# Os coeficientes locais obtidos foram armazenados como um objeto espacial e posteriormente
# convertidos para sf (Simple Features), permitindo análises cartográficas e 
# identificação de padrões espaciais heterogêneos.
# 
# Referências Sugeridas
# Fotheringham, A. S., Brunsdon, C., & Charlton, M. (2002). Geographically 
# Weighted Regression: The Analysis of Spatially Varying Relationships. Wiley.
# 
# Wheeler, D. C., & Páez, A. (2010). Geographically weighted regression. In: 
#   Fischer, M., Getis, A. (Eds.), Handbook of Applied Spatial Analysis. Springer.
# 
# Charlton, M., & Fotheringham, A. S. (2009). GWmodel: Geographically Weighted Models.
# R package documentation.
# 
# Modelagem Espacial com Geographically Weighted Random Forest (GRF)
# Com o objetivo de capturar heterogeneidade espacial não linear na relação entre
# os determinantes contextuais e a mortalidade ajustada por câncer de mama, foi 
# aplicado o algoritmo Geographically Weighted Random Forest (GRF). Essa técnica 
# combina os fundamentos da Random Forest (Breiman, 2001) com a lógica de ponderação
# espacial adaptativa, permitindo que múltiplas florestas sejam ajustadas localmente
# em diferentes regiões geográficas (Georganos et al., 2019).
# 
# Fundamentação da Técnica
# A GRF é uma extensão não paramétrica e espacialmente explícita da Random Forest,
# na qual a predição em cada unidade geográfica é feita com base em uma floresta 
# ajustada somente com observações vizinhas, ponderadas por sua proximidade. 
# A estrutura é análoga à Regressão Geograficamente Ponderada (GWR), porém mais 
# robusta quanto à não linearidade, à presença de colinearidade entre preditores
# e à variabilidade complexa dos efeitos locais.
# 
# Parâmetros de Modelagem
# O modelo GRF foi ajustado utilizando os hiperparâmetros ntree (número de árvores)
# e mtry (número de variáveis sorteadas em cada divisão), previamente otimizados com 
# base em Random Forest global (vide seção anterior). O kernel de ponderação espacial 
# adotado foi do tipo adaptativo, o que implica a definição da largura de banda como o 
# número de vizinhos mais próximos considerados em cada ponto, ao invés de uma distância fixa.
# 
# Otimização da Largura de Banda (Bandwidth)
# A largura de banda — parâmetro essencial para o controle da suavização espacial —
# foi otimizada por meio de validação cruzada estratificada (k-fold, k = 5), com
# base no erro quadrático médio (RMSE). Foram testados valores de bandwidth entre 
# 50 e 300 vizinhos, respeitando o número mínimo necessário de observações em relação
# ao número de preditores, evitando sobreajuste local.
# 
# Em cada iteração de validação, os dados foram divididos em cinco subconjuntos
# aleatórios. O modelo foi treinado em quatro deles e testado no quinto, repetindo 
# o processo para todos os folds e para cada valor de bandwidth. O RMSE médio entre
# os folds foi calculado e armazenado para cada configuração.
# 
# Seleção do Modelo Final
# O valor de bandwidth que apresentou o menor RMSE médio foi selecionado como ótimo.
# Em seguida, o modelo GRF foi ajustado utilizando a base de treinamento completa, 
# com os parâmetros otimizados (ntree, mtry, bw). O modelo final passou a incorporar
# variações locais tanto na estrutura de predição quanto na importância relativa das
# variáveis explicativas.
# 
# Referências Sugeridas
# Breiman, L. (2001). Random forests. Machine Learning, 45(1), 5–32.
# 
# Georganos, S., Grippa, T., Lennert, M., & Vanhuysse, S. (2019). Geographically 
# weighted random forests: a spatial extension of the random forest algorithm to 
# address spatial heterogeneity in remote sensing and population modelling. 
# Geocarto International, 36(2), 121–136.
# 
# Wang, Y., Li, T., Li, Y., & Zhang, C. (2018). Spatial random forest for estimating
# ground-level PM2.5 in the United States. Environmental Pollution, 240, 996–1006.