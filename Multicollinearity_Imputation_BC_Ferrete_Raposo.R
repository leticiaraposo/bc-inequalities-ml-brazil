# =============================================================================
# Título: Spatially Varying Socioeconomic Determinants of Breast Cancer Mortality in Brazil
# Autores: Letícia Raposo, Lívia Ferrete
# Etapa 3: verificação de multicolinearidade, remoção das variáveis com mais de
# 1% de missing e imputação de dados ausentes
# Data: Julho de 2025
# =============================================================================

# --- CARREGAMENTO DE PACOTES ---
library(sp)             # Objetos espaciais (legacy)
library(sf)             # Simple Features - modelagem espacial moderna
library(spgwr)          # GWR tradicional (obsoleto, mas ainda funcional)
library(GWmodel)        # Modelos geograficamente ponderados (recomendado)
library(SpatialML)      # GWR + Random Forest espacial
library(randomForest)   # Random Forest clássico
library(caret)          # Treinamento, validação e tuning de modelos
library(tmap)           # Visualização geográfica (mapas)
library(pdp)            # Partial dependence plots
library(ggplot2)        # Visualização geral
library(dplyr)          # Manipulação de dados
library(readxl)         # Leitura de arquivos Excel
library(mice)           # Imputação multivariada por chained equations



# --- IMPORTAÇÃO DA BASE DE DADOS PROCESSADA ---
dados <- readRDS("base_regioes_imediatas_cancer_mama_periodo_2015_2019.RDS")
dados$mp_taxa_bruta_mortalidade <- NULL # Variável errada



# --- IMPUTAÇÃO DE DADOS AUSENTES COM MICE ---

# Seleção das variáveis explicativas (variáveis 4 a 62)
variaveis_expl <- colnames(dados)[4:62]

# Variáveis utilizadas no modelo, incluindo a resposta (TAXA_AJUSTADA)
colunas_para_modelo <- c("TAXA_AJUSTADA", variaveis_expl)

# Subconjunto da base contendo apenas as variáveis relevantes para a imputação
dados_para_imputar <- dados %>%
  dplyr::select(regiao_geografica_imediata, Y, X, all_of(colunas_para_modelo))

# Diagnóstico inicial de valores ausentes
message("Status dos NAs antes da imputação:")
print(colSums(is.na(dados_para_imputar)))
# Só a variável mp_recursos_destinados_saude_2010 tem 1 NA

# Aplicação do algoritmo MICE com método Predictive Mean Matching (pmm)
imputacoes_mice <- mice(
  dados_para_imputar,
  m = 5,              # número de conjuntos imputados
  method = 'pmm',     # método de imputação robusto para dados contínuos
  maxit = 5,          # número de iterações
  seed = 123          # reprodutibilidade
)

# Extração do primeiro conjunto imputado (poderia-se também combinar múltiplos imputados)
dados_imputado <- complete(imputacoes_mice, 1)

# Diagnóstico final após imputação
message("Status dos NAs após imputação MICE:")
print(colSums(is.na(dados_imputado)))

# --- SALVAMENTO DA BASE FINAL IMPUTADA ---
# saveRDS(dados_imputado, "base_cancer_mama_imputados_ri_09_07.rds")



# =============================================================================
# ANÁLISE EXPLORATÓRIA E CONTROLE DE MULTICOLINEARIDADE
# Correlação com a variável resposta, análise de Pearson/Spearman, VIF
# =============================================================================

# --- DEFINIÇÃO DA VARIÁVEL RESPOSTA ---
variavel_resposta <- "TAXA_AJUSTADA"

# --- IDENTIFICAÇÃO DAS VARIÁVEIS EXPLICATIVAS ---
variaveis_explicativas <- setdiff(colnames(dados_imputado), 
                                  c("regiao_geografica_imediata", "Y", "X", variavel_resposta))

# --- CORRELAÇÃO COM A VARIÁVEL RESPOSTA ---
correlacoes <- sapply(dados_imputado[variaveis_explicativas], 
                      function(x) cor(x, dados_imputado[[variavel_resposta]], use = "pairwise.complete.obs", method = "spearman"))
correlacoes_ordenadas <- sort(correlacoes, decreasing = TRUE)

message("Correlação das variáveis explicativas com 'TAXA_AJUSTADA':")
print(correlacoes_ordenadas)

# --- MATRIZ DE CORRELAÇÃO SPEARMAN E REMOÇÃO DE MULTICOLINEARIDADE ---
matriz_cor <- cor(dados_imputado[, variaveis_explicativas], 
                  method = "spearman", use = "pairwise.complete.obs")

# Identificação de variáveis altamente correlacionadas (cutoff = 0.8)
vars_a_remover_cor <- findCorrelation(matriz_cor, cutoff = 0.8, names = TRUE) # 26 variáveis

# [1] "mp_extremamente_pobres"            "mp_pobres"                         "mp_domicilios_com_geladeira"      
# [4] "mp_idhm"                           "mp_domicilio_desktop_internet"     "mp_automovel"                     
# [7] "mp_renda_per_capita"               "mp_razao_dependencia"              "mp_expectativa_de_vida"           
# [10] "mp_coleta_lixo"                    "mp_eletricidade"                   "mp_mortalidade_infantil"          
# [13] "mp_analfabetismo"                  "mp_raca_branca"                    "mp_medicos_2010"                  
# [16] "mp_pop_planos_seguros_2010"        "mp_superior_completo"              "mp_mais_de_1_a_2_moradores"       
# [19] "mp_pop_agua_encanada"              "mp_criancas_trabalhando"           "mp_urbana"                        
# [22] "mp_idosos"                         "mp_envelhecimento"                 "mp_homens"                        
# [25] "mp_cobertura_ab_2010"              "mp_trinta_e_um_a_sessenta_minutos"

# Remoção do conjunto original imputado
dados2 <- dados_imputado[, !(colnames(dados_imputado) %in% vars_a_remover_cor)]

# --- AVALIAÇÃO DE MULTICOLINEARIDADE VIA VARIANCE INFLATION FACTOR (VIF) ---

# Definição da variável dependente
variavel_dependente <- "TAXA_AJUSTADA"

# Seleção das variáveis numéricas restantes para cálculo de VIF
vars_com_mortalidade <- dados2[, 4:37] %>%
  mutate(across(.cols = where(is.numeric) & !all_of(variavel_dependente),
                .fns = ~as.numeric(scale(.))))  # Padronização (Z-score)

# Modelo linear para avaliação do VIF
library(car)
modelo_vif <- lm(TAXA_AJUSTADA ~ ., data = vars_com_mortalidade)
resultado_vif <- vif(modelo_vif)

# Identificação de variáveis com VIF > 5
vars_alto_vif <- names(resultado_vif[resultado_vif > 5]) # 10 variáveis

# [1] "mp_aposentados_pensionistas"  "mp_casados"                   "mp_ate_0.5_morador"           "mp_ate_5min"                 
# [5] "mp_seis_a_trinta_minutos"     "mp_mais_de_1h_ate_2h"         "mp_taxa_fecundidade"          "mp_nasc_vivos_pre_natal_2010"
# [9] "mp_raca_parda"                "mp_mulheres"

# Remoção das variáveis com VIF elevado do conjunto final
dados2 <- dados2[, !(colnames(dados2) %in% vars_alto_vif)]

# --- LIMPEZA DE OBJETOS INTERMEDIÁRIOS ---
rm(list = c("matriz_cor", "modelo_vif", "vars_com_mortalidade", "resultado_vif",
            "vars_a_remover_cor", "vars_alto_vif", "correlacoes", "correlacoes_ordenadas"))

# --- RESULTADO FINAL ---
message("Base final pronta para modelagem (dados2):")
print(dim(dados2))

# Salvando a base limpa após todos os pré-processamentos
# saveRDS(dados2, "base_cancer_mama_imputados_sem_correlacao_alta_ri.rds")

# Tratamento dos Dados: Imputação de Valores Ausentes, Correlação e Multicolinearidade
# Após a construção da base analítica contendo indicadores socioeconômicos, demográficos
# e de saúde agregados por Regiões Geográficas Imediatas (RGIs), deu-se início ao tratamento
# estatístico dos dados com vistas à preparação para modelagem espacial e inferência.
# 
# 1. Imputação de Dados Ausentes com MICE
# A base analítica original continha valores ausentes (missing values) em apenas 1
# variável explicativa (mp_recursos_destinados_saude_2010). Para lidar com esse problema, 
# optou-se pelo uso do algoritmo Multiple Imputation by Chained Equations (MICE), amplamente recomendado 
# para imputações em contextos multivariados (Rubin, 1987; van Buuren & Groothuis-Oudshoorn, 2011).
# 
# O método selecionado para imputação foi o Predictive Mean Matching (PMM), que preserva
# a distribuição empírica das variáveis contínuas e reduz o risco de imputações implausíveis.
# Foram geradas cinco iterações (m = 5), com número máximo de 5 ciclos de imputação por 
# variável (maxit = 5). O primeiro conjunto imputado foi utilizado como base completa para
# as análises subsequentes, dado que não havia divergência substancial entre os conjuntos imputados.
# 
# 2. Análise de Correlação com a Variável Resposta
# Concluído o processo de imputação, foi realizada uma análise exploratória para investigar
# a correlação entre cada variável explicativa e a variável dependente — TAXA_AJUSTADA,
# que representa a taxa de mortalidade por câncer de mama ajustada por idade (por 100 mil mulheres).
# Foram calculados os coeficientes de correlação de Spearman, e as variáveis foram ordenadas 
# conforme a magnitude da correlação, permitindo a identificação das variáveis potencialmente 
# mais informativas.
# 
# 3. Detecção e Tratamento de Multicolinearidade
# Dado o elevado número de variáveis explicativas e o potencial de redundância informacional,
# foi conduzido um procedimento sistemático para detecção e remoção de multicolinearidade. 
# Este processo ocorreu em duas etapas complementares:
#   
# 3.1. Correlação entre Variáveis Explicativas (Spearman):
# Foi construída uma matriz de correlação de Spearman entre todas as variáveis explanatórias.
# Utilizou-se o critério de exclusão automática de variáveis altamente correlacionadas (|ρ| > 0,80),
# por meio da função findCorrelation() do pacote caret. Esta função remove iterativamente a
# variável com a maior média de correlação absoluta, resultando em um subconjunto menos redundante.
# Um total de 26 variáveis foi eliminado com base nesse critério.
# 
# 3.2. Variance Inflation Factor (VIF):
# Para avaliar a colinearidade multivariada remanescente, foi ajustado um modelo de regressão 
# linear múltipla com a variável dependente TAXA_AJUSTADA. As variáveis explicativas foram
# previamente padronizadas (Z-score) a fim de evitar distorções por diferenças de escala. 
# Calculou-se o Fator de Inflação da Variância (VIF) para cada variável preditora. Aquelas com VIF 
# superior a 5 — valor frequentemente utilizado como limiar de preocupação (Kutner et al., 2005) — 
# foram removidas. Esta etapa resultou na exclusão de mais 10 variáveis.
# 
# 4. Conjunto de Dados Final
# Como resultado do processo de imputação e controle de multicolinearidade, obteve-se um conjunto 
# de dados final (dados2), completo (sem valores ausentes), com variáveis informativas e 
# estatisticamente adequadas à modelagem. Este conjunto será utilizado nas etapas posteriores 
# de regressão geograficamente ponderada (GWR), Random Forest e outras abordagens de análise 
# espacial e preditiva.
# 
# Referências Sugeridas (se desejar incluir em artigo):
# van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate imputation by chained
# equations in R. Journal of Statistical Software, 45(3), 1–67.
# 
# Rubin, D. B. (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.
# 
# Kutner, M. H., Nachtsheim, C. J., Neter, J., & Li, W. (2005). Applied Linear Statistical 
# Models (5th ed.). McGraw-Hill.
# 
# Dormann, C. F., et al. (2013). Collinearity: a review of methods to deal with it and a 
# simulation study evaluating their performance. Ecography, 36(1), 27–46.