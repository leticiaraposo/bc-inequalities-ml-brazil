# =============================================================================
# Título: Spatially Varying Socioeconomic Determinants of Breast Cancer Mortality in Brazil
# Autores: Letícia Raposo
# Instituição: Universidade Federal do Estado do Rio de Janeiro (UNIRIO)
# Etapa 2: agregação dos dados municipais a nível de regiões geográficas 
# imediatas (RGIs)
# Padronização Direta por Idade com População de Referência (Censo 2010)
# Data: Julho de 2025

# Descrição Metodológica do Script:
# O presente script tem por finalidade realizar a análise da 
# mortalidade por câncer de mama em mulheres no Brasil no período
# de 2015 a 2019, com base em dados do Sistema de Informações 
# sobre Mortalidade (SIM/DATASUS) e projeções populacionais do 
# Instituto Brasileiro de Geografia e Estatística (IBGE). 
# O procedimento metodológico adota uma abordagem quantitativa,
# descritiva e comparativa, estruturada conforme as seguintes etapas:
# 
# 1. Aquisição e Processamento de Dados de Óbito (SIM/DATASUS)
# Inicialmente, são extraídos microdados do SIM via pacote microdatasus,
# contemplando óbitos ocorridos entre os anos de 2015 e 2019. 
# O conjunto de variáveis selecionadas inclui o código do município 
# de residência, a data do óbito, a causa básica de morte (CID-10),
# idade e sexo da pessoa falecida.
# Após o pré-processamento, são selecionados exclusivamente os óbitos
# por câncer de mama (CID-10: C50) em mulheres, com posterior 
# classificação por faixa etária padronizada.
# 
# 2. Importação de Dados Populacionais Femininos (IBGE)
# As estimativas populacionais femininas por município e faixa etária
# para os anos de 2015 a 2019 são importadas a partir de planilhas 
# oficiais do IBGE. Adicionalmente, os dados do Censo Demográfico de
# 2010 são utilizados como população padrão de referência, conforme
# recomendado para padronização direta de taxas de mortalidade.
# 
# 3. Cálculo da Taxa de Mortalidade Específica e Ajustada por Idade
# Os óbitos por câncer de mama são agregados por município e 
# faixa etária. Em paralelo, soma-se a população municipal no mesmo
# recorte.
# Com esses dados, calcula-se:
#   
# A taxa de mortalidade específica por faixa etária, definida 
# como o quociente entre o número de óbitos e a população correspondente;
# 
# A taxa de óbitos esperados, estimada com base na aplicação das 
# taxas específicas sobre a população padrão de 2010;
# 
# A taxa ajustada por idade, obtida pela soma dos óbitos esperados 
# e dividida pela população total padrão, expressa por 100 mil habitantes.
# Este procedimento permite a comparabilidade intermunicipal, 
# mitigando distorções decorrentes da estrutura etária heterogênea 
# entre as localidades.
# 
# 4. Agregação por Regiões Geográficas Imediatas (IBGE)
# A partir da vinculação dos códigos municipais aos agrupamentos 
# oficiais de Regiões Geográficas Imediatas (IBGE, 2024), replicam-se
# os cálculos da etapa anterior em nível regional. Isso permite 
# uma análise espacial mais agregada e estatisticamente robusta 
# da mortalidade ajustada.
# 
# 5. Cálculo de Médias Ponderadas de Indicadores Contextuais
# Indicadores municipais diversos (socioeconômicos, estruturais, etc.)
# previamente organizados são integrados à base por código municipal.
# Utiliza-se o peso populacional relativo de cada município na 
# respectiva região para o cálculo de médias ponderadas regionais, 
# garantindo representatividade proporcional.
# 
# 6. Inclusão de Coordenadas Geográficas para Análise Espacial
# São extraídos os centroides geográficos dos polígonos que representam
# as Regiões Geográficas Imediatas a partir de shapefiles do IBGE 
# (RGI/2017). Essas coordenadas são adicionadas à base analítica 
# final com vistas a análises espaciais.
# 
# 7. Armazenamento Final
# Todos os dados analíticos consolidados, incluindo as taxas 
# ajustadas, indicadores contextuais ponderados e coordenadas 
# geográficas, são armazenados em um arquivo .RDS, garantindo 
# reprodutibilidade e integridade da base para análises posteriores.
# =============================================================================

# --- 1. CONFIGURAÇÃO INICIAL ---

# Instalação de pacotes (descomente e execute se ainda não tiver instalado)
# install.packages("microdatasus")
# install.packages("dplyr")
# install.packages("readr")
# install.packages("readxl") # Necessário para ler arquivos .xlsx
# install.packages("stringr") # Necessário para manipulação de strings (CAUSABAS)
# install.packages("tidyr")   # Necessário para pivot_longer e replace_na

# Carregamento dos pacotes
library(microdatasus)
library(dplyr)
library(readr)
library(readxl)
library(stringr)
library(tidyr)

# URLs para as fontes de dados (para referência)
# População Geral (projeções): http://tabnet.datasus.gov.br/cgi/deftohtm.exe?ibge/cnv/popsvs2024br.def
# População Padrão (Censo 2010): http://tabnet.datasus.gov.br/cgi/deftohtm.exe?ibge/cnv/popbr.def



# --- 2. PREPARAÇÃO DOS DADOS DE ÓBITO (SIM) ---

# Baixando os dados brutos do SIM - 2015 a 2019 (5 anos) - não pegando dados da pandemia
dados_sim_raw <- fetch_datasus(
  year_start = 2015,
  year_end = 2019,
  information_system = "SIM-DO",
  vars = c("CODMUNRES", "DTOBITO", "CAUSABAS", "IDADE", "SEXO")
)

# Processar os dados com parse_sim
dados_sim <- process_sim(dados_sim_raw)

# Filtrar óbitos por câncer de mama em mulheres
library(stringr)
obitos_cancer_mama <- dados_sim %>%
  filter(
    str_starts(CAUSABAS, "C50"),
    SEXO == "Feminino"
  )

# Certifique-se de que a coluna está no formato Date
obitos_cancer_mama$DTOBITO <- as.Date(obitos_cancer_mama$DTOBITO)

# Extrai apenas o ano
obitos_cancer_mama$ANO_OBITO <- format(obitos_cancer_mama$DTOBITO, "%Y")

# Criando as faixas etárias
obitos_cancer_mama$IDADEanos <- as.numeric(obitos_cancer_mama$IDADEanos)
obitos_cancer_mama$FAIXA_ETARIA <- cut(obitos_cancer_mama$IDADEanos,
                                       breaks = c(-Inf, 4, 9, 14, 19, 29, 39, 49, 59, 69, 79, Inf),
                                       labels = c("0 a 4 anos", "5 a 9 anos", "10 a 14 anos", "15 a 19 anos",
                                                  "20 a 29 anos", "30 a 39 anos", "40 a 49 anos", "50 a 59 anos",
                                                  "60 a 69 anos", "70 a 79 anos", "80 anos e mais"),
                                       right = TRUE,
                                       include.lowest = TRUE
)

obitos_cancer_mama$DTOBITO <- NULL
obitos_cancer_mama$IDADEminutos <- NULL
obitos_cancer_mama$IDADEhoras <- NULL
obitos_cancer_mama$IDADEdias <- NULL
obitos_cancer_mama$IDADEmeses <- NULL



# --- 3. PREPARAÇÃO DOS DADOS DE POPULAÇÃO (IBGE) ---

# Lendo os dados de populacao feminina

pop_fem_2010 <- read_excel("Dados População Feminina Brasileira/pop_fem_2010.xlsx")
pop_fem_2010$`0 a 4 anos` <- pop_fem_2010$`Menor 1 ano` + pop_fem_2010$`1 a 4 anos`
pop_fem_2010$`Menor 1 ano` <- NULL
pop_fem_2010$`1 a 4 anos` <- NULL

pop_fem_2010 <- pop_fem_2010[,c(1,2,14,3:13)]

# Extrai apenas os 6 primeiros dígitos como código
pop_fem_2010 <- pop_fem_2010 %>%
  mutate(CODMUNRES = substr(Município, 1, 6))
pop_fem_2010$CODMUNRES <- as.character(pop_fem_2010$CODMUNRES)

pop_fem_2015 <- read_excel("Dados População Feminina Brasileira/pop_fem_2015.xlsx")
# Extrai apenas os 6 primeiros dígitos como código
pop_fem_2015 <- pop_fem_2015 %>%
  mutate(CODMUNRES = substr(Município, 1, 6))
pop_fem_2015$CODMUNRES <- as.character(pop_fem_2015$CODMUNRES)

pop_fem_2016 <- read_excel("Dados População Feminina Brasileira/pop_fem_2016.xlsx")
# Extrai apenas os 6 primeiros dígitos como código
pop_fem_2016 <- pop_fem_2016 %>%
  mutate(CODMUNRES = substr(Município, 1, 6))
pop_fem_2016$CODMUNRES <- as.character(pop_fem_2016$CODMUNRES)

pop_fem_2017 <- read_excel("Dados População Feminina Brasileira/pop_fem_2017.xlsx")
# Extrai apenas os 6 primeiros dígitos como código
pop_fem_2017 <- pop_fem_2017 %>%
  mutate(CODMUNRES = substr(Município, 1, 6))
pop_fem_2017$CODMUNRES <- as.character(pop_fem_2017$CODMUNRES)

pop_fem_2018 <- read_excel("Dados População Feminina Brasileira/pop_fem_2018.xlsx")
# Extrai apenas os 6 primeiros dígitos como código
pop_fem_2018 <- pop_fem_2018 %>%
  mutate(CODMUNRES = substr(Município, 1, 6))
pop_fem_2018$CODMUNRES <- as.character(pop_fem_2018$CODMUNRES)

pop_fem_2019 <- read_excel("Dados População Feminina Brasileira/pop_fem_2019.xlsx")
# Extrai apenas os 6 primeiros dígitos como código
pop_fem_2019 <- pop_fem_2019 %>%
  mutate(CODMUNRES = substr(Município, 1, 6))
pop_fem_2019$CODMUNRES <- as.character(pop_fem_2019$CODMUNRES)



# --- 4. CÁLCULO DA TAXA DE MORTALIDADE AJUSTADA POR IDADE ---

# 1. População nacional feminina 2010 por faixa etária (POPULAÇÃO PADRÃO)
# A população padrão para padronização direta é geralmente a soma das populações
# de todas as faixas etárias de uma população de referência (neste caso, Brasil 2010).

# Selecionar apenas as colunas de faixas etárias e somar para obter o total 
# nacional por faixa
pop_padrao_2010_nacional <- pop_fem_2010 %>%
  # Remover colunas que não são faixas etárias ou são totais já calculados
  select(-Ano, -Município, -Total, -CODMUNRES) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>% # Soma a população de cada faixa etária em nível nacional
  # Transforma o dataframe de formato largo para longo (tidy format)
  pivot_longer(cols = everything(), names_to = "FAIXA_ETARIA", values_to = "POP_PADRAO")

# Visualizar a população padrão para verificar
print("População Padrão (Brasil 2010 por Faixa Etária):")
print(pop_padrao_2010_nacional)

# 2. Unificar as populações municipais de 2015 a 2019
# Lista dos dataframes de população para cada ano
# Adapte esta lista se seus dataframes tiverem nomes diferentes (ex: pop_fem_2015, pop_fem_2016, etc.)
list_pop_municipais_anos <- list(pop_fem_2015, pop_fem_2016, 
                                 pop_fem_2017, pop_fem_2018, 
                                 pop_fem_2019)

# Combina todos os dataframes em um só
pop_municipais_2015_2019 <- bind_rows(list_pop_municipais_anos) %>%
  # Transforma o dataframe de formato largo para longo
  pivot_longer(
    cols = `0 a 4 anos`:`80 anos e mais`, # Seleciona as colunas de faixas etárias
    names_to = "FAIXA_ETARIA",
    values_to = "POP_MUN"
  ) %>%
  select(Ano, CODMUNRES, FAIXA_ETARIA, POP_MUN)

# Visualizar uma amostra dos dados de população municipal unificados
print("Amostra das Populações Municipais (2015-2019):")
print(head(pop_municipais_2015_2019))

# 3. Preparar os dados de óbitos por câncer de mama (2015-2019)
obitos_por_municipio_ano_faixa <- obitos_cancer_mama %>%
  filter(ANO_OBITO %in% 2015:2019) %>% # Filtra para os anos de interesse
  group_by(ANO_OBITO, CODMUNRES, FAIXA_ETARIA) %>% # Agrupa por ano, município e faixa etária
  summarise(OBITOS = n(), .groups = "drop") # Conta o número de óbitos em cada grupo

# Visualizar uma amostra dos óbitos preparados
print("Amostra dos Óbitos por Município, Ano e Faixa Etária:")
print(head(obitos_por_municipio_ano_faixa))

# 4. Calcular Taxas de Mortalidade Específicas por Idade por Município e Ano
# Junta os dados de óbitos com os dados de população municipal
pop_municipais_2015_2019$Ano <- as.character(pop_municipais_2015_2019$Ano)

# Agregue os dados de população municipal por município e faixa etária (somando os anos 2015 a 2019):

pop_mun_agr <- pop_municipais_2015_2019 %>%
  group_by(CODMUNRES, FAIXA_ETARIA) %>%
  summarise(POP_MUN = sum(POP_MUN, na.rm = TRUE), .groups = "drop")

# Agregue os óbitos por município e faixa etária (somando 2015–2019):

obitos_mun_agr <- obitos_por_municipio_ano_faixa %>%
  group_by(CODMUNRES, FAIXA_ETARIA) %>%
  summarise(OBITOS = sum(OBITOS, na.rm = TRUE), .groups = "drop")

# Juntar população e óbitos
base_calculo <- pop_mun_agr %>%
  left_join(obitos_mun_agr, by = c("CODMUNRES", "FAIXA_ETARIA")) %>%
  mutate(OBITOS = replace_na(OBITOS, 0)) %>%
  left_join(pop_padrao_2010_nacional, by = "FAIXA_ETARIA")

# Calcular a taxa específica de mortalidade por faixa (por município)

library(dplyr)

base_calculo <- base_calculo %>%
  mutate(
    TX_ESPECIFICA = OBITOS / POP_MUN,
    OBITOS_ESPERADOS = TX_ESPECIFICA * POP_PADRAO
  )

pop_padrao_total <- sum(pop_padrao_2010_nacional$POP_PADRAO)

taxa_ajustada_municipio <- base_calculo %>%
  group_by(CODMUNRES) %>%
  summarise(
    OBITOS_ESPERADOS = sum(OBITOS_ESPERADOS, na.rm = TRUE),
    TAXA_AJUSTADA = (OBITOS_ESPERADOS / pop_padrao_total) * 100000
  )

#---------------------#

library(dplyr)
library(tidyr)

# 1. Corrigir código de município para 6 dígitos
RELATORIO_DTB_BRASIL_2024_MUNICIPIOS <- read_excel("RELATORIO_DTB_BRASIL_2024_MUNICIPIOS.xls")

RELATORIO_DTB_BRASIL_2024_MUNICIPIOS <- RELATORIO_DTB_BRASIL_2024_MUNICIPIOS %>%
  mutate(codigo6 = substr(`Código Município Completo`, 1, 6))

# 2. Associar Regiões Imediatas aos dados de população e óbitos

# População
pop_municipais_com_regiao <- pop_municipais_2015_2019 %>%
  mutate(CODMUNRES = substr(CODMUNRES, 1, 6)) %>%
  left_join(
    RELATORIO_DTB_BRASIL_2024_MUNICIPIOS %>%
      select(codigo6, `Região Geográfica Imediata`),
    by = c("CODMUNRES" = "codigo6")
  )

# Óbitos
obitos_com_regiao <- obitos_por_municipio_ano_faixa %>%
  mutate(CODMUNRES = substr(CODMUNRES, 1, 6)) %>%
  left_join(
    RELATORIO_DTB_BRASIL_2024_MUNICIPIOS %>%
      select(codigo6, `Região Geográfica Imediata`),
    by = c("CODMUNRES" = "codigo6")
  )

# 3. Agregar por Região Imediata e Faixa Etária

pop_regiao <- pop_municipais_com_regiao %>%
  group_by(`Região Geográfica Imediata`, FAIXA_ETARIA) %>%
  summarise(POP_REGIAO = sum(POP_MUN, na.rm = TRUE), .groups = "drop")

obitos_regiao <- obitos_com_regiao %>%
  group_by(`Região Geográfica Imediata`, FAIXA_ETARIA) %>%
  summarise(OBITOS = sum(OBITOS, na.rm = TRUE), .groups = "drop")

# 4. Juntar com a população padrão nacional

base_regional <- pop_regiao %>%
  left_join(obitos_regiao, by = c("Região Geográfica Imediata", "FAIXA_ETARIA")) %>%
  mutate(OBITOS = replace_na(OBITOS, 0)) %>%
  left_join(pop_padrao_2010_nacional, by = "FAIXA_ETARIA")

# 5. Calcular taxa específica e óbitos esperados por faixa

base_regional <- base_regional %>%
  mutate(
    TX_ESPECIFICA = OBITOS / POP_REGIAO,
    OBITOS_ESPERADOS = TX_ESPECIFICA * POP_PADRAO
  )

# 6. Calcular taxa ajustada por idade (por 100 mil habitantes)

pop_padrao_total <- sum(pop_padrao_2010_nacional$POP_PADRAO, na.rm = TRUE)

taxa_ajustada_regiao <- base_regional %>%
  group_by(`Região Geográfica Imediata`) %>%
  summarise(
    OBITOS_ESPERADOS = sum(OBITOS_ESPERADOS, na.rm = TRUE),
    TAXA_AJUSTADA = (OBITOS_ESPERADOS / pop_padrao_total) * 100000,
    .groups = "drop"
  )

# 7. Visualizar resultado
print(taxa_ajustada_regiao)

# Dados coletados pela Livia
library(readxl)
dados <- read_excel("dados_CM_Ferrete_Raposo.xlsx")
# Convertendo a coluna 'codigo' em 'dados' para character para garantir a junção correta
dados$codigo <- as.character(dados$codigo)
head(dados)

dados$taxa_bruta_mortalidade <- NULL # errada!

# Garantir código compatível
RELATORIO_DTB_BRASIL_2024_MUNICIPIOS <- RELATORIO_DTB_BRASIL_2024_MUNICIPIOS %>%
  mutate(codigo6 = substr(`Código Município Completo`, 1, 6))

pop_municipios_total <- pop_municipais_2015_2019 %>%
  group_by(CODMUNRES) %>%
  summarise(Populacao_Total = sum(POP_MUN, na.rm = TRUE), .groups = "drop")

# Corrigir código se necessário
pop_municipios_total <- pop_municipios_total %>%
  mutate(CODMUNRES = substr(CODMUNRES, 1, 6))

dados <- dados %>%
  mutate(codigo = substr(codigo, 1, 6)) %>%
  left_join(pop_municipios_total, by = c("codigo" = "CODMUNRES"))

dados <- dados %>%
  left_join(
    RELATORIO_DTB_BRASIL_2024_MUNICIPIOS %>%
      mutate(codigo6 = substr(`Código Município Completo`, 1, 6)) %>%
      select(codigo6, `Região Geográfica Imediata`),
    by = c("codigo" = "codigo6")
  )

pop_total_regiao <- pop_regiao %>%
  group_by(`Região Geográfica Imediata`) %>%
  summarise(pop_total_regiao = sum(POP_REGIAO, na.rm = TRUE))

# Identificar variáveis a ponderar
variaveis_numericas <- dados %>%
  select(where(is.numeric)) %>%
  colnames()

variaveis_ponderadas <- setdiff(variaveis_numericas, c("Y", "X", "Populacao_Total"))

# Calcular fração da população de cada município na sua região
dados_com_peso <- dados %>%
  left_join(pop_total_regiao, by = "Região Geográfica Imediata") %>%
  mutate(peso_municipal = Populacao_Total / pop_total_regiao)

# Média ponderada por região
media_ponderada_regiao <- dados_com_peso %>%
  group_by(`Região Geográfica Imediata`) %>%
  summarise(
    across(all_of(variaveis_ponderadas), ~ weighted.mean(.x, peso_municipal, na.rm = TRUE), .names = "mp_{.col}"),
    .groups = "drop"
  )

# Unir taxa bruta à ajustada
base_final <- taxa_ajustada_regiao %>%
  left_join(media_ponderada_regiao, by = "Região Geográfica Imediata")

# Gerando os centroides das regioes imediatas
library(dplyr)
library(janitor)

# Padronizar nomes das colunas
relatorio_limpo <- RELATORIO_DTB_BRASIL_2024_MUNICIPIOS %>%
  clean_names() %>%
  mutate(codigo6 = substr(codigo_municipio_completo, 1, 6))

# Garantir que 'dados$codigo' também tem 6 dígitos
dados_corrigido <- dados %>%
  mutate(codigo = substr(codigo, 1, 6)) %>%
  left_join(relatorio_limpo %>% select(codigo6, regiao_geografica_imediata),
            by = c("codigo" = "codigo6"))

# Verifique se a coluna foi criada corretamente:
stopifnot("regiao_geografica_imediata" %in% names(dados_corrigido))

# Agora calcule as coordenadas ponderadas
coordenadas_regiao <- dados_corrigido %>%
  group_by(regiao_geografica_imediata) %>%
  summarise(
    X_regiao = weighted.mean(X, Populacao_Total, na.rm = TRUE),
    Y_regiao = weighted.mean(Y, Populacao_Total, na.rm = TRUE),
    .groups = "drop"
  )


library(sf)

# Supondo que você tenha o shapefile das regiões imediatas:
shapefile_regioes <- st_read("RG2017_rgi_20180911/RG2017_rgi.shp")

sf::sf_use_s2(FALSE)

shapefile_regioes <- shapefile_regioes %>%
  st_zm(drop = TRUE, what = "ZM")  # remove Z e M


# Calcular o centroide de cada polígono
centroides <- shapefile_regioes %>%
  mutate(centroide = st_centroid(geometry)) %>%
  mutate(
    X = st_coordinates(centroide)[, 1],
    Y = st_coordinates(centroide)[, 2]
  ) %>%
  select(rgi, nome_rgi, X, Y)

base_final <- base_final %>%
  left_join(centroides, by = c("Região Geográfica Imediata" = "rgi"))


saveRDS(base_final, "base_regioes_imediatas_cancer_mama_periodo_2015_2019.RDS")
#------------------------#