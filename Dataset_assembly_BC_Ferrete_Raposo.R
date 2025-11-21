# =====================================================================================================
# Título: Spatially Varying Socioeconomic Determinants of Breast Cancer Mortality in Brazil
# Autores: Lívia Ferrete, Letícia Raposo
# Instituição: Universidade Federal do Estado do Rio de Janeiro (UNIRIO)
# Etapa 1: montagem do dataset utilizado no estudo
# Data: Março de 2025
#
# Descrição:
# Este script realiza a montagem do banco de dados utilizado nas análises realizadas em um estudo
# de inferência acerca da influência de fatores socioeconômicos e de saúde na mortalidade
# por neoplasia maligna da mama no Brasil. Foram utilizadas variáveis presentes em 4 grandes
# aplicações - Ipeadata (IPEA), Proadess (Fiocruz), Sidra (IBGE) e Tabnet (Datasus/SUS).
# 
# Saídas:
# - dados_CM_Ferrete_Raposo.xls: planilha de variáveis com os municípios brasileiros  
#   como unidade de observação, salva como pasta de trabalho do Excel.
#
# Observações:
# - para as planilhas com entradas para vários anos, foi preservado separadamente o ano de 2010,
#   por ser o ano de realização do último censo demográfico com dados completos disponíveis e 
#   foi calculada a media dos valores (Anderson et al, 2023) encontrados entre os anos de 2010 a 2019;
# - certificar de que os pacotes do R necessários para cada comando estejam instalados;
# - alterar os caminhos dos diretórios e nomes de documentos como for necessário.
# =====================================================================================================

# 1. Baixando pacotes e carregando bibliotecas ----

# Instalando pacotes
install.packages(c("readxl", "dplyr", "stringr", "openxlsx", 
                        "geobr", "sf", "foreign"))

# Carregando bibliotecas necessárias
library(readxl) # Leitura de arquivos excel
library(dplyr) # Manipulação de dados
library(stringr) # Manipulação de strings
library(openxlsx) # Exportação de dados em excel
library(geobr) # Download de dados geográficos
library(sf) # Para simple features (modelagem espacial)
library(foreign) # Leitura de e escrita em arquivos


# 2. Baixando indicadores ----

## 2.1. Planilhas do Sistema IBGE de Recuperação Automática (SIDRA) ----

# Nas planilhas obtidas pelo Sidra, os códigos dos municípios possuem 7 dígitos.

# Densidade demográfica
densidade_demografica = read_excel("densidade_demografica.xlsx")

# Domicílios chefiados por mulheres
domicilio_chefe_mulher = read_excel("domicilio_mulher_responsavel.xlsx")

# População por raça/cor
raca = read_excel("populacao_residente_raca.xlsx",
                  col_types = c("text", "numeric", "numeric", "numeric", "numeric", "numeric"))

# População casada
casados = read_excel("populacao_casada.xlsx")

# População por sexo
sexo = read_excel("populacao_sexo.xlsx")

# Densidade domiciliar
densidade_domiciliar = read_excel("densidade_domiciliar.xlsx", 
                                   col_types = c("text", "numeric", "numeric", 
                                                 "numeric", "numeric"))

# População cujo trabalho excede 44h semanais
trab_excedente = read_excel("populacao_trab_excedente.xlsx")

# População urbana x rural
urbana_rural = read_excel("populacao_urbana_rural.xlsx",
                          col_types = c("text", "numeric", "numeric"))

# Domicílios com geladeira
geladeira = read_excel("domicilios_com_geladeira.xlsx")

# Domicílios com pelo menos um microcomputador para uso pessoal e acesso a internet
domicilios_desktop_internet = read_excel("domicilios_com_pc_internet.xlsx",
                                         col_types = c("text", "numeric"))

# População com deficiência
pcd = read_excel("populacao_PCD.xlsx")

# População por tempo de deslocamento para o trabalho
deslocamento = read_excel("tempo_deslocamento_trabalho.xlsx",
                          col_types = c("text", "numeric", "numeric", "numeric", "numeric", "numeric"))

# População residente em domicílios com entorno arborizado
domicilio_arborizado = read_excel("domicilio_entorno_arborizado.xlsx",
                                  col_types = c("text", "numeric"))

# População residente em domicílios com entorno pavimentado
domicilio_entorno_pavimentado = read_excel("domicilio_entorno_pavimentado.xlsx")

# População idosa
idosos = read_excel("populacao_idosos.xlsx")

# População evangélica
evangelicos = read_excel("populacao_evangelica.xlsx")

# População imigrante
estrangeiros = read_excel("populacao_estrangeira.xlsx",
                          col_types = c("text", "numeric"))

# Domicíios com pelo menos um automóvel particular
automovel = read_excel("domicilios_com_automovel.xlsx",
                       col_types = c("text", "numeric"))

# Trabalho infantil
trab_infantil = read_excel("trab_infantil.xlsx",
                           col_types = c("text", "numeric"))

# Aposentados ou pensionistas
aposentados_pensionistas = read_excel("aposentados_ou_pensionistas.xlsx")

# População de 25 anos ou mais com ensino médio completo
medio_completo = read_excel("medio_completo.xlsx")

# População de 25 anos ou mais com ensino superior completo
superior_completo = read_excel("superior_completo.xlsx")



## 2.2. Planilhas do Ipeadata ----

# Nas planilhas obtidas pelo Ipeadata, os códigos dos municípios possuem 7 dígitos.

# População residindo em domicílio com água encanada
agua_encanada = read_excel("pop_agua_encanada.xls")

# Índice de Desenvolvimento Humano Municipal (IDHM)
idhm = read_excel("idh_municipal.xls")

# Taxa de fecundidade
fecundidade = read_excel("taxa_fecundidade.xls")

# Renda per capita media
renda_per_capita = read_excel("renda_per_capita.xls")

# População residindo em domicílio com energia elétrica
eletricidade = read_excel("populacao_com_eletricidade.xls")

# Índice de Gini
gini = read_excel("indice_gini.xls")

# Domicílios com rede de esgoto e/ou fossa séptica
esgoto_ou_fossa = read_excel("esgoto_ou_fossa_septica.xls")

# Taxa de analfabetismo (15 anos ou mais)
analfabetismo = read_excel("analfabetismo_15_mais.xls")

# Proporção de pobres
pobres = read_excel("proporcao_pobres.xls")

# Proporção de extremamente pobres
extremamente_pobres = read_excel("proporcao_extremamente_pobres.xls")

# Mortalidade infantil
mortalidade_infantil = read_excel("mortalidade_infantil.xls")

# População residente em domicílios com coleta de lixo
coleta_lixo = read_excel("coleta_lixo.xls")

# Razão de dependência
razao_dependencia = read_excel("razao_dependencia.xls")

# Expectativa de vida
expectativa_de_vida = read_excel("expectativa_de_vida.xls")

# Taxa de envelhecimento
envelhecimento = read_excel("envelhecimento.xls")

# Taxa de homicídio
homicidio = read_excel("taxa_homicidio.xls")
homicidio = homicidio %>%  
  select(-c("2011":"2019")) %>% 
  rename(homicidios_2010 = `2010`)

# Taxa de suicídio
suicidio = read_excel("taxa_suicidio.xls")
suicidio = suicidio %>% 
  select(-c("2011":"2019")) %>% 
  rename(suicidios_2010 = `2010`)



## 2.3. Planilhas do Tabnet ----

# Nas planilhas obtidas pelo Tabnet, os códigos dos municípios possuem 6 dígitos.
# Para estas planilhas, foi necessário extrair apenas o código da célula, visto que o
# Tabnet não cria duas colunas separadas para o nome do município e para o código.

# Taxa de desemprego
desemprego = read_excel("desemprego.xlsx")

# Extraindo o código da coluna com o nome do município
desemprego = desemprego %>% 
  mutate(codigo = str_extract(codigo, "^\\d+")) # remove tudo antes e depois da 
                                                # sequência numérica

# Produto Interno Bruto (PIB) per capita
pib_per_capita = read_excel("pib_per_capita.xlsx")
pib_per_capita = pib_per_capita %>% 
  mutate(codigo = str_extract(codigo, "^\\d+")) 



## 2.4. Planilhas do Proadess ----

# Nas planilhas obtidas pelo Proadess, os códigos dos municípios possuem 6 dígitos.

# Mamógrafos em uso
mamografos = read_excel("mamografos.xlsx",
                        col_types = c("text", "numeric", "numeric", "numeric", "numeric",
                                      "numeric", "numeric", "numeric", "numeric", "numeric",
                                      "numeric"))
mamografos = mamografos %>% 
  select(-c("2011":"2019")) %>% 
  rename(mamografos_2010 = `2010`)

# Nascidos vivos com pré-natal adequado
nasc_vivos_pre_natal = read_excel("./Indicadores/calcular_medias/nasc_vivos_pre_natal.xlsx",
                                  col_types = c("text", "numeric", "numeric", "numeric", "numeric",
                                                "numeric", "numeric", "numeric", "numeric", "numeric",
                                                "numeric"))
nasc_vivos_pre_natal = nasc_vivos_pre_natal %>% 
  select(-c("2011":"2019")) %>% 
  rename(nasc_vivos_pre_natal_2010 = `2010`)

  # População coberta por planos e seguros de saúde
pop_planos_seguros = read_excel("./Indicadores/calcular_medias/cobertura_planos_seguros.xlsx",
                                col_types = c("text", "numeric", "numeric", "numeric", "numeric",
                                              "numeric", "numeric", "numeric", "numeric", "numeric",
                                              "numeric"))
pop_planos_seguros = pop_planos_seguros %>% 
  select(-c("2011":"2019")) %>% 
  rename(pop_planos_seguros_2010 = `2010`)

# Recursos municipais destinados à saúde
recursos_destinados_saude = read_excel("./Indicadores/calcular_medias/recursos_destinados_saude.xlsx",
                                       col_types = c("text", "numeric", "numeric", "numeric", "numeric",
                                                     "numeric", "numeric", "numeric", "numeric", "numeric",
                                                     "numeric"))
recursos_destinados_saude = recursos_destinados_saude %>% 
  select(-c("2011":"2019")) %>% 
  rename(recursos_destinados_saude_2010 = `2010`)

# Total de médicos disponíveis
medicos = read_excel("./Indicadores/calcular_medias/medicos.xlsx",
                     col_types = c("text", "numeric", "numeric", "numeric", "numeric",
                                   "numeric", "numeric", "numeric", "numeric", "numeric",
                                   "numeric"))
medicos = medicos %>% 
  select(-c("2011":"2019")) %>% 
  rename(medicos_2010 = `2010`)

# Total de médicos ginecologistas e obstetras ou mastologistas e cirurgiões da mama disponíveis
medicos_gyn_ob = read_excel("./Indicadores/calcular_medias/gyn_ob.xlsx",
                            col_types = c("text", "numeric", "numeric", "numeric", "numeric",
                                          "numeric", "numeric", "numeric", "numeric", "numeric",
                                          "numeric"))

medicos_gyn_ob = medicos_gyn_ob %>% 
  select(-c("2011":"2019")) %>% 
  rename(medicos_gyn_ob_2010 = `2010`)


# Mortalidade por doenças no aparelho circulatório (estresse, estilo de vida)
mortalidade_ap_circulatorio = read_excel("./Indicadores/calcular_medias/mortalidade_ap_circulatorio.xlsx",
                                         col_types = c("text", "numeric", "numeric", "numeric", "numeric",
                                                       "numeric", "numeric", "numeric", "numeric", "numeric",
                                                       "numeric"))
mortalidade_ap_circulatorio = mortalidade_ap_circulatorio %>% 
  select(-c("2011":"2019")) %>% 
  rename(mortalidade_ap_circulatorio_2010 = `2010`)

# Ultrassons em uso
ultrassons = read_excel("./Indicadores/calcular_medias/ultrassons.xlsx",
                        col_types = c("text", "numeric", "numeric", "numeric", "numeric",
                                      "numeric", "numeric", "numeric", "numeric", "numeric",
                                      "numeric"))
ultrassons = ultrassons %>%
  select(-c("2011":"2019")) %>% 
  rename(ultrassons_2010 = `2010`)
  
# Cobertura da Estratégia Saúde da Família
cobertura_esf = read_excel("./Indicadores/calcular_medias/cobertura_ESF.xlsx",
                           col_types = c("text", "numeric", "numeric", "numeric", "numeric",
                                         "numeric", "numeric", "numeric", "numeric", "numeric",
                                         "numeric"))
cobertura_esf = cobertura_esf %>% 
  select(-c("2011":"2019")) %>% 
  rename(cobertura_esf_2010 = `2010`)

# Total de enfermeiros disponíveis
enfermeiros = read_excel("./Indicadores/calcular_medias/enfermeiros.xlsx",
                         col_types = c("text", "numeric", "numeric", "numeric", "numeric",
                                       "numeric", "numeric", "numeric", "numeric", "numeric",
                                       "numeric"))
enfermeiros = enfermeiros %>% 
  select(-c("2011":"2019")) %>% 
  rename(enfermeiros_2010 = `2010`)

# Leitos hospitalares disponíveis
leitos_hosp = read_excel("./Indicadores/calcular_medias/leitos_hospitalares.xlsx",
                         col_types = c("text", "numeric", "numeric", "numeric", "numeric",
                                       "numeric", "numeric", "numeric", "numeric", "numeric",
                                       "numeric"))
leitos_hosp = leitos_hosp %>% 
  select(-c("2011":"2019")) %>% 
  rename(leitos_hosp_2010 = `2010`)

# Cobertura da Atenção Básica
cobertura_ab = read_excel("./Indicadores/calcular_medias/cobertura_AB.xlsx",
                          col_types = c("text", "numeric", "numeric", "numeric", "numeric",
                                        "numeric", "numeric", "numeric", "numeric", "numeric",
                                        "numeric"))
cobertura_ab = cobertura_ab %>% 
  select(-c("2011":"2019")) %>% 
  rename(cobertura_ab_2010 = `2010`)

# Procedimentos de mamografia
razao_mamografia = read_excel("./Indicadores/calcular_medias/razao_mamografia.xlsx",
                              col_types = c("text", "numeric", "numeric", "numeric", "numeric",
                                            "numeric", "numeric", "numeric", "numeric", "numeric",
                                            "numeric"))
razao_mamografia = razao_mamografia %>% 
  select(-c("2011":"2019")) %>% 
  rename(razao_mamografia_2010 = `2010`)



# 3. Baixando dados geográficos dos municípios ----
geo_municipios = read_municipality(code_muni = "all",
                  year = 2010,
                  simplified = F, # indicado para análises estatísticas espaciais
                  keep_areas_operacionais = F
)

# Extraindo coordenadas
coords = st_coordinates(st_point_on_surface(geo_municipios)) # Pode haver erros na 
                                                             # longitude e latitude



# 4. Padronizando todos os códigos ----

# Função que passa todos os códigos (num, char, float ou outros) para numérico, 
# ao mesmo tempo que padroniza as entradas. Como visto anteriormente, algumas tabelas 
# apresentam códigos com 6 dígitos, enquanto outras apresentam códigos com 7 dígitos;
# sendo assim, preservaremos apenas os 6 primeiros dígitos de todas os códigos.
padronizar_codigo = function(tabela) {
  if ("codigo" %in% colnames(tabela)) {
    
    # converte para caracter, extrai os 6 primeiros dígitos e então converte para numérico
    # (primeira etapa é necessária, pois a função substr() não funciona em variáveis ""num")
    tabela$codigo = as.numeric(substr(as.character(tabela$codigo), 1, 6))
  }
  return(tabela)
}

# Lista todos os objetos no ambiente que são dfs/tabelas
objetos = ls()[sapply(ls(), function(x) is.data.frame(get(x)))]

# Aplicando o código de padronização em todas as tabelas
for (obs in objetos) {
  tabela = get(obs)
  tabela = padronizar_codigo(tabela)      # aplica o código de padronização
  assign(obs, tabela, envir = .GlobalEnv) # atualiza tabela no ambiente
}



# 5. Juntando todas as tabelas ----

# Criando uma lista com todos os bancos necessários no ambiente
lista_tabelas = list(agua_encanada, analfabetismo, aposentados_pensionistas, automovel, casados, 
               cobertura_ab, cobertura_esf, coleta_lixo, densidade_demografica, densidade_domiciliar, desemprego, 
               deslocamento, domicilio_arborizado, domicilio_chefe_mulher, domicilio_entorno_pavimentado,
               domicilios_desktop_internet, eletricidade, enfermeiros, envelhecimento, esgoto_ou_fossa,
               estrangeiros, evangelicos, expectativa_de_vida, extremamente_pobres, fecundidade, geladeira,
               gini, homicidio, idhm, idosos, leitos_hosp, mamografos, medicos, medicos_gyn_ob, medio_completo,               
               mortalidade_ap_circulatorio, mortalidade_infantil, nasc_vivos_pre_natal, pcd, pib_per_capita,              
               pobres, pop_planos_seguros, raca, razao_dependencia, razao_mamografia, recursos_destinados_saude,
               renda_per_capita, sexo, suicidio, superior_completo, trab_excedente, trab_infantil, ultrassons, 
               urbana_rural, taxa_mortalidade)

# Juntando os bancos em um único dataset
dados = Reduce(function(x, y) merge(x, y, by = "codigo", all = T), lista_tabelas)

# Tirando o Brasil do banco, mantendo apenas municípios
dados = dados[-1, ]

# Removendo as linhas cujos municípios não têm "nome" preenchido
# (algumas planilhas do Ipeadata geraram muitos campos vazios)
dados = dados %>% 
  filter(!is.na(nome))

# Inserindo as coordenadas no banco
dados$X = coords[, "X"] # Longitude
dados$Y = coords[, "Y"] # Latitude

# Reorganizando as colunas para vir código > coordenadas > resto
dados = dados %>% 
  select(codigo, nome, Y, X, taxa_bruta_mortalidade, everything())

# Visualizando a tabela
View(dados)

# 6. Salvando o banco ----
write.xlsx(dados, "dados_CM_Ferrete_Raposo.xlsx")