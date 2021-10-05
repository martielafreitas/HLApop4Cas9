###############################################################################
# VIA DE PROCESSAMENTO
###############################################################################

# Unindo os resulatdos dos preditores para as três proteínas 
# em uma única tabela que adiciona o nome que está no nome do
# arquivo, na última coluna da tabela.

setwd("/Users/martielafreitas/Desktop/HLA-Rproject")

# variavel que contem o caminho para os arquivos
file_list <- list.files(path="/Users/martielafreitas/Desktop/HLA-Rproject/CInputs/Preditores/1.Processing")

# inicia uma dataframe vazia
process.hlac <- data.frame()

# onde serao processados os arquivos
setwd("/Users/martielafreitas/Desktop/HLA-Rproject/CInputs/Preditores/1.Processing")

# laco que une todos os arquivos
for(i in 1:length(file_list)){
  #each file will be read in, specify which columns you need read in to avoid any errors
  temp_data <- read.csv(file_list[i], sep=",")
  #clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
  temp_data$Class <- sapply(strsplit(gsub(".csv", "", file_list[i]), " - "), function(x){x[3]})
  #for each iteration, bind the new data to the building dataset
  process.hlac <- rbind(process.hlac, temp_data)
}

# salvando o resultado em um arquivo .csv
write.csv(process.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Processing_Geral.csv")

# Separando cada resultado por organismo de onde vem a proteina
aureus.process.hlac<- subset(process.hlac, Class == "Aureus")
pyogenes.process.hlac <- subset(process.hlac, Class == "Pyogenes")
thermophilus.process.hlac <- subset(process.hlac, Class == "Thermophilus")

# salvando o resultado em um arquivo .csv
write.csv(aureus.process.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Process.csv")
write.csv(pyogenes.process.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Process.csv")
write.csv(thermophilus.process.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Process.csv")


###############################################################################
# LIGAÇÃO
###############################################################################

##-- Binding Results --##

# variavel que contem o caminho para os arquivos
file_list <- list.files(path="/Users/martielafreitas/Desktop/HLA-Rproject/CInputs/Preditores/2.Binding")

#inicia uma dataframe vazia
binding.hlac <- data.frame()

setwd("/Users/martielafreitas/Desktop/HLA-Rproject/CInputs/Preditores/2.Binding")
for(i in 1:length(file_list)){
  #each file will be read in, specify which columns you need read in to avoid any errors
  temp_data <- read.csv(file_list[i], sep=";")
  #clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
  temp_data$Class <- sapply(strsplit(gsub(".csv", "", file_list[i]), " - "), function(x){x[3]})
  #for each iteration, bind the new data to the building dataset
  binding.hlac <- rbind(binding.hlac, temp_data)
}
# salvando resultados gerais de ligacao em um arquivo .csv
write.csv(binding.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Binding_Geral.csv")

# separando os resultados por organismo de onde veio cada proteina
aureus.binding.hlac <- subset(binding.hlac, Class == "Aureus")
pyogenes.binding.hlac <- subset(binding.hlac, Class == "Pyogenes")
thermophilus.binding.hlac <- subset(binding.hlac, Class == "Thermophilus")

# salvando resultados gerais por organismo em um arquivo .csv
write.csv(aureus.binding.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Binding.csv")
write.csv(pyogenes.binding.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Binding.csv")
write.csv(thermophilus.binding.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Binding.csv")



###############################################################################
# IMMUNOGENICIDADE
###############################################################################

##-- Comparing Results for immunogenicity predictions --##

aureus.procbind.hlac <- merge(aureus.binding.hlac, 
                              aureus.process.hlac, 
                              by.x=c("allele", "seq_num", "start", "end", "length", "peptide", "Class"),
                              by.y=c("Allele", "X.", "Start", "End", "Peptide.Length", "Peptide", "Class"))

pyogenes.procbind.hlac <- merge(pyogenes.binding.hlac, 
                                pyogenes.process.hlac, 
                                by.x=c("allele", "seq_num", "start", "end", "length", "peptide", "Class"),
                                by.y=c("Allele", "X.", "Start", "End", "Peptide.Length", "Peptide", "Class"))

thermophilus.procbind.hlac <- merge(thermophilus.binding.hlac, 
                                    thermophilus.process.hlac, 
                                    by.x=c("allele", "seq_num", "start", "end", "length", "peptide", "Class"),
                                    by.y=c("Allele", "X.", "Start", "End", "Peptide.Length", "Peptide", "Class"))

write.csv(aureus.procbind.hlac, file="/Users/martielafreitas/Documents/Rprojects/HLA-C/Results/HLA-C_Aureus_Merged_Process_and_Binding.csv")
write.csv(pyogenes.procbind.hlac, file="/Users/martielafreitas/Documents/Rprojects/HLA-C/Results/HLA-C_Pyogenes_Merged_Process_and_Binding.csv")
write.csv(thermophilus.procbind.hlac, file="/Users/martielafreitas/Documents/Rprojects/HLA-C/Results/HLA-C_Thermophilus_Merged_Process_and_Binding.csv")



###############################################################################
# ANALISE DE RESULTADOS
###############################################################################

## Immunogenicity Results

## Obtaining list of peptides

## Peptídeos totais - independente de população ou frequência, por organismo
aureus.peptides.hlac <- unique(aureus.procbind.hlac$peptide)
pyogenes.peptides.hlac <- unique(pyogenes.procbind.hlac$peptide)
thermophilus.peptides.hlac <- unique(thermophilus.procbind.hlac$peptide)

write.csv(aureus.peptides.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_TotalPeptides.csv")
write.csv(pyogenes.peptides.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_TotalPeptides.csv")
write.csv(thermophilus.peptides.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_TotalPeptides.csv")

# Esses peptídeos passaram pela ferramenta de predição de imunogenicidade.

#- Calling result files -#
file_list <- list.files(path="/Users/martielafreitas/Desktop/HLA-Rproject/CInputs/Preditores/3.Immunogenicity")
immunogenicity.hlac <- data.frame()

setwd("/Users/martielafreitas/Desktop/HLA-Rproject/CInputs/Preditores/3.Immunogenicity")
for(i in 1:length(file_list)){
  #each file will be read in, specify which columns you need read in to avoid any errors
  temp_data <- read.csv(file_list[i], sep=",")
  #clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
  temp_data$Class <- sapply(strsplit(gsub(".csv", "", file_list[i]), " - "), function(x){x[3]})
  #for each iteration, bind the new data to the building dataset
  immunogenicity.hlac <- rbind(immunogenicity.hlac, temp_data)
}
# salvando resultados gerais de imunogenicidade em um arquivo .csv
write.csv(immunogenicity.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Immunogenicity_Geral.csv")


# Separando a lista de HLAs totais por organismo - independente de população ou frequência

aureus.immuno.hlac <- subset(immunogenicity.hlac, Class == "Aureus")
pyogenes.immuno.hlac <- subset(immunogenicity.hlac, Class == "Pyogenes")
thermophilus.immuno.hlac <- subset(immunogenicity.hlac, Class == "Thermophilus")

write.csv(aureus.immuno.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Immunogenicity.csv")
write.csv(pyogenes.immuno.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Immunogenicity.csv")
write.csv(thermophilus.immuno.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Immunogenicity.csv")


###############################################################################
# Filtering final results
###############################################################################

#- Immuno: score >= 0.0: immunogenic -#

setwd("/Users/martielafreitas/Desktop/HLA-Rproject/Results/")

aureus.immuno.positive.hlac <- subset(aureus.immuno.hlac, score >= 0.0)
pyogenes.immuno.positive.hlac <- subset(pyogenes.immuno.hlac, score >= 0.0)
thermophilus.immuno.positive.hlac <- subset(thermophilus.immuno.hlac, score >= 0.0)

#- Merging tables to IC50 -#
colnames(aureus.procbind.hlac)
colnames(aureus.immuno.positive.hlac)

aureus.full.path.hlac <- merge(aureus.procbind.hlac, aureus.immuno.positive.hlac, by.x = "peptide", by.y = "peptide")
pyogenes.full.path.hlac <- merge(pyogenes.procbind.hlac, pyogenes.immuno.positive.hlac, by.x = "peptide", by.y = "peptide")
thermophilus.full.path.hlac <- merge(thermophilus.procbind.hlac, thermophilus.immuno.positive.hlac, by.x = "peptide", by.y = "peptide")

#- IC50 < 500nM  -#
aureus.IC50.hlac <- aureus.full.path.hlac[(aureus.full.path.hlac$MHC.IC50 <= 500.0),]
pyogenes.IC50.hlac <- pyogenes.full.path.hlac[(pyogenes.full.path.hlac$MHC.IC50 <= 500.0),]
thermophilus.IC50.hlac <- thermophilus.full.path.hlac[(thermophilus.full.path.hlac$MHC.IC50 <= 500.0),]

#- Final Table -#
aureus.results.hlac <- aureus.IC50.hlac[,c(1,2,9,10,11,12,13,14,15,8,17)]
pyogenes.results.hlac <- pyogenes.IC50.hlac[,c(1,2,9,10,11,12,13,14,15,8,17)]
thermophilus.results.hlac <- thermophilus.IC50.hlac[,c(1,2,9,10,11,12,13,14,15,8,17)]

# corrigindo cabeçalho antes de salvar! :)
colnames(aureus.results.hlac)
cols <- c("Peptide", "Allele", "Percentile.Rank", "Proteasome.Score", "TAP.Score",
          "MHC.Score", "Processing.Score", "Total.Score", "MHC.IC50", "Binding.Score", "Immunogenicity.Score")
colnames(aureus.results.hlac) <- cols
colnames(pyogenes.results.hlac) <- cols
colnames(thermophilus.results.hlac) <- cols

write.csv(aureus.results.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_ToltalResults.csv")
write.csv(pyogenes.results.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_TotalResults.csv")
write.csv(thermophilus.results.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_TotalResults.csv")

#- Obtaining list of alleles por organismo-#
aureus.final.allele.hlac <- data.frame(unique(aureus.results.hlac$Allele)) #221
pyogenes.final.allele.hlac <- data.frame(unique(pyogenes.results.hlac$Allele)) #225
thermophilus.final.allele.hlac <- data.frame(unique(thermophilus.results.hlac$Allele)) #200

write.csv(aureus.final.allele.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_TotalHLAResults.csv")
write.csv(pyogenes.final.allele.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_TotalHLAResults.csv")
write.csv(thermophilus.final.allele.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_TotalHLAResults.csv")

#- Obtaining list of peptides por organismo -#
aureus.final.peptides.hlac <- data.frame(unique(aureus.results.hlac$Peptide)) #121
pyogenes.final.peptides.hlac <- data.frame(unique(pyogenes.results.hlac$Peptide)) #163
thermophilus.final.peptides.hlac <- data.frame(unique(thermophilus.results.hlac$Peptide)) #189

write.csv(aureus.final.peptides.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_TotalPEPResults.csv")
write.csv(pyogenes.final.peptides.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_TotalPEPResults.csv")
write.csv(thermophilus.final.peptides.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_TotalPEPResults.csv")


common.aurpyo.hlac <- data.frame(merge(aureus.final.peptides.hlac, pyogenes.final.peptides.hlac, 
                                  by.x = colnames(aureus.final.peptides.hlac), by.y = colnames(pyogenes.final.peptides.hlac)))
common.aurthe.hlac <- data.frame(merge(aureus.final.peptides.hlac, thermophilus.final.peptides.hlac, 
                                  by.x = colnames(aureus.final.peptides.hlac), by.y = colnames(thermophilus.final.peptides.hlac)))
common.pyothe.hlac <- data.frame(merge(pyogenes.final.peptides.hlac, thermophilus.final.peptides.hlac, 
                                  by.x = colnames(pyogenes.final.peptides.hlac), by.y = colnames(thermophilus.final.peptides.hlac)))

# Qual o score de imunogenicidade pra os comuns pyo/the? 16 no total!
score.common.pyothe.hlac <- merge(common.pyothe.hlac, immunogenicity.hlac, by.x="unique.pyogenes.results.hlac.Peptide.", by.y = "peptide")
write.csv(score.common.pyothe.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_CommonPeptides_Pyogenes_Thermophilus.csv")


###############################################################################
# Logos
###############################################################################

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Aureus_PeptidesLogo.pdf")
ggplot() + geom_logo( aureus.final.peptides.hlac ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Pyogenes_PeptidesLogo.pdf")
ggplot() + geom_logo( pyogenes.final.peptides.hlac ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Thermophilus_PeptidesLogo.pdf")
ggplot() + geom_logo( thermophilus.final.peptides.hlac ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_PyoThe_PeptidesLogo.pdf")
ggplot() + geom_logo( common.pyothe.hlac ) + theme_logo()
dev.off()


###############################################################################
# Populations 
###############################################################################
setwd("/Users/martielafreitas/Desktop/HLA-Rproject/CInputs/Populacoes")
populacoes.hlac <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/CInputs/Populacoes/HLA-C_GOLD_Geral", sep=";")

# Merged tables
aureus.pop.hlac <-  merge(aureus.results.hlac, populacoes.hlac, by.x = "Allele", by.y = "Allele")
pyogenes.pop.hlac <-  merge(pyogenes.results.hlac, populacoes.hlac, by.x = "Allele", by.y = "Allele")
thermophilus.pop.hlac <-  merge(thermophilus.results.hlac, populacoes.hlac, by.x = "Allele", by.y = "Allele")

write.table(aureus.pop.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Path_Populacao.csv")
write.table(pyogenes.pop.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Path_Populacao.csv")
write.table(thermophilus.pop.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Path_Populacao.csv")

# Unique peptides by Population and allele

aureus.populacao.hlac <- aureus.pop.hlac[,c(1,2,13,14,15,16)]
pyogenes.populacao.hlac <- pyogenes.pop.hlac[,c(1,2,13,14,15,16)]
thermophilus.populacao.hlac <- thermophilus.pop.hlac[,c(1,2,13,14,15,16)]

aureus.ethnicity.hlac <- aureus.populacao.hlac[,c(1,2,6)]
pyogenes.ethnicity.hlac <- pyogenes.populacao.hlac[,c(1,2,6)]
thermophilus.ethnicity.hlac <- thermophilus.populacao.hlac[,c(1,2,6)]

aureus.ethnicity.unq.hlac <- unique(aureus.ethnicity.hlac)
pyogenes.ethnicity.unq.hlac <- unique(pyogenes.ethnicity.hlac)
thermophilus.ethnicity.unq.hlac <-unique(thermophilus.ethnicity.hlac)

write.table(aureus.ethnicity.unq.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Ethnicity_Unique.csv")
write.table(pyogenes.ethnicity.unq.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Ethnicity_Unique.csv")
write.table(thermophilus.ethnicity.unq.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Ethnicity_Unique.csv")

### PEPTIDES

aureus.sankey.hlac <- aureus.populacao.hlac[,c(1,2,6,4)]
pyogenes.sankey.hlac <- pyogenes.populacao.hlac[,c(1,2,6,4)]
thermophilus.sankey.hlac <- thermophilus.populacao.hlac[,c(1,2,6,4)]

write.table(aureus.sankey.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Sankey.csv")
write.table(pyogenes.sankey.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Sankey.csv")
write.table(thermophilus.sankey.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Sankey.csv")

#aureus.sankey <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Sankey_HLA-C.csv", sep = ";")
#pyogenes.sankey <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Sankey_HLA-C.csv", sep = ";")
#thermophilus.sankey <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Sankey_HLA-C.csv", sep = ";")

############################################################################### AUREUS

##################### AUREUS 500 ORIENTAL

aureus.sankey.oriental.hlac <- unique(subset(aureus.sankey.hlac, EthnicOrigin == "Oriental"))
aureus.sankey.oriental.hlac <- aggregate(x = as.numeric(aureus.sankey.oriental.hlac$Frequency),
                                    by = c(list(aureus.sankey.oriental.hlac$Allele, 
                                                aureus.sankey.oriental.hlac$Peptide, 
                                                aureus.sankey.oriental.hlac$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlac <- subset(aureus.sankey.oriental.hlac, x >= 0.10)
# exclude lines wth e-4...n
#aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlac$x),
                            by = c(list(aureus.sankey.top.hlac$Group.1, 
                                        aureus.sankey.top.hlac$Group.2, 
                                        aureus.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(aureus.snk.agg$Group.1),
  target = c(aureus.snk.agg$Group.2),
  value = aureus.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Aureus_Oriental_500_10up.html", selfcontained = TRUE)

rm(sankey, aureus.sankey.top, aureus.snk.agg, links, source, target, value, nodes, name, rede)

##################### AUREUS 500 BLACK

aureus.sankey.black.hlac <- unique(subset(aureus.sankey.hlac, EthnicOrigin == "Black"))
aureus.sankey.black.hlac <- aggregate(x = as.numeric(aureus.sankey.black.hlac$Frequency),
                                 by = c(list(aureus.sankey.black.hlac$Allele, 
                                             aureus.sankey.black.hlac$Peptide, 
                                             aureus.sankey.black.hlac$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlac <- subset(aureus.sankey.black.hlac, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlac$x),
                            by = c(list(aureus.sankey.top.hlac$Group.1, 
                                        aureus.sankey.top.hlac$Group.2, 
                                        aureus.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(aureus.snk.agg$Group.1),
  target = c(aureus.snk.agg$Group.2),
  value = aureus.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Aureus_Black_500_10up.html", selfcontained = TRUE)

rm(sankey, aureus.sankey.top, aureus.snk.agg, links, source, target, value, nodes, name, rede)

##################### AUREUS 500 CAUCASOID

aureus.sankey.caucasoid.hlac <- unique(subset(aureus.sankey.hlac, EthnicOrigin == "Caucasoid"))
aureus.sankey.caucasoid.hlac <- aggregate(x = as.numeric(aureus.sankey.caucasoid.hlac$Frequency),
                                     by = c(list(aureus.sankey.caucasoid.hlac$Allele, 
                                                 aureus.sankey.caucasoid.hlac$Peptide, 
                                                 aureus.sankey.caucasoid.hlac$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlac <- subset(aureus.sankey.caucasoid.hlac, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlac$x),
                            by = c(list(aureus.sankey.top.hlac$Group.1, 
                                        aureus.sankey.top.hlac$Group.2, 
                                        aureus.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(aureus.snk.agg$Group.1),
  target = c(aureus.snk.agg$Group.2),
  value = aureus.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Aureus_Caucasoid_500_10up.html", selfcontained = TRUE)

rm(sankey, aureus.sankey.top, aureus.snk.agg, links, source, target, value, nodes, name, rede)


#### RASCUNHO SE PRECISAR
#common.aur.asbl <- data.frame(merge(aureus.sankey.oriental, aureus.sankey.black, 
#                                    by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))
#common.aur.blca <- data.frame(merge(aureus.sankey.black, aureus.sankey.cauc, 
#                                    by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))
#common.aur.asca <- data.frame(merge(aureus.sankey.oriental, aureus.sankey.cauc, 
#                                    by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))

#common.aur.asbl <- unique(common.aur.asbl[,c(1,2)])
#common.aur.blca <- unique(common.aur.blca[,c(1,2)])
#common.aur.asca <- unique(common.aur.asca[,c(1,2)])

#combined <- rbind(common.aur.asbl, common.aur.blca, common.aur.asca)
#duplicate_rows <- unique(combined[duplicated(combined), ])
#combined_peptides <- unique(duplicate_rows$Group.2)

#write.table(combined_peptides, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Aureus_Peptides_3pop_HLA-C.csv")
##write.table(duplicate_rows, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Aureus_Common_3pop_HLA-C.csv")



################################################################################ PYOGENES

##################### PYOGENES 500 ORIENTAL

pyogenes.sankey.oriental.hlac <- unique(subset(pyogenes.sankey.hlac, EthnicOrigin == "Oriental"))
pyogenes.sankey.oriental.hlac <- aggregate(x = as.numeric(pyogenes.sankey.oriental.hlac$Frequency),
                                      by = c(list(pyogenes.sankey.oriental.hlac$Allele, 
                                                  pyogenes.sankey.oriental.hlac$Peptide, 
                                                  pyogenes.sankey.oriental.hlac$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlac <- subset(pyogenes.sankey.oriental.hlac, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlac$x),
                            by = c(list(pyogenes.sankey.top.hlac$Group.1, 
                                        pyogenes.sankey.top.hlac$Group.2, 
                                        pyogenes.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(pyogenes.snk.agg$Group.1),
  target = c(pyogenes.snk.agg$Group.2),
  value = pyogenes.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Pyogenes_Oriental_500_10up.html", selfcontained = TRUE)

rm(sankey, pyogenes.sankey.top, pyogenes.snk.agg, links, source, target, value, nodes, name, rede)


##################### PYOGENES 500 BLACK
pyogenes.sankey.black.hlac <- unique(subset(pyogenes.sankey.hlac, EthnicOrigin == "Black"))
pyogenes.sankey.black.hlac <- aggregate(x = as.numeric(pyogenes.sankey.black.hlac$Frequency),
                                   by = c(list(pyogenes.sankey.black.hlac$Allele, 
                                               pyogenes.sankey.black.hlac$Peptide, 
                                               pyogenes.sankey.black.hlac$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlac <- subset(pyogenes.sankey.black.hlac, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlac$x),
                              by = c(list(pyogenes.sankey.top.hlac$Group.1, 
                                          pyogenes.sankey.top.hlac$Group.2, 
                                          pyogenes.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(pyogenes.snk.agg$Group.1),
  target = c(pyogenes.snk.agg$Group.2),
  value = pyogenes.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Pyogenes_Black_500_10up.html", selfcontained = TRUE)

rm(sankey, pyogenes.sankey.top, pyogenes.snk.agg, links, source, target, value, nodes, name, rede)


##################### PYOGENES 500 CAUCASOID
pyogenes.sankey.cauc.hlac <- unique(subset(pyogenes.sankey.hlac, EthnicOrigin == "Caucasoid"))
pyogenes.sankey.cauc.hlac <- aggregate(x = as.numeric(pyogenes.sankey.cauc.hlac$Frequency),
                                  by = c(list(pyogenes.sankey.cauc.hlac$Allele, 
                                              pyogenes.sankey.cauc.hlac$Peptide, 
                                              pyogenes.sankey.cauc.hlac$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlac <- subset(pyogenes.sankey.oriental.hlac, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlac$x),
                            by = c(list(pyogenes.sankey.top.hlac$Group.1, 
                                        pyogenes.sankey.top.hlac$Group.2, 
                                        pyogenes.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(pyogenes.snk.agg$Group.1),
  target = c(pyogenes.snk.agg$Group.2),
  value = pyogenes.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Pyogenes_Caucasoid_500_10up.html", selfcontained = TRUE)

rm(sankey, pyogenes.sankey.top, pyogenes.snk.agg, links, source, target, value, nodes, name, rede)


### RASCUNHO SE FOR PRECISO

#common.pyo.asbl <- data.frame(merge(pyogenes.sankey.oriental, pyogenes.sankey.black, 
#                                    by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))
#common.pyo.blca <- data.frame(merge(pyogenes.sankey.black, pyogenes.sankey.cauc, 
#                                    by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))
#common.pyo.asca <- data.frame(merge(pyogenes.sankey.oriental, pyogenes.sankey.cauc, 
#                                    by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))

#common.pyo.asbl <- common.pyo.asbl[,c(1,2)]
#common.pyo.blca <- common.pyo.blca[,c(1,2)]
#common.pyo.asca <- common.pyo.asca[,c(1,2)]

#combined <- rbind(common.pyo.asbl, common.pyo.blca, common.pyo.asca)
#duplicate_rows <- unique(combined[duplicated(combined), ])
#combined_peptides <- unique(duplicate_rows$Group.2)

##write.table(combined_peptides, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Pyogenes_Peptides_3pop_HLA-C.csv")
##write.table(duplicate_rows, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Pyogenes_Common_3pop_HLA-C.csv")


############################################################################### THERMOPHILUS

##################### THERMOPHILUS 500 ORIENTAL

thermophilus.sankey.oriental.hlac <- unique(subset(thermophilus.sankey.hlac, EthnicOrigin == "Oriental"))
thermophilus.sankey.oriental.hlac <- aggregate(x = as.numeric(thermophilus.sankey.oriental.hlac$Frequency),
                                          by = c(list(thermophilus.sankey.oriental.hlac$Allele, 
                                                      thermophilus.sankey.oriental.hlac$Peptide, 
                                                      thermophilus.sankey.oriental.hlac$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlac <- subset(thermophilus.sankey.oriental.hlac, x >= 0.10)
# exclude lines wth e-4...n
#aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlac$x),
                            by = c(list(thermophilus.sankey.top.hlac$Group.1, 
                                        thermophilus.sankey.top.hlac$Group.2, 
                                        thermophilus.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(thermophilus.snk.agg$Group.1),
  target = c(thermophilus.snk.agg$Group.2),
  value = thermophilus.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Thermophilus_Oriental_500_10up.html", selfcontained = TRUE)

rm(sankey, thermophilus.sankey.top, thermophilus.snk.agg, links, source, target, value, nodes, name, rede)

##################### THERMOPHILUS 500 BLACK

thermophilus.sankey.black.hlac <- unique(subset(thermophilus.sankey.hlac, EthnicOrigin == "Black"))
thermophilus.sankey.black.hlac <- aggregate(x = as.numeric(thermophilus.sankey.black.hlac$Frequency),
                                       by = c(list(thermophilus.sankey.black.hlac$Allele, 
                                                   thermophilus.sankey.black.hlac$Peptide, 
                                                   thermophilus.sankey.black.hlac$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlac <- subset(thermophilus.sankey.black.hlac, x >= 0.10)
# exclude lines wth e-4...n
# thermophilus.sankey.10 <- thermophilus.sankey.top[- grep("e-", thermophilus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlac$x),
                            by = c(list(thermophilus.sankey.top.hlac$Group.1, 
                                        thermophilus.sankey.top.hlac$Group.2, 
                                        thermophilus.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(thermophilus.snk.agg$Group.1),
  target = c(thermophilus.snk.agg$Group.2),
  value = thermophilus.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/Thermophilus_Black_500_10up.html", selfcontained = TRUE)

rm(sankey, thermophilus.sankey.top, thermophilus.snk.agg, links, source, target, value, nodes, name, rede)


##################### THERMOPHILUS 500 CAUCASOID

thermophilus.sankey.cauc.hlac <- unique(subset(thermophilus.sankey.hlac, EthnicOrigin == "Caucasoid"))
thermophilus.sankey.cauc.hlac <- aggregate(x = as.numeric(thermophilus.sankey.cauc.hlac$Frequency),
                                      by = c(list(thermophilus.sankey.cauc.hlac$Allele, 
                                                  thermophilus.sankey.cauc.hlac$Peptide, 
                                                  thermophilus.sankey.cauc.hlac$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlac <- subset(thermophilus.sankey.cauc.hlac, x >= 0.10)
# exclude lines wth e-4...n
# thermophilus.sankey.10 <- thermophilus.sankey.top[- grep("e-", thermophilus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlac$x),
                            by = c(list(thermophilus.sankey.top.hlac$Group.1, 
                                        thermophilus.sankey.top.hlac$Group.2, 
                                        thermophilus.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(thermophilus.snk.agg$Group.1),
  target = c(thermophilus.snk.agg$Group.2),
  value = thermophilus.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Thermophilus_Caucasoid_500_10up.html", selfcontained = TRUE)

rm(sankey, thermophilus.sankey.top, thermophilus.snk.agg, links, source, target, value, nodes, name, rede)



### RASCUNHO SE PRECISAR
#common.ther.asbl <- data.frame(merge(thermophilus.sankey.oriental, thermophilus.sankey.black, 
#                                     by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))
#common.ther.blca <- data.frame(merge(thermophilus.sankey.black, thermophilus.sankey.cauc, 
#                                     by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))
#common.ther.asca <- data.frame(merge(thermophilus.sankey.oriental, thermophilus.sankey.cauc, 
#                                     by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))

#common.ther.asbl <- common.ther.asbl[,c(1,2)]
#common.ther.blca <- common.ther.blca[,c(1,2)]
#common.ther.asca <- common.ther.asca[,c(1,2)]

#combined <- rbind(common.ther.asbl, common.ther.blca, common.ther.asca)
#duplicate_rows <- unique(combined[duplicated(combined), ])
#combined_peptides <- unique(duplicate_rows$Group.2)

##write.table(combined_peptides, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Thermophilus_Peptides_3pop_HLA-C.csv")
##write.table(duplicate_rows, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Thermophilus_Common_3pop_HLA-C.csv")


###############################################################################


### ANALISE PARA IC50 < 30

aureus.IC50.30.hlac <- aureus.full.path.hlac[(aureus.full.path.hlac$MHC.IC50 <= 30.0),]
pyogenes.IC50.30.hlac <- pyogenes.full.path.hlac[(pyogenes.full.path.hlac$MHC.IC50 <= 30.0),]
thermophilus.IC50.30.hlac <- thermophilus.full.path.hlac[(thermophilus.full.path.hlac$MHC.IC50 <= 30.0),]

#- Final Table -#

aureus.results.30.hlac <- aureus.IC50.30.hlac[,c(1,2,9,10,11,12,13,14,15,8,17)]
pyogenes.results.30.hlac <- pyogenes.IC50.30.hlac[,c(1,2,9,10,11,12,13,14,15,8,17)]
thermophilus.results.30.hlac <- thermophilus.IC50.30.hlac[,c(1,2,9,10,11,12,13,14,15,8,17)]

# corrigindo cabeçalho antes de salvar! :)
colnames(aureus.results.30.hlac)
cols <- c("Peptide", "Allele", "Percentile.Rank", "Proteasome.Score", "TAP.Score",
          "MHC.Score", "Processing.Score", "Total.Score", "MHC.IC50", "Binding.Score", "Immunogenicity.Score")
colnames(aureus.results.30.hlac) <- cols
colnames(pyogenes.results.30.hlac) <- cols
colnames(thermophilus.results.30.hlac) <- cols

write.csv(aureus.results.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Results_IC50_30.csv")
write.csv(pyogenes.results.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Results_IC50_30.csv")
write.csv(thermophilus.results.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Results_IC50_30.csv")

##-- Summary --##

#- Obtaining list of alleles-#

aureus.final.allele.30.hlac <- data.frame(unique(aureus.results.30.hlac$Allele)) #221
pyogenes.final.allele.30.hlac <- data.frame(unique(pyogenes.results.30.hlac$Allele)) #225
thermophilus.final.allele.30.hlac <- data.frame(unique(thermophilus.results.30.hlac$Allele)) #200

write.csv(aureus.final.allele.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Results_HLA_IC50_30.csv")
write.csv(pyogenes.final.allele.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Results_HLA_IC50_30.csv")
write.csv(thermophilus.final.allele.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Results_HLA_IC50_30.csv")



#- Obtaining list of peptides -#

aureus.final.peptides.30.hlac <- data.frame(unique(aureus.results.30.hlac$Peptide)) #121
pyogenes.final.peptides.30.hlac <- data.frame(unique(pyogenes.results.30.hlac$Peptide)) #163
thermophilus.final.peptides.30.hlac <- data.frame(unique(thermophilus.results.30.hlac$Peptide)) #189


write.csv(aureus.final.peptides.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Results_PEP_IC50_30.csv")
write.csv(pyogenes.final.peptides.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Results_PEP_IC50_30.csv")
write.csv(thermophilus.final.peptides.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Results_PEP_IC50_30.csv")


common.aurpyo.30.hlac <- data.frame(merge(aureus.final.peptides.30.hlac, pyogenes.final.peptides.30.hlac, 
                                     by.x = colnames(aureus.final.peptides.30.hlac), by.y = colnames(pyogenes.final.peptides.30.hlac)))
common.aurthe.30.hlac <- data.frame(merge(aureus.final.peptides.30.hlac, thermophilus.final.peptides.30.hlac, 
                                     by.x = colnames(aureus.final.peptides.30.hlac), by.y = colnames(thermophilus.final.peptides.30.hlac)))
common.pyothe.30.hlac <- data.frame(merge(pyogenes.final.peptides.30.hlac, thermophilus.final.peptides.30.hlac, 
                                     by.x = colnames(pyogenes.final.peptides.30.hlac), by.y = colnames(thermophilus.final.peptides.30.hlac)))


write.csv(common.aurpyo.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Common_HLAwithPyogenes_IC50_30.csv")
write.csv(common.aurthe.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Common_HLAwithThermophilus_IC50_30.csv")
write.csv(common.pyothe.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Common_HLAwithThermophilus_IC50_30.csv")


# Qual o score de imunogenicidade pra os comuns pyo/the? 16 no total!
score.common.pyothe.30.hlac <- merge(common.pyothe.30.hlac, immunogenicity.hlac, by.x="unique.pyogenes.results.30.hlac.Peptide.", by.y = "peptide")
write.csv(score.common.pyothe.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_CommonScore_PEP_PyoThe_IC50_30.csv")


####
## LOGOS IC50 < 30nM
####


pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Aureus_Peptides_IC50_30.pdf")
ggplot() + geom_logo( aureus.final.peptides.30.hlac ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Pyogenes_Peptides_IC50_30.pdf")
ggplot() + geom_logo( pyogenes.final.peptides.30.hlac ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Thermophilus_peptides_IC50_30.pdf")
ggplot() + geom_logo( thermophilus.final.peptides.30.hlac ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_PyoThe_peptides_IC50_30.pdf")
ggplot() + geom_logo( common.pyothe.30.hlac ) + theme_logo()
dev.off()

##--    Populations   --#

setwd("/Users/martielafreitas/Desktop/HLA-Rproject/CInputs/Populacoes")
populacoes.hlac <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/CInputs/Populacoes/HLA-C_GOLD_Geral", sep=";")

# Merged tables
aureus.pop.30.hlac <-  merge(aureus.results.30.hlac, populacoes.hlac, by.x = "Allele", by.y = "Allele")
pyogenes.pop.30.hlac <-  merge(pyogenes.results.30.hlac, populacoes.hlac, by.x = "Allele", by.y = "Allele")
thermophilus.pop.30.hlac <-  merge(thermophilus.results.30.hlac, populacoes.hlac, by.x = "Allele", by.y = "Allele")

write.table(aureus.pop.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Path_Populacao_IC50_30.csv")
write.table(pyogenes.pop.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Path_Populacao_IC50_30.csv")
write.table(thermophilus.pop.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Path_Populacao_IC50_30.csv")

# Unique peptides by Population and allele

aureus.populacao.30.hlac <- aureus.pop.30.hlac[,c(1,2,13,14,15,16)]
pyogenes.populacao.30.hlac <- pyogenes.pop.30.hlac[,c(1,2,13,14,15,16)]
thermophilus.populacao.30.hlac <- thermophilus.pop.30.hlac[,c(1,2,13,14,15,16)]

aureus.ethnicity.30.hlac <- aureus.populacao.30.hlac[,c(1,2,6)]
pyogenes.ethnicity.30.hlac <- pyogenes.populacao.30.hlac[,c(1,2,6)]
thermophilus.ethnicity.30.hlac <- thermophilus.populacao.30.hlac[,c(1,2,6)]

aureus.ethnicity.unq.30.hlac <- unique(aureus.ethnicity.30.hlac)
pyogenes.ethnicity.unq.30.hlac <- unique(pyogenes.ethnicity.30.hlac)
thermophilus.ethnicity.unq.30.hlac <-unique(thermophilus.ethnicity.30.hlac)

write.table(aureus.ethnicity.unq.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Ethnicity_Unique_IC50_30.csv")
write.table(pyogenes.ethnicity.unq.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Ethnicity_Unique_IC50_30.csv")
write.table(thermophilus.ethnicity.unq.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Ethnicity_Unique_IC50_30.csv")

### PEPTIDES

aureus.sankey.30.hlac <- aureus.populacao.30.hlac[,c(1,2,6,4)]
pyogenes.sankey.30.hlac <- pyogenes.populacao.30.hlac[,c(1,2,6,4)]
thermophilus.sankey.30.hlac <- thermophilus.populacao.30.hlac[,c(1,2,6,4)]

write.table(aureus.sankey.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Sankey_PEP_IC50_30.csv")
write.table(pyogenes.sankey.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Sankey_PEP_IC50_30.csv")
write.table(thermophilus.sankey.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Sankey_PEP_IC50_30.csv")

#aureus.sankey.30 <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Sankey_HLA-C.30.csv", sep = ";")
#pyogenes.sankey.30 <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Sankey_HLA-C.30.csv", sep = ";")
#thermophilus.sankey.30 <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Sankey_HLA-C.30.csv", sep = ";")


################################################################################## AUREUS

###################### AUREUS ORIENTAL 30

aureus.sankey.oriental.30.hlac <- unique(subset(aureus.sankey.30.hlac, EthnicOrigin == "Oriental"))
aureus.sankey.oriental.30.hlac <- aggregate(x = as.numeric(aureus.sankey.oriental.30.hlac$Frequency),
                                       by = c(list(aureus.sankey.oriental.30.hlac$Allele, 
                                                   aureus.sankey.oriental.30.hlac$Peptide, 
                                                   aureus.sankey.oriental.30.hlac$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlac <- subset(aureus.sankey.oriental.30.hlac, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlac$x),
                            by = c(list(aureus.sankey.top.hlac$Group.1, 
                                        aureus.sankey.top.hlac$Group.2, 
                                        aureus.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(aureus.snk.agg$Group.1),
  target = c(aureus.snk.agg$Group.2),
  value = aureus.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Aureus_Oriental_30_10up.html", selfcontained = TRUE)

rm(sankey, aureus.sankey.top.hlac, aureus.snk.agg, links, source, target, value, nodes, name, rede)

###################### BLACK 30

aureus.sankey.black.30.hlac <- unique(subset(aureus.sankey.30.hlac, EthnicOrigin == "Black"))
aureus.sankey.black.30.hlac <- aggregate(x = as.numeric(aureus.sankey.black.30.hlac$Frequency),
                                    by = c(list(aureus.sankey.black.30.hlac$Allele, 
                                                aureus.sankey.black.30.hlac$Peptide, 
                                                aureus.sankey.black.30.hlac$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlac <- subset(aureus.sankey.black.30.hlac, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlac$x),
                            by = c(list(aureus.sankey.top.hlac$Group.1, 
                                        aureus.sankey.top.hlac$Group.2, 
                                        aureus.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(aureus.snk.agg$Group.1),
  target = c(aureus.snk.agg$Group.2),
  value = aureus.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Aureus_Black_30_10up.html", selfcontained = TRUE)

rm(sankey, aureus.sankey.top.hlac, aureus.snk.agg, links, source, target, value, nodes, name, rede)


##################### CAUCASOIDE 30

aureus.sankey.cauc.30.hlac <- unique(subset(aureus.sankey.30.hlac, EthnicOrigin == "Caucasoid"))
aureus.sankey.cauc.30.hlac <- aggregate(x = as.numeric(aureus.sankey.cauc.30.hlac$Frequency),
                                   by = c(list(aureus.sankey.cauc.30.hlac$Allele, 
                                               aureus.sankey.cauc.30.hlac$Peptide, 
                                               aureus.sankey.cauc.30.hlac$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlac <- subset(aureus.sankey.cauc.30.hlac, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlac$x),
                            by = c(list(aureus.sankey.top.hlac$Group.1, 
                                        aureus.sankey.top.hlac$Group.2, 
                                        aureus.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(aureus.snk.agg$Group.1),
  target = c(aureus.snk.agg$Group.2),
  value = aureus.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Aureus_Caucasoid_30_10up.html", selfcontained = TRUE)

rm(sankey,aureus.sankey.top.hlac, aureus.snk.agg, links, source, target, value, nodes, name, rede)


### RASCUNHO SE PRECISAR
# common.aur.asbl.30 <- data.frame(merge(aureus.sankey.oriental.30, aureus.sankey.black.30, 
#                                        by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))
# common.aur.blca.30 <- data.frame(merge(aureus.sankey.black.30, aureus.sankey.cauc.30, 
#                                        by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))
# common.aur.asca.30 <- data.frame(merge(aureus.sankey.oriental.30, aureus.sankey.cauc.30, 
#                                        by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))
# 
# common.aur.asbl.30 <- unique(common.aur.asbl.30[,c(1,2)])
# common.aur.blca.30 <- unique(common.aur.blca.30[,c(1,2)])
# common.aur.asca.30 <- unique(common.aur.asca.30[,c(1,2)])
# 
# combined.30 <- rbind(common.aur.asbl.30, common.aur.blca.30, common.aur.asca.30)
# duplicate_rows.30 <- unique(combined.30[duplicated(combined.30), ])
# combined_peptides.30 <- unique(duplicate_rows.30$Group.2)

##write.table(combined_peptides.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Aureus_Peptides_3pop_HLA-C.30.csv")
##write.table(duplicate_rows.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Aureus_Common_3pop_HLA-C.30.csv")


################################################################################ PYOGENES

###################### PYOGENES ORIENTAL 30

pyogenes.sankey.oriental.30.hlac <- unique(subset(pyogenes.sankey.30.hlac, EthnicOrigin == "Oriental"))
pyogenes.sankey.oriental.30.hlac <- aggregate(x = as.numeric(pyogenes.sankey.oriental.30.hlac$Frequency),
                                         by = c(list(pyogenes.sankey.oriental.30.hlac$Allele, 
                                                     pyogenes.sankey.oriental.30.hlac$Peptide, 
                                                     pyogenes.sankey.oriental.30.hlac$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlac <- subset(pyogenes.sankey.oriental.30.hlac, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlac$x),
                            by = c(list(pyogenes.sankey.top.hlac$Group.1, 
                                        pyogenes.sankey.top.hlac$Group.2, 
                                        pyogenes.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(pyogenes.snk.agg$Group.1),
  target = c(pyogenes.snk.agg$Group.2),
  value = pyogenes.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Pyogenes_Oriental_30_10up.html", selfcontained = TRUE)

rm(sankey, pyogenes.sankey.top.hlac, pyogenes.snk.agg, links, source, target, value, nodes, name, rede)

###################### PYOGENES BLACK 30

pyogenes.sankey.black.30.hlac <- unique(subset(pyogenes.sankey.30.hlac, EthnicOrigin == "Black"))
pyogenes.sankey.black.30.hlac <- aggregate(x = as.numeric(pyogenes.sankey.black.30.hlac$Frequency),
                                      by = c(list(pyogenes.sankey.black.30.hlac$Allele, 
                                                  pyogenes.sankey.black.30.hlac$Peptide, 
                                                  pyogenes.sankey.black.30.hlac$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlac <- subset(pyogenes.sankey.black.30.hlac, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlac$x),
                            by = c(list(pyogenes.sankey.top.hlac$Group.1, 
                                        pyogenes.sankey.top.hlac$Group.2, 
                                        pyogenes.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(pyogenes.snk.agg$Group.1),
  target = c(pyogenes.snk.agg$Group.2),
  value = pyogenes.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Pyogenes_Black_30_10up.html", selfcontained = TRUE)

rm(sankey, pyogenes.sankey.top.hlac, pyogenes.snk.agg, links, source, target, value, nodes, name, rede)

###################### PYOGENES CAUCASOID 30

pyogenes.sankey.cauc.30.hlac <- unique(subset(pyogenes.sankey.30.hlac, EthnicOrigin == "Caucasoid"))
pyogenes.sankey.cauc.30.hlac <- aggregate(x = as.numeric(pyogenes.sankey.cauc.30.hlac$Frequency),
                                     by = c(list(pyogenes.sankey.cauc.30.hlac$Allele, 
                                                 pyogenes.sankey.cauc.30.hlac$Peptide, 
                                                 pyogenes.sankey.cauc.30.hlac$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlac <- subset(pyogenes.sankey.cauc.30.hlac, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlac$x),
                            by = c(list(pyogenes.sankey.top.hlac$Group.1, 
                                        pyogenes.sankey.top.hlac$Group.2, 
                                        pyogenes.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(pyogenes.snk.agg$Group.1),
  target = c(pyogenes.snk.agg$Group.2),
  value = pyogenes.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Pyogenes_Caucasoid_30_10up.html", selfcontained = TRUE)

rm(sankey, pyogenes.sankey.top.hlac, pyogenes.snk.agg, links, source, target, value, nodes, name, rede)



### RASUCNHO SE PRECISAR
# common.pyo.asbl.30 <- data.frame(merge(pyogenes.sankey.oriental.30, pyogenes.sankey.black.30, 
#                                        by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))
# common.pyo.blca.30 <- data.frame(merge(pyogenes.sankey.black.30, pyogenes.sankey.cauc.30, 
#                                        by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))
# common.pyo.asca.30 <- data.frame(merge(pyogenes.sankey.oriental.30, pyogenes.sankey.cauc.30, 
#                                        by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))
# 
# common.pyo.asbl.30 <- common.pyo.asbl.30[,c(1,2)]
# common.pyo.blca.30 <- common.pyo.blca.30[,c(1,2)]
# common.pyo.asca.30 <- common.pyo.asca.30[,c(1,2)]
# 
# combined <- rbind(common.pyo.asbl.30, common.pyo.blca.30, common.pyo.asca.30)
# duplicate_rows <- unique(combined.30[duplicated(combined.30), ])
# combined_peptides <- unique(duplicate_rows.30$Group.2)

##write.table(combined_peptides.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Pyogenes_Peptides_3pop_HLA-C.30.csv")
##write.table(duplicate_rows.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Pyogenes_Common_3pop_HLA-C.30.csv")


################################################################################ THERMOPHILUS

###################### THERMOPHILUS 30 ORIENTAL

thermophilus.sankey.oriental.30.hlac <- unique(subset(thermophilus.sankey.30.hlac, EthnicOrigin == "Oriental"))
thermophilus.sankey.oriental.30.hlac <- aggregate(x = as.numeric(thermophilus.sankey.oriental.30.hlac$Frequency),
                                             by = c(list(thermophilus.sankey.oriental.30.hlac$Allele, 
                                                         thermophilus.sankey.oriental.30.hlac$Peptide, 
                                                         thermophilus.sankey.oriental.30.hlac$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlac <- subset(thermophilus.sankey.oriental.30.hlac, x >= 0.10)
# exclude lines wth e-4...n
# thermophilus.sankey.10 <- thermophilus.sankey.top[- grep("e-", thermophilus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlac$x),
                            by = c(list(thermophilus.sankey.top.hlac$Group.1, 
                                        thermophilus.sankey.top.hlac$Group.2, 
                                        thermophilus.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(thermophilus.snk.agg$Group.1),
  target = c(thermophilus.snk.agg$Group.2),
  value = thermophilus.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Thermophilus_Oriental_30_10up.html", selfcontained = TRUE)

rm(sankey, thermophilus.sankey.top.hlac, thermophilus.snk.agg, links, source, target, value, nodes, name, rede)


##################### THERMOPHILUS 30 BLACK

thermophilus.sankey.black.30.hlac <- unique(subset(thermophilus.sankey.30.hlac, EthnicOrigin == "Black"))
thermophilus.sankey.black.30.hlac <- aggregate(x = as.numeric(thermophilus.sankey.black.30.hlac$Frequency),
                                          by = c(list(thermophilus.sankey.black.30.hlac$Allele, 
                                                      thermophilus.sankey.black.30.hlac$Peptide, 
                                                      thermophilus.sankey.black.30.hlac$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlac <- subset(thermophilus.sankey.black.30.hlac, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlac$x),
                            by = c(list(thermophilus.sankey.top.hlac$Group.1, 
                                        thermophilus.sankey.top.hlac$Group.2, 
                                        thermophilus.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(thermophilus.snk.agg$Group.1),
  target = c(thermophilus.snk.agg$Group.2),
  value = thermophilus.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Thermophilus_Black_30_10up.html", selfcontained = TRUE)

rm(sankey, thermophilus.sankey.top.hlac, thermophilus.snk.agg, links, source, target, value, nodes, name, rede)

##################### THERMOPHILUS 30 CAUCASOID

thermophilus.sankey.cauc.30.hlac <- unique(subset(thermophilus.sankey.30.hlac, EthnicOrigin == "Caucasoid"))
thermophilus.sankey.cauc.30.hlac <- aggregate(x = as.numeric(thermophilus.sankey.cauc.30.hlac$Frequency),
                                         by = c(list(thermophilus.sankey.cauc.30.hlac$Allele, 
                                                     thermophilus.sankey.cauc.30.hlac$Peptide, 
                                                     thermophilus.sankey.cauc.30.hlac$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlac <- subset(thermophilus.sankey.cauc.30.hlac, x >= 0.10)
# exclude lines wth e-4...n
# thermophilus.sankey.10 <- thermophilus.sankey.top[- grep("e-", thermophilus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlac$x),
                            by = c(list(thermophilus.sankey.top.hlac$Group.1, 
                                        thermophilus.sankey.top.hlac$Group.2, 
                                        thermophilus.sankey.top.hlac$Group.3)), FUN = sum)
links <-data.frame(
  source = c(thermophilus.snk.agg$Group.1),
  target = c(thermophilus.snk.agg$Group.2),
  value = thermophilus.snk.agg$x)
nodes <- data.frame(
  name=c(as.character(links$source),
         as.character(links$target))
  %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
rede <- sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight = TRUE, nodeWidth = 30)
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-C_Thermophilus_Caucasoid_30_10up.html", selfcontained = TRUE)

rm(sankey,thermophilus.sankey.top.hlac, thermophilus.snk.agg, links, source, target, value, nodes, name, rede)


### RASCUNHO SE PRECISAR

# common.ther.asbl.30 <- data.frame(merge(thermophilus.sankey.oriental.30, thermophilus.sankey.black.30, 
#                                         by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))
# common.ther.blca.30 <- data.frame(merge(thermophilus.sankey.black.30, thermophilus.sankey.cauc.30, 
#                                         by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))
# common.ther.asca.30 <- data.frame(merge(thermophilus.sankey.oriental.30, thermophilus.sankey.cauc.30, 
#                                         by.x = c("Group.1", "Group.2"), by.y = c("Group.1", "Group.2")))
# 
# common.ther.asbl.30 <- common.ther.asbl.30[,c(1,2)]
# common.ther.blca.30 <- common.ther.blca.30[,c(1,2)]
# common.ther.asca.30 <- common.ther.asca.30[,c(1,2)]
# 
# combined.30 <- rbind(common.ther.asbl.30, common.ther.blca.30, common.ther.asca.30)
# duplicate_rows.30 <- unique(combined.30[duplicated(combined.30), ])
# combined_peptides.30 <- unique(duplicate_rows.30$Group.2)

##write.table(combined_peptides.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Thermophilus_Peptides_3pop_HLA-C.30.csv")
##write.table(duplicate_rows.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Thermophilus_Common_3pop_HLA-C.30.csv")

###############################################################################
###############################################################################
###############################################################################

# Comparações HLA-C

# Geral IC50 < 500nM
unq.aureus.pep.hlac <- unique(aureus.IC50.hlac$peptide)
unq.pyogenes.pep.hlac <- unique(pyogenes.IC50.hlac$peptide)
unq.thermophilus.pep.hlac <- unique(thermophilus.IC50.hlac$peptide)

write.table(unq.aureus.pep.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Unique_PEP_IC50_500.csv")
write.table(unq.pyogenes.pep.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Unique_PEP_IC50_500.csv")
write.table(unq.thermophilus.pep.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Unique_PEP_IC50_500.csv")

# Geral IC50 < 30nM
unq.aureus.pep.30.hlac <- unique(aureus.IC50.30.hlac$peptide)
unq.pyogenes.pep.30.hlac <- unique(pyogenes.IC50.30.hlac$peptide)
unq.thermophilu.pep.30.hlac <- unique(thermophilus.IC50.30.hlac$peptide)

write.table(unq.aureus.pep.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Unique_PEP_IC50_30.csv")
write.table(unq.pyogenes.pep.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Unique_PEP_IC50_30.csv")
write.table(unq.thermophilus.pep.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Unique_PEP_IC50_30.csv")

# Por população IC50 < 500nM

## Sp. aureus

### Black

hlac.aureus.sankey.black.500.hlac <- unique(aureus.sankey.black.hlac$Group.1)
hlac.aureus.sankey.black.30.hlac <- unique(aureus.sankey.black.30.hlac$Group.1)
hlac.aureus.sankey.black.500.pep <- unique(aureus.sankey.black.hlac$Group.2)
hlac.aureus.sankey.black.30.pep <- unique(aureus.sankey.black.30.hlac$Group.2)

write.table(hlac.aureus.sankey.black.500.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Unique_HLA_Black_IC50_500.csv")
write.table(hlac.aureus.sankey.black.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Unique_HLA_Black_IC50_30.csv")
write.table(hlac.aureus.sankey.black.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Unique_PEP_Black_IC50_500.csv")
write.table(hlac.aureus.sankey.black.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Unique_PEP_Black_IC50_30.csv")

### Caucasoid

hlac.aureus.sankey.cauc.500.hlac <- unique(aureus.sankey.caucasoid.hlac$Group.1)
hlac.aureus.sankey.cauc.30.hlac <- unique(aureus.sankey.cauc.30.hlac$Group.1)
hlac.aureus.sankey.cauc.500.pep <- unique(aureus.sankey.caucasoid.hlac$Group.2)
hlac.aureus.sankey.cauc.30.pep <- unique(aureus.sankey.cauc.30.hlac$Group.2)

write.table(hlac.aureus.sankey.cauc.500.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Unique_HLA_Cauc_IC50_500.csv")
write.table(hlac.aureus.sankey.cauc.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Unique_HLA_Cauc_IC50_30.csv")
write.table(hlac.aureus.sankey.cauc.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Unique_PEP_Cauc_IC50_500.csv")
write.table(hlac.aureus.sankey.cauc.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Unique_PEP_Cauc_IC50_30.csv")

### Oriental

hlac.aureus.sankey.oriental.500.hlac <- unique(aureus.sankey.oriental.hlac$Group.1)
hlac.aureus.sankey.oriental.30.hlac <- unique(aureus.sankey.oriental.30.hlac$Group.1)
hlac.aureus.sankey.oriental.500.pep <- unique(aureus.sankey.oriental.hlac$Group.2)
hlac.aureus.sankey.oriental.30.pep <- unique(aureus.sankey.oriental.30.hlac$Group.2)

write.table(hlac.aureus.sankey.oriental.500.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Unique_HLA_Orient_IC50_500.csv")
write.table(hlac.aureus.sankey.oriental.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Unique_HLA_Orient_IC50_30.csv")
write.table(hlac.aureus.sankey.oriental.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Unique_PEP_Orient_IC50_500.csv")
write.table(hlac.aureus.sankey.oriental.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Unique_PEP_Orient_IC50_30.csv")

## Sp. pyogenes

### Black

hlac.pyogenes.sankey.black.500.hlac <- unique(pyogenes.sankey.black.hlac$Group.1)
hlac.pyogenes.sankey.black.30.hlac <- unique(pyogenes.sankey.black.30.hlac$Group.1)
hlac.pyogenes.sankey.black.500.pep <- unique(pyogenes.sankey.black.hlac$Group.2)
hlac.pyogenes.sankey.black.30.pep <- unique(pyogenes.sankey.black.30.hlac$Group.2)

write.table(hlac.pyogenes.sankey.black.500.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Unique_HLA_Black_IC50_500.csv")
write.table(hlac.pyogenes.sankey.black.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Unique_HLA_Black_IC50_30.csv")
write.table(hlac.pyogenes.sankey.black.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Unique_PEP_Black_IC50_500.csv")
write.table(hlac.pyogenes.sankey.black.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Unique_PEP_Black_IC50_30.csv")

### Caucasoid

hlac.pyogenes.sankey.cauc.500.hla <- unique(pyogenes.sankey.cauc.hlac$Group.1)
hlac.pyogenes.sankey.cauc.30.hla <- unique(pyogenes.sankey.cauc.30.hlac$Group.1)
hlac.pyogenes.sankey.cauc.500.pep <- unique(pyogenes.sankey.cauc.hlac$Group.2)
hlac.pyogenes.sankey.cauc.30.pep <- unique(pyogenes.sankey.cauc.30.hlac$Group.2)

write.table(hlac.pyogenes.sankey.cauc.500.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Unique_HLA_Cauc_IC50_500.csv")
write.table(hlac.pyogenes.sankey.cauc.30.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Unique_HLA_Cauc_IC50_30.csv")
write.table(hlac.pyogenes.sankey.cauc.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Unique_PEP_Cauc_IC50_500.csv")
write.table(hlac.pyogenes.sankey.cauc.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Unique_PEP_Cauc_IC50_30.csv")

### Oriental

hlac.pyogenes.sankey.oriental.500.hla <- unique(pyogenes.sankey.oriental.hlac$Group.1)
hlac.pyogenes.sankey.oriental.30.hla <- unique(pyogenes.sankey.oriental.30.hlac$Group.1)

hlac.pyogenes.sankey.oriental.500.pep <- unique(pyogenes.sankey.oriental.hlac$Group.2)
hlac.pyogenes.sankey.oriental.30.pep <- unique(pyogenes.sankey.oriental.30.hlac$Group.2)

write.table(hlac.pyogenes.sankey.oriental.500.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Unique_HLA_Orient_IC50_500.csv")
write.table(hlac.pyogenes.sankey.oriental.30.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Unique_HLA_Orient_IC50_30.csv")
write.table(hlac.pyogenes.sankey.oriental.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Unique_PEP_Orient_IC50_500.csv")
write.table(hlac.pyogenes.sankey.oriental.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Unique_PEP_Orient_IC50_30.csv")

## Sp. thermophilus
### Black
hlac.thermophilus.sankey.black.500.hlac <- unique(thermophilus.sankey.black.hlac$Group.1)
hlac.thermophilus.sankey.black.30.hlac <- unique(thermophilus.sankey.black.30.hlac$Group.1)

hlac.thermophilus.sankey.black.500.pep <- unique(thermophilus.sankey.black.hlac$Group.2)
hlac.thermophilus.sankey.black.30.pep <- unique(thermophilus.sankey.black.30.hlac$Group.2)

write.table(hlac.thermophilus.sankey.black.500.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Unique_HLA_Black_IC50_500.csv")
write.table(hlac.thermophilus.sankey.black.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Unique_HLA_Black_IC50_30.csv")
write.table(hlac.thermophilus.sankey.black.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Unique_PEP_Black_IC50_500.csv")
write.table(hlac.thermophilus.sankey.black.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Unique_PEP_Black_IC50_30.csv")

### Caucasoid
hlac.thermophilus.sankey.cauc.500.hlac <- unique(thermophilus.sankey.cauc.hlac$Group.1)
hlac.thermophilus.sankey.cauc.30.hlac <- unique(thermophilus.sankey.cauc.30.hlac$Group.1)

hlac.thermophilus.sankey.cauc.500.pep <- unique(thermophilus.sankey.cauc.hlac$Group.2)
hlac.thermophilus.sankey.cauc.30.pep <- unique(thermophilus.sankey.cauc.30.hlac$Group.2)

write.table(hlac.thermophilus.sankey.cauc.500.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Unique_HLA_Cauc_IC50_500.csv")
write.table(hlac.thermophilus.sankey.cauc.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Unique_HLA_Cauc_IC50_30.csv")
write.table(hlac.thermophilus.sankey.cauc.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Unique_PEP_Cauc_IC50_500.csv")
write.table(hlac.thermophilus.sankey.cauc.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Unique_PEP_Cauc_IC50_30.csv")

### Oriental
hlac.thermophilus.sankey.oriental.500.hla <- unique(thermophilus.sankey.oriental.hlac$Group.1)
hlac.thermophilus.sankey.oriental.30.hla <- unique(thermophilus.sankey.oriental.30.hlac$Group.1)

hlac.thermophilus.sankey.oriental.500.pep <- unique(thermophilus.sankey.oriental.hlac$Group.2)
hlac.thermophilus.sankey.oriental.30.pep <- unique(thermophilus.sankey.oriental.30.hlac$Group.2)

write.table(hlac.thermophilus.sankey.oriental.500.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Unique_HLA_Orient_IC50_500.csv")
write.table(hlac.thermophilus.sankey.oriental.30.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Unique_HLA_Orient_IC50_30.csv")
write.table(hlac.thermophilus.sankey.oriental.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Unique_PEP_Orient_IC50_500.csv")
write.table(hlac.pyogenes.sankey.oriental.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Unique_PEP_Orient_IC50_30.csv")

###############################################################################
###############################################################################
###############################################################################

## For HLA-Crena 500
aureus.arena.hlac <- aureus.ethnicity.unq.hlac[,c(1,2)]
aureus.arena.unq.hlac <- unique(aureus.arena.hlac)
aureus.arena.hla.unq.hlac <- unique(aureus.arena.hlac$Allele)
aureus.arena.pep.unq.hlac <- unique(aureus.arena.hlac$Peptide)

write.table(aureus.arena.unq.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Arena_Unq_HLA-C.csv")
write.table(aureus.arena.hla.unq.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Arena_HlaUnq_HLA-C.csv")
write.table(aureus.arena.pep.unq.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Arena_PepUnq_HLA-C.csv")

pyogenes.arena.hlac <- pyogenes.ethnicity.unq.hlac[,c(1,2)]
pyogenes.arena.unq.hlac <- unique(pyogenes.arena.hlac)
pyogenes.arena.hla.unq.hlac <- unique(pyogenes.arena.hlac$Allele)
pyogenes.arena.pep.unq.hlac <- unique(pyogenes.arena.hlac$Peptide)

write.table(pyogenes.arena.unq.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Arena_Unq_HLA-C.csv")
write.table(pyogenes.arena.hla.unq.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Arena_HlaUnq_HLA-C.csv")
write.table(pyogenes.arena.pep.unq.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Arena_PepUnq_HLA-C.csv")

thermophilus.arena.hlac <- thermophilus.ethnicity.unq.hlac[,c(1,2)]
thermophilus.arena.unq.hlac <- unique(thermophilus.arena.hlac)
thermophilus.arena.hla.unq.hlac <- unique(thermophilus.arena.hlac$Allele)
thermophilus.arena.pep.unq.hlac <- unique(thermophilus.arena.hlac$Peptide)

write.table(thermophilus.arena.unq.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Arena_Unq_HLA-C.csv")
write.table(thermophilus.arena.hla.unq.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Arena_HlaUnq_HLA-C.csv")
write.table(thermophilus.arena.pep.unq.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Arena_PepUnq_HLA-C.csv")

## For HLA-Crena 30
aureus.arena.30.hlac <- aureus.ethnicity.unq.30.hlac[,c(1,2)]
aureus.arena.unq.30.hlac <- unique(aureus.arena.30.hlac)
aureus.arena.hla.unq.30.hlac <- unique(aureus.arena.30.hlac$Allele)
aureus.arena.pep.unq.30.hlac <- unique(aureus.arena.30.hlac$Peptide)

write.table(aureus.arena.unq.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Arena_Unq_30_HLA-C.csv")
write.table(aureus.arena.hla.unq.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Arena_HlaUnq_30_HLA-C.csv")
write.table(aureus.arena.pep.unq.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Aureus_Arena_PepUnq_30_HLA-C.csv")

pyogenes.arena.30.hlac <- pyogenes.ethnicity.unq.30.hlac[,c(1,2)]
pyogenes.arena.unq.30.hlac <- unique(pyogenes.arena.30.hlac)
pyogenes.arena.hla.unq.30.hlac <- unique(pyogenes.arena.30.hlac$Allele)
pyogenes.arena.pep.unq.30.hlac <- unique(pyogenes.arena.30.hlac$Peptide)

write.table(pyogenes.arena.unq.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Arena_Unq_30_HLA-C.csv")
write.table(pyogenes.arena.hla.unq.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Arena_HlaUnq_30_HLA-C.csv")
write.table(pyogenes.arena.pep.unq.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Pyogenes_Arena_PepUnq_30_HLA-C.csv")

thermophilus.arena.30.hlac <- thermophilus.ethnicity.unq.30.hlac[,c(1,2)]
thermophilus.arena.unq.30.hlac <- unique(thermophilus.arena.30.hlac)
thermophilus.arena.hla.unq.30.hlac <- unique(thermophilus.arena.30.hlac$Allele)
thermophilus.arena.pep.unq.30.hlac <- unique(thermophilus.arena.30.hlac$Peptide)

write.table(thermophilus.arena.unq.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Arena_Unq_30_HLA-C.csv")
write.table(thermophilus.arena.hla.unq.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Arena_HlaUnq_30_HLA-C.csv")
write.table(thermophilus.arena.pep.unq.30.hlac, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-C_Thermophilus_Arena_PepUnq_30_HLA-C.csv")

