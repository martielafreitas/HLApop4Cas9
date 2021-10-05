###############################################################################
# VIA DE PROCESSAMENTO
###############################################################################

# Unindo os resulatdos dos preditores para as três proteínas 
# em uma única tabela que adiciona o nome que está no nome do
# arquivo, na última coluna da tabela.

setwd("/Users/martielafreitas/Desktop/HLA-Rproject")

# variavel que contem o caminho para os arquivos
file_list <- list.files(path="/Users/martielafreitas/Desktop/HLA-Rproject/AInputs/Preditores/1.Processing")

# inicia uma dataframe vazia
process.hlaa <- data.frame()

# onde serao processados os arquivos
setwd("/Users/martielafreitas/Desktop/HLA-Rproject/AInputs/Preditores/1.Processing")

# laco que une todos os arquivos
for(i in 1:length(file_list)){
  #each file will be read in, specify which columns you need read in to avoid any errors
  temp_data <- read.csv(file_list[i], sep=",")
  #clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
  temp_data$Class <- sapply(strsplit(gsub(".csv", "", file_list[i]), " - "), function(x){x[3]})
  #for each iteration, bind the new data to the building dataset
  process.hlaa <- rbind(process.hlaa, temp_data)
}

# salvando o resultado em um arquivo .csv
write.csv(process.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Processing_Geral.csv")

# Separando cada resultado por organismo de onde vem a proteina
aureus.process.hlaa<- subset(process.hlaa, Class == "Aureus")
pyogenes.process.hlaa <- subset(process.hlaa, Class == "Pyogenes")
thermophilus.process.hlaa <- subset(process.hlaa, Class == "Thermophilus")

# salvando o resultado em um arquivo .csv
write.csv(aureus.process.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Process.csv")
write.csv(pyogenes.process.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Process.csv")
write.csv(thermophilus.process.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Process.csv")


###############################################################################
# LIGAÇÃO
###############################################################################

##-- Binding Results --##

# variavel que contem o caminho para os arquivos
file_list <- list.files(path="/Users/martielafreitas/Desktop/HLA-Rproject/AInputs/Preditores/2.Binding")

#inicia uma dataframe vazia
binding.hlaa <- data.frame()

setwd("/Users/martielafreitas/Desktop/HLA-Rproject/AInputs/Preditores/2.Binding")
for(i in 1:length(file_list)){
  #each file will be read in, specify which columns you need read in to avoid any errors
  temp_data <- read.csv(file_list[i], sep=";")
  #clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
  temp_data$Class <- sapply(strsplit(gsub(".csv", "", file_list[i]), " - "), function(x){x[3]})
  #for each iteration, bind the new data to the building dataset
  binding.hlaa <- rbind(binding.hlaa, temp_data)
}
# salvando resultados gerais de ligacao em um arquivo .csv
write.csv(binding.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Binding_Geral.csv")

# separando os resultados por organismo de onde veio cada proteina
aureus.binding.hlaa <- subset(binding.hlaa, Class == "Aureus")
pyogenes.binding.hlaa <- subset(binding.hlaa, Class == "Pyogenes")
thermophilus.binding.hlaa <- subset(binding.hlaa, Class == "Thermophilus")

# salvando resultados gerais por organismo em um arquivo .csv
write.csv(aureus.binding.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Binding.csv")
write.csv(pyogenes.binding.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Binding.csv")
write.csv(thermophilus.binding.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Binding.csv")



###############################################################################
# IMMUNOGENICIDADE
###############################################################################

##-- Comparing Results for immunogenicity predictions --##

aureus.procbind.hlaa <- merge(aureus.binding.hlaa, 
                              aureus.process.hlaa, 
                              by.x=c("allele", "seq_num", "start", "end", "length", "peptide", "Class"),
                              by.y=c("Allele", "X.", "Start", "End", "Peptide.Length", "Peptide", "Class"))

pyogenes.procbind.hlaa <- merge(pyogenes.binding.hlaa, 
                                pyogenes.process.hlaa, 
                                by.x=c("allele", "seq_num", "start", "end", "length", "peptide", "Class"),
                                by.y=c("Allele", "X.", "Start", "End", "Peptide.Length", "Peptide", "Class"))

thermophilus.procbind.hlaa <- merge(thermophilus.binding.hlaa, 
                                    thermophilus.process.hlaa, 
                                    by.x=c("allele", "seq_num", "start", "end", "length", "peptide", "Class"),
                                    by.y=c("Allele", "X.", "Start", "End", "Peptide.Length", "Peptide", "Class"))

write.csv(aureus.procbind.hlaa, file="/Users/martielafreitas/Documents/Rprojects/HLA-A/Results/HLA-A_Aureus_Merged_Process_and_Binding.csv")
write.csv(pyogenes.procbind.hlaa, file="/Users/martielafreitas/Documents/Rprojects/HLA-A/Results/HLA-A_Pyogenes_Merged_Process_and_Binding.csv")
write.csv(thermophilus.procbind.hlaa, file="/Users/martielafreitas/Documents/Rprojects/HLA-A/Results/HLA-A_Thermophilus_Merged_Process_and_Binding.csv")



###############################################################################
# ANALISE DE RESULTADOS
###############################################################################

## Immunogenicity Results

## Obtaining list of peptides

## Peptídeos totais - independente de população ou frequência, por organismo
aureus.peptides.hlaa <- unique(aureus.procbind.hlaa$peptide)
pyogenes.peptides.hlaa <- unique(pyogenes.procbind.hlaa$peptide)
thermophilus.peptides.hlaa <- unique(thermophilus.procbind.hlaa$peptide)

write.csv(aureus.peptides.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_TotalPeptides.csv")
write.csv(pyogenes.peptides.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_TotalPeptides.csv")
write.csv(thermophilus.peptides.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_TotalPeptides.csv")

# Esses peptídeos passaram pela ferramenta de predição de imunogenicidade.

#- Calling result files -#
file_list <- list.files(path="/Users/martielafreitas/Desktop/HLA-Rproject/AInputs/Preditores/3.Immunogenicity")
immunogenicity.hlaa <- data.frame()

setwd("/Users/martielafreitas/Desktop/HLA-Rproject/AInputs/Preditores/3.Immunogenicity")
for(i in 1:length(file_list)){
  #each file will be read in, specify which columns you need read in to avoid any errors
  temp_data <- read.csv(file_list[i], sep=",")
  #clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
  temp_data$Class <- sapply(strsplit(gsub(".csv", "", file_list[i]), " - "), function(x){x[3]})
  #for each iteration, bind the new data to the building dataset
  immunogenicity.hlaa <- rbind(immunogenicity.hlaa, temp_data)
}
# salvando resultados gerais de imunogenicidade em um arquivo .csv
write.csv(immunogenicity.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Immunogenicity_Geral.csv")


# Separando a lista de HLAs totais por organismo - independente de população ou frequência

aureus.immuno.hlaa <- subset(immunogenicity.hlaa, Class == "Aureus")
pyogenes.immuno.hlaa <- subset(immunogenicity.hlaa, Class == "Pyogenes")
thermophilus.immuno.hlaa <- subset(immunogenicity.hlaa, Class == "Thermophilus")

write.csv(aureus.immuno.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Immunogenicity.csv")
write.csv(pyogenes.immuno.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Immunogenicity.csv")
write.csv(thermophilus.immuno.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Immunogenicity.csv")


###############################################################################
# Filtering final results
###############################################################################

#- Immuno: score >= 0.0: immunogenic -#

setwd("/Users/martielafreitas/Desktop/HLA-Rproject/Results/")

aureus.immuno.positive.hlaa <- subset(aureus.immuno.hlaa, score >= 0.0)
pyogenes.immuno.positive.hlaa <- subset(pyogenes.immuno.hlaa, score >= 0.0)
thermophilus.immuno.positive.hlaa <- subset(thermophilus.immuno.hlaa, score >= 0.0)

#- Merging tables to IC50 -#
colnames(aureus.procbind.hlaa)
colnames(aureus.immuno.positive.hlaa)

aureus.full.path.hlaa <- merge(aureus.procbind.hlaa, aureus.immuno.positive.hlaa, by.x = "peptide", by.y = "peptide")
pyogenes.full.path.hlaa <- merge(pyogenes.procbind.hlaa, pyogenes.immuno.positive.hlaa, by.x = "peptide", by.y = "peptide")
thermophilus.full.path.hlaa <- merge(thermophilus.procbind.hlaa, thermophilus.immuno.positive.hlaa, by.x = "peptide", by.y = "peptide")

#- IC50 < 500nM  -#
aureus.IC50.hlaa <- aureus.full.path.hlaa[(aureus.full.path.hlaa$MHC.IC50 <= 500.0),]
pyogenes.IC50.hlaa <- pyogenes.full.path.hlaa[(pyogenes.full.path.hlaa$MHC.IC50 <= 500.0),]
thermophilus.IC50.hlaa <- thermophilus.full.path.hlaa[(thermophilus.full.path.hlaa$MHC.IC50 <= 500.0),]

#- Final Table -#
aureus.results.hlaa <- aureus.IC50.hlaa[,c(1,2,9,10,11,12,13,14,15,8,17)]
pyogenes.results.hlaa <- pyogenes.IC50.hlaa[,c(1,2,9,10,11,12,13,14,15,8,17)]
thermophilus.results.hlaa <- thermophilus.IC50.hlaa[,c(1,2,9,10,11,12,13,14,15,8,17)]

# corrigindo cabeçalho antes de salvar! :)
colnames(aureus.results.hlaa)
cols <- c("Peptide", "Allele", "Percentile.Rank", "Proteasome.Score", "TAP.Score",
          "MHC.Score", "Processing.Score", "Total.Score", "MHC.IC50", "Binding.Score", "Immunogenicity.Score")
colnames(aureus.results.hlaa) <- cols
colnames(pyogenes.results.hlaa) <- cols
colnames(thermophilus.results.hlaa) <- cols

write.csv(aureus.results.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_ToltalResults.csv")
write.csv(pyogenes.results.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_TotalResults.csv")
write.csv(thermophilus.results.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_TotalResults.csv")

#- Obtaining list of alleles por organismo-#
aureus.final.allele.hlaa <- data.frame(unique(aureus.results.hlaa$Allele)) #221
pyogenes.final.allele.hlaa <- data.frame(unique(pyogenes.results.hlaa$Allele)) #225
thermophilus.final.allele.hlaa <- data.frame(unique(thermophilus.results.hlaa$Allele)) #200

write.csv(aureus.final.allele.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_TotalHLAResults.csv")
write.csv(pyogenes.final.allele.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_TotalHLAResults.csv")
write.csv(thermophilus.final.allele.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_TotalHLAResults.csv")

#- Obtaining list of peptides por organismo -#
aureus.final.peptides.hlaa <- data.frame(unique(aureus.results.hlaa$Peptide)) #121
pyogenes.final.peptides.hlaa <- data.frame(unique(pyogenes.results.hlaa$Peptide)) #163
thermophilus.final.peptides.hlaa <- data.frame(unique(thermophilus.results.hlaa$Peptide)) #189

write.csv(aureus.final.peptides.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_TotalPEPResults.csv")
write.csv(pyogenes.final.peptides.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_TotalPEPResults.csv")
write.csv(thermophilus.final.peptides.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_TotalPEPResults.csv")


common.aurpyo.hlaa <- data.frame(merge(aureus.final.peptides.hlaa, pyogenes.final.peptides.hlaa, 
                                  by.x = colnames(aureus.final.peptides.hlaa), by.y = colnames(pyogenes.final.peptides.hlaa)))
common.aurthe.hlaa <- data.frame(merge(aureus.final.peptides.hlaa, thermophilus.final.peptides.hlaa, 
                                  by.x = colnames(aureus.final.peptides.hlaa), by.y = colnames(thermophilus.final.peptides.hlaa)))
common.pyothe.hlaa <- data.frame(merge(pyogenes.final.peptides.hlaa, thermophilus.final.peptides.hlaa, 
                                  by.x = colnames(pyogenes.final.peptides.hlaa), by.y = colnames(thermophilus.final.peptides.hlaa)))

# Qual o score de imunogenicidade pra os comuns pyo/the? 16 no total!
score.common.pyothe.hlaa <- merge(common.pyothe.hlaa, immunogenicity.hlaa, by.x="unique.pyogenes.results.hlaa.Peptide.", by.y = "peptide")
write.csv(score.common.pyothe.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_CommonPeptides_Pyogenes_Thermophilus.csv")


###############################################################################
# Logos
###############################################################################

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Aureus_PeptidesLogo.pdf")
ggplot() + geom_logo( aureus.final.peptides.hlaa ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Pyogenes_PeptidesLogo.pdf")
ggplot() + geom_logo( pyogenes.final.peptides.hlaa ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Thermophilus_PeptidesLogo.pdf")
ggplot() + geom_logo( thermophilus.final.peptides.hlaa ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_PyoThe_PeptidesLogo.pdf")
ggplot() + geom_logo( common.pyothe.hlaa ) + theme_logo()
dev.off()


###############################################################################
# Populations 
###############################################################################
setwd("/Users/martielafreitas/Desktop/HLA-Rproject/AInputs/Populacoes")
populacoes.hlaa <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/AInputs/Populacoes/HLA-A_GOLD_Geral", sep=";")

# Merged tables
aureus.pop.hlaa <-  merge(aureus.results.hlaa, populacoes.hlaa, by.x = "Allele", by.y = "Allele")
pyogenes.pop.hlaa <-  merge(pyogenes.results.hlaa, populacoes.hlaa, by.x = "Allele", by.y = "Allele")
thermophilus.pop.hlaa <-  merge(thermophilus.results.hlaa, populacoes.hlaa, by.x = "Allele", by.y = "Allele")

write.table(aureus.pop.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Path_Populacao.csv")
write.table(pyogenes.pop.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Path_Populacao.csv")
write.table(thermophilus.pop.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Path_Populacao.csv")

# Unique peptides by Population and allele

aureus.populacao.hlaa <- aureus.pop.hlaa[,c(1,2,13,14,15,16)]
pyogenes.populacao.hlaa <- pyogenes.pop.hlaa[,c(1,2,13,14,15,16)]
thermophilus.populacao.hlaa <- thermophilus.pop.hlaa[,c(1,2,13,14,15,16)]

aureus.ethnicity.hlaa <- aureus.populacao.hlaa[,c(1,2,6)]
pyogenes.ethnicity.hlaa <- pyogenes.populacao.hlaa[,c(1,2,6)]
thermophilus.ethnicity.hlaa <- thermophilus.populacao.hlaa[,c(1,2,6)]

aureus.ethnicity.unq.hlaa <- unique(aureus.ethnicity.hlaa)
pyogenes.ethnicity.unq.hlaa <- unique(pyogenes.ethnicity.hlaa)
thermophilus.ethnicity.unq.hlaa <-unique(thermophilus.ethnicity.hlaa)

write.table(aureus.ethnicity.unq.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Ethnicity_Unique.csv")
write.table(pyogenes.ethnicity.unq.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Ethnicity_Unique.csv")
write.table(thermophilus.ethnicity.unq.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Ethnicity_Unique.csv")

### PEPTIDES

aureus.sankey.hlaa <- aureus.populacao.hlaa[,c(1,2,6,4)]
pyogenes.sankey.hlaa <- pyogenes.populacao.hlaa[,c(1,2,6,4)]
thermophilus.sankey.hlaa <- thermophilus.populacao.hlaa[,c(1,2,6,4)]

write.table(aureus.sankey.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Sankey.csv")
write.table(pyogenes.sankey.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Sankey.csv")
write.table(thermophilus.sankey.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Sankey.csv")

#aureus.sankey <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Sankey_HLA-A.csv", sep = ";")
#pyogenes.sankey <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Sankey_HLA-A.csv", sep = ";")
#thermophilus.sankey <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Sankey_HLA-A.csv", sep = ";")

############################################################################### AUREUS

##################### AUREUS 500 ORIENTAL

aureus.sankey.oriental.hlaa <- unique(subset(aureus.sankey.hlaa, EthnicOrigin == "Oriental"))
aureus.sankey.oriental.hlaa <- aggregate(x = as.numeric(aureus.sankey.oriental.hlaa$Frequency),
                                    by = c(list(aureus.sankey.oriental.hlaa$Allele, 
                                                aureus.sankey.oriental.hlaa$Peptide, 
                                                aureus.sankey.oriental.hlaa$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlaa <- subset(aureus.sankey.oriental.hlaa, x >= 0.10)
# exclude lines wth e-4...n
#aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlaa$x),
                            by = c(list(aureus.sankey.top.hlaa$Group.1, 
                                        aureus.sankey.top.hlaa$Group.2, 
                                        aureus.sankey.top.hlaa$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Aureus_Oriental_500_10up.html", selfcontained = TRUE)

rm(sankey, aureus.sankey.top, aureus.snk.agg, links, source, target, value, nodes, name, rede)

##################### AUREUS 500 BLACK

aureus.sankey.black.hlaa <- unique(subset(aureus.sankey.hlaa, EthnicOrigin == "Black"))
aureus.sankey.black.hlaa <- aggregate(x = as.numeric(aureus.sankey.black.hlaa$Frequency),
                                 by = c(list(aureus.sankey.black.hlaa$Allele, 
                                             aureus.sankey.black.hlaa$Peptide, 
                                             aureus.sankey.black.hlaa$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlaa <- subset(aureus.sankey.black.hlaa, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlaa$x),
                            by = c(list(aureus.sankey.top.hlaa$Group.1, 
                                        aureus.sankey.top.hlaa$Group.2, 
                                        aureus.sankey.top.hlaa$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Aureus_Black_500_10up.html", selfcontained = TRUE)

rm(sankey, aureus.sankey.top, aureus.snk.agg, links, source, target, value, nodes, name, rede)

##################### AUREUS 500 CAUCASOID

aureus.sankey.caucasoid.hlaa <- unique(subset(aureus.sankey.hlaa, EthnicOrigin == "Caucasoid"))
aureus.sankey.caucasoid.hlaa <- aggregate(x = as.numeric(aureus.sankey.caucasoid.hlaa$Frequency),
                                     by = c(list(aureus.sankey.caucasoid.hlaa$Allele, 
                                                 aureus.sankey.caucasoid.hlaa$Peptide, 
                                                 aureus.sankey.caucasoid.hlaa$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlaa <- subset(aureus.sankey.caucasoid.hlaa, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlaa$x),
                            by = c(list(aureus.sankey.top.hlaa$Group.1, 
                                        aureus.sankey.top.hlaa$Group.2, 
                                        aureus.sankey.top.hlaa$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Aureus_Caucasoid_500_10up.html", selfcontained = TRUE)

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

#write.table(combined_peptides, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Aureus_Peptides_3pop_HLA-A.csv")
##write.table(duplicate_rows, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Aureus_Common_3pop_HLA-A.csv")



################################################################################ PYOGENES

##################### PYOGENES 500 ORIENTAL

pyogenes.sankey.oriental.hlaa <- unique(subset(pyogenes.sankey.hlaa, EthnicOrigin == "Oriental"))
pyogenes.sankey.oriental.hlaa <- aggregate(x = as.numeric(pyogenes.sankey.oriental.hlaa$Frequency),
                                      by = c(list(pyogenes.sankey.oriental.hlaa$Allele, 
                                                  pyogenes.sankey.oriental.hlaa$Peptide, 
                                                  pyogenes.sankey.oriental.hlaa$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlaa <- subset(pyogenes.sankey.oriental.hlaa, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlaa$x),
                            by = c(list(pyogenes.sankey.top.hlaa$Group.1, 
                                        pyogenes.sankey.top.hlaa$Group.2, 
                                        pyogenes.sankey.top.hlaa$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Pyogenes_Oriental_500_10up.html", selfcontained = TRUE)

rm(sankey, pyogenes.sankey.top, pyogenes.snk.agg, links, source, target, value, nodes, name, rede)


##################### PYOGENES 500 BLACK
pyogenes.sankey.black.hlaa <- unique(subset(pyogenes.sankey.hlaa, EthnicOrigin == "Black"))
pyogenes.sankey.black.hlaa <- aggregate(x = as.numeric(pyogenes.sankey.black.hlaa$Frequency),
                                   by = c(list(pyogenes.sankey.black.hlaa$Allele, 
                                               pyogenes.sankey.black.hlaa$Peptide, 
                                               pyogenes.sankey.black.hlaa$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlaa <- subset(pyogenes.sankey.black.hlaa, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlaa$x),
                              by = c(list(pyogenes.sankey.top.hlaa$Group.1, 
                                          pyogenes.sankey.top.hlaa$Group.2, 
                                          pyogenes.sankey.top.hlaa$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Pyogenes_Black_500_10up.html", selfcontained = TRUE)

rm(sankey, pyogenes.sankey.top, pyogenes.snk.agg, links, source, target, value, nodes, name, rede)


##################### PYOGENES 500 CAUCASOID
pyogenes.sankey.cauc.hlaa <- unique(subset(pyogenes.sankey.hlaa, EthnicOrigin == "Caucasoid"))
pyogenes.sankey.cauc.hlaa <- aggregate(x = as.numeric(pyogenes.sankey.cauc.hlaa$Frequency),
                                  by = c(list(pyogenes.sankey.cauc.hlaa$Allele, 
                                              pyogenes.sankey.cauc.hlaa$Peptide, 
                                              pyogenes.sankey.cauc.hlaa$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlaa <- subset(pyogenes.sankey.oriental.hlaa, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlaa$x),
                            by = c(list(pyogenes.sankey.top.hlaa$Group.1, 
                                        pyogenes.sankey.top.hlaa$Group.2, 
                                        pyogenes.sankey.top.hlaa$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Pyogenes_Caucasoid_500_10up.html", selfcontained = TRUE)

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

##write.table(combined_peptides, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Pyogenes_Peptides_3pop_HLA-A.csv")
##write.table(duplicate_rows, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Pyogenes_Common_3pop_HLA-A.csv")


############################################################################### THERMOPHILUS

##################### THERMOPHILUS 500 ORIENTAL

thermophilus.sankey.oriental.hlaa <- unique(subset(thermophilus.sankey.hlaa, EthnicOrigin == "Oriental"))
thermophilus.sankey.oriental.hlaa <- aggregate(x = as.numeric(thermophilus.sankey.oriental.hlaa$Frequency),
                                          by = c(list(thermophilus.sankey.oriental.hlaa$Allele, 
                                                      thermophilus.sankey.oriental.hlaa$Peptide, 
                                                      thermophilus.sankey.oriental.hlaa$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlaa <- subset(thermophilus.sankey.oriental.hlaa, x >= 0.10)
# exclude lines wth e-4...n
#aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlaa$x),
                            by = c(list(thermophilus.sankey.top.hlaa$Group.1, 
                                        thermophilus.sankey.top.hlaa$Group.2, 
                                        thermophilus.sankey.top.hlaa$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Thermophilus_Oriental_500_10up.html", selfcontained = TRUE)

rm(sankey, thermophilus.sankey.top, thermophilus.snk.agg, links, source, target, value, nodes, name, rede)

##################### THERMOPHILUS 500 BLACK

thermophilus.sankey.black.hlaa <- unique(subset(thermophilus.sankey.hlaa, EthnicOrigin == "Black"))
thermophilus.sankey.black.hlaa <- aggregate(x = as.numeric(thermophilus.sankey.black.hlaa$Frequency),
                                       by = c(list(thermophilus.sankey.black.hlaa$Allele, 
                                                   thermophilus.sankey.black.hlaa$Peptide, 
                                                   thermophilus.sankey.black.hlaa$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlaa <- subset(thermophilus.sankey.black.hlaa, x >= 0.10)
# exclude lines wth e-4...n
# thermophilus.sankey.10 <- thermophilus.sankey.top[- grep("e-", thermophilus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlaa$x),
                            by = c(list(thermophilus.sankey.top.hlaa$Group.1, 
                                        thermophilus.sankey.top.hlaa$Group.2, 
                                        thermophilus.sankey.top.hlaa$Group.3)), FUN = sum)
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

thermophilus.sankey.cauc.hlaa <- unique(subset(thermophilus.sankey.hlaa, EthnicOrigin == "Caucasoid"))
thermophilus.sankey.cauc.hlaa <- aggregate(x = as.numeric(thermophilus.sankey.cauc.hlaa$Frequency),
                                      by = c(list(thermophilus.sankey.cauc.hlaa$Allele, 
                                                  thermophilus.sankey.cauc.hlaa$Peptide, 
                                                  thermophilus.sankey.cauc.hlaa$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlaa <- subset(thermophilus.sankey.cauc.hlaa, x >= 0.10)
# exclude lines wth e-4...n
# thermophilus.sankey.10 <- thermophilus.sankey.top[- grep("e-", thermophilus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlaa$x),
                            by = c(list(thermophilus.sankey.top.hlaa$Group.1, 
                                        thermophilus.sankey.top.hlaa$Group.2, 
                                        thermophilus.sankey.top.hlaa$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Thermophilus_Caucasoid_500_10up.html", selfcontained = TRUE)

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

##write.table(combined_peptides, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Thermophilus_Peptides_3pop_HLA-A.csv")
##write.table(duplicate_rows, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Thermophilus_Common_3pop_HLA-A.csv")


###############################################################################


### ANALISE PARA IC50 < 30

aureus.IC50.30.hlaa <- aureus.full.path.hlaa[(aureus.full.path.hlaa$MHC.IC50 <= 30.0),]
pyogenes.IC50.30.hlaa <- pyogenes.full.path.hlaa[(pyogenes.full.path.hlaa$MHC.IC50 <= 30.0),]
thermophilus.IC50.30.hlaa <- thermophilus.full.path.hlaa[(thermophilus.full.path.hlaa$MHC.IC50 <= 30.0),]

#- Final Table -#

aureus.results.30.hlaa <- aureus.IC50.30.hlaa[,c(1,2,9,10,11,12,13,14,15,8,17)]
pyogenes.results.30.hlaa <- pyogenes.IC50.30.hlaa[,c(1,2,9,10,11,12,13,14,15,8,17)]
thermophilus.results.30.hlaa <- thermophilus.IC50.30.hlaa[,c(1,2,9,10,11,12,13,14,15,8,17)]

# corrigindo cabeçalho antes de salvar! :)
colnames(aureus.results.30.hlaa)
cols <- c("Peptide", "Allele", "Percentile.Rank", "Proteasome.Score", "TAP.Score",
          "MHC.Score", "Processing.Score", "Total.Score", "MHC.IC50", "Binding.Score", "Immunogenicity.Score")
colnames(aureus.results.30.hlaa) <- cols
colnames(pyogenes.results.30.hlaa) <- cols
colnames(thermophilus.results.30.hlaa) <- cols

write.csv(aureus.results.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Results_IC50_30.csv")
write.csv(pyogenes.results.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Results_IC50_30.csv")
write.csv(thermophilus.results.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Results_IC50_30.csv")

##-- Summary --##

#- Obtaining list of alleles-#

aureus.final.allele.30.hlaa <- data.frame(unique(aureus.results.30.hlaa$Allele)) #221
pyogenes.final.allele.30.hlaa <- data.frame(unique(pyogenes.results.30.hlaa$Allele)) #225
thermophilus.final.allele.30.hlaa <- data.frame(unique(thermophilus.results.30.hlaa$Allele)) #200

write.csv(aureus.final.allele.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Results_HLA_IC50_30.csv")
write.csv(pyogenes.final.allele.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Results_HLA_IC50_30.csv")
write.csv(thermophilus.final.allele.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Results_HLA_IC50_30.csv")



#- Obtaining list of peptides -#

aureus.final.peptides.30.hlaa <- data.frame(unique(aureus.results.30.hlaa$Peptide)) #121
pyogenes.final.peptides.30.hlaa <- data.frame(unique(pyogenes.results.30.hlaa$Peptide)) #163
thermophilus.final.peptides.30.hlaa <- data.frame(unique(thermophilus.results.30.hlaa$Peptide)) #189


write.csv(aureus.final.peptides.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Results_PEP_IC50_30.csv")
write.csv(pyogenes.final.peptides.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Results_PEP_IC50_30.csv")
write.csv(thermophilus.final.peptides.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Results_PEP_IC50_30.csv")


common.aurpyo.30.hlaa <- data.frame(merge(aureus.final.peptides.30.hlaa, pyogenes.final.peptides.30.hlaa, 
                                     by.x = colnames(aureus.final.peptides.30.hlaa), by.y = colnames(pyogenes.final.peptides.30.hlaa)))
common.aurthe.30.hlaa <- data.frame(merge(aureus.final.peptides.30.hlaa, thermophilus.final.peptides.30.hlaa, 
                                     by.x = colnames(aureus.final.peptides.30.hlaa), by.y = colnames(thermophilus.final.peptides.30.hlaa)))
common.pyothe.30.hlaa <- data.frame(merge(pyogenes.final.peptides.30.hlaa, thermophilus.final.peptides.30.hlaa, 
                                     by.x = colnames(pyogenes.final.peptides.30.hlaa), by.y = colnames(thermophilus.final.peptides.30.hlaa)))


write.csv(common.aurpyo.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Common_HLAwithPyogenes_IC50_30.csv")
write.csv(common.aurthe.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Common_HLAwithThermophilus_IC50_30.csv")
write.csv(common.pyothe.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Common_HLAwithThermophilus_IC50_30.csv")


# Qual o score de imunogenicidade pra os comuns pyo/the? 16 no total!
score.common.pyothe.30.hlaa <- merge(common.pyothe.30.hlaa, immunogenicity.hlaa, by.x="unique.pyogenes.results.30.hlaa.Peptide.", by.y = "peptide")
write.csv(score.common.pyothe.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_CommonScore_PEP_PyoThe_IC50_30.csv")


####
## LOGOS IC50 < 30nM
####


pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Aureus_Peptides_IC50_30.pdf")
ggplot() + geom_logo( aureus.final.peptides.30.hlaa ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Pyogenes_Peptides_IC50_30.pdf")
ggplot() + geom_logo( pyogenes.final.peptides.30.hlaa ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Thermophilus_peptides_IC50_30.pdf")
ggplot() + geom_logo( thermophilus.final.peptides.30.hlaa ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_PyoThe_peptides_IC50_30.pdf")
ggplot() + geom_logo( common.pyothe.30.hlaa ) + theme_logo()
dev.off()

##--    Populations   --#

setwd("/Users/martielafreitas/Desktop/HLA-Rproject/AInputs/Populacoes")
populacoes.hlaa <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/AInputs/Populacoes/HLA-A_GOLD_Geral", sep=";")

# Merged tables
aureus.pop.30.hlaa <-  merge(aureus.results.30.hlaa, populacoes.hlaa, by.x = "Allele", by.y = "Allele")
pyogenes.pop.30.hlaa <-  merge(pyogenes.results.30.hlaa, populacoes.hlaa, by.x = "Allele", by.y = "Allele")
thermophilus.pop.30.hlaa <-  merge(thermophilus.results.30.hlaa, populacoes.hlaa, by.x = "Allele", by.y = "Allele")

write.table(aureus.pop.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Path_Populacao_IC50_30.csv")
write.table(pyogenes.pop.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Path_Populacao_IC50_30.csv")
write.table(thermophilus.pop.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Path_Populacao_IC50_30.csv")

# Unique peptides by Population and allele

aureus.populacao.30.hlaa <- aureus.pop.30.hlaa[,c(1,2,13,14,15,16)]
pyogenes.populacao.30.hlaa <- pyogenes.pop.30.hlaa[,c(1,2,13,14,15,16)]
thermophilus.populacao.30.hlaa <- thermophilus.pop.30.hlaa[,c(1,2,13,14,15,16)]

aureus.ethnicity.30.hlaa <- aureus.populacao.30.hlaa[,c(1,2,6)]
pyogenes.ethnicity.30.hlaa <- pyogenes.populacao.30.hlaa[,c(1,2,6)]
thermophilus.ethnicity.30.hlaa <- thermophilus.populacao.30.hlaa[,c(1,2,6)]

aureus.ethnicity.unq.30.hlaa <- unique(aureus.ethnicity.30.hlaa)
pyogenes.ethnicity.unq.30.hlaa <- unique(pyogenes.ethnicity.30.hlaa)
thermophilus.ethnicity.unq.30.hlaa <-unique(thermophilus.ethnicity.30.hlaa)

write.table(aureus.ethnicity.unq.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Ethnicity_Unique_IC50_30.csv")
write.table(pyogenes.ethnicity.unq.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Ethnicity_Unique_IC50_30.csv")
write.table(thermophilus.ethnicity.unq.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Ethnicity_Unique_IC50_30.csv")

### PEPTIDES

aureus.sankey.30.hlaa <- aureus.populacao.30.hlaa[,c(1,2,6,4)]
pyogenes.sankey.30.hlaa <- pyogenes.populacao.30.hlaa[,c(1,2,6,4)]
thermophilus.sankey.30.hlaa <- thermophilus.populacao.30.hlaa[,c(1,2,6,4)]

write.table(aureus.sankey.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Sankey_PEP_IC50_30.csv")
write.table(pyogenes.sankey.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Sankey_PEP_IC50_30.csv")
write.table(thermophilus.sankey.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Sankey_PEP_IC50_30.csv")

#aureus.sankey.30 <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Sankey_HLA-A.30.csv", sep = ";")
#pyogenes.sankey.30 <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Sankey_HLA-A.30.csv", sep = ";")
#thermophilus.sankey.30 <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Sankey_HLA-A.30.csv", sep = ";")


################################################################################## AUREUS

###################### AUREUS ORIENTAL 30

aureus.sankey.oriental.30.hlaa <- unique(subset(aureus.sankey.30.hlaa, EthnicOrigin == "Oriental"))
aureus.sankey.oriental.30.hlaa <- aggregate(x = as.numeric(aureus.sankey.oriental.30.hlaa$Frequency),
                                       by = c(list(aureus.sankey.oriental.30.hlaa$Allele, 
                                                   aureus.sankey.oriental.30.hlaa$Peptide, 
                                                   aureus.sankey.oriental.30.hlaa$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlaa <- subset(aureus.sankey.oriental.30.hlaa, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlaa$x),
                            by = c(list(aureus.sankey.top.hlaa$Group.1, 
                                        aureus.sankey.top.hlaa$Group.2, 
                                        aureus.sankey.top.hlaa$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Aureus_Oriental_30_10up.html", selfcontained = TRUE)

rm(sankey, aureus.sankey.top.hlaa, aureus.snk.agg, links, source, target, value, nodes, name, rede)

###################### BLACK 30

aureus.sankey.black.30.hlaa <- unique(subset(aureus.sankey.30.hlaa, EthnicOrigin == "Black"))
aureus.sankey.black.30.hlaa <- aggregate(x = as.numeric(aureus.sankey.black.30.hlaa$Frequency),
                                    by = c(list(aureus.sankey.black.30.hlaa$Allele, 
                                                aureus.sankey.black.30.hlaa$Peptide, 
                                                aureus.sankey.black.30.hlaa$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlaa <- subset(aureus.sankey.black.30.hlaa, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlaa$x),
                            by = c(list(aureus.sankey.top.hlaa$Group.1, 
                                        aureus.sankey.top.hlaa$Group.2, 
                                        aureus.sankey.top.hlaa$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Aureus_Black_30_10up.html", selfcontained = TRUE)

rm(sankey, aureus.sankey.top.hlaa, aureus.snk.agg, links, source, target, value, nodes, name, rede)


##################### CAUCASOIDE 30

aureus.sankey.cauc.30.hlaa <- unique(subset(aureus.sankey.30.hlaa, EthnicOrigin == "Caucasoid"))
aureus.sankey.cauc.30.hlaa <- aggregate(x = as.numeric(aureus.sankey.cauc.30.hlaa$Frequency),
                                   by = c(list(aureus.sankey.cauc.30.hlaa$Allele, 
                                               aureus.sankey.cauc.30.hlaa$Peptide, 
                                               aureus.sankey.cauc.30.hlaa$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlaa <- subset(aureus.sankey.cauc.30.hlaa, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlaa$x),
                            by = c(list(aureus.sankey.top.hlaa$Group.1, 
                                        aureus.sankey.top.hlaa$Group.2, 
                                        aureus.sankey.top.hlaa$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Aureus_Caucasoid_30_10up.html", selfcontained = TRUE)

rm(sankey,aureus.sankey.top.hlaa, aureus.snk.agg, links, source, target, value, nodes, name, rede)


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

##write.table(combined_peptides.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Aureus_Peptides_3pop_HLA-A.30.csv")
##write.table(duplicate_rows.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Aureus_Common_3pop_HLA-A.30.csv")


################################################################################ PYOGENES

###################### PYOGENES ORIENTAL 30

pyogenes.sankey.oriental.30.hlaa <- unique(subset(pyogenes.sankey.30.hlaa, EthnicOrigin == "Oriental"))
pyogenes.sankey.oriental.30.hlaa <- aggregate(x = as.numeric(pyogenes.sankey.oriental.30.hlaa$Frequency),
                                         by = c(list(pyogenes.sankey.oriental.30.hlaa$Allele, 
                                                     pyogenes.sankey.oriental.30.hlaa$Peptide, 
                                                     pyogenes.sankey.oriental.30.hlaa$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlaa <- subset(pyogenes.sankey.oriental.30.hlaa, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlaa$x),
                            by = c(list(pyogenes.sankey.top.hlaa$Group.1, 
                                        pyogenes.sankey.top.hlaa$Group.2, 
                                        pyogenes.sankey.top.hlaa$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Pyogenes_Oriental_30_10up.html", selfcontained = TRUE)

rm(sankey, pyogenes.sankey.top.hlaa, pyogenes.snk.agg, links, source, target, value, nodes, name, rede)

###################### PYOGENES BLACK 30

pyogenes.sankey.black.30.hlaa <- unique(subset(pyogenes.sankey.30.hlaa, EthnicOrigin == "Black"))
pyogenes.sankey.black.30.hlaa <- aggregate(x = as.numeric(pyogenes.sankey.black.30.hlaa$Frequency),
                                      by = c(list(pyogenes.sankey.black.30.hlaa$Allele, 
                                                  pyogenes.sankey.black.30.hlaa$Peptide, 
                                                  pyogenes.sankey.black.30.hlaa$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlaa <- subset(pyogenes.sankey.black.30.hlaa, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlaa$x),
                            by = c(list(pyogenes.sankey.top.hlaa$Group.1, 
                                        pyogenes.sankey.top.hlaa$Group.2, 
                                        pyogenes.sankey.top.hlaa$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Pyogenes_Black_30_10up.html", selfcontained = TRUE)

rm(sankey, pyogenes.sankey.top.hlaa, pyogenes.snk.agg, links, source, target, value, nodes, name, rede)

###################### PYOGENES CAUCASOID 30

pyogenes.sankey.cauc.30.hlaa <- unique(subset(pyogenes.sankey.30.hlaa, EthnicOrigin == "Caucasoid"))
pyogenes.sankey.cauc.30.hlaa <- aggregate(x = as.numeric(pyogenes.sankey.cauc.30.hlaa$Frequency),
                                     by = c(list(pyogenes.sankey.cauc.30.hlaa$Allele, 
                                                 pyogenes.sankey.cauc.30.hlaa$Peptide, 
                                                 pyogenes.sankey.cauc.30.hlaa$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlaa <- subset(pyogenes.sankey.cauc.30.hlaa, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlaa$x),
                            by = c(list(pyogenes.sankey.top.hlaa$Group.1, 
                                        pyogenes.sankey.top.hlaa$Group.2, 
                                        pyogenes.sankey.top.hlaa$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Pyogenes_Caucasoid_30_10up.html", selfcontained = TRUE)

rm(sankey, pyogenes.sankey.top.hlaa, pyogenes.snk.agg, links, source, target, value, nodes, name, rede)



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

##write.table(combined_peptides.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Pyogenes_Peptides_3pop_HLA-A.30.csv")
##write.table(duplicate_rows.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Pyogenes_Common_3pop_HLA-A.30.csv")


################################################################################ THERMOPHILUS

###################### THERMOPHILUS 30 ORIENTAL

thermophilus.sankey.oriental.30.hlaa <- unique(subset(thermophilus.sankey.30.hlaa, EthnicOrigin == "Oriental"))
thermophilus.sankey.oriental.30.hlaa <- aggregate(x = as.numeric(thermophilus.sankey.oriental.30.hlaa$Frequency),
                                             by = c(list(thermophilus.sankey.oriental.30.hlaa$Allele, 
                                                         thermophilus.sankey.oriental.30.hlaa$Peptide, 
                                                         thermophilus.sankey.oriental.30.hlaa$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlaa <- subset(thermophilus.sankey.oriental.30.hlaa, x >= 0.10)
# exclude lines wth e-4...n
# thermophilus.sankey.10 <- thermophilus.sankey.top[- grep("e-", thermophilus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlaa$x),
                            by = c(list(thermophilus.sankey.top.hlaa$Group.1, 
                                        thermophilus.sankey.top.hlaa$Group.2, 
                                        thermophilus.sankey.top.hlaa$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Thermophilus_Oriental_30_10up.html", selfcontained = TRUE)

rm(sankey, thermophilus.sankey.top.hlaa, thermophilus.snk.agg, links, source, target, value, nodes, name, rede)


##################### THERMOPHILUS 30 BLACK

thermophilus.sankey.black.30.hlaa <- unique(subset(thermophilus.sankey.30.hlaa, EthnicOrigin == "Black"))
thermophilus.sankey.black.30.hlaa <- aggregate(x = as.numeric(thermophilus.sankey.black.30.hlaa$Frequency),
                                          by = c(list(thermophilus.sankey.black.30.hlaa$Allele, 
                                                      thermophilus.sankey.black.30.hlaa$Peptide, 
                                                      thermophilus.sankey.black.30.hlaa$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlaa <- subset(thermophilus.sankey.black.30.hlaa, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlaa$x),
                            by = c(list(thermophilus.sankey.top.hlaa$Group.1, 
                                        thermophilus.sankey.top.hlaa$Group.2, 
                                        thermophilus.sankey.top.hlaa$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Thermophilus_Black_30_10up.html", selfcontained = TRUE)

rm(sankey, thermophilus.sankey.top.hlaa, thermophilus.snk.agg, links, source, target, value, nodes, name, rede)

##################### THERMOPHILUS 30 CAUCASOID

thermophilus.sankey.cauc.30.hlaa <- unique(subset(thermophilus.sankey.30.hlaa, EthnicOrigin == "Caucasoid"))
thermophilus.sankey.cauc.30.hlaa <- aggregate(x = as.numeric(thermophilus.sankey.cauc.30.hlaa$Frequency),
                                         by = c(list(thermophilus.sankey.cauc.30.hlaa$Allele, 
                                                     thermophilus.sankey.cauc.30.hlaa$Peptide, 
                                                     thermophilus.sankey.cauc.30.hlaa$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlaa <- subset(thermophilus.sankey.cauc.30.hlaa, x >= 0.10)
# exclude lines wth e-4...n
# thermophilus.sankey.10 <- thermophilus.sankey.top[- grep("e-", thermophilus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlaa$x),
                            by = c(list(thermophilus.sankey.top.hlaa$Group.1, 
                                        thermophilus.sankey.top.hlaa$Group.2, 
                                        thermophilus.sankey.top.hlaa$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-A_Thermophilus_Caucasoid_30_10up.html", selfcontained = TRUE)

rm(sankey,thermophilus.sankey.top.hlaa, thermophilus.snk.agg, links, source, target, value, nodes, name, rede)


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

##write.table(combined_peptides.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Thermophilus_Peptides_3pop_HLA-A.30.csv")
##write.table(duplicate_rows.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Thermophilus_Common_3pop_HLA-A.30.csv")

###############################################################################
###############################################################################
###############################################################################

# Comparações HLA-A

# Geral IC50 < 500nM
unq.aureus.pep.hlaa <- unique(aureus.IC50.hlaa$peptide)
unq.pyogenes.pep.hlaa <- unique(pyogenes.IC50.hlaa$peptide)
unq.thermophilus.pep.hlaa <- unique(thermophilus.IC50.hlaa$peptide)

write.table(unq.aureus.pep.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Unique_PEP_IC50_500.csv")
write.table(unq.pyogenes.pep.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Unique_PEP_IC50_500.csv")
write.table(unq.thermophilus.pep.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Unique_PEP_IC50_500.csv")

# Geral IC50 < 30nM
unq.aureus.pep.30.hlaa <- unique(aureus.IC50.30.hlaa$peptide)
unq.pyogenes.pep.30.hlaa <- unique(pyogenes.IC50.30.hlaa$peptide)
unq.thermophilu.pep.30.hlaa <- unique(thermophilus.IC50.30.hlaa$peptide)

write.table(unq.aureus.pep.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Unique_PEP_IC50_30.csv")
write.table(unq.pyogenes.pep.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Unique_PEP_IC50_30.csv")
write.table(unq.thermophilus.pep.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Unique_PEP_IC50_30.csv")

# Por população IC50 < 500nM

## Sp. aureus

### Black

hlaa.aureus.sankey.black.500.hlaa <- unique(aureus.sankey.black.hlaa$Group.1)
hlaa.aureus.sankey.black.30.hlaa <- unique(aureus.sankey.black.30.hlaa$Group.1)
hlaa.aureus.sankey.black.500.pep <- unique(aureus.sankey.black.hlaa$Group.2)
hlaa.aureus.sankey.black.30.pep <- unique(aureus.sankey.black.30.hlaa$Group.2)

write.table(hlaa.aureus.sankey.black.500.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Unique_HLA_Black_IC50_500.csv")
write.table(hlaa.aureus.sankey.black.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Unique_HLA_Black_IC50_30.csv")
write.table(hlaa.aureus.sankey.black.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Unique_PEP_Black_IC50_500.csv")
write.table(hlaa.aureus.sankey.black.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Unique_PEP_Black_IC50_30.csv")

### Caucasoid

hlaa.aureus.sankey.cauc.500.hlaa <- unique(aureus.sankey.caucasoid.hlaa$Group.1)
hlaa.aureus.sankey.cauc.30.hlaa <- unique(aureus.sankey.cauc.30.hlaa$Group.1)
hlaa.aureus.sankey.cauc.500.pep <- unique(aureus.sankey.caucasoid.hlaa$Group.2)
hlaa.aureus.sankey.cauc.30.pep <- unique(aureus.sankey.cauc.30.hlaa$Group.2)

write.table(hlaa.aureus.sankey.cauc.500.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Unique_HLA_Cauc_IC50_500.csv")
write.table(hlaa.aureus.sankey.cauc.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Unique_HLA_Cauc_IC50_30.csv")
write.table(hlaa.aureus.sankey.cauc.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Unique_PEP_Cauc_IC50_500.csv")
write.table(hlaa.aureus.sankey.cauc.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Unique_PEP_Cauc_IC50_30.csv")

### Oriental

hlaa.aureus.sankey.oriental.500.hlaa <- unique(aureus.sankey.oriental.hlaa$Group.1)
hlaa.aureus.sankey.oriental.30.hlaa <- unique(aureus.sankey.oriental.30.hlaa$Group.1)
hlaa.aureus.sankey.oriental.500.pep <- unique(aureus.sankey.oriental.hlaa$Group.2)
hlaa.aureus.sankey.oriental.30.pep <- unique(aureus.sankey.oriental.30.hlaa$Group.2)

write.table(hlaa.aureus.sankey.oriental.500.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Unique_HLA_Orient_IC50_500.csv")
write.table(hlaa.aureus.sankey.oriental.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Unique_HLA_Orient_IC50_30.csv")
write.table(hlaa.aureus.sankey.oriental.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Unique_PEP_Orient_IC50_500.csv")
write.table(hlaa.aureus.sankey.oriental.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Unique_PEP_Orient_IC50_30.csv")

## Sp. pyogenes

### Black

hlaa.pyogenes.sankey.black.500.hlaa <- unique(pyogenes.sankey.black.hlaa$Group.1)
hlaa.pyogenes.sankey.black.30.hlaa <- unique(pyogenes.sankey.black.30.hlaa$Group.1)
hlaa.pyogenes.sankey.black.500.pep <- unique(pyogenes.sankey.black.hlaa$Group.2)
hlaa.pyogenes.sankey.black.30.pep <- unique(pyogenes.sankey.black.30.hlaa$Group.2)

write.table(hlaa.pyogenes.sankey.black.500.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Unique_HLA_Black_IC50_500.csv")
write.table(hlaa.pyogenes.sankey.black.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Unique_HLA_Black_IC50_30.csv")
write.table(hlaa.pyogenes.sankey.black.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Unique_PEP_Black_IC50_500.csv")
write.table(hlaa.pyogenes.sankey.black.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Unique_PEP_Black_IC50_30.csv")

### Caucasoid

hlaa.pyogenes.sankey.cauc.500.hla <- unique(pyogenes.sankey.cauc.hlaa$Group.1)
hlaa.pyogenes.sankey.cauc.30.hla <- unique(pyogenes.sankey.cauc.30.hlaa$Group.1)
hlaa.pyogenes.sankey.cauc.500.pep <- unique(pyogenes.sankey.cauc.hlaa$Group.2)
hlaa.pyogenes.sankey.cauc.30.pep <- unique(pyogenes.sankey.cauc.30.hlaa$Group.2)

write.table(hlaa.pyogenes.sankey.cauc.500.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Unique_HLA_Cauc_IC50_500.csv")
write.table(hlaa.pyogenes.sankey.cauc.30.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Unique_HLA_Cauc_IC50_30.csv")
write.table(hlaa.pyogenes.sankey.cauc.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Unique_PEP_Cauc_IC50_500.csv")
write.table(hlaa.pyogenes.sankey.cauc.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Unique_PEP_Cauc_IC50_30.csv")

### Oriental

hlaa.pyogenes.sankey.oriental.500.hla <- unique(pyogenes.sankey.oriental.hlaa$Group.1)
hlaa.pyogenes.sankey.oriental.30.hla <- unique(pyogenes.sankey.oriental.30.hlaa$Group.1)

hlaa.pyogenes.sankey.oriental.500.pep <- unique(pyogenes.sankey.oriental.hlaa$Group.2)
hlaa.pyogenes.sankey.oriental.30.pep <- unique(pyogenes.sankey.oriental.30.hlaa$Group.2)

write.table(hlaa.pyogenes.sankey.oriental.500.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Unique_HLA_Orient_IC50_500.csv")
write.table(hlaa.pyogenes.sankey.oriental.30.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Unique_HLA_Orient_IC50_30.csv")
write.table(hlaa.pyogenes.sankey.oriental.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Unique_PEP_Orient_IC50_500.csv")
write.table(hlaa.pyogenes.sankey.oriental.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Unique_PEP_Orient_IC50_30.csv")

## Sp. thermophilus
### Black
hlaa.thermophilus.sankey.black.500.hlaa <- unique(thermophilus.sankey.black.hlaa$Group.1)
hlaa.thermophilus.sankey.black.30.hlaa <- unique(thermophilus.sankey.black.30.hlaa$Group.1)

hlaa.thermophilus.sankey.black.500.pep <- unique(thermophilus.sankey.black.hlaa$Group.2)
hlaa.thermophilus.sankey.black.30.pep <- unique(thermophilus.sankey.black.30.hlaa$Group.2)

write.table(hlaa.thermophilus.sankey.black.500.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Unique_HLA_Black_IC50_500.csv")
write.table(hlaa.thermophilus.sankey.black.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Unique_HLA_Black_IC50_30.csv")
write.table(hlaa.thermophilus.sankey.black.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Unique_PEP_Black_IC50_500.csv")
write.table(hlaa.thermophilus.sankey.black.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Unique_PEP_Black_IC50_30.csv")

### Caucasoid
hlaa.thermophilus.sankey.cauc.500.hlaa <- unique(thermophilus.sankey.cauc.hlaa$Group.1)
hlaa.thermophilus.sankey.cauc.30.hlaa <- unique(thermophilus.sankey.cauc.30.hlaa$Group.1)

hlaa.thermophilus.sankey.cauc.500.pep <- unique(thermophilus.sankey.cauc.hlaa$Group.2)
hlaa.thermophilus.sankey.cauc.30.pep <- unique(thermophilus.sankey.cauc.30.hlaa$Group.2)

write.table(hlaa.thermophilus.sankey.cauc.500.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Unique_HLA_Cauc_IC50_500.csv")
write.table(hlaa.thermophilus.sankey.cauc.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Unique_HLA_Cauc_IC50_30.csv")
write.table(hlaa.thermophilus.sankey.cauc.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Unique_PEP_Cauc_IC50_500.csv")
write.table(hlaa.thermophilus.sankey.cauc.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Unique_PEP_Cauc_IC50_30.csv")

### Oriental
hlaa.thermophilus.sankey.oriental.500.hla <- unique(thermophilus.sankey.oriental.hlaa$Group.1)
hlaa.thermophilus.sankey.oriental.30.hla <- unique(thermophilus.sankey.oriental.30.hlaa$Group.1)

hlaa.thermophilus.sankey.oriental.500.pep <- unique(thermophilus.sankey.oriental.hlaa$Group.2)
hlaa.thermophilus.sankey.oriental.30.pep <- unique(thermophilus.sankey.oriental.30.hlaa$Group.2)

write.table(hlaa.thermophilus.sankey.oriental.500.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Unique_HLA_Orient_IC50_500.csv")
write.table(hlaa.thermophilus.sankey.oriental.30.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Unique_HLA_Orient_IC50_30.csv")
write.table(hlaa.thermophilus.sankey.oriental.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Unique_PEP_Orient_IC50_500.csv")
write.table(hlaa.pyogenes.sankey.oriental.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Unique_PEP_Orient_IC50_30.csv")

###############################################################################
###############################################################################
###############################################################################

## For HLA-Arena 500
aureus.arena.hlaa <- aureus.ethnicity.unq.hlaa[,c(1,2)]
aureus.arena.unq.hlaa <- unique(aureus.arena.hlaa)
aureus.arena.hla.unq.hlaa <- unique(aureus.arena.hlaa$Allele)
aureus.arena.pep.unq.hlaa <- unique(aureus.arena.hlaa$Peptide)

write.table(aureus.arena.unq.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Arena_Unq_HLA-A.csv")
write.table(aureus.arena.hla.unq.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Arena_HlaUnq_HLA-A.csv")
write.table(aureus.arena.pep.unq.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Arena_PepUnq_HLA-A.csv")

pyogenes.arena.hlaa <- pyogenes.ethnicity.unq.hlaa[,c(1,2)]
pyogenes.arena.unq.hlaa <- unique(pyogenes.arena.hlaa)
pyogenes.arena.hla.unq.hlaa <- unique(pyogenes.arena.hlaa$Allele)
pyogenes.arena.pep.unq.hlaa <- unique(pyogenes.arena.hlaa$Peptide)

write.table(pyogenes.arena.unq.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Arena_Unq_HLA-A.csv")
write.table(pyogenes.arena.hla.unq.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Arena_HlaUnq_HLA-A.csv")
write.table(pyogenes.arena.pep.unq.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Arena_PepUnq_HLA-A.csv")

thermophilus.arena.hlaa <- thermophilus.ethnicity.unq.hlaa[,c(1,2)]
thermophilus.arena.unq.hlaa <- unique(thermophilus.arena.hlaa)
thermophilus.arena.hla.unq.hlaa <- unique(thermophilus.arena.hlaa$Allele)
thermophilus.arena.pep.unq.hlaa <- unique(thermophilus.arena.hlaa$Peptide)

write.table(thermophilus.arena.unq.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Arena_Unq_HLA-A.csv")
write.table(thermophilus.arena.hla.unq.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Arena_HlaUnq_HLA-A.csv")
write.table(thermophilus.arena.pep.unq.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Arena_PepUnq_HLA-A.csv")

## For HLA-Arena 30
aureus.arena.30.hlaa <- aureus.ethnicity.unq.30.hlaa[,c(1,2)]
aureus.arena.unq.30.hlaa <- unique(aureus.arena.30.hlaa)
aureus.arena.hla.unq.30.hlaa <- unique(aureus.arena.30.hlaa$Allele)
aureus.arena.pep.unq.30.hlaa <- unique(aureus.arena.30.hlaa$Peptide)

write.table(aureus.arena.unq.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Arena_Unq_30_HLA-A.csv")
write.table(aureus.arena.hla.unq.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Arena_HlaUnq_30_HLA-A.csv")
write.table(aureus.arena.pep.unq.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_Arena_PepUnq_30_HLA-A.csv")

pyogenes.arena.30.hlaa <- pyogenes.ethnicity.unq.30.hlaa[,c(1,2)]
pyogenes.arena.unq.30.hlaa <- unique(pyogenes.arena.30.hlaa)
pyogenes.arena.hla.unq.30.hlaa <- unique(pyogenes.arena.30.hlaa$Allele)
pyogenes.arena.pep.unq.30.hlaa <- unique(pyogenes.arena.30.hlaa$Peptide)

write.table(pyogenes.arena.unq.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Arena_Unq_30_HLA-A.csv")
write.table(pyogenes.arena.hla.unq.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Arena_HlaUnq_30_HLA-A.csv")
write.table(pyogenes.arena.pep.unq.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Pyogenes_Arena_PepUnq_30_HLA-A.csv")

thermophilus.arena.30.hlaa <- thermophilus.ethnicity.unq.30.hlaa[,c(1,2)]
thermophilus.arena.unq.30.hlaa <- unique(thermophilus.arena.30.hlaa)
thermophilus.arena.hla.unq.30.hlaa <- unique(thermophilus.arena.30.hlaa$Allele)
thermophilus.arena.pep.unq.30.hlaa <- unique(thermophilus.arena.30.hlaa$Peptide)

write.table(thermophilus.arena.unq.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Arena_Unq_30_HLA-A.csv")
write.table(thermophilus.arena.hla.unq.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Arena_HlaUnq_30_HLA-A.csv")
write.table(thermophilus.arena.pep.unq.30.hlaa, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Thermophilus_Arena_PepUnq_30_HLA-A.csv")

