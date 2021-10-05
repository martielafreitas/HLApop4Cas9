###############################################################################
# VIA DE PROCESSAMENTO
###############################################################################

# Unindo os resulatdos dos preditores para as três proteínas 
# em uma única tabela que adiciona o nome que está no nome do
# arquivo, na última coluna da tabela.

setwd("/Users/martielafreitas/Desktop/HLA-Rproject")

# variavel que contem o caminho para os arquivos
file_list <- list.files(path="/Users/martielafreitas/Desktop/HLA-Rproject/BInputs/Preditores/1.Processing")

# inicia uma dataframe vazia
process.hlab <- data.frame()

# onde serao processados os arquivos
setwd("/Users/martielafreitas/Desktop/HLA-Rproject/BInputs/Preditores/1.Processing")

# laco que une todos os arquivos
for(i in 1:length(file_list)){
  #each file will be read in, specify which columns you need read in to avoid any errors
  temp_data <- read.csv(file_list[i], sep=",")
  #clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
  temp_data$Class <- sapply(strsplit(gsub(".csv", "", file_list[i]), " - "), function(x){x[3]})
  #for each iteration, bind the new data to the building dataset
  process.hlab <- rbind(process.hlab, temp_data)
}

# salvando o resultado em um arquivo .csv
write.csv(process.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Processing_Geral.csv")

# Separando cada resultado por organismo de onde vem a proteina
aureus.process.hlab<- subset(process.hlab, Class == "Aureus")
pyogenes.process.hlab <- subset(process.hlab, Class == "Pyogenes")
thermophilus.process.hlab <- subset(process.hlab, Class == "Thermophilus")

# salvando o resultado em um arquivo .csv
write.csv(aureus.process.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Process.csv")
write.csv(pyogenes.process.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Process.csv")
write.csv(thermophilus.process.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Process.csv")


###############################################################################
# LIGAÇÃO
###############################################################################

##-- Binding Results --##

# variavel que contem o caminho para os arquivos
file_list <- list.files(path="/Users/martielafreitas/Desktop/HLA-Rproject/BInputs/Preditores/2.Binding")

#inicia uma dataframe vazia
binding.hlab <- data.frame()

setwd("/Users/martielafreitas/Desktop/HLA-Rproject/BInputs/Preditores/2.Binding")
for(i in 1:length(file_list)){
  #each file will be read in, specify which columns you need read in to avoid any errors
  temp_data <- read.csv(file_list[i], sep=";")
  #clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
  temp_data$Class <- sapply(strsplit(gsub(".csv", "", file_list[i]), " - "), function(x){x[3]})
  #for each iteration, bind the new data to the building dataset
  binding.hlab <- rbind(binding.hlab, temp_data)
}
# salvando resultados gerais de ligacao em um arquivo .csv
write.csv(binding.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Binding_Geral.csv")

# separando os resultados por organismo de onde veio cada proteina
aureus.binding.hlab <- subset(binding.hlab, Class == "Aureus")
pyogenes.binding.hlab <- subset(binding.hlab, Class == "Pyogenes")
thermophilus.binding.hlab <- subset(binding.hlab, Class == "Thermophilus")

# salvando resultados gerais por organismo em um arquivo .csv
write.csv(aureus.binding.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Binding.csv")
write.csv(pyogenes.binding.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Binding.csv")
write.csv(thermophilus.binding.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Binding.csv")



###############################################################################
# IMMUNOGENICIDADE
###############################################################################

##-- Comparing Results for immunogenicity predictions --##

aureus.procbind.hlab <- merge(aureus.binding.hlab, 
                              aureus.process.hlab, 
                              by.x=c("allele", "seq_num", "start", "end", "length", "peptide", "Class"),
                              by.y=c("Allele", "X.", "Start", "End", "Peptide.Length", "Peptide", "Class"))

pyogenes.procbind.hlab <- merge(pyogenes.binding.hlab, 
                                pyogenes.process.hlab, 
                                by.x=c("allele", "seq_num", "start", "end", "length", "peptide", "Class"),
                                by.y=c("Allele", "X.", "Start", "End", "Peptide.Length", "Peptide", "Class"))

thermophilus.procbind.hlab <- merge(thermophilus.binding.hlab, 
                                    thermophilus.process.hlab, 
                                    by.x=c("allele", "seq_num", "start", "end", "length", "peptide", "Class"),
                                    by.y=c("Allele", "X.", "Start", "End", "Peptide.Length", "Peptide", "Class"))

write.csv(aureus.procbind.hlab, file="/Users/martielafreitas/Documents/Rprojects/HLA-B/Results/HLA-B_Aureus_Merged_Process_and_Binding.csv")
write.csv(pyogenes.procbind.hlab, file="/Users/martielafreitas/Documents/Rprojects/HLA-B/Results/HLA-B_Pyogenes_Merged_Process_and_Binding.csv")
write.csv(thermophilus.procbind.hlab, file="/Users/martielafreitas/Documents/Rprojects/HLA-B/Results/HLA-B_Thermophilus_Merged_Process_and_Binding.csv")



###############################################################################
# ANALISE DE RESULTADOS
###############################################################################

## Immunogenicity Results

## Obtaining list of peptides

## Peptídeos totais - independente de população ou frequência, por organismo
aureus.peptides.hlab <- unique(aureus.procbind.hlab$peptide)
pyogenes.peptides.hlab <- unique(pyogenes.procbind.hlab$peptide)
thermophilus.peptides.hlab <- unique(thermophilus.procbind.hlab$peptide)

write.csv(aureus.peptides.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_TotalPeptides.csv")
write.csv(pyogenes.peptides.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_TotalPeptides.csv")
write.csv(thermophilus.peptides.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_TotalPeptides.csv")

# Esses peptídeos passaram pela ferramenta de predição de imunogenicidade.

#- Calling result files -#
file_list <- list.files(path="/Users/martielafreitas/Desktop/HLA-Rproject/BInputs/Preditores/3.Immunogenicity")
immunogenicity.hlab <- data.frame()

setwd("/Users/martielafreitas/Desktop/HLA-Rproject/BInputs/Preditores/3.Immunogenicity")
for(i in 1:length(file_list)){
  #each file will be read in, specify which columns you need read in to avoid any errors
  temp_data <- read.csv(file_list[i], sep=",")
  #clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
  temp_data$Class <- sapply(strsplit(gsub(".csv", "", file_list[i]), " - "), function(x){x[3]})
  #for each iteration, bind the new data to the building dataset
  immunogenicity.hlab <- rbind(immunogenicity.hlab, temp_data)
}
# salvando resultados gerais de imunogenicidade em um arquivo .csv
write.csv(immunogenicity.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Immunogenicity_Geral.csv")


# Separando a lista de HLAs totais por organismo - independente de população ou frequência

aureus.immuno.hlab <- subset(immunogenicity.hlab, Class == "Aureus")
pyogenes.immuno.hlab <- subset(immunogenicity.hlab, Class == "Pyogenes")
thermophilus.immuno.hlab <- subset(immunogenicity.hlab, Class == "Thermophilus")

write.csv(aureus.immuno.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Immunogenicity.csv")
write.csv(pyogenes.immuno.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Immunogenicity.csv")
write.csv(thermophilus.immuno.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Immunogenicity.csv")


###############################################################################
# Filtering final results
###############################################################################

#- Immuno: score >= 0.0: immunogenic -#

setwd("/Users/martielafreitas/Desktop/HLA-Rproject/Results/")

aureus.immuno.positive.hlab <- subset(aureus.immuno.hlab, score >= 0.0)
pyogenes.immuno.positive.hlab <- subset(pyogenes.immuno.hlab, score >= 0.0)
thermophilus.immuno.positive.hlab <- subset(thermophilus.immuno.hlab, score >= 0.0)

#- Merging tables to IC50 -#
colnames(aureus.procbind.hlab)
colnames(aureus.immuno.positive.hlab)

aureus.full.path.hlab <- merge(aureus.procbind.hlab, aureus.immuno.positive.hlab, by.x = "peptide", by.y = "peptide")
pyogenes.full.path.hlab <- merge(pyogenes.procbind.hlab, pyogenes.immuno.positive.hlab, by.x = "peptide", by.y = "peptide")
thermophilus.full.path.hlab <- merge(thermophilus.procbind.hlab, thermophilus.immuno.positive.hlab, by.x = "peptide", by.y = "peptide")

#- IC50 < 500nM  -#
aureus.IC50.hlab <- aureus.full.path.hlab[(aureus.full.path.hlab$MHC.IC50 <= 500.0),]
pyogenes.IC50.hlab <- pyogenes.full.path.hlab[(pyogenes.full.path.hlab$MHC.IC50 <= 500.0),]
thermophilus.IC50.hlab <- thermophilus.full.path.hlab[(thermophilus.full.path.hlab$MHC.IC50 <= 500.0),]

#- Final Table -#
aureus.results.hlab <- aureus.IC50.hlab[,c(1,2,9,10,11,12,13,14,15,8,17)]
pyogenes.results.hlab <- pyogenes.IC50.hlab[,c(1,2,9,10,11,12,13,14,15,8,17)]
thermophilus.results.hlab <- thermophilus.IC50.hlab[,c(1,2,9,10,11,12,13,14,15,8,17)]

# corrigindo cabeçalho antes de salvar! :)
colnames(aureus.results.hlab)
cols <- c("Peptide", "Allele", "Percentile.Rank", "Proteasome.Score", "TAP.Score",
          "MHC.Score", "Processing.Score", "Total.Score", "MHC.IC50", "Binding.Score", "Immunogenicity.Score")
colnames(aureus.results.hlab) <- cols
colnames(pyogenes.results.hlab) <- cols
colnames(thermophilus.results.hlab) <- cols

write.csv(aureus.results.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_ToltalResults.csv")
write.csv(pyogenes.results.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_TotalResults.csv")
write.csv(thermophilus.results.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_TotalResults.csv")

#- Obtaining list of alleles por organismo-#
aureus.final.allele.hlab <- data.frame(unique(aureus.results.hlab$Allele)) #221
pyogenes.final.allele.hlab <- data.frame(unique(pyogenes.results.hlab$Allele)) #225
thermophilus.final.allele.hlab <- data.frame(unique(thermophilus.results.hlab$Allele)) #200

write.csv(aureus.final.allele.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_TotalHLAResults.csv")
write.csv(pyogenes.final.allele.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_TotalHLAResults.csv")
write.csv(thermophilus.final.allele.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_TotalHLAResults.csv")

#- Obtaining list of peptides por organismo -#
aureus.final.peptides.hlab <- data.frame(unique(aureus.results.hlab$Peptide)) #121
pyogenes.final.peptides.hlab <- data.frame(unique(pyogenes.results.hlab$Peptide)) #163
thermophilus.final.peptides.hlab <- data.frame(unique(thermophilus.results.hlab$Peptide)) #189

write.csv(aureus.final.peptides.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_TotalPEPResults.csv")
write.csv(pyogenes.final.peptides.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_TotalPEPResults.csv")
write.csv(thermophilus.final.peptides.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_TotalPEPResults.csv")


common.aurpyo.hlab <- data.frame(merge(aureus.final.peptides.hlab, pyogenes.final.peptides.hlab, 
                                  by.x = colnames(aureus.final.peptides.hlab), by.y = colnames(pyogenes.final.peptides.hlab)))
common.aurthe.hlab <- data.frame(merge(aureus.final.peptides.hlab, thermophilus.final.peptides.hlab, 
                                  by.x = colnames(aureus.final.peptides.hlab), by.y = colnames(thermophilus.final.peptides.hlab)))
common.pyothe.hlab <- data.frame(merge(pyogenes.final.peptides.hlab, thermophilus.final.peptides.hlab, 
                                  by.x = colnames(pyogenes.final.peptides.hlab), by.y = colnames(thermophilus.final.peptides.hlab)))

# Qual o score de imunogenicidade pra os comuns pyo/the? 16 no total!
score.common.pyothe.hlab <- merge(common.pyothe.hlab, immunogenicity.hlab, by.x="unique.pyogenes.results.hlab.Peptide.", by.y = "peptide")
write.csv(score.common.pyothe.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_CommonPeptides_Pyogenes_Thermophilus.csv")


###############################################################################
# Logos
###############################################################################

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Aureus_PeptidesLogo.pdf")
ggplot() + geom_logo( aureus.final.peptides.hlab ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Pyogenes_PeptidesLogo.pdf")
ggplot() + geom_logo( pyogenes.final.peptides.hlab ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Thermophilus_PeptidesLogo.pdf")
ggplot() + geom_logo( thermophilus.final.peptides.hlab ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_PyoThe_PeptidesLogo.pdf")
ggplot() + geom_logo( common.pyothe.hlab ) + theme_logo()
dev.off()


###############################################################################
# Populations 
###############################################################################
setwd("/Users/martielafreitas/Desktop/HLA-Rproject/BInputs/Populacoes")
populacoes.hlab <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/BInputs/Populacoes/HLA-B_GOLD_Geral", sep=";")

# Merged tables
aureus.pop.hlab <-  merge(aureus.results.hlab, populacoes.hlab, by.x = "Allele", by.y = "Allele")
pyogenes.pop.hlab <-  merge(pyogenes.results.hlab, populacoes.hlab, by.x = "Allele", by.y = "Allele")
thermophilus.pop.hlab <-  merge(thermophilus.results.hlab, populacoes.hlab, by.x = "Allele", by.y = "Allele")

write.table(aureus.pop.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Path_Populacao.csv")
write.table(pyogenes.pop.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Path_Populacao.csv")
write.table(thermophilus.pop.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Path_Populacao.csv")

# Unique peptides by Population and allele

aureus.populacao.hlab <- aureus.pop.hlab[,c(1,2,13,14,15,16)]
pyogenes.populacao.hlab <- pyogenes.pop.hlab[,c(1,2,13,14,15,16)]
thermophilus.populacao.hlab <- thermophilus.pop.hlab[,c(1,2,13,14,15,16)]

aureus.ethnicity.hlab <- aureus.populacao.hlab[,c(1,2,6)]
pyogenes.ethnicity.hlab <- pyogenes.populacao.hlab[,c(1,2,6)]
thermophilus.ethnicity.hlab <- thermophilus.populacao.hlab[,c(1,2,6)]

aureus.ethnicity.unq.hlab <- unique(aureus.ethnicity.hlab)
pyogenes.ethnicity.unq.hlab <- unique(pyogenes.ethnicity.hlab)
thermophilus.ethnicity.unq.hlab <-unique(thermophilus.ethnicity.hlab)

write.table(aureus.ethnicity.unq.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Ethnicity_Unique.csv")
write.table(pyogenes.ethnicity.unq.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Ethnicity_Unique.csv")
write.table(thermophilus.ethnicity.unq.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Ethnicity_Unique.csv")

### PEPTIDES

aureus.sankey.hlab <- aureus.populacao.hlab[,c(1,2,6,4)]
pyogenes.sankey.hlab <- pyogenes.populacao.hlab[,c(1,2,6,4)]
thermophilus.sankey.hlab <- thermophilus.populacao.hlab[,c(1,2,6,4)]

write.table(aureus.sankey.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Sankey.csv")
write.table(pyogenes.sankey.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Sankey.csv")
write.table(thermophilus.sankey.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Sankey.csv")

#aureus.sankey <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Sankey_HLA-B.csv", sep = ";")
#pyogenes.sankey <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Sankey_HLA-B.csv", sep = ";")
#thermophilus.sankey <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Sankey_HLA-B.csv", sep = ";")

############################################################################### AUREUS

##################### AUREUS 500 ORIENTAL

aureus.sankey.oriental.hlab <- unique(subset(aureus.sankey.hlab, EthnicOrigin == "Oriental"))
aureus.sankey.oriental.hlab <- aggregate(x = as.numeric(aureus.sankey.oriental.hlab$Frequency),
                                    by = c(list(aureus.sankey.oriental.hlab$Allele, 
                                                aureus.sankey.oriental.hlab$Peptide, 
                                                aureus.sankey.oriental.hlab$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlab <- subset(aureus.sankey.oriental.hlab, x >= 0.10)
# exclude lines wth e-4...n
#aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlab$x),
                            by = c(list(aureus.sankey.top.hlab$Group.1, 
                                        aureus.sankey.top.hlab$Group.2, 
                                        aureus.sankey.top.hlab$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Aureus_Oriental_500_10up.html", selfcontained = TRUE)

rm(sankey, aureus.sankey.top, aureus.snk.agg, links, source, target, value, nodes, name, rede)

##################### AUREUS 500 BLACK

aureus.sankey.black.hlab <- unique(subset(aureus.sankey.hlab, EthnicOrigin == "Black"))
aureus.sankey.black.hlab <- aggregate(x = as.numeric(aureus.sankey.black.hlab$Frequency),
                                 by = c(list(aureus.sankey.black.hlab$Allele, 
                                             aureus.sankey.black.hlab$Peptide, 
                                             aureus.sankey.black.hlab$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlab <- subset(aureus.sankey.black.hlab, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlab$x),
                            by = c(list(aureus.sankey.top.hlab$Group.1, 
                                        aureus.sankey.top.hlab$Group.2, 
                                        aureus.sankey.top.hlab$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Aureus_Black_500_10up.html", selfcontained = TRUE)

rm(sankey, aureus.sankey.top, aureus.snk.agg, links, source, target, value, nodes, name, rede)

##################### AUREUS 500 CAUCASOID

aureus.sankey.caucasoid.hlab <- unique(subset(aureus.sankey.hlab, EthnicOrigin == "Caucasoid"))
aureus.sankey.caucasoid.hlab <- aggregate(x = as.numeric(aureus.sankey.caucasoid.hlab$Frequency),
                                     by = c(list(aureus.sankey.caucasoid.hlab$Allele, 
                                                 aureus.sankey.caucasoid.hlab$Peptide, 
                                                 aureus.sankey.caucasoid.hlab$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlab <- subset(aureus.sankey.caucasoid.hlab, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlab$x),
                            by = c(list(aureus.sankey.top.hlab$Group.1, 
                                        aureus.sankey.top.hlab$Group.2, 
                                        aureus.sankey.top.hlab$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Aureus_Caucasoid_500_10up.html", selfcontained = TRUE)

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

#write.table(combined_peptides, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Aureus_Peptides_3pop_HLA-B.csv")
##write.table(duplicate_rows, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Aureus_Common_3pop_HLA-B.csv")



################################################################################ PYOGENES

##################### PYOGENES 500 ORIENTAL

pyogenes.sankey.oriental.hlab <- unique(subset(pyogenes.sankey.hlab, EthnicOrigin == "Oriental"))
pyogenes.sankey.oriental.hlab <- aggregate(x = as.numeric(pyogenes.sankey.oriental.hlab$Frequency),
                                      by = c(list(pyogenes.sankey.oriental.hlab$Allele, 
                                                  pyogenes.sankey.oriental.hlab$Peptide, 
                                                  pyogenes.sankey.oriental.hlab$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlab <- subset(pyogenes.sankey.oriental.hlab, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlab$x),
                            by = c(list(pyogenes.sankey.top.hlab$Group.1, 
                                        pyogenes.sankey.top.hlab$Group.2, 
                                        pyogenes.sankey.top.hlab$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Pyogenes_Oriental_500_10up.html", selfcontained = TRUE)

rm(sankey, pyogenes.sankey.top, pyogenes.snk.agg, links, source, target, value, nodes, name, rede)


##################### PYOGENES 500 BLACK
pyogenes.sankey.black.hlab <- unique(subset(pyogenes.sankey.hlab, EthnicOrigin == "Black"))
pyogenes.sankey.black.hlab <- aggregate(x = as.numeric(pyogenes.sankey.black.hlab$Frequency),
                                   by = c(list(pyogenes.sankey.black.hlab$Allele, 
                                               pyogenes.sankey.black.hlab$Peptide, 
                                               pyogenes.sankey.black.hlab$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlab <- subset(pyogenes.sankey.black.hlab, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlab$x),
                              by = c(list(pyogenes.sankey.top.hlab$Group.1, 
                                          pyogenes.sankey.top.hlab$Group.2, 
                                          pyogenes.sankey.top.hlab$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Pyogenes_Black_500_10up.html", selfcontained = TRUE)

rm(sankey, pyogenes.sankey.top, pyogenes.snk.agg, links, source, target, value, nodes, name, rede)


##################### PYOGENES 500 CAUCASOID
pyogenes.sankey.cauc.hlab <- unique(subset(pyogenes.sankey.hlab, EthnicOrigin == "Caucasoid"))
pyogenes.sankey.cauc.hlab <- aggregate(x = as.numeric(pyogenes.sankey.cauc.hlab$Frequency),
                                  by = c(list(pyogenes.sankey.cauc.hlab$Allele, 
                                              pyogenes.sankey.cauc.hlab$Peptide, 
                                              pyogenes.sankey.cauc.hlab$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlab <- subset(pyogenes.sankey.oriental.hlab, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlab$x),
                            by = c(list(pyogenes.sankey.top.hlab$Group.1, 
                                        pyogenes.sankey.top.hlab$Group.2, 
                                        pyogenes.sankey.top.hlab$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Pyogenes_Caucasoid_500_10up.html", selfcontained = TRUE)

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

##write.table(combined_peptides, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Pyogenes_Peptides_3pop_HLA-B.csv")
##write.table(duplicate_rows, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Pyogenes_Common_3pop_HLA-B.csv")


############################################################################### THERMOPHILUS

##################### THERMOPHILUS 500 ORIENTAL

thermophilus.sankey.oriental.hlab <- unique(subset(thermophilus.sankey.hlab, EthnicOrigin == "Oriental"))
thermophilus.sankey.oriental.hlab <- aggregate(x = as.numeric(thermophilus.sankey.oriental.hlab$Frequency),
                                          by = c(list(thermophilus.sankey.oriental.hlab$Allele, 
                                                      thermophilus.sankey.oriental.hlab$Peptide, 
                                                      thermophilus.sankey.oriental.hlab$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlab <- subset(thermophilus.sankey.oriental.hlab, x >= 0.10)
# exclude lines wth e-4...n
#aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlab$x),
                            by = c(list(thermophilus.sankey.top.hlab$Group.1, 
                                        thermophilus.sankey.top.hlab$Group.2, 
                                        thermophilus.sankey.top.hlab$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Thermophilus_Oriental_500_10up.html", selfcontained = TRUE)

rm(sankey, thermophilus.sankey.top, thermophilus.snk.agg, links, source, target, value, nodes, name, rede)

##################### THERMOPHILUS 500 BLACK

thermophilus.sankey.black.hlab <- unique(subset(thermophilus.sankey.hlab, EthnicOrigin == "Black"))
thermophilus.sankey.black.hlab <- aggregate(x = as.numeric(thermophilus.sankey.black.hlab$Frequency),
                                       by = c(list(thermophilus.sankey.black.hlab$Allele, 
                                                   thermophilus.sankey.black.hlab$Peptide, 
                                                   thermophilus.sankey.black.hlab$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlab <- subset(thermophilus.sankey.black.hlab, x >= 0.10)
# exclude lines wth e-4...n
# thermophilus.sankey.10 <- thermophilus.sankey.top[- grep("e-", thermophilus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlab$x),
                            by = c(list(thermophilus.sankey.top.hlab$Group.1, 
                                        thermophilus.sankey.top.hlab$Group.2, 
                                        thermophilus.sankey.top.hlab$Group.3)), FUN = sum)
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

thermophilus.sankey.cauc.hlab <- unique(subset(thermophilus.sankey.hlab, EthnicOrigin == "Caucasoid"))
thermophilus.sankey.cauc.hlab <- aggregate(x = as.numeric(thermophilus.sankey.cauc.hlab$Frequency),
                                      by = c(list(thermophilus.sankey.cauc.hlab$Allele, 
                                                  thermophilus.sankey.cauc.hlab$Peptide, 
                                                  thermophilus.sankey.cauc.hlab$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlab <- subset(thermophilus.sankey.cauc.hlab, x >= 0.10)
# exclude lines wth e-4...n
# thermophilus.sankey.10 <- thermophilus.sankey.top[- grep("e-", thermophilus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlab$x),
                            by = c(list(thermophilus.sankey.top.hlab$Group.1, 
                                        thermophilus.sankey.top.hlab$Group.2, 
                                        thermophilus.sankey.top.hlab$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Thermophilus_Caucasoid_500_10up.html", selfcontained = TRUE)

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

##write.table(combined_peptides, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Thermophilus_Peptides_3pop_HLA-B.csv")
##write.table(duplicate_rows, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Thermophilus_Common_3pop_HLA-B.csv")


###############################################################################


### ANALISE PARA IC50 < 30

aureus.IC50.30.hlab <- aureus.full.path.hlab[(aureus.full.path.hlab$MHC.IC50 <= 30.0),]
pyogenes.IC50.30.hlab <- pyogenes.full.path.hlab[(pyogenes.full.path.hlab$MHC.IC50 <= 30.0),]
thermophilus.IC50.30.hlab <- thermophilus.full.path.hlab[(thermophilus.full.path.hlab$MHC.IC50 <= 30.0),]

#- Final Table -#

aureus.results.30.hlab <- aureus.IC50.30.hlab[,c(1,2,9,10,11,12,13,14,15,8,17)]
pyogenes.results.30.hlab <- pyogenes.IC50.30.hlab[,c(1,2,9,10,11,12,13,14,15,8,17)]
thermophilus.results.30.hlab <- thermophilus.IC50.30.hlab[,c(1,2,9,10,11,12,13,14,15,8,17)]

# corrigindo cabeçalho antes de salvar! :)
colnames(aureus.results.30.hlab)
cols <- c("Peptide", "Allele", "Percentile.Rank", "Proteasome.Score", "TAP.Score",
          "MHC.Score", "Processing.Score", "Total.Score", "MHC.IC50", "Binding.Score", "Immunogenicity.Score")
colnames(aureus.results.30.hlab) <- cols
colnames(pyogenes.results.30.hlab) <- cols
colnames(thermophilus.results.30.hlab) <- cols

write.csv(aureus.results.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Results_IC50_30.csv")
write.csv(pyogenes.results.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Results_IC50_30.csv")
write.csv(thermophilus.results.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Results_IC50_30.csv")

##-- Summary --##

#- Obtaining list of alleles-#

aureus.final.allele.30.hlab <- data.frame(unique(aureus.results.30.hlab$Allele)) #221
pyogenes.final.allele.30.hlab <- data.frame(unique(pyogenes.results.30.hlab$Allele)) #225
thermophilus.final.allele.30.hlab <- data.frame(unique(thermophilus.results.30.hlab$Allele)) #200

write.csv(aureus.final.allele.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Results_HLA_IC50_30.csv")
write.csv(pyogenes.final.allele.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Results_HLA_IC50_30.csv")
write.csv(thermophilus.final.allele.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Results_HLA_IC50_30.csv")



#- Obtaining list of peptides -#

aureus.final.peptides.30.hlab <- data.frame(unique(aureus.results.30.hlab$Peptide)) #121
pyogenes.final.peptides.30.hlab <- data.frame(unique(pyogenes.results.30.hlab$Peptide)) #163
thermophilus.final.peptides.30.hlab <- data.frame(unique(thermophilus.results.30.hlab$Peptide)) #189


write.csv(aureus.final.peptides.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Results_PEP_IC50_30.csv")
write.csv(pyogenes.final.peptides.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Results_PEP_IC50_30.csv")
write.csv(thermophilus.final.peptides.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Results_PEP_IC50_30.csv")


common.aurpyo.30.hlab <- data.frame(merge(aureus.final.peptides.30.hlab, pyogenes.final.peptides.30.hlab, 
                                     by.x = colnames(aureus.final.peptides.30.hlab), by.y = colnames(pyogenes.final.peptides.30.hlab)))
common.aurthe.30.hlab <- data.frame(merge(aureus.final.peptides.30.hlab, thermophilus.final.peptides.30.hlab, 
                                     by.x = colnames(aureus.final.peptides.30.hlab), by.y = colnames(thermophilus.final.peptides.30.hlab)))
common.pyothe.30.hlab <- data.frame(merge(pyogenes.final.peptides.30.hlab, thermophilus.final.peptides.30.hlab, 
                                     by.x = colnames(pyogenes.final.peptides.30.hlab), by.y = colnames(thermophilus.final.peptides.30.hlab)))


write.csv(common.aurpyo.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Common_HLAwithPyogenes_IC50_30.csv")
write.csv(common.aurthe.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Common_HLAwithThermophilus_IC50_30.csv")
write.csv(common.pyothe.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Common_HLAwithThermophilus_IC50_30.csv")


# Qual o score de imunogenicidade pra os comuns pyo/the? 16 no total!
score.common.pyothe.30.hlab <- merge(common.pyothe.30.hlab, immunogenicity.hlab, by.x="unique.pyogenes.results.30.hlab.Peptide.", by.y = "peptide")
write.csv(score.common.pyothe.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_CommonScore_PEP_PyoThe_IC50_30.csv")


####
## LOGOS IC50 < 30nM
####


pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Aureus_Peptides_IC50_30.pdf")
ggplot() + geom_logo( aureus.final.peptides.30.hlab ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Pyogenes_Peptides_IC50_30.pdf")
ggplot() + geom_logo( pyogenes.final.peptides.30.hlab ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Thermophilus_peptides_IC50_30.pdf")
ggplot() + geom_logo( thermophilus.final.peptides.30.hlab ) + theme_logo()
dev.off()

pdf("/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_PyoThe_peptides_IC50_30.pdf")
ggplot() + geom_logo( common.pyothe.30.hlab ) + theme_logo()
dev.off()

##--    Populations   --#

setwd("/Users/martielafreitas/Desktop/HLA-Rproject/BInputs/Populacoes")
populacoes.hlab <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/BInputs/Populacoes/HLA-B_GOLD_Geral", sep=";")

# Merged tables
aureus.pop.30.hlab <-  merge(aureus.results.30.hlab, populacoes.hlab, by.x = "Allele", by.y = "Allele")
pyogenes.pop.30.hlab <-  merge(pyogenes.results.30.hlab, populacoes.hlab, by.x = "Allele", by.y = "Allele")
thermophilus.pop.30.hlab <-  merge(thermophilus.results.30.hlab, populacoes.hlab, by.x = "Allele", by.y = "Allele")

write.table(aureus.pop.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Path_Populacao_IC50_30.csv")
write.table(pyogenes.pop.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Path_Populacao_IC50_30.csv")
write.table(thermophilus.pop.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Path_Populacao_IC50_30.csv")

# Unique peptides by Population and allele

aureus.populacao.30.hlab <- aureus.pop.30.hlab[,c(1,2,13,14,15,16)]
pyogenes.populacao.30.hlab <- pyogenes.pop.30.hlab[,c(1,2,13,14,15,16)]
thermophilus.populacao.30.hlab <- thermophilus.pop.30.hlab[,c(1,2,13,14,15,16)]

aureus.ethnicity.30.hlab <- aureus.populacao.30.hlab[,c(1,2,6)]
pyogenes.ethnicity.30.hlab <- pyogenes.populacao.30.hlab[,c(1,2,6)]
thermophilus.ethnicity.30.hlab <- thermophilus.populacao.30.hlab[,c(1,2,6)]

aureus.ethnicity.unq.30.hlab <- unique(aureus.ethnicity.30.hlab)
pyogenes.ethnicity.unq.30.hlab <- unique(pyogenes.ethnicity.30.hlab)
thermophilus.ethnicity.unq.30.hlab <-unique(thermophilus.ethnicity.30.hlab)

write.table(aureus.ethnicity.unq.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Ethnicity_Unique_IC50_30.csv")
write.table(pyogenes.ethnicity.unq.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Ethnicity_Unique_IC50_30.csv")
write.table(thermophilus.ethnicity.unq.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Ethnicity_Unique_IC50_30.csv")

### PEPTIDES

aureus.sankey.30.hlab <- aureus.populacao.30.hlab[,c(1,2,6,4)]
pyogenes.sankey.30.hlab <- pyogenes.populacao.30.hlab[,c(1,2,6,4)]
thermophilus.sankey.30.hlab <- thermophilus.populacao.30.hlab[,c(1,2,6,4)]

write.table(aureus.sankey.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Sankey_PEP_IC50_30.csv")
write.table(pyogenes.sankey.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Sankey_PEP_IC50_30.csv")
write.table(thermophilus.sankey.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Sankey_PEP_IC50_30.csv")

#aureus.sankey.30 <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Sankey_HLA-B.30.csv", sep = ";")
#pyogenes.sankey.30 <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Sankey_HLA-B.30.csv", sep = ";")
#thermophilus.sankey.30 <- read.csv("/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Sankey_HLA-B.30.csv", sep = ";")


################################################################################## AUREUS

###################### AUREUS ORIENTAL 30

aureus.sankey.oriental.30.hlab <- unique(subset(aureus.sankey.30.hlab, EthnicOrigin == "Oriental"))
aureus.sankey.oriental.30.hlab <- aggregate(x = as.numeric(aureus.sankey.oriental.30.hlab$Frequency),
                                       by = c(list(aureus.sankey.oriental.30.hlab$Allele, 
                                                   aureus.sankey.oriental.30.hlab$Peptide, 
                                                   aureus.sankey.oriental.30.hlab$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlab <- subset(aureus.sankey.oriental.30.hlab, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlab$x),
                            by = c(list(aureus.sankey.top.hlab$Group.1, 
                                        aureus.sankey.top.hlab$Group.2, 
                                        aureus.sankey.top.hlab$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Aureus_Oriental_30_10up.html", selfcontained = TRUE)

rm(sankey, aureus.sankey.top.hlab, aureus.snk.agg, links, source, target, value, nodes, name, rede)

###################### BLACK 30

aureus.sankey.black.30.hlab <- unique(subset(aureus.sankey.30.hlab, EthnicOrigin == "Black"))
aureus.sankey.black.30.hlab <- aggregate(x = as.numeric(aureus.sankey.black.30.hlab$Frequency),
                                    by = c(list(aureus.sankey.black.30.hlab$Allele, 
                                                aureus.sankey.black.30.hlab$Peptide, 
                                                aureus.sankey.black.30.hlab$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlab <- subset(aureus.sankey.black.30.hlab, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlab$x),
                            by = c(list(aureus.sankey.top.hlab$Group.1, 
                                        aureus.sankey.top.hlab$Group.2, 
                                        aureus.sankey.top.hlab$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Aureus_Black_30_10up.html", selfcontained = TRUE)

rm(sankey, aureus.sankey.top.hlab, aureus.snk.agg, links, source, target, value, nodes, name, rede)


##################### CAUCASOIDE 30

aureus.sankey.cauc.30.hlab <- unique(subset(aureus.sankey.30.hlab, EthnicOrigin == "Caucasoid"))
aureus.sankey.cauc.30.hlab <- aggregate(x = as.numeric(aureus.sankey.cauc.30.hlab$Frequency),
                                   by = c(list(aureus.sankey.cauc.30.hlab$Allele, 
                                               aureus.sankey.cauc.30.hlab$Peptide, 
                                               aureus.sankey.cauc.30.hlab$EthnicOrigin)), FUN = sum)

aureus.sankey.top.hlab <- subset(aureus.sankey.cauc.30.hlab, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
aureus.snk.agg <- aggregate(x = as.numeric(aureus.sankey.top.hlab$x),
                            by = c(list(aureus.sankey.top.hlab$Group.1, 
                                        aureus.sankey.top.hlab$Group.2, 
                                        aureus.sankey.top.hlab$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Aureus_Caucasoid_30_10up.html", selfcontained = TRUE)

rm(sankey,aureus.sankey.top.hlab, aureus.snk.agg, links, source, target, value, nodes, name, rede)


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

##write.table(combined_peptides.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Aureus_Peptides_3pop_HLA-B.30.csv")
##write.table(duplicate_rows.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Aureus_Common_3pop_HLA-B.30.csv")


################################################################################ PYOGENES

###################### PYOGENES ORIENTAL 30

pyogenes.sankey.oriental.30.hlab <- unique(subset(pyogenes.sankey.30.hlab, EthnicOrigin == "Oriental"))
pyogenes.sankey.oriental.30.hlab <- aggregate(x = as.numeric(pyogenes.sankey.oriental.30.hlab$Frequency),
                                         by = c(list(pyogenes.sankey.oriental.30.hlab$Allele, 
                                                     pyogenes.sankey.oriental.30.hlab$Peptide, 
                                                     pyogenes.sankey.oriental.30.hlab$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlab <- subset(pyogenes.sankey.oriental.30.hlab, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlab$x),
                            by = c(list(pyogenes.sankey.top.hlab$Group.1, 
                                        pyogenes.sankey.top.hlab$Group.2, 
                                        pyogenes.sankey.top.hlab$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Pyogenes_Oriental_30_10up.html", selfcontained = TRUE)

rm(sankey, pyogenes.sankey.top.hlab, pyogenes.snk.agg, links, source, target, value, nodes, name, rede)

###################### PYOGENES BLACK 30

pyogenes.sankey.black.30.hlab <- unique(subset(pyogenes.sankey.30.hlab, EthnicOrigin == "Black"))
pyogenes.sankey.black.30.hlab <- aggregate(x = as.numeric(pyogenes.sankey.black.30.hlab$Frequency),
                                      by = c(list(pyogenes.sankey.black.30.hlab$Allele, 
                                                  pyogenes.sankey.black.30.hlab$Peptide, 
                                                  pyogenes.sankey.black.30.hlab$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlab <- subset(pyogenes.sankey.black.30.hlab, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlab$x),
                            by = c(list(pyogenes.sankey.top.hlab$Group.1, 
                                        pyogenes.sankey.top.hlab$Group.2, 
                                        pyogenes.sankey.top.hlab$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Pyogenes_Black_30_10up.html", selfcontained = TRUE)

rm(sankey, pyogenes.sankey.top.hlab, pyogenes.snk.agg, links, source, target, value, nodes, name, rede)

###################### PYOGENES CAUCASOID 30

pyogenes.sankey.cauc.30.hlab <- unique(subset(pyogenes.sankey.30.hlab, EthnicOrigin == "Caucasoid"))
pyogenes.sankey.cauc.30.hlab <- aggregate(x = as.numeric(pyogenes.sankey.cauc.30.hlab$Frequency),
                                     by = c(list(pyogenes.sankey.cauc.30.hlab$Allele, 
                                                 pyogenes.sankey.cauc.30.hlab$Peptide, 
                                                 pyogenes.sankey.cauc.30.hlab$EthnicOrigin)), FUN = sum)

pyogenes.sankey.top.hlab <- subset(pyogenes.sankey.cauc.30.hlab, x >= 0.10)
# exclude lines wth e-4...n
# pyogenes.sankey.10 <- pyogenes.sankey.top[- grep("e-", pyogenes.sankey.top$Frequency),]
# aggregate by frequency
pyogenes.snk.agg <- aggregate(x = as.numeric(pyogenes.sankey.top.hlab$x),
                            by = c(list(pyogenes.sankey.top.hlab$Group.1, 
                                        pyogenes.sankey.top.hlab$Group.2, 
                                        pyogenes.sankey.top.hlab$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Pyogenes_Caucasoid_30_10up.html", selfcontained = TRUE)

rm(sankey, pyogenes.sankey.top.hlab, pyogenes.snk.agg, links, source, target, value, nodes, name, rede)



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

##write.table(combined_peptides.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Pyogenes_Peptides_3pop_HLA-B.30.csv")
##write.table(duplicate_rows.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Pyogenes_Common_3pop_HLA-B.30.csv")


################################################################################ THERMOPHILUS

###################### THERMOPHILUS 30 ORIENTAL

thermophilus.sankey.oriental.30.hlab <- unique(subset(thermophilus.sankey.30.hlab, EthnicOrigin == "Oriental"))
thermophilus.sankey.oriental.30.hlab <- aggregate(x = as.numeric(thermophilus.sankey.oriental.30.hlab$Frequency),
                                             by = c(list(thermophilus.sankey.oriental.30.hlab$Allele, 
                                                         thermophilus.sankey.oriental.30.hlab$Peptide, 
                                                         thermophilus.sankey.oriental.30.hlab$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlab <- subset(thermophilus.sankey.oriental.30.hlab, x >= 0.10)
# exclude lines wth e-4...n
# thermophilus.sankey.10 <- thermophilus.sankey.top[- grep("e-", thermophilus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlab$x),
                            by = c(list(thermophilus.sankey.top.hlab$Group.1, 
                                        thermophilus.sankey.top.hlab$Group.2, 
                                        thermophilus.sankey.top.hlab$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Thermophilus_Oriental_30_10up.html", selfcontained = TRUE)

rm(sankey, thermophilus.sankey.top.hlab, thermophilus.snk.agg, links, source, target, value, nodes, name, rede)


##################### THERMOPHILUS 30 BLACK

thermophilus.sankey.black.30.hlab <- unique(subset(thermophilus.sankey.30.hlab, EthnicOrigin == "Black"))
thermophilus.sankey.black.30.hlab <- aggregate(x = as.numeric(thermophilus.sankey.black.30.hlab$Frequency),
                                          by = c(list(thermophilus.sankey.black.30.hlab$Allele, 
                                                      thermophilus.sankey.black.30.hlab$Peptide, 
                                                      thermophilus.sankey.black.30.hlab$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlab <- subset(thermophilus.sankey.black.30.hlab, x >= 0.10)
# exclude lines wth e-4...n
# aureus.sankey.10 <- aureus.sankey.top[- grep("e-", aureus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlab$x),
                            by = c(list(thermophilus.sankey.top.hlab$Group.1, 
                                        thermophilus.sankey.top.hlab$Group.2, 
                                        thermophilus.sankey.top.hlab$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Thermophilus_Black_30_10up.html", selfcontained = TRUE)

rm(sankey, thermophilus.sankey.top.hlab, thermophilus.snk.agg, links, source, target, value, nodes, name, rede)

##################### THERMOPHILUS 30 CAUCASOID

thermophilus.sankey.cauc.30.hlab <- unique(subset(thermophilus.sankey.30.hlab, EthnicOrigin == "Caucasoid"))
thermophilus.sankey.cauc.30.hlab <- aggregate(x = as.numeric(thermophilus.sankey.cauc.30.hlab$Frequency),
                                         by = c(list(thermophilus.sankey.cauc.30.hlab$Allele, 
                                                     thermophilus.sankey.cauc.30.hlab$Peptide, 
                                                     thermophilus.sankey.cauc.30.hlab$EthnicOrigin)), FUN = sum)

thermophilus.sankey.top.hlab <- subset(thermophilus.sankey.cauc.30.hlab, x >= 0.10)
# exclude lines wth e-4...n
# thermophilus.sankey.10 <- thermophilus.sankey.top[- grep("e-", thermophilus.sankey.top$Frequency),]
# aggregate by frequency
thermophilus.snk.agg <- aggregate(x = as.numeric(thermophilus.sankey.top.hlab$x),
                            by = c(list(thermophilus.sankey.top.hlab$Group.1, 
                                        thermophilus.sankey.top.hlab$Group.2, 
                                        thermophilus.sankey.top.hlab$Group.3)), FUN = sum)
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
saveNetwork(rede, "/Users/martielafreitas/Desktop/HLA-Rproject/Images/HLA-B_Thermophilus_Caucasoid_30_10up.html", selfcontained = TRUE)

rm(sankey,thermophilus.sankey.top.hlab, thermophilus.snk.agg, links, source, target, value, nodes, name, rede)


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

##write.table(combined_peptides.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Thermophilus_Peptides_3pop_HLA-B.30.csv")
##write.table(duplicate_rows.30, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/Thermophilus_Common_3pop_HLA-B.30.csv")

###############################################################################
###############################################################################
###############################################################################

# Comparações HLA-B

# Geral IC50 < 500nM
unq.aureus.pep.hlab <- unique(aureus.IC50.hlab$peptide)
unq.pyogenes.pep.hlab <- unique(pyogenes.IC50.hlab$peptide)
unq.thermophilus.pep.hlab <- unique(thermophilus.IC50.hlab$peptide)

write.table(unq.aureus.pep.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Unique_PEP_IC50_500.csv")
write.table(unq.pyogenes.pep.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Unique_PEP_IC50_500.csv")
write.table(unq.thermophilus.pep.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Unique_PEP_IC50_500.csv")

# Geral IC50 < 30nM
unq.aureus.pep.30.hlab <- unique(aureus.IC50.30.hlab$peptide)
unq.pyogenes.pep.30.hlab <- unique(pyogenes.IC50.30.hlab$peptide)
unq.thermophilu.pep.30.hlab <- unique(thermophilus.IC50.30.hlab$peptide)

write.table(unq.aureus.pep.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Unique_PEP_IC50_30.csv")
write.table(unq.pyogenes.pep.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Unique_PEP_IC50_30.csv")
write.table(unq.thermophilus.pep.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Unique_PEP_IC50_30.csv")

# Por população IC50 < 500nM

## Sp. aureus

### Black

hlab.aureus.sankey.black.500.hlab <- unique(aureus.sankey.black.hlab$Group.1)
hlab.aureus.sankey.black.30.hlab <- unique(aureus.sankey.black.30.hlab$Group.1)
hlab.aureus.sankey.black.500.pep <- unique(aureus.sankey.black.hlab$Group.2)
hlab.aureus.sankey.black.30.pep <- unique(aureus.sankey.black.30.hlab$Group.2)

write.table(hlab.aureus.sankey.black.500.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Unique_HLA_Black_IC50_500.csv")
write.table(hlab.aureus.sankey.black.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Unique_HLA_Black_IC50_30.csv")
write.table(hlab.aureus.sankey.black.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Unique_PEP_Black_IC50_500.csv")
write.table(hlab.aureus.sankey.black.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Unique_PEP_Black_IC50_30.csv")

### Caucasoid

hlab.aureus.sankey.cauc.500.hlab <- unique(aureus.sankey.caucasoid.hlab$Group.1)
hlab.aureus.sankey.cauc.30.hlab <- unique(aureus.sankey.cauc.30.hlab$Group.1)
hlab.aureus.sankey.cauc.500.pep <- unique(aureus.sankey.caucasoid.hlab$Group.2)
hlab.aureus.sankey.cauc.30.pep <- unique(aureus.sankey.cauc.30.hlab$Group.2)

write.table(hlab.aureus.sankey.cauc.500.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Unique_HLA_Cauc_IC50_500.csv")
write.table(hlab.aureus.sankey.cauc.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Unique_HLA_Cauc_IC50_30.csv")
write.table(hlab.aureus.sankey.cauc.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Unique_PEP_Cauc_IC50_500.csv")
write.table(hlab.aureus.sankey.cauc.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Unique_PEP_Cauc_IC50_30.csv")

### Oriental

hlab.aureus.sankey.oriental.500.hlab <- unique(aureus.sankey.oriental.hlab$Group.1)
hlab.aureus.sankey.oriental.30.hlab <- unique(aureus.sankey.oriental.30.hlab$Group.1)
hlab.aureus.sankey.oriental.500.pep <- unique(aureus.sankey.oriental.hlab$Group.2)
hlab.aureus.sankey.oriental.30.pep <- unique(aureus.sankey.oriental.30.hlab$Group.2)

write.table(hlab.aureus.sankey.oriental.500.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Unique_HLA_Orient_IC50_500.csv")
write.table(hlab.aureus.sankey.oriental.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Unique_HLA_Orient_IC50_30.csv")
write.table(hlab.aureus.sankey.oriental.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Unique_PEP_Orient_IC50_500.csv")
write.table(hlab.aureus.sankey.oriental.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Unique_PEP_Orient_IC50_30.csv")

## Sp. pyogenes

### Black

hlab.pyogenes.sankey.black.500.hlab <- unique(pyogenes.sankey.black.hlab$Group.1)
hlab.pyogenes.sankey.black.30.hlab <- unique(pyogenes.sankey.black.30.hlab$Group.1)
hlab.pyogenes.sankey.black.500.pep <- unique(pyogenes.sankey.black.hlab$Group.2)
hlab.pyogenes.sankey.black.30.pep <- unique(pyogenes.sankey.black.30.hlab$Group.2)

write.table(hlab.pyogenes.sankey.black.500.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Unique_HLA_Black_IC50_500.csv")
write.table(hlab.pyogenes.sankey.black.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Unique_HLA_Black_IC50_30.csv")
write.table(hlab.pyogenes.sankey.black.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Unique_PEP_Black_IC50_500.csv")
write.table(hlab.pyogenes.sankey.black.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Unique_PEP_Black_IC50_30.csv")

### Caucasoid

hlab.pyogenes.sankey.cauc.500.hla <- unique(pyogenes.sankey.cauc.hlab$Group.1)
hlab.pyogenes.sankey.cauc.30.hla <- unique(pyogenes.sankey.cauc.30.hlab$Group.1)
hlab.pyogenes.sankey.cauc.500.pep <- unique(pyogenes.sankey.cauc.hlab$Group.2)
hlab.pyogenes.sankey.cauc.30.pep <- unique(pyogenes.sankey.cauc.30.hlab$Group.2)

write.table(hlab.pyogenes.sankey.cauc.500.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Unique_HLA_Cauc_IC50_500.csv")
write.table(hlab.pyogenes.sankey.cauc.30.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Unique_HLA_Cauc_IC50_30.csv")
write.table(hlab.pyogenes.sankey.cauc.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Unique_PEP_Cauc_IC50_500.csv")
write.table(hlab.pyogenes.sankey.cauc.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Unique_PEP_Cauc_IC50_30.csv")

### Oriental

hlab.pyogenes.sankey.oriental.500.hla <- unique(pyogenes.sankey.oriental.hlab$Group.1)
hlab.pyogenes.sankey.oriental.30.hla <- unique(pyogenes.sankey.oriental.30.hlab$Group.1)

hlab.pyogenes.sankey.oriental.500.pep <- unique(pyogenes.sankey.oriental.hlab$Group.2)
hlab.pyogenes.sankey.oriental.30.pep <- unique(pyogenes.sankey.oriental.30.hlab$Group.2)

write.table(hlab.pyogenes.sankey.oriental.500.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Unique_HLA_Orient_IC50_500.csv")
write.table(hlab.pyogenes.sankey.oriental.30.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Unique_HLA_Orient_IC50_30.csv")
write.table(hlab.pyogenes.sankey.oriental.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Unique_PEP_Orient_IC50_500.csv")
write.table(hlab.pyogenes.sankey.oriental.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Unique_PEP_Orient_IC50_30.csv")

## Sp. thermophilus
### Black
hlab.thermophilus.sankey.black.500.hlab <- unique(thermophilus.sankey.black.hlab$Group.1)
hlab.thermophilus.sankey.black.30.hlab <- unique(thermophilus.sankey.black.30.hlab$Group.1)

hlab.thermophilus.sankey.black.500.pep <- unique(thermophilus.sankey.black.hlab$Group.2)
hlab.thermophilus.sankey.black.30.pep <- unique(thermophilus.sankey.black.30.hlab$Group.2)

write.table(hlab.thermophilus.sankey.black.500.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Unique_HLA_Black_IC50_500.csv")
write.table(hlab.thermophilus.sankey.black.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Unique_HLA_Black_IC50_30.csv")
write.table(hlab.thermophilus.sankey.black.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Unique_PEP_Black_IC50_500.csv")
write.table(hlab.thermophilus.sankey.black.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Unique_PEP_Black_IC50_30.csv")

### Caucasoid
hlab.thermophilus.sankey.cauc.500.hlab <- unique(thermophilus.sankey.cauc.hlab$Group.1)
hlab.thermophilus.sankey.cauc.30.hlab <- unique(thermophilus.sankey.cauc.30.hlab$Group.1)

hlab.thermophilus.sankey.cauc.500.pep <- unique(thermophilus.sankey.cauc.hlab$Group.2)
hlab.thermophilus.sankey.cauc.30.pep <- unique(thermophilus.sankey.cauc.30.hlab$Group.2)

write.table(hlab.thermophilus.sankey.cauc.500.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Unique_HLA_Cauc_IC50_500.csv")
write.table(hlab.thermophilus.sankey.cauc.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Unique_HLA_Cauc_IC50_30.csv")
write.table(hlab.thermophilus.sankey.cauc.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Unique_PEP_Cauc_IC50_500.csv")
write.table(hlab.thermophilus.sankey.cauc.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Unique_PEP_Cauc_IC50_30.csv")

### Oriental
hlab.thermophilus.sankey.oriental.500.hla <- unique(thermophilus.sankey.oriental.hlab$Group.1)
hlab.thermophilus.sankey.oriental.30.hla <- unique(thermophilus.sankey.oriental.30.hlab$Group.1)

hlab.thermophilus.sankey.oriental.500.pep <- unique(thermophilus.sankey.oriental.hlab$Group.2)
hlab.thermophilus.sankey.oriental.30.pep <- unique(thermophilus.sankey.oriental.30.hlab$Group.2)

write.table(hlab.thermophilus.sankey.oriental.500.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Unique_HLA_Orient_IC50_500.csv")
write.table(hlab.thermophilus.sankey.oriental.30.hla, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Unique_HLA_Orient_IC50_30.csv")
write.table(hlab.thermophilus.sankey.oriental.500.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Unique_PEP_Orient_IC50_500.csv")
write.table(hlab.pyogenes.sankey.oriental.30.pep, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Unique_PEP_Orient_IC50_30.csv")

###############################################################################
###############################################################################
###############################################################################

## For HLA-Brena 500
aureus.arena.hlab <- aureus.ethnicity.unq.hlab[,c(1,2)]
aureus.arena.unq.hlab <- unique(aureus.arena.hlab)
aureus.arena.hla.unq.hlab <- unique(aureus.arena.hlab$Allele)
aureus.arena.pep.unq.hlab <- unique(aureus.arena.hlab$Peptide)

write.table(aureus.arena.unq.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Arena_Unq_HLA-B.csv")
write.table(aureus.arena.hla.unq.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Arena_HlaUnq_HLA-B.csv")
write.table(aureus.arena.pep.unq.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Arena_PepUnq_HLA-B.csv")

pyogenes.arena.hlab <- pyogenes.ethnicity.unq.hlab[,c(1,2)]
pyogenes.arena.unq.hlab <- unique(pyogenes.arena.hlab)
pyogenes.arena.hla.unq.hlab <- unique(pyogenes.arena.hlab$Allele)
pyogenes.arena.pep.unq.hlab <- unique(pyogenes.arena.hlab$Peptide)

write.table(pyogenes.arena.unq.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Arena_Unq_HLA-B.csv")
write.table(pyogenes.arena.hla.unq.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Arena_HlaUnq_HLA-B.csv")
write.table(pyogenes.arena.pep.unq.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Arena_PepUnq_HLA-B.csv")

thermophilus.arena.hlab <- thermophilus.ethnicity.unq.hlab[,c(1,2)]
thermophilus.arena.unq.hlab <- unique(thermophilus.arena.hlab)
thermophilus.arena.hla.unq.hlab <- unique(thermophilus.arena.hlab$Allele)
thermophilus.arena.pep.unq.hlab <- unique(thermophilus.arena.hlab$Peptide)

write.table(thermophilus.arena.unq.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Arena_Unq_HLA-B.csv")
write.table(thermophilus.arena.hla.unq.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Arena_HlaUnq_HLA-B.csv")
write.table(thermophilus.arena.pep.unq.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Arena_PepUnq_HLA-B.csv")

## For HLA-Brena 30
aureus.arena.30.hlab <- aureus.ethnicity.unq.30.hlab[,c(1,2)]
aureus.arena.unq.30.hlab <- unique(aureus.arena.30.hlab)
aureus.arena.hla.unq.30.hlab <- unique(aureus.arena.30.hlab$Allele)
aureus.arena.pep.unq.30.hlab <- unique(aureus.arena.30.hlab$Peptide)

write.table(aureus.arena.unq.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Arena_Unq_30_HLA-B.csv")
write.table(aureus.arena.hla.unq.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Arena_HlaUnq_30_HLA-B.csv")
write.table(aureus.arena.pep.unq.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Aureus_Arena_PepUnq_30_HLA-B.csv")

pyogenes.arena.30.hlab <- pyogenes.ethnicity.unq.30.hlab[,c(1,2)]
pyogenes.arena.unq.30.hlab <- unique(pyogenes.arena.30.hlab)
pyogenes.arena.hla.unq.30.hlab <- unique(pyogenes.arena.30.hlab$Allele)
pyogenes.arena.pep.unq.30.hlab <- unique(pyogenes.arena.30.hlab$Peptide)

write.table(pyogenes.arena.unq.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Arena_Unq_30_HLA-B.csv")
write.table(pyogenes.arena.hla.unq.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Arena_HlaUnq_30_HLA-B.csv")
write.table(pyogenes.arena.pep.unq.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Pyogenes_Arena_PepUnq_30_HLA-B.csv")

thermophilus.arena.30.hlab <- thermophilus.ethnicity.unq.30.hlab[,c(1,2)]
thermophilus.arena.unq.30.hlab <- unique(thermophilus.arena.30.hlab)
thermophilus.arena.hla.unq.30.hlab <- unique(thermophilus.arena.30.hlab$Allele)
thermophilus.arena.pep.unq.30.hlab <- unique(thermophilus.arena.30.hlab$Peptide)

write.table(thermophilus.arena.unq.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Arena_Unq_30_HLA-B.csv")
write.table(thermophilus.arena.hla.unq.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Arena_HlaUnq_30_HLA-B.csv")
write.table(thermophilus.arena.pep.unq.30.hlab, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-B_Thermophilus_Arena_PepUnq_30_HLA-B.csv")

