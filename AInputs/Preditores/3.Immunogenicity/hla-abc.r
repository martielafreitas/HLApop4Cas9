
# Comuns para Grafico da Proteína (Cobrinha)

# Peptideos únicos abaixo de 30nM, AUREUS
AA <- unique(aureus.arena.30.hlaa$Peptide)
AB <- unique(aureus.arena.30.hlab$Peptide)
AC <- unique(aureus.arena.30.hlac$Peptide)

write.csv(AA, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_IC50_30_AA.csv")
write.csv(AB, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_IC50_30_AB.csv")
write.csv(AC, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_IC50_30_AC.csv")

# Peptideos únicos abaixo de 30nM, PYOGENES
PA <- unique(pyogenes.arena.30.hlaa$Peptide)
PB <- unique(pyogenes.arena.30.hlab$Peptide)
PC <- unique(pyogenes.arena.30.hlac$Peptide)

write.csv(PA, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_IC50_30_PA.csv")
write.csv(PB, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_IC50_30_PB.csv")
write.csv(PC, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_IC50_30_PC.csv")

# Peptideos únicos abaixo de 30nM, THERMOPHILUS
TA <- unique(thermophilus.arena.30.hlaa$Peptide)
TB <- unique(thermophilus.arena.30.hlab$Peptide)
TC <- unique(thermophilus.arena.30.hlac$Peptide)

write.csv(TA, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_IC50_30_TA.csv")
write.csv(TB, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_IC50_30_TB.csv")
write.csv(TC, file="/Users/martielafreitas/Desktop/HLA-Rproject/Results/HLA-A_Aureus_IC50_30_TC.csv")


AAB <-  merge(as.data.frame(AA), as.data.frame(AB), by.x = "AA", by.y = "AB")
AAC <-  merge(as.data.frame(AA), as.data.frame(AC), by.x = "AA", by.y = "AC")
ACB <-  merge(as.data.frame(AB), as.data.frame(AC), by.x = "AB", by.y = "AC")

PABC <-  merge(as.data.frame(PA), as.data.frame(PB), as.data.frame(PC), by.x = "PA", by.y = "PB", by.z="PC")
TABC <-  merge(as.data.frame(TA), as.data.frame(TB), as.data.frame(TC), by.x = "TA", by.y = "TB", by.z="TC")

AABC <-  merge(as.data.frame(PA), as.data.frame(PB), as.data.frame(PC), by.x = "PA", by.y = "PB", by.z="PC")
AABC <-  merge(as.data.frame(AA), as.data.frame(AB), as.data.frame(AC), by.x = "AA", by.y = "AB", by.z="AC")
PABC <-  merge(as.data.frame(PA), as.data.frame(PB), as.data.frame(PC), by.x = "PA", by.y = "PB", by.z="PC")
TABC <-  merge(as.data.frame(TA), as.data.frame(TB), as.data.frame(TC), by.x = "TA", by.y = "TB", by.z="TC")

###############################################################################
###############################################################################
###############################################################################

aureus.pep.hlac <- read.csv("/Users/martielafreitas/Documents/Rprojects/HLA-C/Results/Aureus_Peptides_3pop_HLA-C.csv", sep = " ")
pyogenes.pep <- read.csv("/Users/martielafreitas/Documents/Rprojects/HLA-C/Results/Pyogenes_Peptides_3pop_HLA-C.csv", sep = " ")
thermophilus.pep <- read.csv("/Users/martielafreitas/Documents/Rprojects/HLA-C/Results/Thermophilus_Peptides_3pop_HLA-C.csv", sep = " ")

common_pyo_thermo <- data.frame(merge(pyogenes.pep, thermophilus.pep, by = "x"))
common_pyo_thermo


aureus.path <- read.csv("/Users/martielafreitas/Documents/Rprojects/HLA-C/Results/Aureus_Path_Populacao_HLA-C.csv", sep = " ")
pyogenes.path <- read.csv("/Users/martielafreitas/Documents/Rprojects/HLA-C/Results/Pyogenes_Path_Populacao_HLA-C.csv", sep = " ")
thermophilus.path <- read.csv("/Users/martielafreitas/Documents/Rprojects/HLA-C/Results/Thermophilus_Path_Populacao_HLA-C.csv", sep = " ")

aureus.pep
pyogenes.pep
thermophilus.pep

unique(aureus.path$Peptide)
unique(pyogenes.path$Peptide)
unique(thermophilus.path$Peptide)

pyogenes.pep.path <-pyogenes.path[(pyogenes.path$Peptide %in% c('GLDIGTNSV', 'GLYETRIDL', 'GTNSVGWAV', 'HAHDAYLNA', 'ILTFRIPYY', 'KILTFRIPY', 'KQRTFDNGS', 'LTFRIPYYV', 'LVETRQITK', 'RIPYYVGPL', 'RLKRTARRR', 'RTARRRYTR', 'TARRRYTRR')),]
pyogenes.pep.path.menor <- unique(pyogenes.pep.path[,c(2,8,9,10,11,16)])

pyogenes.pep.path.menor <- pyogenes.pep.path.menor[(pyogenes.pep.path.menor$EthnicOrigin %in% c('Caucasoid', 'Black', 'Oriental')),]


pyogenes.path[(pyogenes.path$Peptide %in% c('GLDIGTNSV', 'GLYETRIDL', 'GTNSVGWAV', 'HAHDAYLNA', 'ILTFRIPYY', 'KILTFRIPY', 'KQRTFDNGS', 'LTFRIPYYV', 'LVETRQITK', 'RIPYYVGPL', 'RLKRTARRR', 'RTARRRYTR', 'TARRRYTRR')),]

thermophilus.pep.path <-thermophilus.path[(thermophilus.path$Peptide %in% c('GLDIGTNSV', 'GLYETRIDL', 'GTNSVGWAV', 'HAHDAYLNA', 'ILTFRIPYY', 'KILTFRIPY', 'KQRTFDNGS', 'LTFRIPYYV', 'LVETRQITK', 'RIPYYVGPL', 'RLKRTARRR', 'RTARRRYTR', 'TARRRYTRR')),]
thermophilus.pep.path <-thermophilus.path[(thermophilus.path$Peptide %in% c('GLDIGTNSV', 'GLYETRIDL', 'GTNSVGWAV', 'HAHDAYLNA', 'ILTFRIPYY', 'KILTFRIPY', 'KQRTFDNGS', 'LTFRIPYYV', 'LVETRQITK', 'RIPYYVGPL', 'RLKRTARRR', 'RTARRRYTR', 'TARRRYTRR')),]
