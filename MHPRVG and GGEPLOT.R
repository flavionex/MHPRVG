
# Setting Work Directory --------------------------------------------------


setwd("C:\\Users\\arantefl\\Desktop\\Rice_data_2021")
getwd()
Dados<- read.csv("C:\\Users\\arantefl\\Desktop\\Rice_data_2021\\VCU SC CL.csv")#***Mudar ao Analisar

#Dados<-Dados[complete.cases(Dados) & Dados$Plot.Discarded!="Yes",]

#View(Dados)

###
# Observacoes
####
# 1. AMB substituido por Book.Name
# 2. DAP substituido por Yield.kg.ha
# 3. PROG substituido por Material.Name
# 4. REP substituido por RepNo

str(Dados)
summary(Dados)

# Analise exploratoria ----------------------------------------------------

resp<- Dados$Yield.kg.ha ###***Mudar ao analisar caractere diferente

#View(Dados) 
data<-resp
expAnalysis<- function(data, digits=2){
  NAs <- sum(is.na(data))
  data<- data[!is.na(data)]
  exp<- round(c(mean(data), sd(data), min(data), max(data), length(data)),digits)
  exp<- round(c(exp,(exp[2]/exp[1])*100),digits)
  if(any(NAs))
    exp<- c(exp,NAs)
  else exp<- c(exp,0)
  names(exp)<- c("Mean", "sd","Min", "Max","Length", "CV%", "NAs")
  return(exp)
  }

expAnalysis(resp)

# Analise Exploratoria por Local ------------------------------------------
envMeans<- aggregate(resp, by=list(Dados$Book.Name), FUN=expAnalysis)
envMeans

# Estimativa de Parametros Geneticos --------------------------------------


str(Dados)

Dados$Book.Name<- as.factor(Dados$Book.Name)
Dados$Material.Name<- as.factor(Dados$Material.Name)
Dados$RepNo<- as.factor(Dados$RepNo)

str(Dados)

library(lme4)

### parcela Linear

nREP<- nlevels(Dados$RepNo) #numero de blocos
#nPP<- 1 # numero de plantas por parcela caso houver
nAMB<- nlevels(Dados$Book.Name) #Numero de ambientes

m1<- lmer(resp ~ Book.Name + Book.Name:RepNo + (1|Material.Name) + (1|Book.Name:Material.Name) , data=Dados) #modelo de intercpto aleatorio

summary(m1)

vc<- as.data.frame(VarCorr(m1),comp="Variance")
vc<- subset(vc,select = c("grp","vcov"))
vc
Media<- mean(resp,na.rm=T)
std<- sd(resp,na.rm=T)


vGxE = vc[1,2]
vProg = vc[2,2]
vRes = vc[3,2]
#teria um VParc da parcela

vPhen = vGxE + vProg + vRes #aqui somaria o vParc
h2a = 4*vProg/vPhen
# h2m = vProg/(vProg+(vGxE/nAMB)+(vParc/nREP)+(vRes/(nAMB+nREP+nPP)))  #se tivesse o VParc
h2m = vProg/(vProg+(vGxE/nAMB)+(vRes/(nAMB+nREP)))
rgg = sqrt(h2m)
c2ge = vGxE/vPhen
rgloc = h2a / (h2a+c2ge) #correlacao genotipica entre ambientes. Quanto mais alto menos o ambiente esta interferindo no genotipo acima de 0.7 já é alto
CVgi = (sqrt(4*vProg)/Media)*100
CVg = (sqrt(vProg)/Media)*100
CVe = (sqrt(vRes)/Media)*100

Parametros<- c( "vProg", "vGxE", "vRes", "vPhen", "h2a", "h2m", "rgg", "c2ge", "rgloc", "CVgi", "CVg", "CVe", "Media", "std" )
Valores<- round(c( vProg, vGxE, vRes, vPhen, h2a, h2m, rgg, c2ge, rgloc, CVgi, CVg, CVe, Media, std),3)
Paramdf<- data.frame(Parametros,Valores)
Paramdf

########################
# BLUP's Todos os locais
########################
library(tidyverse) # para usar a funcao rownames_to_column

#essa funcao nao precisaria ser criada.

BLUP<- function(data){
  Bl<- rownames_to_column(data)
  names(Bl)<- c("Genotipo", "g")
  Bl$"g+u"<- Bl$g + mean(resp,na.rm=T)
  arrange(Bl, desc(Bl$"g+u"))
}

Blup_ord<- BLUP(ranef(m1)$Material.Name)
Blup_ord


################################
# BLUPS - Por local
###############################

BLUPge<- ranef(m1)$"Book.Name:Material.Name"
BLUPge<- rownames_to_column(BLUPge)
BLUPge<- separate(BLUPge, rowname, into = c("Book.Name", "Genotipo"),sep=":")
colnames(BLUPge)<- c("Ambiente","Genotipo", "ge")

Book.Name<- levels(Dados$Book.Name)
BLUPge_ord<- list()

for(i in Book.Name){
  Blge<- BLUPge[BLUPge$Ambiente == i,]  #extrai somente o ambiente i
  Blge<- merge(Blge, Blup_ord, by="Genotipo", )  #merge uni dataframes de acordo com a coluna # *** tirei o "all=T" desta linha; Blge<- merge(Blge, Blup_ord, by="Genotipo", all=T) 
  Blege<- drop_na(Blge)                                 #caso exista Na ele tira
  u<- mean(resp[Dados$Book.Name==i], na.rm=T)                 #media do ambiente i
  Blge$"g+ge"<- Blge$g + Blge$ge                   
  Blge$"g+ge+u"<- Blge$"g+ge"+u
  Blge<- Blge[,-c(3,4)]                            #tirar a coluna 3 e 4                       
  Blge<- arrange(Blge,desc(Blge$"g+ge"))           #colocando em order decrescente
  BLUPge_ord[[i]]<- Blge                           #Adicionando na lista que criamos.
}

BLUPge_ord

##############################
# Estabilidade de valores geneticos (MHVG)
##############################
library(plyr)
library(dplyr)

MHVG_data<- ldply(BLUPge_ord, data.frame)[,-1]

MHVG<- function(data){
  MHVGm<- aggregate(1/data$g.ge.u, list(data$Genotipo), FUN = sum, na.rm = TRUE)
  colnames(MHVGm)<- c("Genotipo", "MHVGm")
  n<- table(data$Genotipo)
  n<- data.frame(n)
  MHVGm$MHVG <- n$Freq/MHVGm$MHVGm
  MHVGm<- MHVGm[,-2]
  arrange(MHVGm, desc(MHVGm$MHVG))
}

MHVG_ord<- MHVG(MHVG_data)

MHVG_ord

#############################
# Adaptabilidade de Valores Geneticos (PRVG)
############################

PRVG_data<- list()

for(i in Book.Name){
  u<- mean(resp[Dados$Book.Name==i], na.rm=T)
  dt<- BLUPge_ord[[i]]
  dt$"VGij/VGj"<- dt$"g+ge+u"/u
  PRVG_data[[i]]<- dt
}

PRVG_data<- ldply(PRVG_data,data.frame)[,-1]

PRVG<- function(data){
  PRVG<- aggregate(data$"VGij.VGj", list(data$Genotipo), FUN = mean, na.rm = T)
  colnames(PRVG)<- c("Genotipo","PRVG")
  PRVG$PRVG<- PRVG$PRVG*(mean(resp, na.rm=T))
  arrange(PRVG, desc(PRVG))
}

PRVG_ord<- PRVG(PRVG_data)

PRVG_ord

################################
# Estabilidade e adaptabilidade de valores geneticos (MHPRVG)
################################

MHPRVG<- function(data){
  data$"1/PRVG"<-1/data$"VGij.VGj"
  MHPRVG<- aggregate(data$"1/PRVG", list(data$Genotipo), FUN = sum, na.rm = T)
  colnames(MHPRVG)<- c("Genotipo", "MHPRVG")
  n<- table(data[,1])
  n<- data.frame(n)
  MHPRVG$MHPRVG<- n$Freq/MHPRVG$MHPRVG
  MHPRVG$MHPRVG<- MHPRVG$MHPRVG*(mean(resp, na.rm=T))
  arrange(MHPRVG, desc(MHPRVG))
}

MHPRVG_ord<- MHPRVG(PRVG_data)

MHPRVG_ord

#######################
# All Results from the Methodology Improve the results stopped here
######################



#####################################
# Output de resultados
###################################

getwd()

dir.create("Analises")
setwd("Analises")

arquivo<- "TESTING_OUTPUT_VCU SC CL"#***Mudar ao rodar a analise

if(file.exists(arquivo)){
  file.remove(arquivo)
}

sink(arquivo, append = TRUE)

cat("\n--------------------------------------------------------------------------------------------\n")
cat("                                       Analise Genotipo x Ambiente \n")
cat("--------------------------------------------------------------------------------------------\n\n")
cat("---------------------------------------Analise exploratoria---------------------------------\n\n")
print(envMeans)
cat("\n")
cat("--------------------------------------Parametros Geneticos----------------------------------\n\n")
print(Paramdf)
cat("\n")
cat("--------------------------------------BLUP - Todos os locais--------------------------------\n\n")
print(Blup_ord)
cat("\n")
cat("--------------------------------------BLUP - Por local--------------------------------------\n\n")
print(BLUPge_ord)
cat("\n")
cat("---------------------------------------------MHVG-------------------------------------------\n\n")
print(MHVG_ord)
cat("\n")
cat("-------------------------------------------PRVG---------------------------------------------\n\n")
print(PRVG_ord)
cat("\n")
cat("-------------------------------------------MHPRVG-------------------------------------------\n\n")
print(MHPRVG_ord)

sink()

##############################
# GGE Biplot
##############################

#library(reshape2)

#meanDados<- aggregate(subset(Dados, select = "Yield.kg.ha"), list(Dados$Material.Name,Dados$Book.Name), FUN = mean,  na.rm=T)

#GE_matrix<- as.matrix(acast(meanDados, Group.1~Group.2, value.var = "Yield.kg.ha"))

#install.packages("GGEBiplots")
#library(GGEBiplots)

#rowGE<- GGEModel(GE_matrix, SVP= "row")
#colGE<- GGEModel(GE_matrix, SVP= "column")

#Which won foca mais nos ambientes baiscamente atribuido; Linhas verdes formam os mega ambientes. Genotipos , mais longe da origem ou seja nos vertices, estap mais longe da media.. isso pode ser bom ou ruim.genotipos no vertice do Q1 sao ruins.
#Which_won<- WhichWon(colGE, sizeGen=3, colEnv="dark red", sizeEnv=4, colGen="black", colSegment= "forest green")


#Eixo x em vermelhor eixo do ambiente media (EAM)  a seta nos diz a direcao que um genotipo tem que estar pra ser mais produtivo
# a linha vertical representa a media... genotipos a esquerda dela, menos produtivo ... genotipos a direita , superiores a media
# O vetor trasejado indica a estabilidade do genotipo. Quanto mais longo, mas instavel ele é.
#Mean_stab<- MeanStability(rowGE, axis_expand = 1.03, sizeGen=4, colEnv= "dark red", sizeEnv = 4, colGen= "black")

# bi plot baseado nas difrencas dos ambientes. PLot discriminativo. o genotipo nao importa tanto. 
#Disc_rep<- DiscRep(colGE, axis_expand = 1.03, sizeGen = 4, colEnv= "dark red", sizeEnv=4, colGen= "black")

# seta : genotipos mais proximos do centro sao os genotipos ideiais.
#Rank_gen<- RankGen(rowGE, axis_expand = 1.03, sizeGen = 4, colEnv= "dark red", sizeEnv=4, colGen= "black")

#dir.create("GGEBiplot")

