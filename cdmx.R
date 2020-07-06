#### Database management ####
library(dummies);library(readr); library(tidyverse); library(survival); library(mediation); library(ggpubr); library(rms); library(CAinterprTools)
library(survminer); library(haven); library(rsq); library(ResourceSelection); library(ggsci);library(timereg); library(coxme); library(FactoMineR)
library(pROC);library(sf); library(rgdal); library(ggsci); library(ggmap); library(scales); library(jtools); library(cowplot); library(factoextra)
library(ggstance); library(flextable); library(simPH); library(ggthemes); library(lme4); library(lmerTest); library(prismatic); library(viridis)

setwd("C:/Users/HP-PC/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/COVID-CDMX")
setwd("C:/Users/Usuario Lab. Datos/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/COVID-CDMX")
cdmx<-read.csv("base-covid-sinave.csv", na.strings = c("", " ", "SE IGNORA","<NA>"),stringsAsFactors = FALSE)
cdmx$covid[cdmx$resultado_definitivo=="SARS-CoV-2"]<-1;cdmx$covid[cdmx$resultado_definitivo=="NEGATIVO"]<-0
cdmx[,c(32:61)][cdmx[,c(32:61)]=="NO"]<-0;cdmx[,c(32:61)][cdmx[,c(32:61)]=="SI"]<-1
cdmx[,c(32:61)][is.na(cdmx[,c(32:61)])]<-0
cdmx[,c(32:61)]<-sapply(cdmx[,c(32:61)], as.numeric)
cdmx$intubado[cdmx$intubado=="NO" | is.na(cdmx$intubado)]<-0;cdmx$intubado[cdmx$intubado=="SI"]<-1
cdmx$uci[cdmx$unidad_cuidados_intensivos=="NO" | is.na(cdmx$unidad_cuidados_intensivos)]<-0;cdmx$uci[cdmx$unidad_cuidados_intensivos=="SI"]<-1
cdmx$tipo_paciente[cdmx$tipo_paciente=="AMBULATORIO" | is.na(cdmx$tipo_paciente)]<-0;cdmx$tipo_paciente[cdmx$tipo_paciente=="HOSPITALIZADO"]<-1
cdmx$contacto_infeccion_viral[cdmx$contacto_infeccion_viral=="NO" | is.na(cdmx$contacto_infeccion_viral)]<-0;cdmx$contacto_infeccion_viral[cdmx$contacto_infeccion_viral=="SI"]<-1
cdmx$tipo_paciente<-as.numeric(cdmx$tipo_paciente)
cdmx$diagnostico_clinico_neumonia[cdmx$diagnostico_clinico_neumonia=="NO" | is.na(cdmx$diagnostico_clinico_neumonia)]<-0;cdmx$diagnostico_clinico_neumonia[cdmx$diagnostico_clinico_neumonia=="SI"]<-1
cdmx$fecha_defuncion<-as.character(cdmx$fecha_defuncion);cdmx$fecha_defuncion[cdmx$fecha_defuncion==""]<-NA
cdmx$numero_sintomas<-cdmx$fiebre+cdmx$tos+cdmx$disnea+cdmx$irritabilidad+cdmx$odinofagia+cdmx$diarrea+
  cdmx$dolor_toracico+cdmx$calofrios+cdmx$cefalea+cdmx$mialgias+cdmx$artralgias+cdmx$ataque_al_estado_general+
  cdmx$rinorrea+cdmx$polipnea+cdmx$vomito+cdmx$conjuntivitis+cdmx$cianosis+cdmx$inicio_subito_sintomas
cdmx$sintomas_resp<-cdmx$fiebre+cdmx$tos+cdmx$disnea+cdmx$cianosis+cdmx$polipnea
cdmx$asint[cdmx$sintomas_resp>0]<-0;cdmx$asint[cdmx$sintomas_resp==0]<-1
cdmx$any_sy[(cdmx$irritabilidad+cdmx$odinofagia+cdmx$diarrea+
               cdmx$dolor_toracico+cdmx$calofrios+cdmx$cefalea+cdmx$mialgias+cdmx$artralgias+cdmx$ataque_al_estado_general+
               cdmx$rinorrea+cdmx$conjuntivitis)>0]<-1;cdmx$any_sy[is.na(cdmx$any_sy)]<-0
cdmx$only_resp[cdmx$sintomas_resp>0 & cdmx$any_sy==0]<-1;cdmx$only_resp[!(cdmx$sintomas_resp>0 & cdmx$any_sy==0)]<-0
cdmx$comorb<-cdmx$diabetes+cdmx$hipertension+cdmx$enfermedad_cardiaca+cdmx$insuficiencia_renal_cronica+
  cdmx$asma+cdmx$VIH_SIDA+cdmx$epoc+cdmx$obesidad+cdmx$tabaquismo+cdmx$otra_condicion
cdmx$asintomaticos[cdmx$numero_sintomas>0]<-0;cdmx$asintomaticos[cdmx$numero_sintomas==0]<-1
cdmx$FU_time<-as.Date(cdmx$fecha_defuncion)-as.Date(cdmx$fecha_inicio_sintomas)
cdmx$fecha_inicio_sintomas<-as.Date(cdmx$fecha_inicio_sintomas)
cdmx$FU_time[is.na(cdmx$FU_time)]<-(as.Date(cdmx$fecha_ingreso)-as.Date(cdmx$fecha_inicio_sintomas));cdmx$FU_time<-as.numeric(cdmx$FU_time)
cdmx$FU_admission<-as.numeric(as.Date(cdmx$fecha_ingreso)-as.Date(cdmx$fecha_inicio_sintomas))
cdmx$sintomas<-cdmx$asint+cdmx$asintomaticos
cdmx$sintomas<-factor(cdmx$sintomas, labels = c("Symptomatic","Non-respiratory", "Asymptomatic"))
cdmx$sintomas2<-cdmx$asint*2+cdmx$asintomaticos+cdmx$only_resp
cdmx$sintomas2<-factor(cdmx$sintomas2, labels = c("Both","Only respiratory","Only non-respiratory", "Asymptomatic"))
cdmx$sintomas<-relevel(cdmx$sintomas, ref="Asymptomatic")
cdmx$sintomas2<-relevel(cdmx$sintomas2, ref="Asymptomatic")
cdmx$edad40[cdmx$edad<40]<-1;cdmx$edad40[!cdmx$edad<40]<-0
cdmx$edad60[cdmx$edad>60]<-1;cdmx$edad60[!cdmx$edad>60]<-0
cdmx$edad_cat[cdmx$edad<40]<-1;cdmx$edad_cat[cdmx$edad>=40 & cdmx$edad<60]<-2;cdmx$edad_cat[cdmx$edad>=60]<-3
cdmx$edad_cat<-factor(cdmx$edad_cat, labels = c("<40y", "40-60y", ">60y"))
cdmx$edad_cat<-relevel(cdmx$edad_cat, ref = "<40y")
d1<-dummy(cdmx$evolucion_caso)
cdmx<-cbind(cdmx, d1)
d2<-dummy(cdmx$ocupacion)
cdmx<-cbind(cdmx, d2)
cdmx$pers_sal<-cdmx$MEDICOS+cdmx$ENFERMERAS+cdmx$DENTISTAS+cdmx$`OTROS TRABAJADORES DE LA SALUD`+
  cdmx$LABORATORISTAS;cdmx$pers_sal<-na.tools::na.replace(cdmx$pers_sal,0)
cdmx$irc_has<-cdmx$hipertension+cdmx$insuficiencia_renal_cronica*2


cdmx1<- cdmx %>% filter(resultado_definitivo=="SARS-CoV-2")

#### Sin sintomas respiratorios ####

## Desenlaces sin sintomas respiratorios ##
table(cdmx1$asint, cdmx1$tipo_paciente) %>% prop.table(1)
table(cdmx1$asint, cdmx1$diagnostico_clinico_neumonia) %>% prop.table(1)
table(cdmx1$asint, cdmx1$intubado) %>% prop.table(1)
table(cdmx1$asint, cdmx1$uci) %>% prop.table(1)
table(cdmx1$asint, cdmx1$DEFUNCION) %>% prop.table(1)
table(cdmx1$asint, cdmx1$tipo_paciente, cdmx1$DEFUNCION)
table(cdmx1$asint, cdmx1$evolucion_caso)

tapply(cdmx1$FU_admission, cdmx1$sintomas, quantile)
dunn.test::dunn.test(cdmx1$FU_admission, cdmx1$sintomas, kw = T)

##Modelo hospitalizacion
resp<- cdmx1%>%filter(asint==1)

m1<-coxph(Surv(FU_admission, tipo_paciente)~vomito+conjuntivitis+dolor_abdominal+rinorrea+enfermedad_cardiaca+inmunosupresivo+diabetes+
            hipertension+contacto_infeccion_viral+pers_sal+frailty(cve_entidad_unidad_medica), data=resp)
summary(m1)

## Modelo mortalidad asintomáticos ##
m1<-coxph(Surv(FU_time, DEFUNCION)~edad+diabetes+epoc+insuficiencia_renal_cronica+contacto_infeccion_viral+any_sy+vomito+frailty(cve_entidad_unidad_medica), data=resp)
summary(m1)

table(cdmx1$asintomaticos)%>%prop.table()
#### Asintomaticos ####

table(cdmx1$asintomaticos, cdmx1$tipo_paciente) %>% prop.table(1)
table(cdmx1$asintomaticos, cdmx1$uci) %>% prop.table(1)
table(cdmx1$asintomaticos, cdmx1$intubado) %>% prop.table(1)
table(cdmx1$asintomaticos, cdmx1$DEFUNCION) %>% prop.table(1)
table(cdmx1$asintomaticos, cdmx1$tipo_paciente, cdmx1$DEFUNCION)
table(cdmx1$asintomaticos, cdmx1$evolucion_caso)

## Modelo hospitalizaci?n ##
resp<- cdmx1%>%filter(asintomaticos==1)
m1<-coxph(Surv(FU_admission, tipo_paciente)~edad+enfermedad_cardiaca+pers_sal+contacto_infeccion_viral+
            insuficiencia_renal_cronica+diabetes+frailty(cve_entidad_unidad_medica), data=resp)
summary(m1)

## Modelo mortalidad asintomáticos ##
m1<-coxph(Surv(FU_time, DEFUNCION)~edad60+contacto_infeccion_viral+insuficiencia_renal_cronica+frailty(cve_entidad_unidad_medica), data=resp)
summary(m1)
cox.zph(m1)

table(cdmx1$asintomaticos, cdmx1$evolucion_caso)
#### Con sintomas no respiratorios ####
cdmx1<- cdmx %>% filter(resultado_definitivo=="SARS-CoV-2")
d3<-dummy(cdmx1$sintomas)
cdmx1<-cbind(cdmx1, d3)
table(cdmx1$`Non-respiratory`)%>%prop.table()
## Desenlaces sin sintomas respiratorios ##
table(cdmx1$`sintomasNon-respiratory`, cdmx1$tipo_paciente) %>% prop.table(1)
table(cdmx1$`sintomasNon-respiratory`, cdmx1$unidad_cuidados_intensivos) %>% prop.table(1)
table(cdmx1$`sintomasNon-respiratory`, cdmx1$intubado) %>% prop.table(1)
table(cdmx1$`sintomasNon-respiratory`, cdmx1$DEFUNCION) %>% prop.table(1)
table(cdmx1$`sintomasNon-respiratory`, cdmx1$tipo_paciente, cdmx1$DEFUNCION)
table(cdmx1$`sintomasNon-respiratory`, cdmx1$evolucion_caso)

## Modelo hospitalizaci?n ##
resp<- cdmx1%>%filter(sintomas=="Non-respiratory")
m1<-coxph(Surv(FU_admission, tipo_paciente)~contacto_infeccion_viral+vomito+dolor_abdominal+conjuntivitis+inmunosupresivo+diabetes+hipertension+frailty(cve_entidad_unidad_medica), data=resp)
summary(m1)
cox.zph(m1)

## Modelo mortalidad asintom?ticos ##
m1<-coxph(Surv(FU_time, DEFUNCION)~contacto_infeccion_viral+edad+diabetes+epoc+insuficiencia_renal_cronica+strata(sexo)+frailty(cve_entidad_unidad_medica), data=resp)
summary(m1)
cox.zph(m1)

table(cdmx1$`sintomasNon-respiratory`, cdmx1$evolucion_caso)
#### Outcomes ####
cdmx1<- cdmx %>% filter(resultado_definitivo=="SARS-CoV-2")
m1<-coxph(Surv(FU_time, DEFUNCION)~sintomas+comorb+edad+sexo+frailty(cve_entidad_unidad_medica), data=cdmx1)
summary(m1)

m1<-coxph(Surv(FU_time, DEFUNCION)~sintomas2+comorb+edad+sexo+frailty(cve_entidad_unidad_medica), data=cdmx1)
summary(m1)

m1<-coxph(Surv(FU_admission, tipo_paciente)~sintomas2+comorb+edad+sexo+frailty(cve_entidad_unidad_medica), data=cdmx1)
summary(m1)

m1<-coxph(Surv(FU_admission, tipo_paciente)~sintomas+comorb+edad+sexo+frailty(cve_entidad_unidad_medica), data=cdmx1)
summary(m1)

#### Figura 1 ####
sint<-cdmx1 %>% group_by(rango_de_edad, sintomas, .drop=T)%>%
  summarise(n=n())%>%drop_na()%>%spread(sintomas, n)%>%
  mutate(freq1=Symptomatic/49507, 
         freq2=`Non-respiratory`/3168 ,
         freq3=Asymptomatic/2211)%>%
  dplyr::select(rango_de_edad, freq1, freq2, freq3)
edad<-c(rep(1:11,3))
edad<-factor(edad, labels = c("0-05","06-15", "16-20",  "21-30",  "31-40",  "41-50",  "51-60",
                              "61-70",  "71-80",  "81-90",  "91-100"))
freq<-c(sint$freq1, sint$freq2, sint$freq3)
group<-c(rep("RS+NRS",11), rep("NRS",11),rep("Asymptomatic",11))
sint2<-data.frame(edad, freq, group)
class(sint2$edad)
fig1a<-ggplot(sint2, aes(y=freq, x=edad,group=group, fill=group)) +
  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8))+
  xlab("Age range (years)")+ylab("Frequency(%)")+
  scale_y_continuous(labels = scales::percent_format())+
  scale_fill_jama()+labs(fill="Group",ylab="Frequency(%)")+theme_classic()+
  theme(axis.text.x = element_text(angle=90), legend.position = "top")

sint<-cdmx1 %>% filter(asintomaticos==0)%>%group_by(asint)%>%
  summarise(diarr=sum(diarrea,na.rm=T), cef=sum(cefalea, na.rm=T), mial=sum(mialgias, na.rm=T), art=sum(artralgias, na.rm=T),n=n(),
            at=sum(ataque_al_estado_general, na.rm=T), rin=sum(rinorrea, na.rm=T), da=sum(dolor_abdominal, na.rm=T), dt=sum(dolor_toracico, na.rm=T),
            vom=sum(vomito, na.rm=T), conj=sum(conjuntivitis, na.rm=T), irri=sum(irritabilidad, na.rm=T), sudden=sum(inicio_subito_sintomas, na.rm=T)) %>% 
  mutate(freq1 = diarr /n,freq2 = cef / n,freq3 = mial / n, freq4=art/n, freq5=at/n, freq6=rin/n, freq7=da/n, freq8=dt/n, 
         freq9=vom/n,freq10=conj/n, freq11=irri/n, freq12=sudden/n)%>%drop_na()
sint1<-sint[,c(1,15:26)]
names(sint1)<-c("Status", "Diarrhea", "Headache", "Myalgias", "Arthralgias", "General malaise", "Rinorrhea", "Abdominal pain", "Chest pain",
                "Vomiting", "Conjunctivitis", "Irritability", "Sudden onset")
freq1<-c(unlist(sint1[1,-1]));freq2<-c(unlist(sint1[2,-1]))
freq<-c(freq1, freq2)
names<-c(rep(names(sint1)[-1],2))
group<-c(rep("RS+NRS",12), rep("NRS",12)); group<-factor(group, label=c("NRS","RS+NRS"))
sint2<-data.frame(names, freq, group)
fig1b<-ggplot(sint2%>% arrange(), aes(y=freq, x=group,fill=group)) +
  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8))+
  ylab("Frequency (%)")+facet_wrap(~names, scales = "free")+
  scale_y_continuous(labels = scales::percent_format())+xlab("")+
  geom_text(position = position_dodge(width= 0.8),aes(label=round(freq*100,2)), vjust=1.6,color="white", size=3.5)+
  scale_fill_jama()+labs(fill="Symptoms")+theme_pubr()+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

fig1<-ggarrange(fig1a, fig1b, labels = c("A", "B"))

ggsave(file = "Figure1_asint.jpg", 
       fig1,
       width = 40, 
       height = 15,
       units=c("cm"),
       dpi = 600,
       limitsize = FALSE)

#### Figure 2 ####
m1<-glmer((sintomas=="Symptomatic")~edad_cat+pers_sal+diabetes+contacto_infeccion_viral+obesidad+insuficiencia_renal_cronica+asma+sexo+epoc+hipertension+tabaquismo+inmunosupresivo+(1|cve_entidad_unidad_medica), family="binomial",data=cdmx1)
m2<-glmer((sintomas=="Non-respiratory")~edad_cat+pers_sal+diabetes+contacto_infeccion_viral+insuficiencia_renal_cronica+obesidad+asma+edad_cat+sexo+epoc+hipertension+tabaquismo+inmunosupresivo+(1|cve_entidad_unidad_medica), family="binomial",data=cdmx1[cdmx1$asint==1,]);summary(m1)
fig2a<-plot_summs(m1,m2, exp=T, model.names = c("RS+NRS vs WoRS", "NRS vs. Asymptomatic"), legend.title = "Symptoms", 
           coefs = c("40-60y"="edad_cat40-60y",">60y"="edad_cat>60y","Diabetes"="diabetes","Obesity"="obesidad","Asthma"="asma",
                     "Male sex"="sexoMASCULINO", "COPD"="epoc", "Hypertension"="hipertension",
                     "Smoking"="tabaquismo", "Immunosupression"="inmunosupresivo", "CKD"="insuficiencia_renal_cronica",
                     "HCW"="pers_sal", "Previous viral exposure"="contacto_infeccion_viral1"))+
  xlab("Mixed effect model OR (95%CI)")
mod1_km<-survfit(Surv(FU_time, DEFUNCION) ~ sintomas, data = cdmx1)
fig1<-ggsurvplot(mod1_km, data = cdmx1, size = 1,palette = "bw",conf.int = F,
                 risk.table = T,pval = TRUE,ggtheme = theme_classic(),xlab="Time from evaluation or symptom onset (Days)",
                 ylab="Survival probability",
                 legend.labs = c("Asymptomatic",
                                 "RS+NRS",
                                 "NRS"),
                 xlim = c(0,30),
                 ylim= c(0,1.0),
                 break.y.by= c(0.1),
                 break.x.by= c(5),
                 pval.coord = c(0, 0.25))+theme_survminer(base_size = 9,
                                                          base_family = "Arial")

fig2b<-ggarrange(fig1$plot, fig1$table, heights = c(2, 0.7),
                 ncol = 1, nrow = 2)

fig2<-ggarrange(fig2a, fig2b, labels = c("A", "B"))

ggsave(file = "Figure2_asint.jpg", 
       fig2,
       width = 35, 
       height = 12,
       units=c("cm"),
       dpi = 600,
       limitsize = FALSE)

#### Table 1 ####
table(cdmx1$sintomas)

table(cdmx1$sintomas, cdmx1$sexo)
table(cdmx1$sintomas, cdmx1$sexo) %>% prop.table(1)

tapply(cdmx1$edad, cdmx1$sintomas, mean)
tapply(cdmx1$edad, cdmx1$sintomas, sd)

table(cdmx1$sintomas, cdmx1$DEFUNCION)
table(cdmx1$sintomas, cdmx1$DEFUNCION) %>% prop.table(1)

table(cdmx1$sintomas, cdmx1$tipo_paciente)
table(cdmx1$sintomas, cdmx1$tipo_paciente) %>% prop.table(1)

table(cdmx1$sintomas, cdmx1$uci)
table(cdmx1$sintomas, cdmx1$uci) %>% prop.table(1)

table(cdmx1$sintomas, cdmx1$intubado)
table(cdmx1$sintomas, cdmx1$intubado) %>% prop.table(1)

table(cdmx1$sintomas, cdmx1$contacto_infeccion_viral)
table(cdmx1$sintomas, cdmx1$contacto_infeccion_viral) %>% prop.table(1)

table(cdmx1$sintomas, cdmx1$diabetes)
table(cdmx1$sintomas, cdmx1$diabetes) %>% prop.table(1)

table(cdmx1$sintomas, cdmx1$epoc)
table(cdmx1$sintomas, cdmx1$epoc) %>% prop.table(1)

table(cdmx1$sintomas, cdmx1$asma)
table(cdmx1$sintomas, cdmx1$asma) %>% prop.table(1)
chisq.test(cdmx1$sintomas, cdmx1$asma)

table(cdmx1$sintomas, cdmx1$VIH_SIDA)
table(cdmx1$sintomas, cdmx1$VIH_SIDA) %>% prop.table(1)
chisq.test(cdmx1$sintomas, cdmx1$VIH_SIDA)

table(cdmx1$sintomas, cdmx1$obesidad)
table(cdmx1$sintomas, cdmx1$obesidad) %>% prop.table(1)
chisq.test(cdmx1$sintomas, cdmx1$obesidad)

table(cdmx1$sintomas, cdmx1$tabaquismo)
table(cdmx1$sintomas, cdmx1$tabaquismo) %>% prop.table(1)
chisq.test(cdmx1$sintomas, cdmx1$tabaquismo)

table(cdmx1$sintomas, cdmx1$enfermedad_cardiaca)
table(cdmx1$sintomas, cdmx1$enfermedad_cardiaca) %>% prop.table(1)
chisq.test(cdmx1$sintomas, cdmx1$enfermedad_cardiaca)

table(cdmx1$sintomas, cdmx1$hipertension)
table(cdmx1$sintomas, cdmx1$hipertension) %>% prop.table(1)
chisq.test(cdmx1$sintomas, cdmx1$hipertension)

table(cdmx1$sintomas, cdmx1$inmunosupresivo)
table(cdmx1$sintomas, cdmx1$inmunosupresivo) %>% prop.table(1)
chisq.test(cdmx1$sintomas, cdmx1$inmunosupresivo)

table(cdmx1$sintomas, cdmx1$insuficiencia_renal_cronica)
table(cdmx1$sintomas, cdmx1$insuficiencia_renal_cronica) %>% prop.table(1)
chisq.test(cdmx1$sintomas, cdmx1$insuficiencia_renal_cronica)

table(cdmx1$sintomas, cdmx1$pers_sal)
table(cdmx1$sintomas, cdmx1$pers_sal) %>% prop.table(1)
chisq.test(cdmx1$sintomas, cdmx1$pers_sal)



#### Correspondence Analysis ####

cdmx2<-cdmx1[,c(32:49)]
names(cdmx2)<-c("Fever", "Cough", "Odynophagia", "Dyspnea", "Irritability", "Diarrhea", "Chest pain",
                "Shivering", "Headache", "Myalgias", "Arthalgias", "General malaise", "Rinorrhea",
                "Polypnea", "Vomiting", "Abdominal pain", "Conjunctivitis", "Cyanosis")

res.ca<-CA(cdmx2, ncp = 5)
fviz_screeplot(res.ca, addlabels = TRUE)
res <- caCluster(cdmx2, opt.part=TRUE, which="cols", dim=4)

