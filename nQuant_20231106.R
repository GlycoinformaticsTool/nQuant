#-------------------------nQuant for single running by label-free quantification-------------------------------------------#
rm(list=ls(all=TRUE))
install.packages("readxl")
install.packages("dplyr")
library("readxl")
library("dplyr")
setwd("path")
MS1<-read.table("ms1.txt",header=FALSE,fill = TRUE,sep = " ")
MS1<-as.data.frame(MS1)
MS1<-MS1[,-c(3:5)]
MS1<-MS1[-c(1:11),]
colnames(MS1)<-c("Mz","Intensity")
rownames(MS1)<-c(1:nrow(MS1))
RetTimeIndex<-grep("RetTime", MS1$Mz)
NglycanList<-read_excel("Nglycan_Identification_Result.xlsx")
NglycanList<-as.data.frame(NglycanList)
nQuant<-c("Nglycan","Code","Type","TotalIntensity")
for (a in 1:nrow(NglycanList)){
  print(a)
  TotalIntensity<-0
  AllMz<-vector()
  for (b in 1:5){                                                       #Range of charge states 
    Adducts<-c(rep(c(1.00783, 18.03437, 22.98977,38.96371), times=b))   #Different adducts 
    Sum<-combn(Adducts, b, sum)
    Sum<-unique(Sum)
    MplusA<-NglycanList$Mass[a]+Sum
    AllMz1<-MplusA/b
    AllMz1<-as.vector(AllMz1)
    AllMz<-c(AllMz,AllMz1)
    }
    for (c in RetTimeIndex){
      t<-gsub("I\tRetTime	", "", MS1$Mz[c])
      t<-as.numeric(t)
      if (between(t,NglycanList$Start1[a]*60,NglycanList$End1[a]*60) | between(t,NglycanList$Start2[a]*60,NglycanList$End2[a]*60)) {
         MzList<-c("Mz","Intensity")
         w<-which(RetTimeIndex==c)
         Row1<-c+3
         Row2<-RetTimeIndex[w+1]-3
         for (d in Row1:Row2){
            Mz1<-MS1[d,]
            MzList<-rbind(MzList,Mz1)
            }
            colnames(MzList)<-MzList[1,]
            MzList<-MzList[-1,]
            MzList<-as.data.frame(MzList)
            MzList$Mz<-as.numeric(MzList$Mz)
            MzList$Intensity<-as.numeric(MzList$Intensity)
            Error1<-MzList$Mz-NglycanList$'m/z'[a]
            Error1<-Error1/NglycanList$'m/z'[a]*1000000
            Error1<-as.vector(Error1)
            Min1<-which.min(abs(Error1))
            if (abs(Error1[Min1])<=3){
              for (e in AllMz){
                Error2<-MzList$Mz-e
                Error2<-Error2/e*1000000
                Error2<-as.vector(Error2)
                Min2<-which.min(abs(Error2))
                if (abs(Error2[Min2])<=3){
                  TotalIntensity<-TotalIntensity+MzList$Intensity[Min2]
                  Result<-paste(as.character(a),as.character(MzList$Intensity[Min2]),as.character(TotalIntensity), sep = "_")
                  print(Result)
                  }
                }
              }
            }
          }   
         C<-NglycanList$Neu5Gc[a]*19+NglycanList$Neu5Ac[a]*18+NglycanList$HexNAc[a]*13+NglycanList$Hex[a]*11+NglycanList$Fuc[a]*10+NglycanList$`Red-HexNAc`[a]*14-(NglycanList$Neu5Gc[a]+NglycanList$Neu5Ac[a]+NglycanList$HexNAc[a]+NglycanList$Hex[a]+NglycanList$Fuc[a])*2
         H<-NglycanList$Neu5Gc[a]*35+NglycanList$Neu5Ac[a]*33+NglycanList$HexNAc[a]*25+NglycanList$Hex[a]*22+NglycanList$Fuc[a]*20+NglycanList$`Red-HexNAc`[a]*29-(NglycanList$Neu5Gc[a]+NglycanList$Neu5Ac[a]+NglycanList$HexNAc[a]+NglycanList$Hex[a]+NglycanList$Fuc[a])*6
         N<-NglycanList$Neu5Gc[a]+NglycanList$Neu5Ac[a]+NglycanList$HexNAc[a]+NglycanList$`Red-HexNAc`[a]
         O<-NglycanList$Neu5Gc[a]*10+NglycanList$Neu5Ac[a]*9+NglycanList$HexNAc[a]*6+NglycanList$Hex[a]*6+NglycanList$Fuc[a]*5+NglycanList$`Red-HexNAc`[a]*6-(NglycanList$Neu5Gc[a]+NglycanList$Neu5Ac[a]+NglycanList$HexNAc[a]+NglycanList$Hex[a]+NglycanList$Fuc[a])
         MonoisotopicAbundance<-0.9894^C*0.999855^H*0.996205^N*0.99757^O   
         TotalIntensity<-TotalIntensity/MonoisotopicAbundance             #The novel label-free quantification algorithm
         nQuant1<-c(NglycanList$Nglycan[a],NglycanList$Code[a],NglycanList$Type[a],TotalIntensity)
         nQuant<-rbind(nQuant,nQuant1)
         }
nQuant<-as.data.frame(nQuant)
colnames(nQuant)<-nQuant[1,]
nQuant<-nQuant[-1,]
write.table(nQuant,"nQuant.csv",sep=",",row.names=FALSE,col.names = TRUE)



#-------------------------nQuant for multiple running by label-free quantification-------------------------------------------#
rm(list=ls(all=TRUE))
install.packages("readxl")
install.packages("dplyr")
library("readxl")
library("dplyr")               
setwd("path")
MS1_1<-read.table("ms1_1.txt",header=FALSE,fill = TRUE,sep = " ")
MS1_2<-read.table("ms1_2.txt",header=FALSE,fill = TRUE,sep = " ")
MS1_3<-read.table("ms1_3.txt",header=FALSE,fill = TRUE,sep = " ")
MS1_1<-as.data.frame(MS1_1)
MS1_2<-as.data.frame(MS1_2)
MS1_3<-as.data.frame(MS1_3)
MS1_1<-MS1_1[,-c(3:5)]
MS1_2<-MS1_2[,-c(3:5)]
MS1_3<-MS1_3[,-c(3:5)]
colnames(MS1_1)<-c("Mz1","Intensity1")
colnames(MS1_2)<-c("Mz2","Intensity2")
colnames(MS1_3)<-c("Mz3","Intensity3")
MS1_1_n<-nrow(MS1_3)-nrow(MS1_1)
MS1_2_n<-nrow(MS1_3)-nrow(MS1_2)
MS1_1_add<-cbind(rep(0,times=MS1_1_n),rep(0,times=MS1_1_n))
colnames(MS1_1_add)<-c("Mz1","Intensity1")
MS1_1<-rbind(MS1_1,MS1_1_add)
MS1_2_add<-cbind(rep(0,times=MS1_2_n),rep(0,times=MS1_2_n))
colnames(MS1_2_add)<-c("Mz2","Intensity2")
MS1_2<-rbind(MS1_2,MS1_2_add)
MS1<-cbind(MS1_1,MS1_2,MS1_3)
MS1<-MS1[-c(1:11),]
rownames(MS1)<-c(1:nrow(MS1))
RetTimeIndex_MS1_1 <- grep("RetTime", MS1$Mz1)
RetTimeIndex_MS1_2 <- grep("RetTime", MS1$Mz2)
RetTimeIndex_MS1_3 <- grep("RetTime", MS1$Mz3)
NglycanList<-read_excel("Nglycan_Identification_Result.xlsx")
NglycanList<-as.data.frame(NglycanList)
nQuant<-c("Nglycan","Code","Type","TotalIntensity_MS1_1","TotalIntensity_MS1_2","TotalIntensity_MS1_3")
SampleCount<-c(1,2,3)                                                  #Count of running
RetTimeIndex<-list(RetTimeIndex_MS1_1,RetTimeIndex_MS1_2,RetTimeIndex_MS1_3)
for (a in 1:nrow(NglycanList)){
  print(a)
  TotalIntensity<-c(0,0,0)
  AllMz<-vector()
  for (b in 1:5){                                                       #Range of charge states
    Adducts<-c(rep(c(1.00783, 18.03437, 22.98977,38.96371), times=b))   #Different adducts 
    Sum<-combn(Adducts, b, sum)
    Sum<-unique(Sum)
    MplusA<-NglycanList$Mass[a]+Sum
    AllMz1<-MplusA/b
    AllMz1<-as.vector(AllMz1)
    AllMz<-c(AllMz,AllMz1)
    }
    AllMz<-AllMz[AllMz>=450 & AllMz<=2000]
    for (c in SampleCount){
      Mzncol<-c*2-1
      Intensityncol<-c*2
      for (d in RetTimeIndex[[c]]){
        t<-gsub("I\tRetTime	", "", MS1[d,Mzncol])
        t<-as.numeric(t)
        if (between(t,NglycanList$Start1[a]*60,NglycanList$End1[a]*60) | between(t,NglycanList$Start2[a]*60,NglycanList$End2[a]*60)) {
          MzList<-c("Mz","Intensity")
          w<-which(RetTimeIndex[[c]]==d)
          Row1<-d+3
          Row2<-RetTimeIndex[[c]][w+1]-3
          for (e in Row1:Row2){
            Mz1<-MS1[e,c(Mzncol,Intensityncol)]
            MzList<-rbind(MzList,Mz1)
            }
          colnames(MzList)<-MzList[1,]
          MzList<-MzList[-1,]
          MzList<-as.data.frame(MzList)
          MzList$Mz<-as.numeric(MzList$Mz)
          MzList$Intensity<-as.numeric(MzList$Intensity)
          Error1<-MzList$Mz-NglycanList$'m/z'[a]
          Error1<-Error1/NglycanList$'m/z'[a]*1000000
          Error1<-as.vector(Error1)
          Min1<-which.min(abs(Error1))
          if (abs(Error1[Min1])<=3){
            for (f in AllMz){
              Error2<-MzList$Mz-f
              Error2<-Error2/f*1000000
              Error2<-as.vector(Error2)
              Min2<-which.min(abs(Error2))
              if (abs(Error2[Min2])<=3){
                TotalIntensity[c]<-TotalIntensity[c]+MzList$Intensity[Min2]
                Result<-paste(NglycanList$Nglycan[a],as.character(MzList$Intensity[Min2]),as.character(TotalIntensity[c]), sep = "_")
                print(Result)
               }
             }
           }
         }
       }   
    C<-NglycanList$Neu5Gc[a]*19+NglycanList$Neu5Ac[a]*18+NglycanList$HexNAc[a]*13+NglycanList$Hex[a]*11+NglycanList$Fuc[a]*10+NglycanList$`Red-HexNAc`[a]*14-(NglycanList$Neu5Gc[a]+NglycanList$Neu5Ac[a]+NglycanList$HexNAc[a]+NglycanList$Hex[a]+NglycanList$Fuc[a])*2
    H<-NglycanList$Neu5Gc[a]*35+NglycanList$Neu5Ac[a]*33+NglycanList$HexNAc[a]*25+NglycanList$Hex[a]*22+NglycanList$Fuc[a]*20+NglycanList$`Red-HexNAc`[a]*29-(NglycanList$Neu5Gc[a]+NglycanList$Neu5Ac[a]+NglycanList$HexNAc[a]+NglycanList$Hex[a]+NglycanList$Fuc[a])*6
    N<-NglycanList$Neu5Gc[a]+NglycanList$Neu5Ac[a]+NglycanList$HexNAc[a]+NglycanList$`Red-HexNAc`[a]
    O<-NglycanList$Neu5Gc[a]*10+NglycanList$Neu5Ac[a]*9+NglycanList$HexNAc[a]*6+NglycanList$Hex[a]*6+NglycanList$Fuc[a]*5+NglycanList$`Red-HexNAc`[a]*6-(NglycanList$Neu5Gc[a]+NglycanList$Neu5Ac[a]+NglycanList$HexNAc[a]+NglycanList$Hex[a]+NglycanList$Fuc[a])
    MonoisotopicAbundance<-0.9894^C*0.999855^H*0.996205^N*0.99757^O             
    TotalIntensity[c]<-TotalIntensity[c]/MonoisotopicAbundance                  ##The novel label-free quantification algorithm
  }
  nQuant1<-c(NglycanList$Nglycan[a],NglycanList$Code[a],NglycanList$Type[a],TotalIntensity[1],TotalIntensity[2],TotalIntensity[3])
  nQuant<-rbind(nQuant,nQuant1)
}
nQuant<-as.data.frame(nQuant)
colnames(nQuant)<-nQuant[1,]
nQuant<-nQuant[-1,]
write.table(nQuant,"nQuant.csv",sep=",",row.names=FALSE,col.names = TRUE)



#-------------------------nQuant for multiple running by isotopic labelling quantification-------------------------------------------#
rm(list=ls(all=TRUE))
install.packages("readxl")
install.packages("dplyr")
library("readxl")
library("dplyr") 
setwd("path")
MS1_1<-read.table("ms1_1.txt",header=FALSE,fill = TRUE,sep = " ")
MS1_2<-read.table("ms1_2.txt",header=FALSE,fill = TRUE,sep = " ")
MS1_3<-read.table("ms1_3.txt",header=FALSE,fill = TRUE,sep = " ")
MS1_1<-as.data.frame(MS1_1)
MS1_2<-as.data.frame(MS1_2)
MS1_3<-as.data.frame(MS1_3)
MS1_1<-MS1_1[,-c(3:5)]
MS1_2<-MS1_2[,-c(3:5)]
MS1_3<-MS1_3[,-c(3:5)]
colnames(MS1_1)<-c("Mz1","Intensity1")
colnames(MS1_2)<-c("Mz2","Intensity2")
colnames(MS1_3)<-c("Mz3","Intensity3")
MS1_2_n<-nrow(MS1_1)-nrow(MS1_2)
MS1_3_n<-nrow(MS1_1)-nrow(MS1_3)
MS1_2_add<-cbind(rep(0,times=MS1_2_n),rep(0,times=MS1_2_n))
colnames(MS1_2_add)<-c("Mz2","Intensity2")
MS1_2<-rbind(MS1_2,MS1_2_add)
MS1_3_add<-cbind(rep(0,times=MS1_3_n),rep(0,times=MS1_3_n))
colnames(MS1_3_add)<-c("Mz3","Intensity3")
MS1_3<-rbind(MS1_3,MS1_3_add)
MS1<-cbind(MS1_1,MS1_2,MS1_3)
MS1<-MS1[-c(1:11),]
rownames(MS1)<-c(1:nrow(MS1))
RetTimeIndex_MS1_1 <- grep("RetTime", MS1$Mz1)
RetTimeIndex_MS1_2 <- grep("RetTime", MS1$Mz2)
RetTimeIndex_MS1_3 <- grep("RetTime", MS1$Mz3)
RetTimeIndex<-list(RetTimeIndex_MS1_1,RetTimeIndex_MS1_2,RetTimeIndex_MS1_3)
NglycanList<-read_excel("Nglycan_Identification_Result.xlsx")
NglycanList<-as.data.frame(NglycanList)
SampleCount<-c(1,2,3)                                                   #Count of running
nQuant<-c("Nglycan","Code","Type","12C_1","12C_2","12C_3","13C_1","13C_2","13C_3")
for (a in 1:nrow(NglycanList)){
  TotalIntensity<-c(0,0,0,0,0,0)
  for (b in SampleCount){
    Mzncol<-b*2-1
    Intensityncol<-b*2
    for (c in RetTimeIndex[[b]]){
      t<-gsub("I\tRetTime	", "", MS1[c,Mzncol])
      t<-as.numeric(t)
      if (between(t,NglycanList$Start1[a]*60,NglycanList$End1[a]*60) | between(t,NglycanList$Start2[a]*60,NglycanList$End2[a]*60)) {
        MzList<-c("Mz","Intensity")
        w<-which(RetTimeIndex[[b]]==c)
        Row1<-c+3
        Row2<-RetTimeIndex[[b]][w+1]-3
        for (d in Row1:Row2){
          MzList1<-c(MS1[d,Mzncol],MS1[d,Intensityncol])
          MzList<-rbind(MzList,MzList1)
          }
          colnames(MzList)<-MzList[1,]
          MzList<-MzList[-1,]
          MzList<-as.data.frame(MzList)
          MzList$Mz<-as.numeric(MzList$Mz)
          MzList$Intensity<-as.numeric(MzList$Intensity)
          Err_12C<-MzList$Mz-NglycanList$`12C_m/z`[a]
          Err_12C<-Err_12C/NglycanList$`12C_m/z`[a]*1000000
          Err_12C<-as.vector(Err_12C)
          Min1<-which.min(abs(Err_12C))
          Err_13C<-MzList$Mz-NglycanList$`13C_m/z`[a]
          Err_13C<-Err_13C/NglycanList$`13C_m/z`[a]*1000000
          Err_13C<-as.vector(Err_13C)
          Min2<-which.min(abs(Err_13C))
          if (abs(Err_12C[Min1])<=5 & abs(Err_13C[Min2])<=5){
            C<-NglycanList$Neu5Gc[a]*19+NglycanList$Neu5Ac[a]*18+NglycanList$HexNAc[a]*13+NglycanList$Hex[a]*11+NglycanList$Fuc[a]*10+NglycanList$`Red-HexNAc`[a]*14-(NglycanList$Neu5Gc[a]+NglycanList$Neu5Ac[a]+NglycanList$HexNAc[a]+NglycanList$Hex[a]+NglycanList$Fuc[a])*2
            H<-NglycanList$Neu5Gc[a]*35+NglycanList$Neu5Ac[a]*33+NglycanList$HexNAc[a]*25+NglycanList$Hex[a]*22+NglycanList$Fuc[a]*20+NglycanList$`Red-HexNAc`[a]*29-(NglycanList$Neu5Gc[a]+NglycanList$Neu5Ac[a]+NglycanList$HexNAc[a]+NglycanList$Hex[a]+NglycanList$Fuc[a])*6
            N<-NglycanList$Neu5Gc[a]+NglycanList$Neu5Ac[a]+NglycanList$HexNAc[a]+NglycanList$`Red-HexNAc`[a]
            O<-NglycanList$Neu5Gc[a]*10+NglycanList$Neu5Ac[a]*9+NglycanList$HexNAc[a]*6+NglycanList$Hex[a]*6+NglycanList$Fuc[a]*5+NglycanList$`Red-HexNAc`[a]*6-(NglycanList$Neu5Gc[a]+NglycanList$Neu5Ac[a]+NglycanList$HexNAc[a]+NglycanList$Hex[a]+NglycanList$Fuc[a])
            MonoisotopicAbundance_C12<-0.9894^C*0.999855^H*0.996205^N*0.99757^O   
            Intensity_C12<-MzList$Intensity[Min1]/MonoisotopicAbundance_C12              
            C13<-NglycanList$Neu5Gc[a]*8+NglycanList$Neu5Ac[a]*7+NglycanList$HexNAc[a]*5+NglycanList$Hex[a]*5+NglycanList$Fuc[a]*4+NglycanList$`Red-HexNAc`[a]*6-(NglycanList$Neu5Gc[a]+NglycanList$Neu5Ac[a]+NglycanList$HexNAc[a]+NglycanList$Hex[a]+NglycanList$Fuc[a])*2
            C12<-C-C13
            MonoisotopicAbundance_C13<-0.9894^C12*0.999855^H*0.996205^N*0.99757^O*0.991763^C13
            Intensity_C13<-MzList$Intensity[Min2]/MonoisotopicAbundance_C13                    #The calibrated isotopic labelling quantification algorithm
            TotalIntensity[b]<-TotalIntensity[b]+Intensity_C12
            TotalIntensity[b+3]<-TotalIntensity[b+3]+Intensity_C13
            Result<-paste(NglycanList$Nglycan[a],as.character(MzList$Intensity[Min1]),as.character(MzList$Intensity[Min2]),sep = "_")
            print(Result)
            }
          }   
        }
      }
      nQuant1<-c(NglycanList$Nglycan[a],NglycanList$Code[a],NglycanList$Type[a],TotalIntensity[1],TotalIntensity[2],TotalIntensity[3],TotalIntensity[4],TotalIntensity[5],TotalIntensity[6])
      nQuant<-rbind(nQuant,nQuant1)
     }
nQuant<- as.data.frame(nQuant)
colnames(nQuant)<-nQuant[1,]
nQuant<-nQuant[-1,]
write.table(nQuant,"nQuant.csv",sep=",",row.names=FALSE,col.names = TRUE)





