#Functional gene clustering using DAVID algorithms---Hao Ma, Guangtu Gao, Gregory M Weber---2017

#The input file has two columns with hearder names of 'Gene_identifier' and 'Annotation_information', which can be saved as Text (Tab delimited) file") 

Path_file<-readline("Please enter the PATH and FILE seperated by \\:")
# For example, copy and paste: c:\Users\haoma\My Documents\Egg_quality_ paper_CRIS_Project\Cluster analysis\Gene_GO_test.txt

cv<-as.numeric(strsplit(readline("Enter cut-off kappa values seperated by space : "), " ")[[1]])
# For example, type: 0.3 0.35 0.4 0.5

# Execute above line by line, then select the following script and run. kappa values and results will be saved under the folder provided above. The result file can be opened with Excel-> Select Text file->Delimited->Tab->General.

D<-read.table(Path_file,head=TRUE, colClasses =c("character","character"))

D<-D[order(D$Gene_identifier),]

Total_N<-nrow(D)

D_unique<-D[!duplicated(D$Gene_identifier),]   #Extract unique gene name

N<-nrow(D_unique)

Kappa_matrix<-matrix(0,N+1, N+1) #create a matrix for kappa values

for (i in 1:N){
  
  Kappa_matrix[i+1,1]=D_unique[i,1]
  
  Kappa_matrix[1,i+1]=D_unique[i,1]
}

Cat_unique<-D[!duplicated(D$Annotation_information),]

M<-nrow(Cat_unique)

paste("Total number of rows in your input file:", Total_N)

paste ("Total number of unique gene identifiers:", N)

paste ("Total number of unique gene annotation information:", M)

Gene_Cat_matrix<-matrix(0, N+1, M+1)

for (i in 1:N){
  
  Gene_Cat_matrix[i+1,1]=Kappa_matrix[i+1,1]
  
}

for (j in 1:M){
  
  Gene_Cat_matrix[1,j+1]=Cat_unique[j,2]
  
}
rm(D_unique, Cat_unique)

paste ("Generating two-way table:") 

for (i in 1:N){
  
  D1<-subset(D,Gene_identifier==Gene_Cat_matrix[i+1,1]) 
  
  N1<-nrow(D1)
  
  for (j in 1:M){
    
    for (k in 1:N1){
      
      if (Gene_Cat_matrix[i+1,1]==D1[k,1] && Gene_Cat_matrix[1,j+1]==D1[k,2]){ 
        
        Gene_Cat_matrix[i+1,j+1]=1}
      
    }        
    
    print (c(i+1,j+1))
    
  }}

rm(D1)

paste ("Calculating kappa values:")

for(i in 2:(N+1)){
  
  for (j in 2:(N+1)){ 
    
    if (i==j){
      
      Kappa_matrix[i,j]=1}
    
    if(i<j){
      
      a<-0;b<-0;c<-0;d<-0 
      
      sum1 <- sum(as.numeric(Gene_Cat_matrix[i,2:(M+1)]))
      
      sum2 <- sum(as.numeric(Gene_Cat_matrix[j,2:(M+1)]))
      
      for (k in 2:(M+1)){
        
        if (Gene_Cat_matrix[i,k]==1 && Gene_Cat_matrix[j,k]==1){
          
          a<-a+1}
      }                                
      b <- sum1-a; c <- sum2-a; d <- M+a-sum1-sum2                                            
      
      Kappa_matrix[i,j]<-((a+d)*M-(a+b)*(a+c)-(c+d)*(b+d))/(M^2-(a+b)*(a+c)-(c+d)*(b+d))
      
      Kappa_matrix[j,i]=Kappa_matrix[i,j]
      
      print (c(i,j))
      
    }             
    
  }
  
}

Path_file1 <- substr(Path_file,1,nchar(Path_file)-4)

Path_file2 <- paste(c(Path_file1,"_Kappa_value",".txt"),collapse = "")

write.table(Kappa_matrix,file=Path_file2, sep="\t", col.names = F, row.names = F)



for (aa in cv){
  
  D <- read.table(Path_file2, head=TRUE);
  
  Path_file1 <- substr(Path_file2,1,nchar(Path_file)-4)
  
  Path_file3 <- paste(c(Path_file1,"_cut_off_",aa,"_cluster.txt"),collapse = "")
  
  Row<-nrow(D)
  
  Col<-ncol(D)
  
  for (x in 1:Row){                 
    
    for (y in 2:Col){            
      
      if (D[x,y]>=aa){        
        
        D[x,y]<-1}
      
      else {D[x,y]<-0}
      
    }
    
  }
  
  Group<-cbind(D,0)                 
  
  names = c(colnames(D), "Seed")
  
  colnames(Group)<- names
  
  for (i in 1:Row){                  
    
    j<-0
    
    S<-(sum(D[i,2:Col])-1)          
    
    if (S<3){
      
      Group[i,(Col+1)]<-0
      
    }
    else{
      a<-0
      
      Temp<-matrix(0,Row,Col)
      
      for (j in 2:Col){
        
        if((j-1)!=i){
          
          if(D[i,j]==1){
            
            Temp[,j]<-D[,j]
          }
          
        }
        
      }
      
      for (k in 2:Col){
        
        if ((k-1)!=i){
          
          if(D[i,k]==1){ 
            
            a<-a+sum(Temp[(k-1),])
          }
          
        }
      }
      
      if(2*a>=S*(S+1)){
        
        Group[i,(Col+1)]<-1
        
      } else{Group[i,(Col+1)]<-0}
      
    }
    
    print (c(i,j))
    
  }
  
  print("Complete selected seed genes")
  
  D <- Group
  
  Row<-nrow(D)
  
  Col<-ncol(D)
  
  Group<-D
  
  for(i in 1:(Row-1)){
    
    if(Group[i,Col]==1){
      
      for(j in (i+1):Row){
        
        if(Group[j,Col]==1 & Group[i,Col]==1){
          
          a<-0;b<-0;c<-0;s1<-0;s2<-0
          
          s1<-sum(Group[i,2:(Col-1)])
          
          s2<-sum(Group[j,2:(Col-1)])
          
          for(k in 2:(Col-1)){
            
            if(Group[i,k]==1 & Group[j,k]==1){
              
              a<-a+1
              
            }
            
          }
          
          if (2*a>s1){
            
            for (l in 2:(Col-1)){
              
              if (Group[i,l]==1 & Group[j,l]==0){
                
                Group[j,l]<-1
              }
              
            }
            
            Group[i, Col]<-0          
            
          }
          
          else if(2*a>s2){
            
            for (m in 2:(Col-1)){
              
              if (Group[i,m]==0 & Group[j,m]==1){
                
                Group[i,m]<-1
                
              }
              
            }
            Group[j, Col]<-0
            
          }
        }
        print(c(i,j))
      }
    }
  }
  
  names(Group)[Col] <- "Cluster"
  
  names(Group)[1] <- "Seed_gene"
  
  Group<-subset(Group,Cluster==1) 
  
  Row<-nrow(Group)
  
  Col<-ncol(Group)
  
  for (i in 1:Row){
    
    Group[i,Col]<-i
    
    for (j in 2:(Col-1)){
      
      if (Group[i,j]==1){
        
        Group[i,j]<-colnames(Group)[j]
      }
    }
    
  }
  
  write.table(Group, file=Path_file3)
  
}