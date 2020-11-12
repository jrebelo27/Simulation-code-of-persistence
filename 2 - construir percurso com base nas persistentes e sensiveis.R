for(cada_plasmideo in c("R1")){
  
  for (cada_densidade in c("L", "I", "H")){
    
    for (tau_zero in c(20, 30, 50,60,70,80,90,100,110,120,130, 150, 200, 250, 300, 350, 400)){
      
      for(k in c(0.02, 0.025, 0.03,0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080,0.09, 0.095, 0.1, 0.15, 0.2)){
        
        for (beta in c( -1.1, -1.2,-1.5,-1.7, -1.8, -1.9, -2.0, -2.1, -2.2, -2.3, -2.4, -2.5, -2.7, -2.9, -3.1, -3.3, -3.5)){
          
          
          print(paste0("C:/Users/Utilizador/Desktop/Novo com modelo Simsek/dados 2020 sem acordar persistentes ",cada_plasmideo, "/", cada_densidade, "/tau_zero", tau_zero, "-k", k, "-beta", beta))
          setwd(paste0("C:/Users/Utilizador/Desktop/Novo com modelo Simsek/dados 2020 sem acordar persistentes ",cada_plasmideo, "/", cada_densidade, "/tau_zero", tau_zero, "-k", k, "-beta", beta))
          
          for(cada_inicio in c(0.2, 2)){
            
            cada_fim <- cada_inicio*10
            
            for(cada_frequencia in c(1,2,3)){
              
              nome <- paste0("_", cada_plasmideo, "_", cada_densidade, "_tau_zero", tau_zero, "_k", k, "_beta", beta, "_raio_", cada_inicio, "_", cada_fim, "_L", cada_frequencia, ".csv")
              
              
              dados_persistentes <- read.csv(file = paste0("quantidade_persistentes_sobreviventes_por_ger", nome), header = TRUE, sep = ",", stringsAsFactors = FALSE, dec=".")
              dados_sensiveis <- read.csv(file = paste0("quantidade_sensiveis_sobreviventes_por_ger", nome), header = TRUE, sep = ",", stringsAsFactors = FALSE, dec=".")
              
              percurso_final <- c()
              
              for (cada_linha_dados in 1:nrow(dados_sensiveis)){
                
                percurso_sensiveis <- c()
                percurso_persistentes <- c()
                
                for (cada_coluna_dados in 2:ncol(dados_sensiveis)){
                  
                  if(cada_coluna_dados == 2){
                    
                    percurso_sensiveis <- c(percurso_sensiveis, dados_sensiveis[cada_linha_dados, cada_coluna_dados])
                    percurso_persistentes <- c(percurso_persistentes, dados_persistentes[cada_linha_dados, cada_coluna_dados])
                    
                  }else{
                    
                    percurso_sensiveis <- c(percurso_sensiveis, ((percurso_sensiveis[length(percurso_sensiveis)]*2) + dados_sensiveis[cada_linha_dados, cada_coluna_dados]))
                    percurso_persistentes <- c(percurso_persistentes, ((percurso_persistentes[length(percurso_persistentes)]*2) + dados_persistentes[cada_linha_dados, cada_coluna_dados]))
                    
                  }
                }
                
                percurso_total <- percurso_sensiveis + percurso_persistentes
                percurso_final <- rbind(percurso_final, c(dados_sensiveis[cada_linha_dados,1], percurso_total))
              }
              
              nomes_colunas <- colnames(dados_sensiveis)
              percurso_final <- as.data.frame(percurso_final)
              colnames(percurso_final) <- nomes_colunas
              
              write.table(percurso_final,
                          file = paste0("percurso_total_sensiveis_persistentes", nome),
                          sep = ",",
                          row.names = FALSE,
                          col.names = TRUE)
            }
          }
        }
      }
    }
  }
}

