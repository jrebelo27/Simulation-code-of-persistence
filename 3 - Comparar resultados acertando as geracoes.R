###############################################################################################
#Cuidado com as linhas que n?o interessam consoante o plasmideo (os zeros)
###############################################################################################

for (cada_plasmideo in c("R1")){
  
  if (cada_plasmideo == "R1"){
  
    lista_densidades = c("I", "H")
    
  }else{
    
    lista_densidades = c("L", "I", "H")
    
  }
  
  for (cada_densidade in lista_densidades){
    
    
    if(cada_densidade == "L"){
      densidade_nome <- "Low"
      linhas_dados <- c(1)
      
    }
    
    if(cada_densidade == "I"){
      densidade_nome <- "Intermediate" 
      linhas_dados <- c(1, 2, 3)
    }
    
    if(cada_densidade == "H"){
      densidade_nome <- "High" 
      linhas_dados <- c(1, 2, 3)
    }
    
    
    for (erro_inicial in c(4.01)){
      
      valores_maiores_1 <- 0
      valores_menores_1 <- 0
      valores_iguais_1 <- 0
      
      valores_maiores_2 <- 0
      valores_menores_2 <- 0
      valores_iguais_2 <- 0
      
      valores_maiores_3 <- 0
      valores_menores_3 <- 0
      valores_iguais_3 <- 0
      
      linhas_colunas <- c()
      
      
      for (cada_inicio in c(0.2, 2)){
        
        cada_fim <- cada_inicio * 10
        
        for (tau_zero in c(20, 30, 50,60,70,80,90,100,110,120,130, 150, 200, 250, 300, 350, 400)){
          
          for(k in c(0.02, 0.025, 0.03,0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080,0.09, 0.095, 0.1, 0.15, 0.2)){
            
            for (beta in c( -1.1, -1.2,-1.5,-1.7, -1.8, -1.9, -2.0, -2.1, -2.2, -2.3, -2.4, -2.5, -2.7, -2.9, -3.1, -3.3, -3.5)){
              
              
              for(cada_linha_dados in linhas_dados){
                
                print(paste0("C:/Users/Utilizador/Desktop/Novo com modelo Simsek/dados 2020 sem acordar persistentes ",cada_plasmideo, "/", cada_densidade, "/tau_zero", tau_zero, "-k", k, "-beta", beta))
                setwd(paste0("C:/Users/Utilizador/Desktop/Novo com modelo Simsek/dados 2020 sem acordar persistentes ",cada_plasmideo, "/", cada_densidade, "/tau_zero", tau_zero, "-k", k, "-beta", beta))
                
                
                nome <- paste0("_", cada_plasmideo, "_", cada_densidade, "_tau_zero", tau_zero, "_k", k, "_beta", beta, "_raio_", cada_inicio, "_", cada_fim, "_L", cada_linha_dados, ".csv")
                nome_entrada <- paste0("percurso_total_sensiveis_persistentes", nome)
                dados_simulacoes <- read.csv(file = nome_entrada, header = TRUE, sep = ",", stringsAsFactors = FALSE, dec=".")
                
                
                #print(nome)
                
                #print(dados_simulacoes[5,5])
                
                setwd("C:/Users/Utilizador/Desktop/Novo com modelo Simsek/dados originais")
                
                
                nome_entrada_sensiveis <- paste0("dados_", cada_plasmideo, "_", cada_densidade, "_medias.csv")
                dados_originais <- read.csv(file = nome_entrada_sensiveis, header = TRUE, sep = ",", stringsAsFactors = FALSE, dec=".")
                dados_originais[,4] <- as.numeric(dados_originais[,4])
                
                
                sensiveis_iolanda <- dados_originais[cada_linha_dados, 4]
                produtoras_finais_iolanda <- dados_originais[cada_linha_dados, 2]
                produtoras_inicias <- dados_originais[cada_linha_dados, 1]
                
                
                nr_geracoes_iolanda <- round(log2(produtoras_finais_iolanda/produtoras_inicias), 6)
                nr_geracoes_nossas <- (round(log2(produtoras_finais_iolanda/produtoras_inicias), 0))+1
                
                
                
                produtoras_finais_simulacoes <- produtoras_inicias*(2**nr_geracoes_nossas)
                
                
                diferenca_geracoes <- nr_geracoes_iolanda - nr_geracoes_nossas
                
                raio_espalhamento <- cada_inicio
                quantas_vezes_explica <- 0
                
                #############Atencao, isto e para nao sobrepor o 2 dos dois ficheiros (ultimo do 0.2-2 e primeiro do 2-20)
                
                if (cada_inicio == 0.2){
                  linha_em_que_comeca <- 1
                }else{
                  linha_em_que_comeca <- 4
                  raio_espalhamento <- raio_espalhamento + cada_inicio
                }
                
                
                
                for (cada_linha in linha_em_que_comeca:nrow(dados_simulacoes)){
                  
                  
                  valor_a_comparar <- dados_simulacoes[cada_linha, ncol(dados_simulacoes)]
                  
                 
                  
                  
                  valor_a_comparar <- round(valor_a_comparar * (2 ** diferenca_geracoes), 0)
                 
                  
                  erro <- erro_inicial
                  
                  limite_inferior <- sensiveis_iolanda/erro
                  limite_superior <- sensiveis_iolanda*erro
                  
                  if(cada_linha_dados == 1){
                    
                    if(valor_a_comparar > sensiveis_iolanda){
                      valores_maiores_1 <- valores_maiores_1 + 1
                    }else if(valor_a_comparar < sensiveis_iolanda){
                      valores_menores_1 <- valores_menores_1 + 1
                    }else{
                      valores_iguais_1 <- valores_iguais_1 + 1
                    }
                    
                  }
                  
                  if(cada_linha_dados == 2){
                    
                    if(valor_a_comparar > sensiveis_iolanda){
                      valores_maiores_2 <- valores_maiores_2 + 1
                    }else if(valor_a_comparar < sensiveis_iolanda){
                      valores_menores_2 <- valores_menores_2 + 1
                    }else{
                      valores_iguais_2 <- valores_iguais_2 + 1
                    }
                    
                  }
                  
                  if(cada_linha_dados == 3){
                    
                    if(valor_a_comparar > sensiveis_iolanda){
                      valores_maiores_3 <- valores_maiores_3 + 1
                    }else if(valor_a_comparar < sensiveis_iolanda){
                      valores_menores_3 <- valores_menores_3 + 1
                    }else{
                      valores_iguais_3 <- valores_iguais_3 + 1
                    }
                    
                  }
                  
                  
                  
                  
                  if (valor_a_comparar >= limite_inferior & valor_a_comparar <= limite_superior) {
                    
                    quantas_vezes_explica <- quantas_vezes_explica + 1
                  }
                  
                  
                  if (cada_linha %% 3 == 0){
                    
                    
                    
                    if(quantas_vezes_explica > 0){
                      
                      
                      
                      linhas_colunas <- rbind(linhas_colunas, c(raio_espalhamento, tau_zero, k, beta, quantas_vezes_explica, cada_linha_dados))
                    }
                    
                    raio_espalhamento <- raio_espalhamento + cada_inicio
                    
                    
                    
                    
                    
                    quantas_vezes_explica <- 0
                    
                  }
                }
              }
            }
          }
          
          
          if(length(linhas_colunas > 0)){  
            for (i in 1:nrow(linhas_colunas)){
              
              if(linhas_colunas[i,6] == 1){
                
                linhas_colunas[i,6] <- "1R:99S"
                
              }else if(linhas_colunas[i,6] == 2){
                
                linhas_colunas[i,6] <- "50R:50S"
                
              }else if (linhas_colunas[i,6] == 3){
                
                linhas_colunas[i,6] <- "99R:1S"
              }
              
            }  
            
            
            
            colnames(linhas_colunas) <- c("raio", "tau_zero", "k", "beta", "quantas_vezes_explica", "linha")
            
            setwd(paste0("C:/Users/Utilizador/Desktop/Analises artigo Joao/Resultados/", cada_plasmideo))
            
            write.table(linhas_colunas,
                        file = paste0("resultados_ajuste_ger_exp_", cada_plasmideo, "_", cada_densidade, "_erro_", erro_inicial, ".csv"),
                        sep = ",",
                        row.names = FALSE,
                        col.names = TRUE)
          }
        }
      }
      
      
      if (erro_inicial == 2){
        
        #print(cada_densidade)
        
        #print(paste0("valores1: ", valores_menores_1, " - ", valores_iguais_1,  " - ", valores_maiores_1))
        #print(paste0("valores2: ", valores_menores_2,  " - ",  valores_iguais_2, " - ",  valores_maiores_2))
        #print(paste0("valores3: ", valores_menores_3,  " - ", valores_iguais_3, " - ",  valores_maiores_3))
      }
      
    }
  }
}



