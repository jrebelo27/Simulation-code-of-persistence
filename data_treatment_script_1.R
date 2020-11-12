#########################################################################################################
#Calculate the number of susceptible bacteria in each generation, considering duplications
#For power-law data, use k1 and beta. For exponential data, use k1 and k2
#########################################################################################################

for(each_plasmid in c("R1")){
  for (each_density in c("L", "I", "H")){
    for (tau_zero in c(20,30,50,60,70,80,90,100,110,120,130,150,200,250,300,350,400)){
      for(k1 in c(0.015,0.020,0.025,0.030,0.040,0.045,0.050,0.055,0.060,0.065,0.070,0.075,0.080,0.090,0.095,0.1,0.2)){
        for (beta in c(-1.1,-1.2,-1.5,-1.7,-1.8,-1.9,-2.0,-2.1,-2.2,-2.3,-2.4,-2.5,-2.7,-2.9,-3.1,-3.3,-3.5)){
        #for (k2 in c(0.001, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050)){

          setwd(paste0("~/data", each_plasmid, "/", each_density, "/tau_zero", tau_zero, "-k1", k1, "-beta", beta))
          
          for(each_begin in c(0.2, 2)){
            each_end <- each_begin*10
            for(each_frequency in c(1,2,3)){
              name <- paste0("_", each_plasmid, "_", each_density, "_tau_zero", tau_zero, "_k1", k1, "_beta", beta, "_radius_", each_begin, "_", each_end, "_L", each_frequency, ".csv")
              persistent_data <- read.csv(file = paste0("amount_surviving_persistents_generation", name), header = TRUE, sep = ",", stringsAsFactors = FALSE, dec=".")
              susceptible_data <- read.csv(file = paste0("amount_surviving_susceptibles_generation", name), header = TRUE, sep = ",", stringsAsFactors = FALSE, dec=".")
              
              final_path <- c()
              
              for (each_line_data in 1:nrow(susceptible_data)){
                susceptibles_path <- c()
                persistents_path <- c()
                
                for (each_collumn_data in 2:ncol(susceptible_data)){
                  if(each_collumn_data == 2){
                    susceptibles_path <- c(susceptibles_path, susceptible_data[each_line_data, each_collumn_data])
                    persistents_path <- c(persistents_path, persistent_data[each_line_data, each_collumn_data])
                  }else{
                    susceptibles_path <- c(susceptibles_path, ((susceptibles_path[length(susceptibles_path)]*2) + susceptible_data[each_line_data, each_collumn_data]))
                    persistents_path <- c(persistents_path, ((persistents_path[length(persistents_path)]*2) + persistent_data[each_line_data, each_collumn_data]))
                  }
                }
                total_path <- susceptibles_path + persistents_path
                final_path <- rbind(final_path, c(susceptible_data[each_line_data,1], total_path))
              }
              collumn_names <- colnames(susceptible_data)
              final_path <- as.data.frame(final_path)
              colnames(final_path) <- collumn_names
              
              write.table(final_path,
                          file = paste0("total_path_susceptibles_persistents", name),
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

