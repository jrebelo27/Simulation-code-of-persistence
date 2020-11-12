###################################################################################
#Persistent population decay according to a power law
#Persister cells leave the dormant state as soon as the medium becomes detoxified
###################################################################################

#Function of the detoxification radius of each resistant bacterium
resistant_detoxification_radius <- function(generation){
  return(initial_diffusion_rate*sqrt(generation))}

#Exponential decay of the non-persistent population
function_f_1 <- function(generation){
  return(A_1*exp(-k*(generation)))
}

#Power-law decay of the persistent population
function_f_2 <- function(generation){
  return(A_2 * generation^beta)
}

#Function that determines the percentage of susceptible bacteria in the dormant state at a certain time
function_g <- function(generation){
  time <- (generation*time_generation)
  if(time < tau_zero){
    return(1-integrate(function_f_1, lower= 0, upper = time)$value)
  }else{
    return(1- integrate(function_f_1, lower= 0, upper = tau_zero)$value - integrate(function_f_2, lower= tau_zero, upper = time)$value)
  }
}

#Bacteria duplicate every thirty minutes
time_generation <- 30

#Iterate over all conditions and all parameters
for (each_plasmid in c("R1")){
  dir.create(paste0("~/data",each_plasmid))
  for (each_density in c("L", "I", "H")){
    dir.create(paste0("~/data", each_plasmid, "/", each_density))
    for (tau_zero in c(20,30,50,60,70,80,90,100,110,120,130,150,200,250,300,350,400)){
      for(k in c(0.015,0.020,0.025,0.030,0.040,0.045,0.050,0.055,0.060,0.065,0.070,0.075,0.080,0.090,0.095,0.1,0.2)){
        for (beta in c(-1.1,-1.2,-1.5,-1.7,-1.8,-1.9,-2.0,-2.1,-2.2,-2.3,-2.4,-2.5,-2.7,-2.9,-3.1,-3.3,-3.5)){
          
          #Calculate the values of A1 and A2 to ensure that the two functions have the same value when time is equal to tau zero 
          p <- exp(k*tau_zero) - 1 - (tau_zero*k/(1+beta))
          A_1 <- (k*exp(k*tau_zero))/p
          A_2 <- (k/(tau_zero**beta))/p
          
          #Calculate integrals to save time
          integrals <- c()
          for (i in 1:31){
            if(function_g(i) > 0){
              integrals <- c(integrals, function_g(i))
            }else{
              integrals <- c(integrals, 0)
            }
          }
          
          plasmid <- each_plasmid
          density <- each_density
          dir.create(paste0("~/data", each_plasmid, "/", each_density, "/tau_zero", as.character(tau_zero), "-k", as.character(k), "-beta", as.character(beta)))
          setwd("~/original data")
          
          data_excel <- read.csv(file = paste0("data_", plasmid, "_", density, "_averages.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE, dec=",")
          data_excel[,1] <- as.numeric(data_excel[,1])
          data_excel[,2] <- as.numeric(data_excel[,2])
          data_excel[,3] <- as.numeric(data_excel[,3])
          data_excel[,4] <- as.numeric(data_excel[,4])
          
          colnames(data_excel) <- c("initial_resistant", "final_resistant", "initial_susceptible", "final_susceptible")
          
          #Iterate to test diferent diffusion rates
          for (each_initial_diffusion_rate in c(0.2, 2)){
            radius_value <- 0
            incrementing_value <- each_initial_diffusion_rate
            
            for (each_row in 1:3){
              number_repetitions <- 1
              number_conditions <- 10
              initial_diffusion_rate <- radius_value
              incrementing_radius <- incrementing_value
              surviving_susceptibles_generation_final_path <- c()
              surviving_persistents_generation_final_path <- c()
              persistents_detoxified_area_generation_final_path <- c()
              susceptibles_asleep_generation_final_path <- c()
              
              actual_row <- -1
              for (each_condition in 1:number_conditions){
                initial_diffusion_rate <- initial_diffusion_rate + incrementing_radius
                
                for (each_repetition in 1:3){
                  actual_row <- actual_row + 1
                  setwd("~/original data")
                  #Import the dataframes of the distances between each susceptible cell and the closer resistant one
                  new_dataframe <- read.csv(file = paste0("distances_", plasmid, "_", density, "_L", each_row, "_R", each_repetition, ".csv"), header = TRUE, sep = ";", stringsAsFactors = FALSE, dec=",") 
                  assign(paste0(plasmid, "_L", each_row, "_R", each_repetition), new_dataframe)
                  
                  for (each_equal_repetition in 1:number_repetitions){
                    all_distances <- new_dataframe[,1]
                    all_distances <- as.data.frame(all_distances)
                    number_initial_susceptible <- round(data_excel$initial_susceptible[[each_row]])
                    number_generations <- round(log2(data_excel$final_resistant[[each_row]]/data_excel$initial_resistant[[each_row]]))
                    surviving_susceptibles_generation <- c()
                    surviving_persistents_generation <- c()
                    persistents_detoxified_area_generation <- c()
                    susceptibles_asleep_generation <- c()
                    
                    #Iterate as many genarations as the ones completed by resistant cells
                    for (each_generation in 0:(number_generations)){
                      resistant_detoxification_radius_bacterium <- resistant_detoxification_radius(each_generation)
                      awake_susceptibles <- c()
                      
                      #Calculate the number of susceptible cells that leave the dormant state in this cycle
                      if (nrow(all_distances) >= 1){
                        amount_susceptibles_alive <- round(integrals[each_generation+1]*number_initial_susceptible,0)
                        if(amount_susceptibles_alive < nrow(all_distances)){
                          amount_susceptibles_to_awake <- nrow(all_distances) - amount_susceptibles_alive
                          awake_susceptibles <- (1:amount_susceptibles_to_awake)
                        }
                      }
                      
                      awake_susceptibles_distances <- c()
                      
                      if (length(awake_susceptibles) >= 1 && nrow(all_distances)>0){
                        awake_susceptibles_distances <- all_distances[awake_susceptibles, 1]
                        all_distances <- as.data.frame(all_distances[-awake_susceptibles,])
                      }
                      
                      #Susceptibles in a detoxified area leave the dormant state
                      if ((each_generation*time_generation-time_generation) > tau_zero){
                        if(nrow(all_distances) > 0){
                          detoxified_persistents <- which(all_distances[,1] <= resistant_detoxification_radius_bacterium)
                          awake_persistents_distances <- all_distances[detoxified_persistents,]
                          awake_persistents_distances <- as.data.frame(awake_persistents_distances)
                          if(length(detoxified_persistents) > 0){
                            all_distances <- all_distances[-detoxified_persistents,]
                            all_distances <- as.data.frame(all_distances)
                          }
                        }else{
                          awake_persistents_distances <- c()
                          awake_persistents_distances <- as.data.frame(awake_persistents_distances)
                        }
                      }
                      
                      if((each_generation*time_generation-time_generation) > tau_zero){
                        persistents_detoxified_area <- nrow(awake_persistents_distances)
                      }else{
                        persistents_detoxified_area <- 0
                      }
                      
                      surviving_susceptibles_this_generation <- 0
                      surviving_persistents_this_generation <- 0
                      
                      #Intervals before tau zero (all susceptible bacteria are considered non-persistent)
                      if((each_generation*time_generation) <= (tau_zero%/%time_generation * time_generation)){
                        #Verify how many susceptible bacteria survive
                        if (length(awake_susceptibles_distances) != 0) {
                          surviving_susceptibles_this_generation <- sum(awake_susceptibles_distances < resistant_detoxification_radius_bacterium)
                        }
                      #Intervals after tau zero (all susceptible bacteria are considered persistent)  
                      }else if((each_generation*time_generation) > (tau_zero%/%time_generation * time_generation + time_generation)){
                        #Verify how many susceptible bacteria survive
                        if (length(awake_susceptibles_distances) != 0) {
                          surviving_persistents_this_generation <- sum(awake_susceptibles_distances <= resistant_detoxification_radius_bacterium)
                        }
                        if (nrow(awake_persistents_distances) != 0) {
                          surviving_persistents_this_generation  <- surviving_persistents_this_generation + persistents_detoxified_area
                        }
                      #Intervals after tau zero (susceptible bacteria can be considered persistent or non-persistent) 
                      }else{
                        #Verify how many susceptible bacteria survive
                        if (length(awake_susceptibles_distances) != 0) {
                          surviving_susceptibles_this_generation_all <- sum(awake_susceptibles_distances <= resistant_detoxification_radius_bacterium)
                        }
                        #Verify how many susceptible are persistent and how many are non-persistent
                        total_weight <- integrate(function_f_1, lower= 0, upper = tau_zero)$value + integrate(function_f_2, lower= tau_zero, upper = (each_generation*time_generation))$value
                        susceptibles_percentage <- integrate(function_f_1, lower= 0, upper = tau_zero)$value/total_weight
                        surviving_susceptibles_this_generation <- round(surviving_susceptibles_this_generation_all*susceptibles_percentage,0)
                        surviving_persistents_this_generation <- surviving_susceptibles_this_generation_all - surviving_susceptibles_this_generation
                      }
                      
                      surviving_susceptibles_generation <- c(surviving_susceptibles_generation, surviving_susceptibles_this_generation)
                      surviving_persistents_generation <- c(surviving_persistents_generation, surviving_persistents_this_generation)
                      persistents_detoxified_area_generation <- c(persistents_detoxified_area_generation, persistents_detoxified_area)
                      susceptibles_asleep_generation <- c(susceptibles_asleep_generation, nrow(all_distances))
                      
                      if (each_generation == number_generations && (nrow(all_distances) != 0)){
                        surviving_persistents_generation[length(surviving_persistents_generation)] <- (surviving_persistents_generation[length(surviving_persistents_generation)]+nrow(all_distances))
                      }
                    }
                    
                    if (each_equal_repetition == 1){
                      surviving_susceptibles_generation_repetition_path <- surviving_susceptibles_generation
                      surviving_persistents_generation_repetition_path <- surviving_persistents_generation
                      persistents_detoxified_area_generation_repetition_path <- persistents_detoxified_area_generation
                      susceptibles_asleep_generation_repetition_path <- susceptibles_asleep_generation
                      
                    }else{
                      surviving_susceptibles_generation_repetition_path <- rbind(surviving_susceptibles_generation_repetition_path, surviving_susceptibles_generation)
                      surviving_persistents_generation_repetition_path <- rbind(surviving_persistents_generation_repetition_path, surviving_persistents_generation)
                      persistents_detoxified_area_generation_repetition_path <- rbind(persistents_detoxified_area_generation_repetition_path, persistents_detoxified_area_generation)
                      susceptibles_asleep_generation_repetition_path <- rbind(susceptibles_asleep_generation_repetition_path, susceptibles_asleep_generation)
                    }
                  }
                  
                  if(actual_row == 0){
                    surviving_susceptibles_generation_provisional_path <- c()
                    surviving_persistents_generation_provisional_path <- c()
                    persistents_detoxified_area_generation_provisional_path <- c()
                    susceptibles_asleep_generation_provisional_path <- c()
                  }
                  
                  if (actual_row%%10 == 0 && actual_row != 0){
                    surviving_susceptibles_generation_final_path <- unname(surviving_susceptibles_generation_final_path)
                    surviving_susceptibles_generation_provisional_path <- unname(surviving_susceptibles_generation_provisional_path)
                    surviving_persistents_generation_final_path <- unname(surviving_persistents_generation_final_path)
                    surviving_persistents_generation_provisional_path <- unname(surviving_persistents_generation_provisional_path)
                    persistents_detoxified_area_generation_final_path <- unname(persistents_detoxified_area_generation_final_path)
                    persistents_detoxified_area_generation_provisional_path <- unname(persistents_detoxified_area_generation_provisional_path)
                    susceptibles_asleep_generation_final_path <- unname(susceptibles_asleep_generation_final_path)
                    susceptibles_asleep_generation_provisional_path <- unname(susceptibles_asleep_generation_provisional_path)
                    
                    if (actual_row > 10){
                      column_names <- c(1:ncol(surviving_susceptibles_generation_final_path))
                      colnames(surviving_susceptibles_generation_final_path) <- column_names
                      colnames(surviving_susceptibles_generation_provisional_path) <- column_names
                      colnames(surviving_persistents_generation_final_path) <- column_names
                      colnames(surviving_persistents_generation_provisional_path) <- column_names
                      colnames(persistents_detoxified_area_generation_final_path) <- column_names
                      colnames(persistents_detoxified_area_generation_provisional_path) <- column_names
                      colnames(susceptibles_asleep_generation_final_path) <- column_names
                      colnames(susceptibles_asleep_generation_provisional_path) <- column_names
                    }
                    
                    surviving_susceptibles_generation_final_path <- rbind(surviving_susceptibles_generation_final_path, surviving_susceptibles_generation_provisional_path)
                    surviving_susceptibles_generation_final_path <- as.data.frame(surviving_susceptibles_generation_final_path)
                    surviving_susceptibles_generation_provisional_path <- c(surviving_susceptibles_generation_repetition_path)
                    surviving_persistents_generation_final_path <- rbind(surviving_persistents_generation_final_path, surviving_persistents_generation_provisional_path)
                    surviving_persistents_generation_final_path <- as.data.frame(surviving_persistents_generation_final_path)
                    surviving_persistents_generation_provisional_path <- c(surviving_persistents_generation_repetition_path)
                    persistents_detoxified_area_generation_final_path <- rbind(persistents_detoxified_area_generation_final_path, persistents_detoxified_area_generation_provisional_path)
                    persistents_detoxified_area_generation_final_path <- as.data.frame(persistents_detoxified_area_generation_final_path)
                    persistents_detoxified_area_generation_provisional_path <- c(persistents_detoxified_area_generation_repetition_path)
                    susceptibles_asleep_generation_final_path <- rbind(susceptibles_asleep_generation_final_path, susceptibles_asleep_generation_provisional_path)
                    susceptibles_asleep_generation_final_path <- as.data.frame(susceptibles_asleep_generation_final_path)
                    susceptibles_asleep_generation_provisional_path <- c(susceptibles_asleep_generation_repetition_path)
                    
                  }else{
                    surviving_susceptibles_generation_provisional_path <- rbind(surviving_susceptibles_generation_provisional_path, surviving_susceptibles_generation_repetition_path)
                    surviving_susceptibles_generation_provisional_path <- as.data.frame(surviving_susceptibles_generation_provisional_path)
                    surviving_persistents_generation_provisional_path <- rbind(surviving_persistents_generation_provisional_path, surviving_persistents_generation_repetition_path)
                    surviving_persistents_generation_provisional_path <- as.data.frame(surviving_persistents_generation_provisional_path)
                    persistents_detoxified_area_generation_provisional_path <- rbind(persistents_detoxified_area_generation_provisional_path, persistents_detoxified_area_generation_repetition_path)
                    persistents_detoxified_area_generation_provisional_path <- as.data.frame(persistents_detoxified_area_generation_provisional_path)
                    susceptibles_asleep_generation_provisional_path <- rbind(susceptibles_asleep_generation_provisional_path, susceptibles_asleep_generation_repetition_path)
                    susceptibles_asleep_generation_provisional_path <- as.data.frame(susceptibles_asleep_generation_provisional_path)
                  }
                }
              }
              
              surviving_susceptibles_generation_final_path <- as.data.frame(surviving_susceptibles_generation_final_path)
              surviving_persistents_generation_final_path <- as.data.frame(surviving_persistents_generation_final_path)
              persistents_detoxified_area_generation_final_path <- as.data.frame(persistents_detoxified_area_generation_final_path)
              susceptibles_asleep_generation_final_path <- as.data.frame(susceptibles_asleep_generation_final_path)
              out_repetitions <- c(1:3)
              in_repetitions <- c(1:number_repetitions)
              conditions <- seq((initial_diffusion_rate - incrementing_radius*(number_conditions-1)), incrementing_radius*(number_conditions), incrementing_radius)
              names <- c()
              
              for (each_condition in conditions){
                for (each_in_repetitions in in_repetitions){
                  for(each_out_repetitions in out_repetitions){
                    names <- c(names, paste0("R", each_in_repetitions , "_R", each_out_repetitions, "_radius_", each_condition))
                  }
                }
              }
              
              column_names <- c(1:ncol(surviving_susceptibles_generation_final_path))
              colnames(surviving_susceptibles_generation_final_path) <- column_names
              colnames(surviving_susceptibles_generation_provisional_path) <- column_names
              colnames(surviving_persistents_generation_final_path) <- column_names
              colnames(surviving_persistents_generation_provisional_path) <- column_names
              colnames(persistents_detoxified_area_generation_final_path) <- column_names
              colnames(persistents_detoxified_area_generation_provisional_path) <- column_names
              colnames(susceptibles_asleep_generation_final_path) <- column_names
              colnames(susceptibles_asleep_generation_provisional_path) <- column_names
              surviving_susceptibles_generation_final_path <- rbind(surviving_susceptibles_generation_final_path, surviving_susceptibles_generation_provisional_path)
              surviving_persistents_generation_final_path <- rbind(surviving_persistents_generation_final_path, surviving_persistents_generation_provisional_path)
              persistents_detoxified_area_generation_final_path <- rbind(persistents_detoxified_area_generation_final_path, persistents_detoxified_area_generation_provisional_path)
              susceptibles_asleep_generation_final_path <- rbind(susceptibles_asleep_generation_final_path, susceptibles_asleep_generation_provisional_path)
              surviving_susceptibles_generation_final_path <- cbind(names, surviving_susceptibles_generation_final_path, row.names = NULL)
              surviving_persistents_generation_final_path <- cbind(names, surviving_persistents_generation_final_path, row.names = NULL)
              persistents_detoxified_area_generation_final_path  <- cbind(names, persistents_detoxified_area_generation_final_path, row.names = NULL)
              susceptibles_asleep_generation_final_path  <- cbind(names, susceptibles_asleep_generation_final_path, row.names = NULL)
              column_names <- c("Condicao", paste0("G", c(1:(ncol(surviving_susceptibles_generation_final_path)-1))))
              colnames(surviving_susceptibles_generation_final_path) <- column_names
              colnames(surviving_persistents_generation_final_path) <- column_names
              colnames(persistents_detoxified_area_generation_final_path) <- column_names
              colnames(susceptibles_asleep_generation_final_path) <- column_names
              rownames(surviving_susceptibles_generation_final_path) <- NULL
              rownames(surviving_persistents_generation_final_path) <- NULL
              rownames(persistents_detoxified_area_generation_final_path) <- NULL
              rownames(susceptibles_asleep_generation_final_path) <- NULL
              surviving_susceptibles_generation_final_path <- apply(surviving_susceptibles_generation_final_path, 2, as.character)
              surviving_persistents_generation_final_path <- apply(surviving_persistents_generation_final_path, 2, as.character)
              persistents_detoxified_area_generation_final_path <- apply(persistents_detoxified_area_generation_final_path, 2, as.character)
              susceptibles_asleep_generation_final_path <- apply(susceptibles_asleep_generation_final_path, 2, as.character)
              setwd(paste0("~/data", each_plasmid, "/", each_density, "/tau_zero", as.character(tau_zero), "-k", as.character(k), "-beta", as.character(beta)))
             
              write.table(surviving_susceptibles_generation_final_path,
                          file = paste0("amount_surviving_susceptibles_generation_", plasmid, "_", each_density, "_tau_zero", as.character(tau_zero), "_k", as.character(k), "_beta", as.character(beta), "_radius_", incrementing_radius, "_", incrementing_radius*number_conditions, "_L", each_row, ".csv"),
                          sep = ",",
                          row.names = FALSE,
                          col.names = TRUE)
              write.table(surviving_persistents_generation_final_path,
                          file = paste0("amount_surviving_persistents_generation_", plasmid, "_", each_density, "_tau_zero", as.character(tau_zero), "_k", as.character(k), "_beta", as.character(beta), "_radius_", incrementing_radius, "_", incrementing_radius*number_conditions, "_L", each_row, ".csv"),
                          sep = ",",
                          row.names = FALSE,
                          col.names = TRUE)
              write.table(persistents_detoxified_area_generation_final_path ,
                          file = paste0("amount_susceptibles_detoxified_area_generation_", plasmid, "_", each_density, "_tau_zero", as.character(tau_zero), "_k", as.character(k), "_beta", as.character(beta), "_radius_", incrementing_radius, "_", incrementing_radius*number_conditions, "_L", each_row, ".csv"),
                          sep = ",",
                          row.names = FALSE,
                          col.names = TRUE)
              write.table(susceptibles_asleep_generation_final_path ,
                          file = paste0("susceptibles_asleep_generation_", plasmid, "_", each_density, "_tau_zero", as.character(tau_zero), "_k", as.character(k), "_beta", as.character(beta), "_radius_", incrementing_radius, "_", incrementing_radius*number_conditions, "_L", each_row, ".csv"),
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



