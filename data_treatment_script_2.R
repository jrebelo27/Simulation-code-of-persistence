##########################################################################
#Compare the experimental and computational results
#For power-law data, use k1 and beta. For exponential data, use k1 and k2
##########################################################################

for (each_plasmid in c("R1")){

  for (each_density in c("I", "H")){
    if(each_density == "I"){
      density_name <- "Intermediate" 
      data_lines <- c(1, 2, 3)
    }
    if(each_density == "H"){
      density_name <- "High" 
      data_lines <- c(1, 2, 3)
    }
    
    for (initial_error in c(2, 4)){
      lines_collumns <- c()
      
      for (each_begin in c(0.2, 2)){
        each_end <- each_begin * 10
        
        for (tau_zero in c(20,30,50,60,70,80,90,100,110,120,130,150,200,250,300,350,400)){
          for(k1 in c(0.015,0.020,0.025,0.030,0.040,0.045,0.050,0.055,0.060,0.065,0.070,0.075,0.080,0.090,0.095,0.1,0.2)){
            for (beta in c(-1.1,-1.2,-1.5,-1.7,-1.8,-1.9,-2.0,-2.1,-2.2,-2.3,-2.4,-2.5,-2.7,-2.9,-3.1,-3.3,-3.5)){
            #for (k2 in c(0.001, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050)){
              
              for(each_line_data in data_lines){
                setwd(paste0("~/data", each_plasmid, "/", each_density, "/tau_zero", tau_zero, "-k", k, "-beta", beta))
                
                name <- paste0("_", each_plasmid, "_", each_density, "_tau_zero", tau_zero, "_k", k, "_beta", beta, "_raio_", each_begin, "_", each_end, "_L", each_line_data, ".csv")
                entry_name <- paste0("total_path_susceptibles_persistents", name)
                simulations_data <- read.csv(file = entry_name, header = TRUE, sep = ",", stringsAsFactors = FALSE, dec=".")

                setwd("~/original data")
                
                entry_name_susceptibles <- paste0("data_", each_plasmid, "_", each_density, "_averages.csv")
                original_data <- read.csv(file = entry_name_susceptibles, header = TRUE, sep = ",", stringsAsFactors = FALSE, dec=".")
                original_data[,4] <- as.numeric(original_data[,4])
                
                experimental_susceptibles <- original_data[each_line_data, 4]
                experimental_final_resistants <- original_data[each_line_data, 2]
                initial_resistants <- original_data[each_line_data, 1]
                number_generations_experiments <- round(log2(experimental_final_resistants/initial_resistants), 6)
                number_generations_simulations <- (round(log2(experimental_final_resistants/initial_resistants), 0))+1
                generations_difference <- number_generations_experiments - number_generations_simulations
                diffusion_radius <- each_begin
                times_that_explained <- 0

                if (each_begin == 0.2){
                  start_line <- 1
                }else{
                  start_line <- 4
                  diffusion_radius <- diffusion_radius + each_begin
                }
                
                for (each_line in start_line:nrow(simulations_data)){
                  value_to_compare <- simulations_data[each_line, ncol(simulations_data)]
                  value_to_compare <- round(value_to_compare * (2 ** generations_difference), 0)
                  error <- initial_error
                  inferior_limit <- experimental_susceptibles/error
                  upper_limit <- experimental_susceptibles*error
                  
                  if (value_to_compare >= inferior_limit & value_to_compare <= upper_limit){
                    times_that_explained <- times_that_explained + 1
                  }
                  if (each_line %% 3 == 0){
                    if(times_that_explained > 0){
                      lines_collumns <- rbind(lines_collumns, c(diffusion_radius, tau_zero, k, beta, times_that_explained, each_line_data))
                    }
                    diffusion_radius <- diffusion_radius + each_begin
                    times_that_explained <- 0
                  }
                }
              }
            }
          }
          if(length(lines_collumns > 0)){  
            for (i in 1:nrow(lines_collumns)){
              if(lines_collumns[i,6] == 1){
                lines_collumns[i,6] <- "1R:99S"
              }else if(lines_collumns[i,6] == 2){
                lines_collumns[i,6] <- "50R:50S"
              }else if (lines_collumns[i,6] == 3){
                lines_collumns[i,6] <- "99R:1S"
              }
            }  
            colnames(lines_collumns) <- c("radius", "tau_zero", "k", "beta", "times_that_explained", "line")
            
            setwd(paste0("~/Results/", each_plasmid))
            
            write.table(lines_collumns,
                        file = paste0("results_adjust_ger_exp", each_plasmid, "_", each_density, "_error_", initial_error, ".csv"),
                        sep = ",",
                        row.names = FALSE,
                        col.names = TRUE)
          }
        }
      }
    }
  }
}