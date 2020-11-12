############################################################################################################
#Construction of the dataframes with the distances between susceptible cells and the nearest resistant cell
#These dataframes will be used in the simulations
############################################################################################################

library(dplyr)

increment <- 200
susceptibles_dataframes <- 500

#Number of diferent dataframes
number_repetitions <- 3
plate_radius <- 45000
plate_area <- (pi*(plate_radius^2))
excel_data <- read.csv(file = "R-L.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE, dec=",")
data_average <- c()

#Treat the data to use the mean values
for (each_line in seq(1,9,3)){
  initial_resistants <- 0
  final_resistants <- 0
  initial_susceptibles <- 0
  final_susceptibles <- 0
  
  for (each_line_2 in each_line:(each_line + 2)){
    initial_resistants <- initial_resistants + excel_data[each_line_2, 1]
    final_resistants <- final_resistants + excel_data[each_line_2, 2]
    initial_susceptibles <- initial_susceptibles + excel_data[each_line_2, 3]
    final_susceptibles <- final_susceptibles + excel_data[each_line_2, 4]
  }
  data_average <- rbind(data_average, c(round(initial_resistants/3), round(final_resistants/3), round(initial_susceptibles/3), round(final_susceptibles/3)))
}

data_average <- as.data.frame(data_average)
excel_data <- data_average
colnames(excel_data) <- c("initial_resistants", "final_resistants", "initial_susceptibles", "final_susceptibles")
write.table(excel_data,
            file = paste0("data_R1_L_average.csv"),
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)

#Function to generate random coordinates within a given radius
random_coordinates <- function (R){
  r = R * sqrt(rbeta(1,1,1))
  theta = rbeta(1,1,1) * 2 * pi
  return(list(c(round(r * cos(theta), 0), round(r * sin(theta), 0))))
}

for (each_line in 1:3){
  resistants_initial_number <- excel_data[each_line, 1]
  susceptible_initial_number <- excel_data[each_line, 3]
  bacterium_diameter <- plate_radius*2 * sqrt(1/(susceptible_initial_number+resistants_initial_number))
  bacterium_radius <- round(bacterium_diameter/2, 0)
  
  for (each_repetition in 1:number_repetitions){
    #Create smaller dataframes for producers to optimize code speed
    for (each_name in seq(-plate_radius+increment,plate_radius-increment,increment)){
      new_dataframe <- data.frame(matrix(vector(), 0, 2,
                                          dimnames=list(c(), c("X", "Y"))),
                                   stringsAsFactors=F)
      assign(paste("resistants_", each_name, sep = ''), new_dataframe)
    }
    
    for (each_initial_resistant in 1:resistants_initial_number){
      coordinates_p <- random_coordinates(plate_radius)
      coordinate_x <- coordinates_p[[1]][1]
      coordinate_y <- coordinates_p[[1]][2]
      if(coordinate_x == 45000){
        name_x <- coordinate_x-increment
      }else if (coordinate_x == -45000){
        name_x <- coordinate_x+increment
      }else{
        if (coordinate_x == 0){
          name_x <- 0
        }else{
          name_x <- coordinate_x - (abs(coordinate_x)%%increment)*(coordinate_x/abs(coordinate_x))}
      }
      
      used_dataframe <- get(paste("resistants_", name_x, sep=""))
      coordinates <- c(coordinate_x, coordinate_y)
      used_dataframe <- rbind(used_dataframe, coordinates)
      colnames(used_dataframe) <- c("X", "Y")
      assign(paste("resistants_", name_x, sep = ''), used_dataframe)
    }
    
    #Define the positions of the sensitive
    susceptibles_by_dataframe <- susceptible_initial_number%/%susceptibles_dataframes
    remaining_susceptibles <- susceptible_initial_number%%susceptibles_dataframes
    
    for (each_dataframe in 1:susceptibles_dataframes){
      new_dataframe <- data.frame(matrix(vector(), 0, 2,
                                          dimnames=list(c(), c("X", "Y"))),
                                   stringsAsFactors=F)
      if (each_dataframe != susceptibles_dataframes){
        provisional_positions <- data.frame(matrix(vector(), 0, 2,
                                                  dimnames=list(c(), c("X", "Y"))),
                                           stringsAsFactors=F)
        provisional_positions_other <- data.frame(matrix(vector(), 0, 2,
                                                         dimnames=list(c(), c("X", "Y"))),
                                                  stringsAsFactors=F)
        line <- 1
        line_other <-1
        
        for (each_susceptibles in 1:susceptibles_by_dataframe){
          coordinates_s <- random_coordinates(plate_radius)
          provisional_positions[line, 1] <- coordinates_s[[1]][1]
          provisional_positions[line, 2] <- coordinates_s[[1]][2]
          if (line%%300 == 0 && line!=0 ){
            provisional_positions_other <- rbind(provisional_positions_other, provisional_positions)
            provisional_positions <- data.frame(matrix(vector(), 0, 2,
                                                      dimnames=list(c(), c("X", "Y"))),
                                               stringsAsFactors=F)
            line <- 0
          }
          if (line_other%%10000 == 0 && line_other!=0 ){
            print(paste("each_susceptibles_inicial: ", each_susceptibles))
            new_dataframe <- rbind(new_dataframe, provisional_positions_other)
            provisional_positions_other <- data.frame(matrix(vector(), 0, 2,
                                                             dimnames=list(c(), c("X", "Y"))),
                                                      stringsAsFactors=F)
            line_other <- 0
          }
          line <- line +1
          line_other <- line_other + 1
        }
        provisional_positions_other <- rbind(provisional_positions_other, provisional_positions)
        new_dataframe <- rbind(new_dataframe, provisional_positions_other)
        assign(paste("sensiveis_", each_dataframe, sep = ''), new_dataframe)
      }else{
        provisional_positions <- data.frame(matrix(vector(), 0, 2,
                                                  dimnames=list(c(), c("X", "Y"))),
                                           stringsAsFactors=F)
        provisional_positions_other <- data.frame(matrix(vector(), 0, 2,
                                                         dimnames=list(c(), c("X", "Y"))),
                                                  stringsAsFactors=F)
        line <- 1
        line_other <-1
        
        for (each_susceptibles in 1:(susceptibles_by_dataframe + remaining_susceptibles)){
          coordinates_s <- random_coordinates(plate_radius)
          provisional_positions[line, 1] <- coordinates_s[[1]][1]
          provisional_positions[line, 2] <- coordinates_s[[1]][2]
          if (line%%300 == 0 && line!=0 ){
            provisional_positions_other <- rbind(provisional_positions_other, provisional_positions)
            provisional_positions <- data.frame(matrix(vector(), 0, 2,
                                                      dimnames=list(c(), c("X", "Y"))),
                                               stringsAsFactors=F)
            line <- 0
          }
          if (line_other%%10000 == 0 && line_other!=0 ){
            new_dataframe <- rbind(new_dataframe, provisional_positions_other)
            provisional_positions_other <- data.frame(matrix(vector(), 0, 2,
                                                             dimnames=list(c(), c("X", "Y"))),
                                                      stringsAsFactors=F)
            line_other <- 0
          }
          line <- line +1
          line_other <- line_other + 1
        }
        
        provisional_positions_other <- rbind(provisional_positions_other, provisional_positions)
        new_dataframe <- rbind(new_dataframe, provisional_positions_other)
        assign(paste("susceptibles_", each_dataframe, sep = ''), new_dataframe)
      }
    }
    
    provisional_distances <- c()
    distances <- c()
    
    for (each_dataframe in 1:susceptibles_dataframes){
      used_susceptibles_dataframe <- (mget(paste("susceptibles_", c(each_dataframe), sep="")))
      used_susceptibles_dataframe <- bind_rows(used_susceptibles_dataframe)
      
      for (each_susceptibles in 1:nrow(used_susceptibles_dataframe)){
        x_susceptible <- used_susceptibles_dataframe[each_susceptibles, 1]
        y_susceptible <- used_susceptibles_dataframe[each_susceptibles, 2]
        if (x_susceptible != 0){
          x_susceptible <- x_susceptible - (abs(x_susceptible)%%increment)*(x_susceptible/abs(x_susceptible))
        }
        if (x_susceptible == plate_radius){
          x_susceptible <- plate_radius - increment
        }
        if (x_susceptible == -plate_radius){
          x_susceptible <- -plate_radius + increment
        }
        dataframes_compared_resistants <- data.frame(matrix(vector(), 0, 2,
                                                            dimnames=list(c(), c("X", "Y"))),
                                                     stringsAsFactors=F)
        i <- 1
        
        while (nrow(dataframes_compared_resistants) == 0){
          x_inferior <- x_susceptible - i*increment
          x_superior <- x_susceptible + i*increment
          
          if (x_superior >= plate_radius){
            x_superior <- plate_radius - increment
          }
          if (x_inferior <= -plate_radius){
            x_inferior <- -plate_radius + increment
          }
          
          numbers_dataframes <- seq(x_inferior, x_superior, increment)
          dataframes_compared_resistants <- (mget(paste("resistants_", numbers_dataframes, sep="")))
          dataframes_compared_resistants <- bind_rows(dataframes_compared_resistants)
          colnames(dataframes_compared_resistants) <- c("X", "Y")
          y_inferior <- y_susceptible - i*increment
          y_superior <- y_susceptible + i*increment
          dataframes_compared_resistants <- dataframes_compared_resistants[dataframes_compared_resistants$Y <= y_superior & dataframes_compared_resistants$Y >= y_inferior, ]
          i <- i + 1
        }
        
        resistant_distance <- plate_radius*2
        
        for (each_compared_resistant in 1:nrow(dataframes_compared_resistants)){
          distance <- sqrt((dataframes_compared_resistants[each_compared_resistant,1] - used_susceptibles_dataframe[each_susceptibles,1])^2 + 
                              (dataframes_compared_resistants[each_compared_resistant,2] - used_susceptibles_dataframe[each_susceptibles,2])^2)
          if (distance <= resistant_distance) {
            #Update the distance     
            resistant_distance <- distance}
        }
        
        resistant_distance <- round(resistant_distance, 0)
        
        if (each_susceptibles%%300 == 0){
          provisional_distances <- rbind(provisional_distances, c(resistant_distance))
          distances <- rbind(distances, provisional_distances)
          provisional_distances <- c()
        }else{
          provisional_distances <- rbind(provisional_distances, c(resistant_distance))
        }
      }
    }
    
    distances <- rbind(distances, provisional_distances)
    provisional_distances <- c()
    distances <- as.data.frame(distances, stringsAsFactors = F)
    colnames(distances) <- c("Closer_resistant")
    distances$Closer_resistant <- as.integer(distances$Closer_resistant)
    write.table(distances,
                file = paste0("distances_pBR322_L_L",each_line, "_R",each_repetition,".csv"),
                sep = ",",
                row.names = FALSE,
                col.names = TRUE)
  }
}
    
    