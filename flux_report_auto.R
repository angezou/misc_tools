# Step 1: load appropriate R libraries

library(readxl)
library(dplyr)
library(tidyr)
library(openxlsx)
library(plotrix)

# Step 2: Run lines 11 to 38 (this loads in the function flux_calc_auto) which allows for automated generation of N₂O-N (mg m^-2 h^-1), CH₄-C (mg m^-2 h^-1), and CO₂-N (mg m^-2 h^-1) values

# Note: In order for the function to run properly, the input file format MUST follow the format found in columns A to J of the tab "Beauregard" in the example input file: input_PALCCA_test.xlsx
# This includes column names, white spaces, and special characters

flux_calc_auto <- function(input_path, input_sheet){
  
  input_df <- read_excel(input_path)
  
  input_df$`N₂O-N (μg/L)` <- (input_df$`N₂O (ppm)`*44)/(0.082*(input_df$`Temperature (°C)` + 273.15))
  input_df$`CH₄-C (μg/L)` <- (input_df$`CH₄ (ppm)`*16)/(0.082*(input_df$`Temperature (°C)` + 273.15))
  input_df$`CO₂-C (μg/L)`<- (input_df$`CO₂ (ppm)`*44)/(0.082*(input_df$`Temperature (°C)` + 273.15))

  slope_df <- input_df %>% 
    drop_na(Treatment) %>%
    mutate(group = (row_number() - 1) %/% 4) %>% # make groups of 4 rows
    # compute slope per group
    group_by(group) %>%
    mutate(`N₂O-N (µg L⁻¹ min⁻¹)` = coef(lm(`N₂O-N (μg/L)` ~ `Time points (min)`))[2]) %>%
    mutate(`CH₄-C (µg L⁻¹ min⁻¹)` = coef(lm(`CH₄-C (μg/L)` ~ `Time points (min)`))[2]) %>%
    mutate(`CO₂-C (µg L⁻¹ min⁻¹)` = coef(lm(`CO₂-C (μg/L)` ~ `Time points (min)`))[2]) %>%
    ungroup()
  
  flux_df <- slope_df %>%
    mutate(
      `N₂O-N (mg m-2 h-1)` = `N₂O-N (µg L⁻¹ min⁻¹)`*`Chamber volume (L)`/`Chamber area (m²)`*60/1000,
      `CH₄-C (mg m-2 h-1)` = `CH₄-C (µg L⁻¹ min⁻¹)`*`Chamber volume (L)`/`Chamber area (m²)`*60/1000,
      `CO₂-N (mg m-2 h-1)` = `CO₂-C (µg L⁻¹ min⁻¹)`*`Chamber volume (L)`/`Chamber area (m²)`*60/1000
    )
  
  return(subset(flux_df, select = -c(group)))
  
}

# Step 3: Read in input file that contains flux values that need to be converted, additionally specify the list of excel tabs that contain flux values needing to be converted, the names of the newly updated excel spreadsheets can be used to replace the placeholder values "B1", "A2", "C1", and the names of input tabs can be used to replace the placeholder values "Beauregard"

# Note: the number of tabs can be adjusted according to the input file, if there are fewer than 3 input tabs, remove the additional values in the "sheets" list, if there are more than 3 input tabs, feel free to add more lines to the "sheets" list, just don't forget to add a comma to the pre-existing list before adding a new value

input_path = "soil/flux_auto/Input_PALCCA_test.xlsx"

# Step 4: Generate the converted flux values using the flux_calc_auto function on each excel tab listed

sheets <- list(
  "B1" = flux_calc_auto(input_path, "Beauregard"),
  "A2" = flux_calc_auto(input_path, "Beauregard"),
  "C1" = flux_calc_auto(input_path, "Beauregard")
)

# Step 5: write to excel output, the file path can be subject to change depending on user preferences 

write.xlsx(sheets, "soil/flux_auto/output_test.xlsx") 


# Example:

input_path = "soil/flux_auto/Ridgetown - Field 1 (AD7407) - Total flux results calculation - Summer 2025_test.xls"

sheets <- list(
  #"B1" = flux_calc_auto(input_path, "Beauregard"),
  #"A2" = flux_calc_auto(input_path, "Beauregard"),
  "Ridge_AD7407" = flux_calc_auto(input_path, "Ridge. AD7407")
)

write.xlsx(sheets, "soil/flux_auto/Ridge_AD7407_output_test.xlsx") 


# Step 6: Run lines 62 to 105 to load in the function flux_val_gen, which automates the interpoloation of flux values 

#Note: in order for the function to function properly, the input file format MUST follow the format found in the example input file "cumul_test_input.txt, this includes column names, special characters, and white spaces

flux_val_gen <- function(input_path){
  
  cumul_input <- read_excel(input_path)
  
  ncol_ci <- ncol(cumul_input)
  
  dates <- as.Date(as.numeric(names(cumul_input[5:ncol_ci])), origin = "1899-12-30")
  
  day_diff <- diff(dates)
  day_diff <- as.numeric(day_diff)
  
  new_names <- paste0("rate_", seq_along(day_diff))
  
  rc_df <- cumul_input[1:4]
  
  rc_df[new_names] <- lapply(seq_along(day_diff), function(i) {
    -(cumul_input[[i+4]]*24 - cumul_input[[i+5]]*24) / day_diff[i]
  })
  
  ncol_rc <- ncol(rc_df)
  
  li_df <- cumul_input[1:4]
  
  for (i in seq_along(rc_df[5:ncol_rc-1])) {
    
    lin_interp1 <- cumul_input[5:ncol_ci][[i]] * 24
    
    li_df[[paste0("day_",i)]] <- lin_interp1
    
    new_col <- lin_interp1 + rc_df[5:ncol_rc][[i]]
    
    li_df[[ paste0("day_", i,"_2") ]] <- new_col
    
    for (j in 3:day_diff[i]) {
      new_col <- new_col + rc_df[5:ncol_rc][[i]]
      li_df[[ paste0("day_", i,"_",j) ]] <- new_col
    }
    
  }
  
  li_df[[paste0("day_",length(day_diff)+1)]] <- cumul_input[[ncol_ci]] * 24
  
  li_df$`Sum N₂O-N (mg m ⁻²)` <- rowSums(li_df[grep("^day_", names(li_df))])
  
  return(li_df)
}

# Step 7: generates sum of interpolated flux values using the flux_val_gen function and stores them in the li_df variable

li_df <- flux_val_gen("soil/flux_auto/cumul_test_input- Champs 1 - 2025.xlsx")

# Step 8: Generates the summary file (average and standard deviation) of N2O flux for each treatment

sum_li <- li_df %>%
  group_by(Gas, Treatments) %>%
  summarise(
    avg = mean(`Sum N₂O-N (mg m ⁻²)`, na.rm =TRUE)/100,
    std_err = std.error(`Sum N₂O-N (mg m ⁻²)`, na.rm = TRUE)/100
  )

# Step 9: Write to excel output files, change file path according to user needs

write.xlsx(li_df,"soil/flux_auto/li_df.xlsx")
write.xlsx(sum_li,"soil/flux_auto/sum_li.xlsx")



