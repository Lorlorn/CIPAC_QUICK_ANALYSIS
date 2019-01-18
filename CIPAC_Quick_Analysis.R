#' Title: CIPAC Quick Analysis
#' Author: Lawrence Chen 
#' Date:  20190110
#' Version: v1.0.0 (Faggocity)
#' Purpose: This script is an enhanced one for performing CIPAC evaluations
#' Guide: 
#'  1.) Please DO Load 'Tool Function' into the R Environment.
#'  2.) Load your weak data tables, bro.
library(stringr)
library(tidyr)
library(dplyr)
library(metRology) 
library(outliers) 
library(reshape2)
#### Tool Functions ####
## Horwitz RSDR ##
RSDR_Hor <- function(C)
{
  if(log10(C) >= 0) {
    warning("RSDR_Hor()::Input 'C' > 1, normally the concentration would not exceed 100%.")
  }
  return(2^(1-0.5*log10(C)));
}

## Analysis One Sample ##
#' m_lab: labels, could be factors.
#' m_response: responses, concentrations.
#' unit_ratio: unit ratio, 1/(unit(solute)/unite(solvent)), eg. 1 / (g/kg) = 1000
#' Utility: function will eliminate NA from any position of input along with the other.
GetSampleSummary <- function(m_lab, m_response, unit_ratio = 1000)
{
  # check inputs & eliminate NAs
  if (length(m_lab) != length(m_response))  {
    warning("GetSampleSummary()::Input vectors length INequivalent, STUPID!");
    return(NA);
  }
  posi = unique(c(which(is.na(m_lab)), which(is.na(m_response))));
  if(length(posi) != 0) {
    m_lab = m_lab[-posi];
    m_response = m_response[-posi];
  }
  m_lab = as.character(m_lab) %>% as.numeric;
  
  # define output obj.
  output.items <- c("Xm", "L", "Sr", "SR", "RSDr", "RSDR", "RSDR(Hor)", "Hor. Rat.");
  output <- as.data.frame(matrix(NA, ncol = length(output.items)));
  colnames(output) = output.items;
  
  # compute universal statistics
  L = length(unique(m_lab)); # L
  N = length(m_lab); # Total Number
  total.mean = mean(m_response); # Yo
  
  # compute group-specific statistics
  lab.K = tapply(m_response, m_lab, length); # K
  lab.mean = tapply(m_response, m_lab, mean); # Group Mean
  lab.var = tapply(m_response, m_lab, var);
  lab.sd = sqrt(lab.var);
  
  # compute more...
  KH = 1 / mean(1 / lab.K); # K Harmonic
  dev = lab.mean - total.mean; # Dev(group, total mean)
  MSU = KH*sum(dev^2) / (L-1); #BQI equation 11: MSamong using Unweighted SSamong
  Kminus1 = lab.K - 1;
  MSE = sum(lab.var * Kminus1) / sum(Kminus1);
  varA = max( 0, (MSU - MSE) / KH);
  
  # compute output region -> easier for dummies...
  sr = sqrt(MSE) # Rpeatability
  sR = sqrt((MSU / KH) + ((KH - 1) * MSE / KH)) # Reproducibility
  RSDr = sr / total.mean * 100;
  RSDR = sR / total.mean * 100;
  pRSDR = RSDR_Hor(total.mean / unit_ratio);
  hor_ratio = RSDR / pRSDR;
  r = varA / (sR^2)
  # gain output
  output$Xm = total.mean;
  output$L = L;
  output$Sr = sr;
  output$SR = sR;
  output$RSDr = RSDr;
  output$RSDR = RSDR;
  output$`RSDR(Hor)` = pRSDR;
  output$`Hor. Rat.` = hor_ratio;
  
  return(t(output));
}

## Column Extreme Value Pin Point ##
ReturnVecExtremePosition <- function(x)
{
  profile <- (x - mean(x, na.rm = T))^2;
  return(which.max(profile))
}
#### Working Region ####
# Load dataset into environment, bro.
## Experiment Settings
sample_number <- 5; # Number of Sample, one sample a column.
day <- 2; # Number of experiment performing day by each sample
ori_matrix <- matrix(scan(), ncol = sample_number * day + 1, byrow = T);

## Data Tidyr
lab <- ori_matrix[, 1] %>% as.factor; # fist column is default to be 'labels' for grouping factors, s.a. Lab. Num..
data_by_sample <- matrix(c(ori_matrix[, 2:ncol(ori_matrix)]), ncol = sample_number);
lab_tags <- levels(lab);
sample_tags <- paste("Sample", 1:(ncol(data_by_sample)), sep = "_") %>% as.factor;
colnames(data_by_sample) <- sample_tags;

## Compute CIPAC Summary for each
final_output <- c();
for(i in 1:ncol(data_by_sample)) {
  final_output <- cbind(final_output, GetSampleSummary(rep(lab, day), data_by_sample[, i]));
}
colnames(final_output) <- sample_tags;
print(final_output);

## Outlier Detection
#' Grubb's Test
#' Compute Set Mean
#' Compute one-sided outliers
# Cleaning Data
data_by_sample_lab <- cbind(Lab = rep(lab, day), data_by_sample) %>% as.data.frame;
data_by_sample_lab$Lab <- data_by_sample_lab$Lab %>% as.factor
working_data <- melt(data_by_sample_lab, id.vars = "Lab") %>% group_by(variable, Lab) %>% summarise(Mean = mean(value));

## Perform Grubb's Test
cat("\nPerforming Grubb's Test:\n")
output_data <- list(SampleAveragePerLab = matrix(NA, ncol = sample_number, nrow = length(lab_tags)),
                    Straggler = matrix(1, ncol = sample_number, nrow = length(lab_tags)),
                    Outlier = matrix(1, ncol = sample_number, nrow = length(lab_tags)));
output_list <- list(Straggler = c(),
                    Outlier = c());

grubbs_two_sided_result <- c();
for(i in 1:length(sample_tags))  {
  output_data$SampleAveragePerLab[, i] = working_data[which(working_data$variable == sample_tags[i]),]$Mean;
  temp_result = grubbs.test(working_data[which(working_data$variable == sample_tags[i]),]$Mean, type =10) %>% unlist;
  grubbs_two_sided_result <- cbind(grubbs_two_sided_result, temp_result);
  
  p = as.numeric(temp_result[4]);
  extreme_tag = temp_result[3] %>% strsplit(split = " ")  %>% unlist %>% .[1];
  if (p < 0.05)  {
    extreme_row = ifelse( "highest" == extreme_tag, 
                          which.max(output_data$SampleAveragePerLab[, i]), 
                          which.min(output_data$SampleAveragePerLab[, i]));
    output_data$Straggler[extreme_row, i] = NA;
    m_staggler = data.frame(Lab = lab_tags[extreme_row], Sample = sample_tags[i]);
    output_list$Straggler = rbind(output_list$Straggler, m_staggler);
    if (p < 0.01) {
      output_data$Outlier[extreme_row, i] = NA;
      m_outlier = data.frame(Lab = lab_tags[extreme_row], Sample = sample_tags[i]);
      output_list$Outlier = rbind(output_list$Outlier, m_outlier);
    }
  }
  cat(as.character(sample_tags[i]), ": p-value = ", round(p, 5), "| Ha = ", temp_result[3], "\n");
}
print(output_list);

## Final Ouput 2
# Table formation
final_output_2 <- list(Without_Outliers = c(),
                       Without_Stragglers = c());
working_data <- list(Straggler = data_by_sample,
                      Outlier = data_by_sample);
n_stragglers <- nrow(output_list$Straggler);
n_outliers <- nrow(output_list$Outlier);

# Compute Stragglers/Outliers
if (n_stragglers != 0)  {
  cat("Number of Stragglers:", n_stragglers);
  for (i in 1:n_stragglers) {
    row_factor = which(data_by_sample_lab$Lab == as.character(output_list$Straggler$Lab[i]));
    col_factor = which(sample_tags == output_list$Straggler$Sample[i]);
    working_data$Straggler[row_factor, col_factor] <- NA;
  }
}
if (n_outliers != 0)  {
  cat("\nNumber of Outliers:", n_outliers);  
  for (i in 1:n_outliers) {
    row_factor = which(data_by_sample_lab$Lab == as.character(output_list$Outlier$Lab[i]));
    col_factor = which(sample_tags == output_list$Outlier$Sample[i]);
    working_data$Outlier[row_factor, col_factor] <- NA;
  }
}

for(i in 1:ncol(data_by_sample)) {
  final_output_2$Without_Stragglers <- cbind(final_output_2$Without_Stragglers, GetSampleSummary(rep(lab, day), working_data$Straggler[, i]));
  final_output_2$Without_Outliers <- cbind(final_output_2$Without_Outliers, GetSampleSummary(rep(lab, day), working_data$Outlier[, i]));
}
colnames(final_output_2$Without_Outliers) <- sample_tags;
colnames(final_output_2$Without_Stragglers) <- sample_tags;
cat("\n")
print(final_output_2)

#### UNUSE REGION ####
#' Mandel's h for each day
#' abs critical values : [95%, 99%] = [1.84, 2.27]

#' Mandel's k for lab mean
#' abs critical values : [95%, 99%] = [1.51, 1.76]
# mandel_h_test_obj <- list(
#   Ori_Data = ori_matrix[, 2:ncol(ori_matrix)],
#   CI_95 = list(
#     Template = matrix(1, nrow = nrow(ori_matrix), ncol = (ncol(ori_matrix) - 1)),
#     Outliers = rep(0, each = (ncol(ori_matrix) - 1))
#   ),
#   CI_99 = list(
#     Template = matrix(1, nrow = nrow(ori_matrix), ncol = (ncol(ori_matrix) - 1)),
#     Outliers = rep(0, each = (ncol(ori_matrix) - 1))
#   ))
# mandel_h_test_results <- apply(mandel_h_test_obj$Ori_Data, 2, function(x){ mandel.h(x) %>% unlist });
# 
# for(i in 1:ncol(mandel_h_test_results)) {
#   for(j in 1:nrow(mandel_h_test_results)) {
#     if(is.nan(mandel_h_test_results[j, i]))  {
#       next;
#     }
#     if (abs(mandel_h_test_results[j, i]) > 1.84)  {
#       mandel_h_test_obj$CI_95$Template[j, i] = NA;
#       if (abs(mandel_h_test_results[j, i]) > 2.27) {
#         mandel_h_test_obj$CI_99$Template[j, i] = NA;
#       } 
#     } 
#   }
# }

# Plot




