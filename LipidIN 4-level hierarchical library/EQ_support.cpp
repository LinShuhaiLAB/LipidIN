#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAX_LINE_LENGTH 1024
#define MAX_TITLE_LENGTH 256

typedef struct {
  int index;
  double mz;
  double rt;
  double intensity;
  char title[MAX_TITLE_LENGTH];
  double score1;
  double score2;
  char oxfa_fragment[MAX_TITLE_LENGTH];
} Record;

// [[Rcpp::export]]
void extract_oxfa_fragment(const std::string& title, std::string& oxfa_fragment) {
  size_t start = title.find("_DES_");
  if (start != std::string::npos) {
    oxfa_fragment = title.substr(start + 5);
  } else {
    oxfa_fragment.clear();
  }
}


// [[Rcpp::export]]
void extract_levels(const std::string& title, int& level1_3, int& level2_4) {
  level1_3 = 0;
  level2_4 = 0;
  size_t level_pos = title.find("LEVEL");
  
  while (level_pos != std::string::npos) {
    int level;
    sscanf(title.c_str() + level_pos, "LEVEL%d", &level);
    
    if (level == 1 || level == 3) {
      if (level > level1_3) level1_3 = level;
    } else if (level == 2 || level == 4) {
      if (level > level2_4) level2_4 = level;
    }
    level_pos = title.find("LEVEL", level_pos + 5);
  }
}

// [[Rcpp::export]]
bool is_valid_record(const int level1_3, const int level2_4) {
  return (level1_3 > 0 && level2_4 > 0);
}

// [[Rcpp::export]]
void process_file(const std::string& input_file) {
  FILE *fp_in = fopen(input_file.c_str(), "r");
  if (!fp_in) {
    Rcpp::stop("Error opening input file");
  }
  
  // Generate output file name
  std::string output_file = input_file;
  size_t dot = output_file.rfind('.');
  if (dot != std::string::npos) output_file.erase(dot);
  output_file += "_processed.csv";
  
  FILE *fp_out = fopen(output_file.c_str(), "w");
  if (!fp_out) {
    fclose(fp_in);
    Rcpp::stop("Error opening output file");
  }
  
  char line[MAX_LINE_LENGTH];
  fgets(line, MAX_LINE_LENGTH, fp_in);  // Read the header
  fprintf(fp_out, "%s", line);  // Write the header to the output file
  
  Record prev_record;
  bool first_record = true;
  
  double max_score1_level1_3 = 0.0, max_score2_level1_3 = 0.0;
  double max_score1_level2_4 = 0.0, max_score2_level2_4 = 0.0;
  int valid_level1_3 = 0, valid_level2_4 = 0;
  
  while (fgets(line, MAX_LINE_LENGTH, fp_in)) {
    Record current_record;
    sscanf(line, "%d,%lf,%lf,%lf,%[^,],%lf,%lf",
           &current_record.index, &current_record.mz, &current_record.rt, &current_record.intensity,
           current_record.title, &current_record.score1, &current_record.score2);
    
    // Skip records where score1 is less than 0.5
    if (current_record.score1 < 0.5) {
      continue;
    }
    
    std::string oxfa_fragment;
    extract_oxfa_fragment(current_record.title, oxfa_fragment);
    strncpy(current_record.oxfa_fragment, oxfa_fragment.c_str(), MAX_TITLE_LENGTH);
    
    int level1_3_curr = 0, level2_4_curr = 0;
    extract_levels(current_record.title, level1_3_curr, level2_4_curr);
    
    if (first_record || (prev_record.index == current_record.index &&
        prev_record.mz == current_record.mz &&
        prev_record.rt == current_record.rt &&
        strcmp(prev_record.oxfa_fragment, current_record.oxfa_fragment) == 0)) {
      
      // Update state
      if (level1_3_curr > 0) {
        valid_level1_3 = 1;
        max_score1_level1_3 = std::max(max_score1_level1_3, current_record.score1);
        max_score2_level1_3 = std::max(max_score2_level1_3, current_record.score2);
      }
      if (level2_4_curr > 0) {
        valid_level2_4 = 1;
        max_score1_level2_4 = std::max(max_score1_level2_4, current_record.score1);
        max_score2_level2_4 = std::max(max_score2_level2_4, current_record.score2);
      }
      
      first_record = false;
    } else {
      // Check and write the previous record
      if (is_valid_record(valid_level1_3, valid_level2_4)) {
        prev_record.score1 = (max_score1_level1_3 + max_score1_level2_4) / 2.0;
        prev_record.score2 = (max_score2_level1_3 + max_score2_level2_4) / 2.0;
        
        fprintf(fp_out, "%d,%.4lf,%.4lf,%.4lf,%s,%.4lf,%.4lf\n",
                prev_record.index, prev_record.mz, prev_record.rt, prev_record.intensity,
                prev_record.oxfa_fragment, prev_record.score1, prev_record.score2);
      }
      
      // Reset state to process new record
      prev_record = current_record;
      valid_level1_3 = (level1_3_curr > 0) ? 1 : 0;
      valid_level2_4 = (level2_4_curr > 0) ? 1 : 0;
      max_score1_level1_3 = (level1_3_curr > 0) ? current_record.score1 : 0.0;
      max_score2_level1_3 = (level1_3_curr > 0) ? current_record.score2 : 0.0;
      max_score1_level2_4 = (level2_4_curr > 0) ? current_record.score1 : 0.0;
      max_score2_level2_4 = (level2_4_curr > 0) ? current_record.score2 : 0.0;
    }
  }
  
  // Process the last record
  if (is_valid_record(valid_level1_3, valid_level2_4)) {
    prev_record.score1 = (max_score1_level1_3 + max_score1_level2_4) / 2.0;
    prev_record.score2 = (max_score2_level1_3 + max_score2_level2_4) / 2.0;
    
    fprintf(fp_out, "%d,%.4lf,%.4lf,%.4lf,%s,%.4lf,%.4lf\n",
            prev_record.index, prev_record.mz, prev_record.rt, prev_record.intensity,
            prev_record.oxfa_fragment, prev_record.score1, prev_record.score2);
  }
  
  fclose(fp_in);
  fclose(fp_out);
}