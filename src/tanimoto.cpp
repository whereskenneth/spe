#include <math.h>
#include <vector>
#include "tanimoto.h"
#include <iostream>


//################################################

int count_bits(int n) {
  unsigned int c;
  for (c=0;n;c++)
    n &= n- 1;
  return c;
}

//###############################################

//###############################################

float compute_distance(std::vector<float> a, std::vector<float> b) {
  float square_sum=0.0;
  for (int i=0; i< a.size(); i++) {
    square_sum+=pow(a[i]-b[i],2);
  }
  return sqrt(square_sum);

}

//###############################################

//###############################################

std::vector<std::vector<float> > compute_tanimoto_matrix(std::vector<std::vector<int> > fp_int_array) {

  std::vector<std::vector<float> > tani_matrix;
  std::vector<float> tani_temp_array;
  int a,b,diffs,ones,diff_count,ones_count;
  float tani_sum;


  for (int i=0; i < fp_int_array.size(); i++) {
    for (int j=i+1; j < fp_int_array.size(); j++) {
      tani_sum=0.0;      
      diff_count=0;
      ones_count=0;
      for (int k=0; k < fp_int_array[i].size(); k++) {
        a=fp_int_array[i][k];
        b=fp_int_array[j][k];
        ones=0;
        diffs=0;
        diffs=a^b;
        ones=a&b;
        diff_count+=count_bits(diffs);
        ones_count+=count_bits(ones);
      }
      tani_sum+=(diff_count)/(diff_count+ones_count+ 2.2250738585072014e-308);    
      tani_temp_array.push_back(tani_sum);
    }  
  tani_matrix.push_back(tani_temp_array);
  tani_temp_array.clear();
  }
  
  return tani_matrix;
}

//##################################################################

//##################################################################

float compute_tanimoto(std::vector<int> a, std::vector<int> b) {
  //computes tanimoto distance in bitspace   
  
    float tani_sum=0.0;
    int ones=0;
    int diffs=0;
    int diff_count=0;
    int ones_count=0;  
    int a_int;
    int b_int;

    for (int k=0; k < a.size(); k++) {
      a_int=a[k];
      b_int=b[k];
      ones=0;
      diffs=0;
      diffs=a_int^b_int;
      ones=a_int&b_int;
      diff_count+=count_bits(diffs);
      ones_count+=count_bits(ones);
    }
    
    tani_sum=(diff_count)/(diff_count+ones_count+2.2250738585072014e-308);


  return tani_sum;

}


float compute_tanimoto(float a, float b) {
  // Computes tanimoto distance in continuous space

  return 1-a*b/(pow(a,2) + pow(b,2) - a*b + 2.2250738585072014e-308);

}




std::vector<double> compute_average_tanimoto(std::vector<std::vector<std::vector<int> > > fp_int_array) { 

  double tanimoto_sum;
  std::vector<double> average_tanimoto;

  for (int i =0; i < fp_int_array.size(); i++) {

    tanimoto_sum=0.0;
    for (int j=0; j < fp_int_array[i].size(); j++ ) {
      
      for (int k=j+1; k < fp_int_array[i].size(); k++) {
        tanimoto_sum+=1-compute_tanimoto( fp_int_array[i][j],fp_int_array[i][k] );
      }
    }
    average_tanimoto.push_back(tanimoto_sum/(0.5*pow(fp_int_array[i].size(),2) - fp_int_array[i].size()) );
  }
  
  return average_tanimoto;

}

std::vector<double> compute_average_tanimoto(std::vector<std::vector<float> > scalar_array) {

  double tanimoto_sum;
  std::vector<double> average_tanimoto;


  for (int i=0; i < scalar_array.size(); i++) {

    tanimoto_sum=0.0;
    for (int j=0; j < scalar_array[i].size(); j++ ) {

      for (int k=j+1; k < scalar_array[i].size(); k++ ) {
        tanimoto_sum+=1-compute_tanimoto(scalar_array[i][j], scalar_array[i][k]);           
      }

    }
    //std::cout << tanimoto_sum << std::endl;
    average_tanimoto.push_back(tanimoto_sum/(0.5*pow(scalar_array[i].size(),2)-scalar_array[i].size()));        
  }

  return average_tanimoto;

}


std::vector<std::vector<double> > compute_statistics(std::vector<std::vector<float> > scalar_array) { 


  std::vector<std::vector<double> > stats_array;
  //0 will be average
  //1 will be variance
  //2 will be Standard deviation
  std::vector<double> tmp_stats;

  double sum; 
  double sum_squared;
  double mean;
  double variance;

  for (int i=0; i < scalar_array.size(); i++) {

    sum=0.0;
    sum_squared=0.0;
    for (int j=0; j < scalar_array[i].size(); j++ ) {

      sum+=scalar_array[i][j];
      sum_squared+=pow(scalar_array[i][j],2);

    }
    mean=sum/scalar_array[i].size();
    variance=(sum_squared - sum*mean)/(scalar_array[i].size() - 1 );
  
    tmp_stats.push_back(mean);
    tmp_stats.push_back(variance);
    tmp_stats.push_back(sqrt(variance));

    stats_array.push_back(tmp_stats);
    tmp_stats.clear();

  }

  return stats_array;


}
  




