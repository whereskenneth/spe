#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <string>
#include <time.h>
#include <math.h>
#include <iterator>
#include <omp.h>
#include "spe_functions.h"
#include "tanimoto.h"





int modified_spe_embed(std::vector<std::vector<float> >* spe_coords,std::vector<std::vector<std::vector<int> > > fingerprint_array,std::vector<float> binary_weight_array,std::vector<std::vector<float> > scalar_array,std::vector<float> scalar_weight_array, unsigned int num_mols, int dimension, unsigned int total_iterations, float cutoff_length, float initial_learning_rate, float embedding_scale) {

//  srand( time(NULL) );
	srand( 1 );

  double tanimoto_distance;
  double euclidean_distance;
  double tanimoto_partial;
  double tanimoto_sum;
  float scale_factor;
  float learning_rate=initial_learning_rate;
  
  unsigned int mol1_index;
  unsigned int mol2_index;
  unsigned int temp;
  unsigned int d;
  unsigned int ilindex;
  unsigned int bin_desc_num;
  unsigned int scalar_desc_num;

  std::vector<float> spe_temp_coords;

  double sum_of_weights=0.0;

  for (int i=0; i < scalar_weight_array.size(); i++) { 
    sum_of_weights+=scalar_weight_array[i];
  }
  for (int i=0; i < binary_weight_array.size(); i++) {
    sum_of_weights+=binary_weight_array[i]; 
  }
    



  for (int i=0; i<num_mols; i++) {
    for (int j=0; j<dimension; j++) {
      spe_temp_coords.push_back((rand()%1000)/500.0-1.0);
    }
    (*spe_coords).push_back(spe_temp_coords);
    spe_temp_coords.clear();
  }
  // Initialize embedded coordinates

  for (int olindex=0; olindex<total_iterations; olindex++) {

    if ((olindex%(total_iterations/25))==0 && olindex%(total_iterations/10) != 0) { 
      printf(".");
      std::cout.flush();
    }
    if (olindex%(total_iterations/10)==0) {
      printf("%1.0f%%",100.0*(olindex)/(total_iterations));
      std::cout.flush();
    }
          
    mol1_index=rand()%num_mols;


    #pragma omp parallel for private(tanimoto_distance) private(euclidean_distance) private(mol2_index) private(scale_factor) private(d) private (tanimoto_sum) private(bin_desc_num) private(tanimoto_partial) private(scalar_desc_num)
    for (int ilindex=mol1_index+1; ilindex< mol1_index+num_mols; ilindex++) {
       
      mol2_index=ilindex%num_mols;      

      tanimoto_sum=0.0;
      for (bin_desc_num=0; bin_desc_num < fingerprint_array.size(); bin_desc_num++) { 

        tanimoto_partial = compute_tanimoto(fingerprint_array[bin_desc_num][mol1_index],fingerprint_array[bin_desc_num][mol2_index]);
        tanimoto_sum+=binary_weight_array[bin_desc_num]*pow(tanimoto_partial, 2);
      }
      for (scalar_desc_num=0; scalar_desc_num < scalar_array.size(); scalar_desc_num++) {
      
        tanimoto_partial = compute_tanimoto(scalar_array[scalar_desc_num][mol1_index], scalar_array[scalar_desc_num][mol2_index]);
        tanimoto_sum+=scalar_weight_array[scalar_desc_num]*pow(tanimoto_partial,2 );

      }


      tanimoto_distance = sqrt( tanimoto_sum / sum_of_weights );
      //std::cout << tanimoto_distance << std::endl;
      euclidean_distance = compute_distance((*spe_coords)[mol1_index],(*spe_coords)[mol2_index]);

      //std::cout << "Tanimoto: " << tanimoto_distance << "  Euclid: " << euclidean_distance << std::endl;

      if (tanimoto_distance <= cutoff_length) {

        scale_factor=0.5*learning_rate*(tanimoto_distance - euclidean_distance)/(euclidean_distance + 1.00e-10);	

      }

      else if (tanimoto_distance > cutoff_length && euclidean_distance < tanimoto_distance ) {

        scale_factor=0.5*learning_rate*(tanimoto_distance - euclidean_distance)/(euclidean_distance + 1.00e-10);
        //std::cout << scale_factor << std::endl;

      }
      else {
        
        scale_factor=0.0;

      }


      for (d=0; d<dimension; d++) {
        (*spe_coords)[mol1_index][d]+=scale_factor*((*spe_coords)[mol1_index][d]-(*spe_coords)[mol2_index][d]);
        (*spe_coords)[mol2_index][d]+=scale_factor*((*spe_coords)[mol2_index][d]-(*spe_coords)[mol1_index][d]);
      }
    }
  
    learning_rate-=initial_learning_rate/total_iterations;
  }


  // Scale the coordinates
  
	for (int scale_all_counter=0; scale_all_counter < num_mols; scale_all_counter++) {

		for (d=0; d<dimension; d++) {

			(*spe_coords)[scale_all_counter][d]*=embedding_scale;

		}

	}
	

  return 1;
}


int static_spe_embed( std::vector<std::vector<float> >* spe_coords, 
											std::vector<std::vector<float> > static_spe_coords, 
											std::vector<std::vector<std::vector<int> > > fingerprint_array, 
											std::vector<std::vector<std::vector<int> > > static_fingerprint_array,
											std::vector<float> binary_weight_array, 
											std::vector<std::vector<float> > scalar_array, 
											std::vector<std::vector<float> > static_scalar_array,
											std::vector<float> scalar_weight_array, 
											unsigned int num_mols, 
											int dimension, 
											unsigned int total_iterations, 
											float cutoff_length, 
											float initial_learning_rate) {

	srand( 1 );

  double tanimoto_distance;
  double euclidean_distance;
  double tanimoto_partial;
  double tanimoto_sum;
  float learning_rate=initial_learning_rate;
	float scale_factor;
  
  unsigned int mol1_index;
  unsigned int mol2_index;
  unsigned int temp;
  unsigned int d;
  unsigned int ilindex;
  unsigned int bin_desc_num;
  unsigned int scalar_desc_num;
	unsigned int static_num_mols = static_spe_coords.size();

  std::vector<float> spe_temp_coords;

  double sum_of_weights=0.0;

  for (int i=0; i < scalar_weight_array.size(); i++) { 
    sum_of_weights+=scalar_weight_array[i];
  }
  for (int i=0; i < binary_weight_array.size(); i++) {
    sum_of_weights+=binary_weight_array[i]; 
  }
    

	


  for (int i=0; i<num_mols; i++) {
    for (int j=0; j<dimension; j++) {
      spe_temp_coords.push_back((rand()%1000)/500.0-1.0);
    }
    (*spe_coords).push_back(spe_temp_coords);
    spe_temp_coords.clear();
  }



	for (int olindex=0; olindex<total_iterations; olindex++) {

    if ((olindex%(total_iterations/25))==0 && olindex%(total_iterations/10) != 0) { 
      printf(".");
      std::cout.flush();
    }
    if (olindex%(total_iterations/10)==0) {
      printf("%1.0f%%",100.0*(olindex)/(total_iterations));
      std::cout.flush();
    }
          
    mol1_index=rand()%num_mols;


    #pragma omp parallel for private(tanimoto_distance) private(euclidean_distance) private(mol2_index) private(scale_factor) private(d) private (tanimoto_sum) private(bin_desc_num) private(tanimoto_partial) private(scalar_desc_num) 
    for (int ilindex=mol1_index+1; ilindex < mol1_index+num_mols+static_num_mols; ilindex++) {
       

      mol2_index=ilindex%(num_mols+static_num_mols);      


      tanimoto_sum=0.0;
      for (bin_desc_num=0; bin_desc_num < fingerprint_array.size(); bin_desc_num++) { 

				if ( mol2_index >= num_mols ) {

				//	std::cout << "Marker x Mol2 =" << mol2_index << std::endl;
				//	std::cout << mol2_index-num_mols << std::endl;
					//std::cout << static_fingerprint_array[bin_desc_num][mol2_index-num_mols][0] << std::endl;
					tanimoto_partial = compute_tanimoto(fingerprint_array[bin_desc_num][mol1_index],static_fingerprint_array[bin_desc_num][mol2_index-num_mols]);


				}
				else {

					
	        tanimoto_partial = compute_tanimoto(fingerprint_array[bin_desc_num][mol1_index],fingerprint_array[bin_desc_num][mol2_index]);
				}

        tanimoto_sum+=binary_weight_array[bin_desc_num]*pow(tanimoto_partial, 2);

      }



      for (scalar_desc_num=0; scalar_desc_num < scalar_array.size(); scalar_desc_num++) {
      
				if ( mol2_index >= num_mols ) {
		
				//	std::cout << "Marker y" << std::endl;
					tanimoto_partial = compute_tanimoto(scalar_array[scalar_desc_num][mol1_index], static_scalar_array[scalar_desc_num][mol2_index-num_mols]);
				}

				else {
				
        	tanimoto_partial = compute_tanimoto(scalar_array[scalar_desc_num][mol1_index], scalar_array[scalar_desc_num][mol2_index]);
				
				}

        tanimoto_sum+=scalar_weight_array[scalar_desc_num]*pow(tanimoto_partial,2 );

      }



	    tanimoto_distance = sqrt( tanimoto_sum / sum_of_weights );
			//std::cout << mol2_index << " " << num_mols << std::endl;
			if (mol2_index >= num_mols) {
				//std::cout << "Marker A" << std::endl;
	  	  euclidean_distance = compute_distance((*spe_coords)[mol1_index],static_spe_coords[mol2_index-num_mols]);
			
			}
			else {
				//std::cout << "Marker B" << std::endl;
				euclidean_distance = compute_distance((*spe_coords)[mol1_index],(*spe_coords)[mol2_index]);
			}



			scale_factor=0.0;
      if (tanimoto_distance <= cutoff_length) {
        scale_factor=0.5*learning_rate*(tanimoto_distance - euclidean_distance)/(euclidean_distance + 1.00e-10);	
      }

      else if (tanimoto_distance > cutoff_length && euclidean_distance < tanimoto_distance ) {
        scale_factor=0.5*learning_rate*(tanimoto_distance - euclidean_distance)/(euclidean_distance + 1.00e-10);

      }
		

      for (d=0; d<dimension; d++) {

				if ( mol2_index >= num_mols ) {

				//	std::cout << "Marker 1" << std::endl;
        	(*spe_coords)[mol1_index][d]+=scale_factor*((*spe_coords)[mol1_index][d]-static_spe_coords[mol2_index-num_mols][d]);
				//	std::cout << "Marker 2" << std::endl;
				}

				else {


					//std::cout << "Marker 3" << std::endl;
			  	(*spe_coords)[mol1_index][d]+=scale_factor*((*spe_coords)[mol1_index][d]-(*spe_coords)[mol2_index][d]);
				//	std::cout << "Marker 4" << std::endl;
        	(*spe_coords)[mol2_index][d]+=scale_factor*((*spe_coords)[mol2_index][d]-(*spe_coords)[mol1_index][d]);
				//	std::cout << "Marker 5" << std::endl;
					//std::cout << "Spe mol1 " << (*spe_coords)[mol1_index][d] << " Spe mol2 " << (*spe_coords)[mol2_index][d] << std::endl;
					//std::cout << " mol 1 " << mol1_index << " mol 2 " << mol2_index << std::endl;
				}


      }
    }
  
    learning_rate-=initial_learning_rate/total_iterations;
  }


  // Scale the coordinates
  /*
	for (int scale_all_counter=0; scale_all_counter < num_mols; scale_all_counter++) {
		for (d=0; d<dimension; d++) {
			(*spe_coords)[scale_all_counter][d]*=embedding_scale;
		}
	}
	*/

  return 1;
}





















/*

std::vector<std::vector<float> > modified_spe_embed_mols(std::vector<float> data_array, int dimension, int outer_loop_length, float cutoff_length, float inner_loop_ratio, float initial_learning_rate ) {

  srand( time(NULL) );
  unsigned int num_mols=data_array.size();
  unsigned int total_iterations=outer_loop_length*inner_loop_ratio;
  double tanimoto_distance;
  double euclidean_distance;
  float scale_factor;
  float learning_rate=initial_learning_rate;

  std::vector<std::vector<float> > spe_coords;
  std::vector<float> spe_temp_coords;

  unsigned int mol1_index;
  unsigned int mol2_index;
  unsigned int temp;
  unsigned int d;
  unsigned int ilindex;

  for (int i=0; i<num_mols; i++) {
    for (int j=0; j<dimension; j++) {
      spe_temp_coords.push_back((rand()%10000)/10000.0-0.5);
    }
    spe_coords.push_back(spe_temp_coords);
    spe_temp_coords.clear();
  }
  // Initialize embedded coordinates

  for (int olindex=0; olindex<total_iterations; olindex++) {

    if ((olindex%(total_iterations/25))==0 && olindex%(total_iterations/10) != 0) { 
      printf(".");
      std::cout.flush();
    }
    if (olindex%(total_iterations/10)==0) {
      printf("%1.0f%%",100.0*(olindex)/(total_iterations));
      std::cout.flush();
    }
          
    mol1_index=rand()%num_mols;

    //mol2_index=rand()%num_mols;
    #pragma omp parallel for private(tanimoto_distance) private(euclidean_distance) private(mol2_index) private(scale_factor) private(d)
    for (int ilindex=mol1_index+1; ilindex< mol1_index+num_mols; ilindex++) {
       
      mol2_index=ilindex%num_mols;      


      tanimoto_distance = compute_tanimoto(data_array[mol1_index],data_array[mol2_index]);
      euclidean_distance = compute_distance(spe_coords[mol1_index],spe_coords[mol2_index]);

      if (tanimoto_distance <= cutoff_length || (tanimoto_distance > cutoff_length && euclidean_distance < tanimoto_distance )) {

        scale_factor=0.5*learning_rate*(tanimoto_distance - euclidean_distance)/(euclidean_distance + 1.00e-10);	

        for (d=0; d<dimension; d++) {
          spe_coords[mol1_index][d]+=scale_factor*(spe_coords[mol1_index][d]-spe_coords[mol2_index][d]);
          spe_coords[mol2_index][d]+=scale_factor*(spe_coords[mol2_index][d]-spe_coords[mol1_index][d]);
        }
      }
    }
    learning_rate-=initial_learning_rate/total_iterations;
  }
  return spe_coords;
}












std::vector<std::vector<float> > modified_spe_embed_mols(std::vector<std::vector<int> > fp_int_array, int dimension, int outer_loop_length, float cutoff_length, float inner_loop_ratio, float initial_learning_rate ) {

  srand( time(NULL) );
  unsigned int num_mols=fp_int_array.size();
  unsigned int total_iterations=outer_loop_length*inner_loop_ratio;
  double tanimoto_distance;
  double euclidean_distance;
  float scale_factor;
  float learning_rate=initial_learning_rate;

  std::vector<std::vector<float> > spe_coords;
  std::vector<float> spe_temp_coords;

  unsigned int mol1_index;
  unsigned int mol2_index;
  unsigned int temp;
  unsigned int d;
  unsigned int ilindex;

  for (int i=0; i<num_mols; i++) {
    for (int j=0; j<dimension; j++) {
      spe_temp_coords.push_back((rand()%10000)/10000.0-0.5);
    }
    spe_coords.push_back(spe_temp_coords);
    spe_temp_coords.clear();
  }
  // Initialize embedded coordinates

  for (int olindex=0; olindex<total_iterations; olindex++) {

    if ((olindex%(total_iterations/25))==0 && olindex%(total_iterations/10) != 0) { 
      printf(".");
      std::cout.flush();
    }
    if (olindex%(total_iterations/10)==0) {
      printf("%1.0f%%",100.0*(olindex)/(total_iterations));
      std::cout.flush();
    }
          
    mol1_index=rand()%num_mols;

    //mol2_index=rand()%num_mols;
    #pragma omp parallel for private(tanimoto_distance) private(euclidean_distance) private(mol2_index) private(scale_factor) private(d)
    for (int ilindex=mol1_index+1; ilindex< mol1_index+num_mols; ilindex++) {
       
      mol2_index=ilindex%num_mols;      


      tanimoto_distance = compute_tanimoto(fp_int_array[mol1_index],fp_int_array[mol2_index]);
      euclidean_distance = compute_distance(spe_coords[mol1_index],spe_coords[mol2_index]);

      if (tanimoto_distance <= cutoff_length || (tanimoto_distance > cutoff_length && euclidean_distance < tanimoto_distance )) {

        scale_factor=0.5*learning_rate*(tanimoto_distance - euclidean_distance)/(euclidean_distance + 1.00e-10);	

        for (d=0; d<dimension; d++) {
          spe_coords[mol1_index][d]+=scale_factor*(spe_coords[mol1_index][d]-spe_coords[mol2_index][d]);
          spe_coords[mol2_index][d]+=scale_factor*(spe_coords[mol2_index][d]-spe_coords[mol1_index][d]);
        }
      }
    }
    learning_rate-=initial_learning_rate/total_iterations;
  }
  return spe_coords;
}






//############################################################################################################



// Deprecated
std::vector<std::vector<float> > spe_embed_mols(std::vector<std::vector<int> > fp_int_array, int dimension, int outer_loop_length, float cutoff_length, float inner_loop_ratio, float initial_learning_rate ) {

  srand( time(NULL) );
  unsigned int num_mols=fp_int_array.size();
  unsigned int inner_loop_length=inner_loop_ratio*num_mols;
  double tanimoto_distance;
  double euclidean_distance;
  float scale_factor;
  float learning_rate=initial_learning_rate;

  std::vector<std::vector<float> > spe_coords;
  std::vector<float> spe_temp_coords;

  unsigned int mol1_index;
  unsigned int mol2_index;
  unsigned int temp;
  unsigned int d;
  unsigned int ilindex;

  for (int i=0; i<num_mols; i++) {
    for (int j=0; j<dimension; j++) {
      spe_temp_coords.push_back((rand()%10000)/10000.0);
    }
    spe_coords.push_back(spe_temp_coords);
    spe_temp_coords.clear();
  }
  // Initialize embedded coordinates

  for (int olindex=0; olindex<outer_loop_length; olindex++) {

    if ((olindex%(outer_loop_length/25))==0 && olindex%(outer_loop_length/10) != 0) { 
      printf(".");
      std::cout.flush();
    }
    if (olindex%(outer_loop_length/10)==0) {
      printf("%1.0f%%",100.0*(olindex)/(outer_loop_length));
      std::cout.flush();
    }
    #pragma omp parallel for private(tanimoto_distance) private(euclidean_distance) private(mol1_index) private(mol2_index)  private(temp) private(scale_factor) private(ilindex) private(d)
    for (int ilindex=0; ilindex<inner_loop_length; ilindex++) {
      mol1_index=rand()%num_mols;
      mol2_index=rand()%num_mols;

      while (mol1_index == mol2_index) { 
        mol2_index=rand()%num_mols;
      }
      if (mol2_index < mol1_index) {
        temp=mol2_index;
        mol2_index=mol1_index;  
        mol1_index=temp;
      }              
      tanimoto_distance = compute_tanimoto(fp_int_array[mol1_index],fp_int_array[mol2_index]);
      euclidean_distance = compute_distance(spe_coords[mol1_index],spe_coords[mol2_index]);

      if (tanimoto_distance > cutoff_length  && euclidean_distance >= tanimoto_distance ) {;}
      else {
        scale_factor=0.5*learning_rate*(tanimoto_distance - euclidean_distance)/(euclidean_distance + 1.00e-100);
	if (scale_factor > 1000) { 
          std::cout << "Embedding is exploding!" << std::endl;
          std::cout << "Scale Factor: " << scale_factor << std::endl;  
          std::cout << "tanimoto: " << tanimoto_distance << std::endl;
          std::cout << "euclidean: " << euclidean_distance << std::endl;
          scale_factor = 0.0;	
        }
        for (d=0; d<dimension; d++) {
          spe_coords[mol1_index][d]+=scale_factor*(spe_coords[mol1_index][d]-spe_coords[mol2_index][d]);
          spe_coords[mol2_index][d]+=scale_factor*(spe_coords[mol2_index][d]-spe_coords[mol1_index][d]);
        }
      }
    }
    learning_rate-=initial_learning_rate/outer_loop_length;
  }
  return spe_coords;
}




//################################################################################################################




// Deprecated
std::vector<std::vector<float> > spe_embed_mols(std::vector<std::vector<float> > tanimoto_matrix, int dimension=2) {

  srand( time(NULL) );
  int num_mols=tanimoto_matrix.size();
  int outer_loop_length=10000;
  int inner_loop_length=0.4*num_mols;
  float cutoff_length=0.4;
  float initial_learning_rate=1.0;
  float tanimoto_distance;
  float euclidean_distance;
  float scale_factor;
  float learning_rate=initial_learning_rate;

  std::vector<std::vector<float> > spe_coords;
  std::vector<float> spe_temp_coords;

  int mol1_index;
  int mol2_index;




  for (int i=0; i<num_mols; i++) {
    for (int j=0; j<dimension; j++) {
      spe_temp_coords.push_back((rand()%100)/100.0);
    }
    spe_coords.push_back(spe_temp_coords);
    spe_temp_coords.clear();
  }
  // Initialize embedded coordinates

  for (int olindex=0; olindex<outer_loop_length; olindex++) {

    #pragma omp parallel for  
    for (int ilindex=0; ilindex<inner_loop_length; ilindex++) {

      mol1_index=rand()%num_mols;
      mol2_index=rand()%num_mols;

      while (mol1_index == mol2_index) { 
        mol2_index=rand()%num_mols;
      }

      if (mol2_index < mol1_index) {
        int temp=mol2_index;
        mol2_index=mol1_index;  
        mol1_index=temp;
      }        

      
      tanimoto_distance = tanimoto_matrix[mol1_index][mol2_index];
      euclidean_distance = compute_distance(spe_coords[mol1_index],spe_coords[mol2_index]);


 
      if (tanimoto_distance > cutoff_length && euclidean_distance >= tanimoto_distance) {;}
      else {
        
        scale_factor=0.5*learning_rate*(tanimoto_distance - euclidean_distance)/(euclidean_distance + 2.2250738585072014e-308);


        for (int d=0; d<dimension; d++) {


          spe_coords[mol1_index][d]+=scale_factor*(spe_coords[mol1_index][d]-spe_coords[mol2_index][d]);
          spe_coords[mol2_index][d]+=scale_factor*(spe_coords[mol2_index][d]-spe_coords[mol1_index][d]);

        }

      }

    
    }


    learning_rate-=initial_learning_rate/outer_loop_length;

  }

  return spe_coords;

}
*/
