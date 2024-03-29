#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <stdio.h>
#include <time.h>

#include <boost/program_options.hpp>
namespace po=boost::program_options;

#include <openbabel/mol.h>

#include "read_write_files.h"
#include "spe_functions.h"
#include "tanimoto.h"

int main(int argc, char* argv[]){


  //srand( 1 );

  std::vector<std::string> input_filenames;
  std::string static_input_filename;
  std::string output_filename;
  std::vector<std::string> TagArgs;
  std::vector<std::string> SDTags;
  std::vector<float> SDWeights;
  
  
  int molcount=0;
  bool use_input_filename;
  float embedding_scale=1.0;


  po::options_description desc("Options:"); 
  desc.add_options()
    ("help,h","Produces this message")
    ("input,i",po::value< std::vector<std::string> >()->multitoken(), "input file")
    ("static,s",po::value< std::string >(), "input file with static SPE coordinates")
    ("output,o",po::value< std::string >(), "output file")
    ("dim,d",po::value<int>()->default_value(2),"Number of Dimensions to embed data into")
    ("cutoff,r",po::value<float>()->default_value(0.4),"Cutoff Radius for embedding - Use -1 for infinite cutoff")
    ("tag,t",po::value< std::vector<std::string> >()->multitoken(), "SD tags that correspond to data used for the embedding algorithm")
	  ("iterations,n",po::value<int>()->default_value(5000),"Number of total iterations of algorithm")
    ("initial-learning-rate",po::value<float>()->default_value(1.0),"Initial learning rate for SPE algorithm")
    ("embedding-scale",po::value<float>()->default_value(1.0),"Used to scale the embedded coordinates up or down")
    ("stats,T","Run statistics on descriptors");

  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc), vm);
  po::notify(vm);

  int dimension;
  float cutoff;
  unsigned int total_iterations;
  float initial_learning_rate;
  bool run_spe=true;
  bool run_stats=false;
  bool run_static=false;

  if (vm.count("help")) { 
      std::cout << desc << "\n";
			printf("Usage examples:       ./spe -h\n");
      printf("                         ^\n");
		  printf("                         ----------- Prints this message\n\n");
			printf("                     ./spe -i jarp.sdf -o test_spe.sdf\n");
      printf("                        ^\n");
			printf("                        ----------- Takes the example jarp.sdf file, embeds\n");
			printf("                                    the CF's into 2 dimensions\n\n");
			printf("                     ./spe -i jarp.sdf -d 3\n");
      printf("                        ^\n");
			printf("                        ------------Embeds jarp.sdf into 3 dimensions using\n");
			printf("                                    fingerprints, output jarp_spe-output.sdf\n\n");
			printf("                     ./spe -i jarp1.sdf jarp2.sdf -o test_spe.sdf\n");
      printf("                        ^\n");
			printf("                        -----------Embeds both files jarp1 and jarp2 into 2d\n");
			printf("                                    space, and combines their output.\n\n");
      return 1;
  }

  if(vm.count("input")) {     
    input_filenames=vm["input"].as< std::vector<std::string> >();
  }
  else {
    std::cout <<desc << std::endl;      
    std::cout << "No input files" << std::endl;
    return 0;
  }

  if (vm.count("static")) {
    static_input_filename=vm["static"].as<std::string>();
    run_static=true;
    run_spe=false;
  }
    


  if (vm.count("output")) { 
    output_filename=vm["output"].as<std::string>();
  }
  else {

    use_input_filename=true;    
  }
  
  if (vm.count("dim")) {
    dimension=vm["dim"].as<int>();
    
    if (dimension > 3) {
      dimension=3;
      std::cout << "Dimensions must be 2, or 3. Setting dimension to 3" << std::endl;
    }
    if (dimension < 2) {
      dimension=2;
      std::cout << "Dimensions must be 2, or 3. Setting dimension to 2" << std::endl;
    }
  }

  if (vm.count("tag")) {
    TagArgs=vm["tag"].as<std::vector<std::string> >();

    //std::cout << TagArgs << std::endl;
  }
  else {

    TagArgs.push_back("CF");
    std::cout << std::endl << "Using default CF tags for embedding..." << std::endl;

  }

  if (vm.count("stats")) {
    
    run_spe=false;
    run_stats=true;

  }
  
  if (vm.count("embedding-scale")) { 

    if (vm["embedding-scale"].as<float>() <= 0.0) {
			
			embedding_scale=1.0;
		}
		
		else {
			embedding_scale=vm["embedding-scale"].as<float>();
		}

	}

    


  cutoff=vm["cutoff"].as<float>();
  total_iterations=vm["iterations"].as<int>();
  initial_learning_rate=vm["initial-learning-rate"].as<float>();


  std::vector<std::vector<std::vector<int> > > fingerprint_array;
	std::vector<std::vector<std::vector<int> > > static_fingerprint_array;
  std::vector<std::vector<float> > scalar_array;  
	std::vector<std::vector<float> > static_scalar_array;
  std::vector<float> binary_weight_array;
  std::vector<float> scalar_weight_array;
  std::vector<std::vector<float> > spe_coordinates;
  std::vector<std::vector<float> > static_spe_coordinates;
  std::vector<std::string> binary_tags;
  std::vector<std::string> scalar_tags;

  

  populate_data_string(&SDTags,&SDWeights,TagArgs);   

  std::cout << "Populating data arrays...   ";
  std::cout.flush();
  clock_t tick = clock();



  

  if (populate_data_array(SDTags,SDWeights, &binary_tags, &scalar_tags, &fingerprint_array, &static_fingerprint_array,&binary_weight_array,&scalar_array,&static_scalar_array,&scalar_weight_array,input_filenames, &molcount, dimension, static_input_filename, &static_spe_coordinates)) {
    std::cout << "Finished " << molcount << " mols in ";
    std::cout.flush();
    printf("%1.3f seconds.\n", (clock()-tick)/1000000.);

  }
  else {
    std::cout << std::endl << "Files are missing data - Exiting" << std::endl; 
    std::cout << "Found ..." << std::endl;
    for (int i=0; i < SDTags.size(); i++) {
      std::cout << SDTags[i] << std::endl;
    }
    if (run_spe) { std::cout << "Can't embed with missing data" << std::endl; return 0;}
 
  }
    
  if (run_spe) {
    std::cout << "Embedding into " << dimension << " dimensional space." << std::endl;
    tick=clock();
    
		//std::cout << total_iterations << ": Total iterations" << std::endl;

    if ( modified_spe_embed(&spe_coordinates,fingerprint_array,binary_weight_array,scalar_array,scalar_weight_array, molcount, dimension, total_iterations, cutoff, initial_learning_rate, embedding_scale) ) { 

      std::cout << "  Finshed!" << std::endl;
      std::cout << "Finished in " << (clock()-tick)/60000000. << " minutes. " << std::endl;
     
      if (use_input_filename) {
        add_spe_coordinates_to_mol(spe_coordinates ,input_filenames);  
      }
      else {
        add_spe_coordinates_to_mol(spe_coordinates ,input_filenames,output_filename);
      }      
    }
    else {
      std::cout << "Error in embedding, exiting" << std::endl;
    }
  }


  if (run_static) {

    std::cout << "Embedding libraries onto " << static_input_filename << std::endl;

    tick=clock();

		std::cout << "Dimension of static_fingerprint_array: " << static_fingerprint_array.size() << " " << static_fingerprint_array[0].size() << " " << static_fingerprint_array[0][0].size() << std::endl;

		std::cout << "Dimension of fingerprint_array: " << fingerprint_array.size() << " " << fingerprint_array[0].size() << " " << fingerprint_array[0][0].size() << std::endl;

    if (static_spe_embed(&spe_coordinates, static_spe_coordinates, fingerprint_array, static_fingerprint_array, binary_weight_array, scalar_array, static_scalar_array, scalar_weight_array, molcount, dimension, total_iterations, cutoff, initial_learning_rate));

  //  for (int i =0; i < static_spe_coordinates.size(); i++ ) {

  //    std::cout << static_spe_coordinates[i][0] << " " << static_spe_coordinates[i][1] << std::endl;
	//		std::cout << "Number of static mols: " << static_spe_coordinates.size() << std::endl;
	//		std::cout << "Number of dynamics mols: " << spe_coordinates.size() << std::endl
  //  }

		std::cout << "Finished in " << (clock() - tick)/60000000. << " minutes. " << std::endl;


		//std::cout << "Writing data to file " << output_filename << std::endl;

		add_spe_coordinates_to_mol(spe_coordinates,static_input_filename,input_filenames,output_filename);

		std::cout << "Finished!" << std::endl;


  }

  









  if (run_stats) { 

    std::vector<double> average_tanimoto_binary;
    std::vector<double> average_tanimoto_scalar;
    std::vector<std::vector<double> > scalar_statistics;        

    average_tanimoto_binary=compute_average_tanimoto(fingerprint_array);
    average_tanimoto_scalar=compute_average_tanimoto(scalar_array);
    
    scalar_statistics = compute_statistics(scalar_array);

    for (int i=0; i < binary_tags.size(); i++ ) {

      std::cout << "Average tanimoto in " << binary_tags[i] << ": " << average_tanimoto_binary[i] << std::endl;

    }
    std::cout<<std::endl;

    for (int i=0; i < scalar_tags.size(); i++ ) {

      std::cout << "Average tanimoto in " << scalar_tags[i] << ": " << average_tanimoto_scalar[i] << std::endl;
      std::cout << "Mean value of " << scalar_tags[i] << ": " << scalar_statistics[i][0] << std::endl;
      std::cout << "Variance in " << scalar_tags[i] << ": " << scalar_statistics[i][1] << std::endl;
      std::cout << "Std. Dev. in " << scalar_tags[i] << ": " << scalar_statistics[i][2] << std::endl;
      std::cout << std::endl;
    }
  }






/*
  int randmol=rand()%molcount;
  for (int i=0; i< scalar_array.size(); i++) {
  
    std::cout << " Data tag: " << scalar_tags[i];
    std::cout << " Number of mols: " << scalar_array[i].size(); 
    std::cout << " Weight: " << scalar_weight_array[i] << std::endl;  
    std::cout << " Example: " << scalar_array[i][randmol] << std::endl;

  }

  std::cout << "\n\n";
*/
  





/*
  std::cout << "Populating mol array...   ";
  std::cout.flush();
  molcount=get_mols_in_array(input_filenames,&molarray);
  std::cout << "Finished!" << std::endl;


    get_data_in_array(input_filenames,&scalar_data_array, &fp_int_array, sd_tags_for_embedding);

    //spe_coordinates=modified_spe_embed_mols( fp_int_array, dimension, outer_loop_number, cutoff, inner_loop_ratio, initial_learning_rate);



   // spe_coordinates=modified_spe_embed_mols( scalar_data_array, dimension, outer_loop_number, cutoff, inner_loop_ratio, initial_learning_rate);



  //tanimoto_matrix=compute_tanimoto_matrix( fp_int_array );

  //spe_coordinates=spe_embed_mols( tanimoto_matrix, dimension);


*/



  return 0;


}
