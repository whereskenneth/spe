








int modified_spe_embed( std::vector<std::vector<float> >* spe_coordinates,
											  std::vector<std::vector<std::vector<int> > > fingerprint_array,
												std::vector<float> binary_weight_array,
												std::vector<std::vector<float> > scalar_array,
												std::vector<float> scalar_weight_array, 
												unsigned int num_mols, 
												int dimension, 
												unsigned int total_iterations, 
												float cutoff_length, 
												float initial_learning_rate, 
												float embedding_scale );


int static_spe_embed( std::vector<std::vector<float> >* spe_coordinates, 
										  std::vector<std::vector<float> > static_spe_coordinates,
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
											float initial_learning_rate);

//


//std::vector<std::vector<float> > spe_embed_mols(std::vector<std::vector<float> > tanimoto_matrix, int dimension); // Not used

//std::vector<std::vector<float> > spe_embed_mols(std::vector<std::vector<int> > fp_int_array,int dimension,int outer_loop_length,float cutoff_length,float inner_loop_ratio,float initial_learning_rate );



//std::vector<std::vector<float> > modified_spe_embed_mols(std::vector<std::vector<int> > fp_int_array, int dimension, int outer_loop_length, float cutoff_length, float inner_loop_ratio, float initial_learning_rate );

//std::vector<std::vector<float> > modified_spe_embed_mols(std::vector<float> data_array, int dimension, int outer_loop_length, float cutoff_length, float inner_loop_ratio, float initial_learning_rate );
