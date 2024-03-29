

std::vector<std::string> split_string(std::string s, 
																			char target);

int test_for_float(std::string TagArg);



void populate_data_string(std::vector<std::string>* SDTags, 
													std::vector<float>* SDWeights,  
													std::vector<std::string> TagArgs);

int populate_data_array(std::vector<std::string> SDTags, 
												std::vector<float> SDWeights, 
												std::vector<std::string>* binary_tags, 
												std::vector<std::string>* scalar_tags, 
												std::vector<std::vector<std::vector<int> > >* fingerprint_array, 	
												std::vector<std::vector<std::vector<int> > >* static_fingerprint_array, 
												std::vector<float>* binary_weight_array, 
												std::vector<std::vector<float> >* scalar_array, 
												std::vector<std::vector<float> >* static_scalar_array, 
												std::vector<float>* scalar_weight_array, 
												std::vector<std::string> input_filenames, 
												int* molcount, 
												int dimension, 
												std::string static_input_filename, 
												std::vector<std::vector<float> >* static_spe_coordinates);




int add_spe_coordinates_to_mol( std::vector<std::vector<float> > spe_coords, 
																std::string static_input_filename, 
																std::vector<std::string> filename, 
																std::string output_filename);

int add_spe_coordinates_to_mol( std::vector<std::vector<float> > spe_coords,
																std::vector<std::string> filenames, 
																std::string output_filename);

int add_spe_coordinates_to_mol(std::vector<std::vector<float> > spe_coords,
															 std::vector<std::string> filenames);


