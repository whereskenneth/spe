

// Ken Tussey 
// Univ. of Illinois




int count_bits(int n);



float compute_distance(std::vector<float> a, std::vector<float> b);

float compute_tanimoto(std::vector<int> a, std::vector<int> b);

float compute_tanimoto(float a, float b);

std::vector<std::vector<float> > compute_tanimoto_matrix(std::vector<std::vector<int> > fp_int_array);

std::vector<double> compute_average_tanimoto(std::vector<std::vector<std::vector<int> > > fp_int_array);

std::vector<double> compute_average_tanimoto(std::vector<std::vector<float> > scalar_array);

std::vector<std::vector<double> > compute_statistics(std::vector<std::vector<float> > scalar_array);
