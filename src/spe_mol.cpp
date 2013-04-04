#include <iostream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <string>
#include <time.h>
#include <math.h>
#include <iterator>
#include <omp.h>

#include <boost/program_options.hpp>
namespace po=boost::program_options;

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/descriptor.h>
#include <openbabel/generic.h>

int count_bits(int n) {
  unsigned int c;
  for (c=0;n;c++)
    n &= n- 1;
  return c;
}


std::vector<std::string> split(std::string s, char target)
{
    // Post: This function is the same as the split() function in python. 
    //       A string is split at every occurence of target & each section 
    //       is stored in a vector element.
    
    std::vector <std::string> result;
    std::string tempStr;
    
    //if (!islower(target)) {
    target = tolower(target);
    //}
    
    for (int i=0; i<s.length(); i++) 
    {
        s[i] = tolower(s[i]);
        if (s[i] == target) {
           std::string res = tempStr;
           result.push_back(res);
           tempStr = "";
        }
        else tempStr += s[i];
    }
    
    if (tempStr.length() > 0) {
         std::string res = tempStr;
         result.push_back(res);
    }
    
    return result; 
}


//##################################################################



std::vector<std::vector<int> > get_fps_in_array(std::vector<std::string> filenames) {

//            Gets chemical fingerprints from file, and stores them in an array of molecules each with
//            an array of integers. Vector types are used to avoid hardcoding of fingerprint length into
//            code


  OpenBabel::OBConversion obconversion;
  obconversion.SetInFormat("sdf"); // FIXME Change to input format agnostic
  OpenBabel::OBMol mol;
  int counter=0;
  int fp_integer;
  std::vector<std::vector<int> > fp_int_array;

  std::vector<int> fp_temp_array;
  std::vector<std::string> fp_string_array;    

  for (int fileno=0; fileno<filenames.size();fileno++) {

    bool notatend = obconversion.ReadFile(&mol, filenames[fileno]);

    while (notatend) {

      std::vector<OpenBabel::OBGenericData*>::iterator k;
      std::vector<OpenBabel::OBGenericData*> vdata = mol.GetData();
      std::vector<std::string>::iterator j;
      
      
      for (k=vdata.begin();k!=vdata.end();++k) {
        if ((*k)->GetAttribute() == "CF" ) {
          std::istringstream fp_string (((OpenBabel::OBPairData*)(*k))->GetValue(),std::istringstream::in);
          
          fp_string_array=split(fp_string.str(),'\t');

          for (j=fp_string_array.begin();j!=fp_string_array.end();++j) {
            fp_temp_array.push_back(atoi((*j).c_str()));       
              //fp_int_array.push_back(atoi((*j).c_str()));    
          }
        }
      }

      fp_int_array.push_back(fp_temp_array);
      fp_temp_array.clear();
      counter++;
      notatend=obconversion.Read(&mol);
    }
  }

  return fp_int_array;

}


float compute_tanimoto(std::vector<int> a, std::vector<int> b) {
   
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





//#############################################################################


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


float compute_distance(std::vector<float> a, std::vector<float> b) {
  float square_sum=0.0;
  for (int i=0; i< a.size(); i++) {
    square_sum+=pow(a[i]-b[i],2);
  }
  return sqrt(square_sum);

}



std::vector<std::vector<float> > spe_embed_mols(std::vector<std::vector<int> > fp_int_array, int dimension, int outer_loop_length, float cutoff_length, float inner_loop_ratio, float initial_learning_rate ) {

  srand( time(NULL) );
  int num_mols=fp_int_array.size();
  int inner_loop_length=inner_loop_ratio*num_mols;
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


  std::cout << std::endl;
  for (int olindex=0; olindex<outer_loop_length; olindex++) {

    if ((olindex%(outer_loop_length/25))==0 && olindex%(outer_loop_length/10) != 0) { 
      printf(".");
      std::cout.flush();
    }
    if (olindex%(outer_loop_length/10)==0) {
      printf("%1.0f%%",100.0*(olindex)/(outer_loop_length));
      std::cout.flush();
    }
  
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

      
      tanimoto_distance = compute_tanimoto(fp_int_array[mol1_index],fp_int_array[mol2_index]);
      euclidean_distance = compute_distance(spe_coords[mol1_index],spe_coords[mol2_index]);


 
      if (tanimoto_distance > cutoff_length  && euclidean_distance >= tanimoto_distance) {;}
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


int add_spe_coordinates_to_mol(std::vector<std::vector<float> > spe_coords ,std::vector<std::string> filenames, std::string output_filename) {

  int molcounter;

  OpenBabel::OBConversion obconversion;
  OpenBabel::OBMol mol;

  obconversion.SetOutFormat("SDF");
  std::ofstream ofs(output_filename.c_str());
  obconversion.SetOutStream(&ofs);
  
  for (int fileno=0; fileno<filenames.size(); fileno++) {

    OpenBabel::OBFormat* inFormat = obconversion.FormatFromExt(filenames[fileno]);

    obconversion.SetInFormat(inFormat);

    bool notatend = obconversion.ReadFile(&mol, filenames[fileno]);    
    molcounter=0;
    while (notatend) {

      OpenBabel::OBPairData *SpeX = new OpenBabel::OBPairData;
      SpeX->SetAttribute("SPE.x");
      std::ostringstream spex;
  
      OpenBabel::OBPairData *SpeY = new OpenBabel::OBPairData;
      SpeY->SetAttribute("SPE.y"); 
      std::ostringstream spey;

      OpenBabel::OBPairData *SpeZ = new OpenBabel::OBPairData;

      std::ostringstream spez; 
      spex << spe_coords[molcounter][0];
      SpeX->SetValue(spex.str());
      spey << spe_coords[molcounter][1];
      SpeY->SetValue(spey.str());
      if (spe_coords[0].size() == 3) {
        SpeZ->SetAttribute("SPE.z");
        spez << spe_coords[molcounter][2];
        SpeZ->SetValue(spez.str());
      }
  
      mol.SetData(SpeX);
      mol.SetData(SpeY);
      mol.SetData(SpeZ);

      std::vector<OpenBabel::OBGenericData*>::iterator k;
      std::vector<OpenBabel::OBGenericData*> vdata = mol.GetData();
      std::vector<std::string>::iterator j;
  
      bool has_clid=false;

      for (k=vdata.begin();k!=vdata.end();++k) {
        if ((*k)->GetAttribute() == "clid" ) {has_clid=true;}
     
      }
      
      if (!has_clid) {
      
        OpenBabel::OBPairData *clid = new OpenBabel::OBPairData;
      
        clid->SetAttribute("clid");
        clid->SetValue("-1");
        mol.SetData(clid);

      }




      obconversion.Write(&mol);
      notatend=obconversion.Read(&mol);
  
      molcounter++;
    }
  
  }

  return 0;

}



int main(int argc, char* argv[]){

  std::vector<std::string> input_filenames;
  std::string output_filename;
  std::vector<std::string> sd_tags_for_embedding;
     

  po::options_description desc("Options:"); 
  desc.add_options()
    ("help,h","Produces this message")
    ("input,i",po::value< std::vector<std::string> >()->multitoken(), "input file")
    ("output,o",po::value< std::string >(), "output file")
    ("dim,d",po::value<int>()->default_value(2),"Number of Dimensions to embed data into")
    ("cutoff,r",po::value<float>()->default_value(0.4),"Cutoff Radius for embedding - Use -1 for infinite cutoff")
    ("tag,t",po::value< std::vector<std::string> >()->multitoken(), "SD tags that correspond to data used for the embedding algorithm")
    ("inner-loop-ratio",po::value<float>()->default_value(0.2),"Scaling factor for number of inner loops multiplied by the number of molecules to run")
    ("outer-loop",po::value<int>()->default_value(10000),"Number of outer loops to run")
    ("initial-learning-rate",po::value<float>()->default_value(1.0),"Initial learning rate for SPE algorithm")

  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc), vm);
  po::notify(vm);

  int dimension;
  int outer_loop_number;
  float cutoff;
  float inner_loop_ratio;
  float initial_learning_rate;



  

  if (vm.count("help")) { 
      std::cout << desc << "\n";
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


  if (vm.count("output")) { 
    output_filename=vm["output"].as<std::string>();
  }
  else {
    std::cout << desc << std::endl;
    std::cout << "No output files" << std::endl;
    return 0;
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
    sd_tags_for_embedding=vm["tag"].as<std::vector<std::string> >();
  }
  else {
    std::string cf_string="CF";
    sd_tags_for_embedding.push_back(cf_string);
    std::cout << "Using default 'CF' as tag for embedding" << std::endl;
  }


  outer_loop_number=vm["outer-loop"].as<int>();
  cutoff=vm["cutoff"].as<float>();
  inner_loop_ratio=vm["inner-loop-ratio"].as<float>();
  initial_learning_rate=vm["initial-learning-rate"].as<float>();


  std::vector<std::vector<int> > fp_int_array;
  std::vector<std::vector<float> > tanimoto_matrix;
  std::vector<std::vector<float> > spe_coordinates;


  std::cout << "Populating fingerprint integer array...    " << std::endl;
  fp_int_array=get_fps_in_array(input_filenames);    
  std::cout << "Finished!" << std::endl;


  //tanimoto_matrix=compute_tanimoto_matrix( fp_int_array );

  //spe_coordinates=spe_embed_mols( tanimoto_matrix, dimension);


  std::cout << "Embedding molecules into "<<dimension;
  std::cout << " dimensional space..." << std::endl;

  spe_coordinates=spe_embed_mols( fp_int_array, dimension, outer_loop_number, cutoff, inner_loop_ratio, initial_learning_rate);

  std::cout << "Finished!" << std::endl;

  std::cout << "Writing coordinates to " << output_filename;
  std::cout << "..." << std::endl;

  add_spe_coordinates_to_mol(spe_coordinates,input_filenames,output_filename);


  std::cout << "Finished!\n\n" << std::endl;

  return 1;


}
