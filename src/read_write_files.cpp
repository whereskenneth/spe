#include <iostream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <string>

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/descriptor.h>
#include <openbabel/generic.h>
#include "read_write_files.h"





std::vector<std::string> split_string(std::string s, 
																			char target)
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



int test_for_float(std::string TagArg) {

  std::istringstream argstream(TagArg);

  float test_holder;

  if (argstream >> test_holder) { return 1;}
  else { return 0;}

}




void populate_data_string(std::vector<std::string>* SDTags, 
													std::vector<float>* SDWeights, 
													std::vector<std::string> TagArgs)  {

  int argno=0;

  while (argno < TagArgs.size()) {

    if (test_for_float(TagArgs[argno])) {
      //std::cout << "Skipping "<< TagArgs[argno] << std::endl;
      argno++;
      continue;
    }
    else {
      //std::cout << "Found descriptor " << TagArgs[argno] << std::endl; 
      (*SDTags).push_back(TagArgs[argno]);      
      //std::cout << "Found weight ";
      //std::cout.flush(); 
      if (argno+1 < TagArgs.size() ) {

        if (test_for_float(TagArgs[argno+1]) ) {
          (*SDWeights).push_back( atof( TagArgs[argno+1].c_str() ) );
          //std::cout << TagArgs[argno+1] << std::endl;
          argno++;
        }
        else{
         //std::cout << 1.0 << std::endl;
         (*SDWeights).push_back(1.0);
        }
      }
      else {
        //std::cout << 1.0 << std::endl;
        (*SDWeights).push_back(1.0);
      }

      
    }
    argno++;

  }        
}



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
												std::vector<std::vector<float> >* static_spe_coordinates) {


  OpenBabel::OBConversion obconversion;
  OpenBabel::OBMol mol;
  

  std::vector<std::string> datastring_array;

  std::vector<std::vector<int> > tmp_binary_array;
  std::vector<int> tmp_fp_array;
  std::vector<float> tmp_scalar_array;
  std::vector<float> tmp_static_spe_coordinates;  
  std::vector<std::vector<int> > tmp_split_array;
 	std::vector<float> tmp_scalar_split_array;



  bool has_data=false;
  //Do First loop through first molecule to see how long our vectors will be

  OpenBabel::OBFormat* inFormat = obconversion.FormatFromExt(input_filenames[0]);
  
  obconversion.SetInFormat(inFormat);


  bool notatend = obconversion.ReadFile(&mol, input_filenames[0]);
	unsigned int attribute_bookmark=0;
	unsigned int attribute_counter=0;
  std::vector<OpenBabel::OBGenericData*>::iterator k;
  std::vector<OpenBabel::OBGenericData*> vdata = mol.GetData();

	

	
	//This block designates data as either a scalar or vector value

  for (int i=0; i<SDTags.size(); i++) {

    has_data=false;
		attribute_counter=0;

		
    for (k=vdata.begin();k!=vdata.end();++k) {

      if ((*k)->GetAttribute() == SDTags[i] ) {
        has_data=true;
        std::istringstream datastring (((OpenBabel::OBPairData*)(*k))->GetValue(),std::istringstream::in);
          
        datastring_array=split_string(datastring.str(),'\t');
       
        if (datastring_array.size() > 1) {
          (*binary_weight_array).push_back(SDWeights[i]);
          (*binary_tags).push_back(SDTags[i]);
        }
        else {
          (*scalar_weight_array).push_back(SDWeights[i]);
          (*scalar_tags).push_back(SDTags[i]);
        }
      }
			
    }
  }

  
//  if ((*binary_weight_array).size()) {std::cout << "# binary descriptors: " << (*binary_weight_array).size() << std::endl;}

//  if ((*scalar_weight_array).size()) {std::cout << "# scalar descriptors: " << (*scalar_weight_array).size() << std::endl;}
  




  if (!static_input_filename.empty()) {

    OpenBabel::OBFormat* inFormat = obconversion.FormatFromExt(static_input_filename);    
    obconversion.SetInFormat(inFormat);


		

    bool notatend=obconversion.ReadFile(&mol, static_input_filename);
  
    std::vector< std::string> SPE_data_handles;    

    SPE_data_handles.push_back("SPE.x");
    SPE_data_handles.push_back("SPE.y");

    if (dimension == 3) { SPE_data_handles.push_back("SPE.z"); }

    while (notatend) {

      std::vector<OpenBabel::OBGenericData*>::iterator k;
      std::vector<OpenBabel::OBGenericData*> vdata = mol.GetData();

      
      for (int i=0; i<dimension; i++) {

        has_data=false;
        for (k=vdata.begin();k != vdata.end();++k) {

          if ( (*k)->GetAttribute() == SPE_data_handles[i] ) {

            has_data=true;
            std::istringstream datastring (((OpenBabel::OBPairData*)(*k))->GetValue(),std::istringstream::in);            
            tmp_static_spe_coordinates.push_back( atof( datastring.str().c_str() ) ); 
          }        
        }

        if (!has_data) {

          std::cout << "Missing " << SPE_data_handles[i] << " data. Exiting" << std::endl;
          return 0;

        }
      }
		
			for (int i=0; i<SDTags.size(); i++) {
        has_data=false;
      
      
        for (k=vdata.begin();k!=vdata.end();++k) {

          if ((*k)->GetAttribute() == SDTags[i] ) {
            has_data=true;
            std::istringstream datastring (((OpenBabel::OBPairData*)(*k))->GetValue(),std::istringstream::in);            
            datastring_array=split_string(datastring.str(),'\t');		

            if (datastring_array.size() > 1) {
              for (int j=0; j< datastring_array.size(); j++) {
                tmp_fp_array.push_back( atoi(datastring_array[j].c_str()) );
              }
              tmp_binary_array.push_back(tmp_fp_array); 
              tmp_fp_array.clear();    
            }    
            else {
             // std::cout << atof(datastring_array[0].c_str()) << std::endl;
              tmp_scalar_array.push_back( atof(datastring_array[0].c_str()) );
            }
          } //if
        } //for vdata

        if (!has_data) { 
          std::cout << std::endl << "Missing - " << SDTags[i] << " data from static input file. " << std::endl;          
          return 0; 
        }
      } // for sd tags
	
      (*static_spe_coordinates).push_back(tmp_static_spe_coordinates);
      tmp_static_spe_coordinates.clear();

      notatend=obconversion.Read(&mol);  
	  }	// while loop end




	  for (int i=0; i < (*binary_weight_array).size(); i++ ) {
  	  for (int j=0; j < tmp_binary_array.size(); j++) {
  		  tmp_split_array.push_back(tmp_binary_array[i+j*(*binary_weight_array).size()]);
   	  }
   	  (*static_fingerprint_array).push_back(tmp_split_array);
 		  tmp_split_array.clear();
 	  }
    

 	  for (int i=0; i < (*scalar_weight_array).size(); i++) {
  	  for (int j=0; j < tmp_scalar_array.size(); j+=(*scalar_weight_array).size()) {
   	    tmp_scalar_split_array.push_back(tmp_scalar_array[i+j]);
   	  }
      (*static_scalar_array).push_back(tmp_scalar_split_array);
      tmp_scalar_split_array.clear();
    }

  } // static input file loop end


  for (int fileno=0; fileno<input_filenames.size();fileno++) {

    
    OpenBabel::OBFormat* inFormat = obconversion.FormatFromExt(input_filenames[fileno]);
  
    obconversion.SetInFormat(inFormat); 

    bool notatend = obconversion.ReadFile(&mol, input_filenames[fileno]);


    while (notatend) {

      std::vector<OpenBabel::OBGenericData*>::iterator k;
      std::vector<OpenBabel::OBGenericData*> vdata = mol.GetData();
      for (int i=0; i<SDTags.size(); i++) {
        has_data=false;
      
      
        for (k=vdata.begin();k!=vdata.end();++k) {

          if ((*k)->GetAttribute() == SDTags[i] ) {
            has_data=true;
            std::istringstream datastring (((OpenBabel::OBPairData*)(*k))->GetValue(),std::istringstream::in);            
            datastring_array=split_string(datastring.str(),'\t');
            if (datastring_array.size() > 1) {
              for (int j=0; j< datastring_array.size(); j++) {
                tmp_fp_array.push_back( atoi(datastring_array[j].c_str()) );
              }
              tmp_binary_array.push_back(tmp_fp_array); 
              tmp_fp_array.clear();    
            }    
            else {
             // std::cout << atof(datastring_array[0].c_str()) << std::endl;
              tmp_scalar_array.push_back( atof(datastring_array[0].c_str()) );
            }
          } //if
        } //for vdata

        if (!has_data) { 
          std::cout << std::endl << "Missing - " << SDTags[i] << " data. " << std::endl;          
          return 0; 
        }
      } // for sd tags

      notatend=obconversion.Read(&mol);
      (*molcount)++;
    } // While mols in file loop close
  } //Filename Loop close


  tmp_split_array.clear();
  
  for (int i=0; i < (*binary_weight_array).size(); i++ ) {
    for (int j=0; j < tmp_binary_array.size(); j++) {
      tmp_split_array.push_back(tmp_binary_array[i+j*(*binary_weight_array).size()]);
    }
    (*fingerprint_array).push_back(tmp_split_array);
    tmp_split_array.clear();
  }
    
  tmp_scalar_split_array.clear();
  for (int i=0; i < (*scalar_weight_array).size(); i++) {
    for (int j=0; j < tmp_scalar_array.size(); j+=(*scalar_weight_array).size()) {
      tmp_scalar_split_array.push_back(tmp_scalar_array[i+j]);
    }
    (*scalar_array).push_back(tmp_scalar_split_array);
    tmp_scalar_split_array.clear();
  }


} // Function close








int add_spe_coordinates_to_mol( std::vector<std::vector<float> > spe_coords ,
																std::string static_input_filename,
																std::vector<std::string> filenames, 
																std::string output_filename) {

  int molcounter;

  OpenBabel::OBConversion obconversion;
  OpenBabel::OBMol mol;

  obconversion.SetOutFormat("SDF");
  std::ofstream ofs(output_filename.c_str());
  obconversion.SetOutStream(&ofs);
  
  std::cout << "Writing data to file " << output_filename << "... ";
  std::cout.flush();

	// Write static coords first


	OpenBabel::OBFormat* inFormat = obconversion.FormatFromExt(static_input_filename);
	obconversion.SetInFormat(inFormat);

	bool notatend = obconversion.ReadFile(&mol, static_input_filename);

	while (notatend) {

		obconversion.Write(&mol);
		notatend=obconversion.Read(&mol);
	
	}



  for (int fileno=0; fileno<filenames.size(); fileno++) {

    inFormat = obconversion.FormatFromExt(filenames[fileno]);

    obconversion.SetInFormat(inFormat);

    notatend = obconversion.ReadFile(&mol, filenames[fileno]);    
    molcounter=0;
    while (notatend) {

      OpenBabel::OBPairData *SpeX = new OpenBabel::OBPairData;
      SpeX->SetAttribute("SPE.x");
      std::ostringstream spex;
  
      OpenBabel::OBPairData *SpeY = new OpenBabel::OBPairData;
      SpeY->SetAttribute("SPE.y"); 
      std::ostringstream spey;


      spex << spe_coords[molcounter][0];
      SpeX->SetValue(spex.str());
      spey << spe_coords[molcounter][1];
      SpeY->SetValue(spey.str());
      mol.SetData(SpeX);
      mol.SetData(SpeY);


      if (spe_coords[0].size() == 3) {

        OpenBabel::OBPairData *SpeZ = new OpenBabel::OBPairData;
        std::ostringstream spez; 
        SpeZ->SetAttribute("SPE.z");
        spez << spe_coords[molcounter][2];
        SpeZ->SetValue(spez.str());
        mol.SetData(SpeZ);
      }
  
      obconversion.Write(&mol);
      notatend=obconversion.Read(&mol);
  
      molcounter++;
    }
  
  }
  std::cout << "Finished!" << std::endl;


  return 0;

}










int add_spe_coordinates_to_mol( std::vector<std::vector<float> > spe_coords ,
																std::vector<std::string> filenames, 
																std::string output_filename) {

  int molcounter;

  OpenBabel::OBConversion obconversion;
  OpenBabel::OBMol mol;

  obconversion.SetOutFormat("SDF");
  std::ofstream ofs(output_filename.c_str());
  obconversion.SetOutStream(&ofs);
  
  std::cout << "Writing data to file " << output_filename << "... ";
  std::cout.flush();


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


      spex << spe_coords[molcounter][0];
      SpeX->SetValue(spex.str());
      spey << spe_coords[molcounter][1];
      SpeY->SetValue(spey.str());
      mol.SetData(SpeX);
      mol.SetData(SpeY);


      if (spe_coords[0].size() == 3) {

        OpenBabel::OBPairData *SpeZ = new OpenBabel::OBPairData;
        std::ostringstream spez; 
        SpeZ->SetAttribute("SPE.z");
        spez << spe_coords[molcounter][2];
        SpeZ->SetValue(spez.str());
        mol.SetData(SpeZ);
      }
  
      obconversion.Write(&mol);
      notatend=obconversion.Read(&mol);
  
      molcounter++;
    }
  
  }
  std::cout << "Finished!" << std::endl;


  return 0;

}


int add_spe_coordinates_to_mol( std::vector<std::vector<float> > spe_coords ,
																std::vector<std::string> filenames) {

  int molcounter;

  OpenBabel::OBConversion obconversion;
  OpenBabel::OBMol mol;

	std::vector<std::string> output_filenames;
  std::vector<std::string> split_input;
  std::string output_name;
  

  obconversion.SetOutFormat("SDF");


  
  for (int fileno=0; fileno<filenames.size(); fileno++) {
    
    split_input=split_string(filenames[fileno],'.');
    output_name=split_input[0];
    for (int i=1; i<split_input.size()-1;i++) {
      output_name+="."+split_input[i];
    }
      
  
		output_filenames.push_back(output_name+".spe.sdf");  	
    std::cout << "Writing data to file " << output_filenames[fileno] << "... ";
    std::cout.flush();	
  	std::ofstream ofs(output_filenames[fileno].c_str());
    OpenBabel::OBFormat* inFormat = obconversion.FormatFromExt(filenames[fileno]);

    obconversion.SetInFormat(inFormat);
    obconversion.SetOutStream(&ofs);
    bool notatend = obconversion.ReadFile(&mol, filenames[fileno]);    
    molcounter=0;
    while (notatend) {

      OpenBabel::OBPairData *SpeX = new OpenBabel::OBPairData;
      SpeX->SetAttribute("SPE.x");
      std::ostringstream spex;
  
      OpenBabel::OBPairData *SpeY = new OpenBabel::OBPairData;
      SpeY->SetAttribute("SPE.y"); 
      std::ostringstream spey;


      spex << spe_coords[molcounter][0];
      SpeX->SetValue(spex.str());
      spey << spe_coords[molcounter][1];
      SpeY->SetValue(spey.str());
      mol.SetData(SpeX);
      mol.SetData(SpeY);


      if (spe_coords[0].size() == 3) {
        OpenBabel::OBPairData *SpeZ = new OpenBabel::OBPairData;
        std::ostringstream spez; 
        SpeZ->SetAttribute("SPE.z");
        spez << spe_coords[molcounter][2];
        SpeZ->SetValue(spez.str());
	mol.SetData(SpeZ);
      }
  

      obconversion.Write(&mol);
      notatend=obconversion.Read(&mol);
  
      molcounter++;
    }
  std::cout<<"Finished!" << std::endl;
  }
  return molcounter;
}





