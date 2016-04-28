#ifndef LIMITE_CRANE_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define LIMITE_CRANE_H


//fichier itk
#include "itkimage.h"
#include "itkImageDuplicator.h"
#include "itkImageBase.h"
#include "gradient.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageIOBase.h"
#include "itkNiftiImageIO.h"
#include "itkIndex.h"
#include <iostream>
#include <cmath>
#include "itkArray.h"

using namespace std;

typedef itk::Image<double,3> ImageType;
typedef itk::ImageFileWriter<ImageType >  WriterType;

class Limite_crane{
public:
 Limite_crane(ImageType::Pointer imageitk_us);
 void findlimite();


private:
	ImageType::Pointer roi_crane;
	bool limite;
	int epaisseur;
	int coord;
	int position;
	int dist_limite[200];
		
	

};

#endif