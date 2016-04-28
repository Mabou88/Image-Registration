#ifndef ROI_CRANE_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define ROI_CRANE_H


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

class Roi_crane{
public:
 Roi_crane(ImageType::Pointer imageitk_us);
 void limiteposterieur();
 void limitedroit( int z);
 void limitegauche( int z);
 void limiteinferieuravant();
 void limiteinferieurarriere();
 void setlimitepost_value(int v);
 void bordureexterieur();
 void zonegrise();
 //void modifiederriere();
 void sauvegardeimage();
 ImageType::Pointer getcraneROI();


private:
	ImageType::Pointer roi_crane;
	int limitepost_value;
	int max(int,int);
	//bool limite;
	//int epaisseur;
	//int coord;
	//int limite_z[200];
	//int limite_x[200];
	//int position;
	

};

#endif