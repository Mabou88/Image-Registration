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
#include "itkLinearInterpolateImageFunction.h"

using namespace std;

typedef itk::Image<double,3> ImageType;
typedef itk::Image<unsigned char, 3> BinaryImageType;
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
 void maskcrane(ImageType::Pointer mask_us);
 void calculvolumecrane();
 //void modifiederriere();
 void sauvegardeimage();
 ImageType::Pointer getcraneROI();

 void setdim();
 


private:
	ImageType::Pointer roi_crane;
	int limitepost_value;
	int max(int,int);
	int dim_x;
	int dim_y;
	int dim_z;
	//bool limite;
	//int epaisseur;
	//int coord;
	//int limite_z[200];
	//int limite_x[200];
	//int position;
	

};

#endif