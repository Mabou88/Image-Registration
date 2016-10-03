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

typedef itk::Image<short,3> ImageType;
typedef itk::Image<unsigned char, 3> BinaryImageType;
typedef itk::ImageFileWriter<ImageType >  WriterType;

class Roi_crane{
public:
 Roi_crane(ImageType::Pointer imageitk_us,int type);
 void limiteposterieur();
 void limitedroit( int z);
 void limitegauche( int z);
 void limiteinferieuravant();
 void limiteinferieurarriere();
 void setlimitepost_value(int v);
 void bordureexterieur();
 void zonegrise();
 void maskcrane(BinaryImageType::Pointer mask_us);
 void calculvolumecrane();
 //void modifiederriere();
 void sauvegardeimage();
 ImageType::Pointer getcraneROI();
 void makeirmus();
 void makeusirm();

 void setdim();
 void setlimitexy(int z);
 void limiteparcercle(int z);
 void regionexterieur(int z);
 void regioninferieur(int z);
 void regionssuperieur();
 int getlargueurcerveau();
 ImageType::IndexType getcentrecerveau();

 
 


private:
	ImageType::Pointer roi_crane;
	int limitepost_value;
	int max(int,int);
	int dim_x;
	int dim_y;
	int dim_z;
	int min_x;
	int max_x;
	int min_y;
	int max_y;
	int lim_min_x;
	int lim_max_x;
	int lim_min_y;
	int lim_max_y;
	
	ImageType::IndexType centre;
	int largueurcerveau;
	int type_;
	//bool limite;
	//int epaisseur;
	//int coord;
	//int limite_z[200];
	//int limite_x[200];
	//int position;
	

};

#endif