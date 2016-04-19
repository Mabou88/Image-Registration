#ifndef LC2_H
#define LC2_H

#include "itkimage.h"
#include "mitkImage.h"
#include "mitkBaseGeometry.h"
#include <mitkImageToItk.h>
#include <mitkITKImageImport.h>

typedef itk::Image<short,3> ImageType;

//mes fichiers
#include "Matrice.h"
#include "Normalisationimage.h"
#include "gradient.h"

class Lc2{
public:
	Lc2(mitk::BaseGeometry::Pointer geo_irm,mitk::BaseGeometry::Pointer geo_us,mitk::Image::Pointer image_irm,mitk::Image::Pointer image_us);
	double getlc2();
private:
	double lc2final;
	mitk::BaseGeometry::Pointer geo_irm_;
	mitk::BaseGeometry::Pointer geo_us_;
	mitk::Image::Pointer image_irm_;
	mitk::Image::Pointer image_us_;
};
#endif