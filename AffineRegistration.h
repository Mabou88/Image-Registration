#ifndef AFFINEREGISTRATION_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define AFFINEREGISTRATION_H


//fichier itk
#include "itkimage.h"
#include "itkImageFileReader.h"
#include "itkImageIOBase.h"
#include "itkNiftiImageIO.h"
#include "itkNrrdImageIO.h"


//fichier registration
#include "itkAmoebaOptimizer.h"
#include "itkCastImageFilter.h"
#include "itkEllipseSpatialObject.h"
#include "itkImage.h"
#include "itkImageRegion.h"
#include "itkIndex.h"
#include "itkSize.h"
#include "itkImageRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkGradientDescentOptimizer.h"
#include <itkQuaternionRigidTransformGradientDescentOptimizer.h>
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkAffineTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkMutualInformationImageToImageMetric.h"

using namespace std;

typedef itk::Image<short,3> ImageType;

void AffineRegistration(ImageType::Pointer itkimageirm ,ImageType::Pointer itkimageus);

#endif