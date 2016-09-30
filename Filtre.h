#ifndef FILTRE_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define FILTRE_H

#include "gradient.h"

//fichier itk
#include "itkimage.h"
#include "itkImageDuplicator.h"
#include "itkImageBase.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageIOBase.h"
#include "itkNiftiImageIO.h"
#include "itkIndex.h"
#include <iostream>
#include <cmath>
#include "itkArray.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCropImageFilter.h"
#include "itkResampleImageFilter.h"

#include "itkFFTWForwardFFTImageFilter.h"
#include "itkFFTWInverseFFTImageFilter.h"
#include "itkVnlForwardFFTImageFilter.h"
#include "itkVnlInverseFFTImageFilter.h"

#include "itkMaskImageFilter.h"
#include "itkConvolutionImageFilter.h"
#include "itkImageRegionIterator.h"

#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkMultiplyImageFilter.h"

#include "itkDiscreteGaussianImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkLaplacianImageFilter.h"
#include "itkZeroCrossingBasedEdgeDetectionImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"


#include "vnl/vnl_vector.txx"
#include "vnl/vnl_c_vector.txx"
#include "vnl/algo/vnl_fft_base.txx"
#include "vnl/algo/vnl_fft_prime_factors.txx"
//#include "vnl/algo/vnl_fft_base.h"
//#include "vnl/algo/vnl_fft_prime_factors.h"
//#include "vnl/algo/vnl_fft.cxx"
//#include "vnl/algo/vnl_fft.h"


using namespace std;

typedef itk::Image<short,3> ImageType;
typedef itk::Image<float,3> ImageTypeLap;
typedef itk::Image<unsigned char, 3> BinaryImageType;

typedef itk::ImageFileWriter<ImageType >  WriterType;
typedef itk::ImageFileWriter<ImageTypeLap> WriterLapType;

//typedef itk::ImageFileWriter<LabelMapType >  LabelWriterType;
typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::ImageFileReader<BinaryImageType> BinaryReaderType;
typedef itk::ImageRegionConstIterator<BinaryImageType> ImageBinaryConstIteratorType;
typedef itk::ImageRegionConstIterator<ImageType> ImageConstIteratorType;
typedef itk::ImageRegionIterator<ImageType> ImageIteratorType;
typedef itk::CropImageFilter<ImageType, ImageType> CropFilter;

typedef itk::ConvolutionImageFilter<ImageType,ImageType,ImageType> FilterType;
typedef itk::AddImageFilter<ImageType,ImageType,ImageType> AddFilter;
typedef itk::SubtractImageFilter<ImageType,ImageType,ImageType> SubtractFilter;
typedef itk::MultiplyImageFilter<ImageType,ImageType,ImageType> MultiplyFilter;

typedef itk::DiscreteGaussianImageFilter<ImageType, ImageTypeLap >  GaussianType;
typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType >  GaussianType2;
typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
typedef itk::RescaleIntensityImageFilter<ImageType,ImageType> RescaleFilterType;

typedef itk::LaplacianImageFilter<ImageTypeLap,ImageTypeLap > LaplacianFilter;
typedef itk::ZeroCrossingImageFilter< ImageTypeLap,  ImageTypeLap>     ZeroCrossingFilter;
typedef itk::ThresholdImageFilter <ImageType> ThresholdFilter;
typedef itk::BinaryThresholdImageFilter <ImageType,ImageType> BinaryThresholdFilter;


//typedef itk::FFTWForwardFFTImageFilter < ImageType > FFTFilterType;
//typedef itk::VnlForwardFFTImageFilter< ImageType >  VnlFFTFilterType;
//typedef FFTFilterType::OutputImageType ComplexImageType;

class Filtre {
public:
	Filtre(ImageType::Pointer imagetofiltre);
	void filtregaussian();
	void filtrepassehaut();
	void fft_transform();
	void createkernel(unsigned int largueur);
	void selectROI();
	void sauvegardeimage();
	void imageaddition();
	void threshold();
	ImageType::Pointer getresult();

private:
	ImageType::Pointer imagefiltre;
	ImageTypeLap::Pointer imagelaplacian;
	ImageType::Pointer kernel;
	ImageType::Pointer imageresult;
	ImageType::Pointer filtre;
	//ComplexImageType::Pointer imagefiltreFFT;
};
#endif FILTRE_H