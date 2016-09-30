#ifndef SIMILARITYMEASURE_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define SIMILARITYMEASURE_H


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
#include "itkSimilarityIndexImageFilter.h"
#include "itkLabelOverlapMeasuresImageFilter.h"
#include "itkLabelMap.h"
#include "itkLabelObject.h"
#include "itkBinaryImageToLabelMapFilter.h"

using namespace std;

typedef itk::Image<short,3> ImageType;
typedef itk::Image<unsigned char, 3> BinaryImageType;
typedef itk::ImageFileWriter<ImageType >  WriterType;
//typedef itk::ImageFileWriter<LabelMapType >  LabelWriterType;
typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::ImageFileReader<BinaryImageType> BinaryReaderType;
typedef itk::SimilarityIndexImageFilter <BinaryImageType,BinaryImageType> SimilarityIndex;
typedef itk::NumericTraits<ImageType>::RealType RealType;
typedef itk::LabelObject<short,3> LabelObjectType;
typedef itk::LabelMap<LabelObjectType> LabelMapType;
typedef itk::LabelOverlapMeasuresImageFilter <LabelMapType> OverlapMeasurement;
typedef itk::BinaryImageToLabelMapFilter < BinaryImageType,LabelMapType> LabelMapFilter;
typedef itk::ImageToImageFilter <LabelMapType,LabelMapType> ImageFilter;
typedef itk::ImageRegionConstIterator<BinaryImageType> ImageBinaryConstIteratorType;

class Similaritymeasure {
public:
	Similaritymeasure();
	float getindex();
	//compute the overlap of the 2 segmentation and calcul dice coefficient
	void computeoverlap();
private:
	BinaryImageType::Pointer image_us_seg;
	BinaryImageType::Pointer image_irm_seg;
	float similarityindex;
};
#endif SIMILARITYMEASURE_H