//
//  main.cpp
//  LC2_M
// Calcul de la metrique LC2
//  Created by Maxime Gérard on 17/12/15.
//  Copyright © 2015 Maxime Gérard. All rights reserved.
//

#include <iostream>
#include <string>
#include <fstream>
#include <chrono>

#include "gradient.h"
#include "LC2_function.hpp"
#include "Roi_crane.h"
#include "AffineRegistration.h"
#include "Similaritymeasure.h"
#include "Filtre.h"


#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegistrationMethod.h"
#include "itkAmoebaOptimizer.h"
#include "itkNiftiImageIO.h"
#include "itkPNGImageIO.h"
#include "itkMetaImageIO.h"
#include "itkImageIterator.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRigid3Dtransform.h"
#include "itkAmoebaOptimizer.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkEuler3DTransform.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include"itkLC2ImageToImageMetric.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkShrinkImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkCommand.h"
#include "itkAffineTransform.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkCastImageFilter.h"
#include "itkCheckerBoardImageFilter.h"
#include "itkMultiResolutionPyramidImageFilter.h"
#include "itkBsplineTransform.h"
#include "itkAdaptiveHistogramEqualizationImageFilter.h"
#include "itkCenteredEuler3DTransform.h"
#include <cmath>


#include <dlib\optimization.h>
#include <dlib\matrix.h>

//#include <optimization.h>
//#include <matrix.h>

using namespace std;
using namespace dlib;




//images
typedef itk::Image<short,3> ImageType;
typedef itk::Image<unsigned char, 3> BinaryImageType;
typedef itk::Image<double,2> Image2DType;
typedef itk::Vector<int,3> VectorGeo;

//IO
typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::ImageFileReader<Image2DType> Reader2DType;
typedef itk::ImageFileReader<BinaryImageType> BinaryReaderType;
typedef itk::ImageFileWriter<ImageType> WriterType;

//iteration over images
typedef itk::ImageRegionIterator<ImageType> IteratorType;

//recalage
typedef itk::AmoebaOptimizer AmoebaOptimizerType;
typedef itk::LevenbergMarquardtOptimizer LMOptimizerType;
typedef itk::LinearInterpolateImageFunction<ImageType,double> InterpolatorType;
typedef itk::ImageRegistrationMethod<ImageType, ImageType> RegistrationType;
typedef itk::Euler3DTransform<double> EulerTransformType;
typedef itk::CenteredEuler3DTransform <double> EulerCenterTransformType;
typedef itk::LC2ImageToImageMetric<ImageType, ImageType> LC2MetricType;
typedef itk::TranslationTransform<double,3> TranslationType;
typedef itk::AffineTransform<double,3> AffineTransformType;
typedef itk::BSplineTransform<double,3,3> BSplineTransformType;
typedef BSplineTransformType::ParametersType     BSParametersType;
typedef itk::AdaptiveHistogramEqualizationImageFilter<ImageType> HistoEqualizerType;
typedef itk::TranslationTransform<double,3> TranslationTransformType;

//pour le mapping de l'image mobile registree ac tsf determinee par registration framework
typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
typedef itk::ResampleImageFilter<MaskType, MaskType> ResamplerBinaryType;

//pour la mise a l'echelle de la bande passante des intensites
typedef itk::RescaleIntensityImageFilter<ImageType,ImageType> RescaleFilterType;
typedef itk::MinimumMaximumImageCalculator<ImageType> MinMaxCalculatorType;
typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkFilterType;

//MultiResolutionRegistration


using namespace std;


//classe pour suivre le recalage
class CommandIterationUpdate : public itk::Command
{
public :
    typedef CommandIterationUpdate Self;
    typedef itk::Command SuperClass;
    typedef itk::SmartPointer<Self> Pointer;
    itkNewMacro(Self);
    
protected:
    CommandIterationUpdate()
    {
        m_IterationNumber =0;
    }
public:
    typedef itk::AmoebaOptimizer OptimizerType;
    typedef const OptimizerType * OptimizerPointer;
    
    void Execute(itk::Object *caller, const itk::EventObject &event) ITK_OVERRIDE
    {
        Execute( (const itk::Object *)caller, event);
    }
    
    void Execute(const itk::Object * object,
                 const itk::EventObject & event) ITK_OVERRIDE
    {
        OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
        if( ! itk::IterationEvent().CheckEvent( &event ) )
        {
            return;
        }
        std::cout << m_IterationNumber++ << "   "<<endl;
        std::cout << optimizer->GetCachedValue() << "   "<<endl;
        std::cout << optimizer->GetCachedCurrentPosition() << std::endl;
    }
private:
    unsigned long m_IterationNumber;
};




int main(int argc, const char * argv[]) {
    
    //to determine time of computation
    std::srand(time(NULL));
    std::time_t tbegin,tend;
    double texec = 0;
    tbegin = std::time(NULL);
   
    //lecture des images
    
    string filenameUS;
    string filenameIRM;
    string outputPath;
    string filenameMaskLiver;
    
    bool useLiverMask=false;
    
    /*********************
     * RUNTIME ARG
     ******************/
    
    cout<<"runtime arguments acquisition"<<endl;
    
    for(int i = 0; i < argc; ++i)
    {
        //input image US
        if(strcmp(argv[i], "-iUS")==0)
        {
            i++;
            filenameUS = argv[i];
        }
        
        //input image IRM
        if(strcmp(argv[i], "-iIRM")==0)
        {
            i++;
            filenameIRM = argv[i];
        }
        
        if(strcmp(argv[i], "-iMaskLiver")==0)
        {
            i++;
            filenameMaskLiver= argv[i];
            useLiverMask=true;
            cout<<"Use of mask image of liver"<<endl;
        }
        
        if(strcmp(argv[i], "-o")==0)
        {
            i++;
            outputPath = argv[i];
        }
        
        
    }
    
    /**********************
     * VERIFICATION INPUTS
     ********************/
    
    cout<<"input error handling"<<endl;
    
    if(filenameUS == "")
    {
        cerr<<"Input US file not provided"<<endl;
        return EXIT_FAILURE;
        
    }
    
    if(filenameIRM == "")
    {
        cerr<<"Input MRI file not provided"<<endl;
        return EXIT_FAILURE;
    }
    
    
    if(outputPath == "")
    {
        cerr<<"output path not provided"<<endl;
        return EXIT_FAILURE;
    }
    

    
    /**********************
     * US READING
     *********************/
    
    cout<<"Reading images"<<endl;
    
    ImageType::Pointer image_US = ImageType::New();
    ImageType::Pointer image_IRM = ImageType::New();
    
    
    ReaderType::Pointer reader1 = ReaderType::New();
    itk::NiftiImageIO::Pointer m_io = itk::NiftiImageIO::New();
    reader1->SetImageIO(m_io);
    reader1->SetFileName(filenameUS);
    try {
        reader1->Update();
    } catch (itk::ExceptionObject &e) {
        cout<<"Error while reading US image"<<endl;
        cout<<e<<endl;
        EXIT_FAILURE;
    }
    
    image_US = reader1->GetOutput();
    cout<<"test lecture US"<<endl;
    cout<<"dimensions US : "<<image_US->GetLargestPossibleRegion().GetSize()<<endl;
    
    
    //min max image US
    MinMaxCalculatorType::Pointer minMaxUS = MinMaxCalculatorType::New();
    minMaxUS->SetImage(image_US);
    minMaxUS->Compute();
    
    cout<<"intensity range US image : "<<"[ "<<minMaxUS->GetMinimum()<<","<<minMaxUS->GetMaximum()<<" ]"<<endl;
    
    /**********************
     * IRM READING
     *********************/
    
    ReaderType::Pointer reader2 = ReaderType::New();
    itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
    reader2->SetImageIO(io);
    reader2->SetFileName(filenameIRM);
    try {
        reader2->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"Error while reading IRM image"<<endl;
        cerr<<e<<endl;
        EXIT_FAILURE;
    }
    
    image_IRM = reader2->GetOutput();
    cout<<"test lecture IRM"<<endl;
    cout<<"dimensions IRM : "<<image_IRM->GetLargestPossibleRegion().GetSize()<<endl;
    
    //min max image IRM
    MinMaxCalculatorType::Pointer minMaxIRM = MinMaxCalculatorType::New();
    minMaxIRM->SetImage(image_IRM);
    minMaxIRM->Compute();
    
    cout<<"initial intensity range IRM image : "<<"[ "<<minMaxIRM->GetMinimum()<<","<<minMaxIRM->GetMaximum()<<" ]"<<endl;
    
    cout<<"done reading images"<<endl;
    
    /********************
     * Liver mask reading
     **********************/
    
    //gotta declare it outside if loop to be visible within the entire main function
    BinaryImageType::Pointer LiverMask = BinaryImageType::New();
    
    if(useLiverMask)
    {
        cout<<"Reading MRI liver mask"<<endl;
        BinaryReaderType::Pointer BinReader = BinaryReaderType::New();
        BinReader->SetImageIO(io);
        BinReader->SetFileName(filenameMaskLiver);
        try {
            BinReader->Update();
        } catch (itk::ExceptionObject &e) {
            cout<<"Error while reading liver mask"<<endl;
            cout<<e<<endl;
        }
        
       LiverMask =BinReader->GetOutput();
    }

    
    
	
	//*************************************************************************
  // Registration
  //*************************************************************************
	//AffineRegistration(image_IRM,image_US);
	//Note mettre le rapport de spacing direct et enlever le shrink de 2 de IRM
    /*******************
     * DOWNSAMPLING US
     ******************/
    //RESOLUTION ADAPTATION US
    ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
    shrinkFilter->SetInput(image_US);
    
    //recuperation des spacing des deux images pour savoir quel facteur utiliser
    ImageType::SpacingType spacingUS = image_US->GetSpacing();
    ImageType::SpacingType spacingIRM = image_IRM->GetSpacing();
   
	//calcul du spacing factor
    //int shrinkX = int(spacingIRM[0]/spacingUS[0]);
    //int shrinkY = int(spacingIRM[1]/spacingUS[1]);
    //int shrinkZ = int(spacingIRM[2]/spacingUS[2]);

	int shrinkX = 1;
    int shrinkY = 1;
    int shrinkZ = 1;
    
    cout<<"shrinking factors : "<<shrinkX<<", "<<shrinkY<<", "<<shrinkZ<<endl;

    //porc 6 : (3,3,2)/ porc 1 (5,4,3)
    shrinkFilter->SetShrinkFactor(0, shrinkX);
    shrinkFilter->SetShrinkFactor(1, shrinkY);
    shrinkFilter->SetShrinkFactor(2, shrinkZ);
    try {
        shrinkFilter->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error while downsampling US image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
	
	
    
   ImageType::Pointer US_shrunk = shrinkFilter->GetOutput();
	// ImageType::Pointer US_shrunk =image_US;
//        //verification ecriture de l'image
//            WriterType::Pointer writer6 = WriterType::New();
//            string out6 = outputPath+"/ShrunkUS.nii.gz";
//            //itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
//            writer6->SetImageIO(io);
//            writer6->SetInput(US_shrunk);
//            writer6->SetFileName(out6);
//            try {
//                writer6->Update();
//            } catch (itk::ExceptionObject &e) {
//                cerr<<"error while writing rescaled image"<<endl;
//                cerr<<e<<endl;
//                return EXIT_FAILURE;
//            }
//
    cout<<"done writing shrunk US"<<endl;
    cout<<"taille US shrunk : "<<US_shrunk->GetLargestPossibleRegion().GetSize()<<endl;
    

    /*******************
     * DOWNSAMPLING IRM
     *******************/
    
    //APPROCHE MUTLIRES
    
    ShrinkFilterType::Pointer shrinkerMRI = ShrinkFilterType::New();
    shrinkerMRI->SetInput(image_IRM);
    shrinkerMRI->SetShrinkFactor(0, 1); //changer à 1 ,(2 pour MRI avant)
    shrinkerMRI->SetShrinkFactor(1, 1);
	shrinkerMRI->SetShrinkFactor(2, 1);
    try {
        shrinkerMRI->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error while shrinking MRI image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    ImageType::Pointer MRI_shrunk = shrinkerMRI->GetOutput();
    
    cout<<"taille MRI shrunk : "<<MRI_shrunk->GetLargestPossibleRegion().GetSize()<<endl;
    
    /****************
     * RESCALING IRM
     ***************/
    
    //MAKE SUR MRI INTENSITIES ARE WITHIN [0,255] -> why not put the US data in the MRI range ?
    cout<<"rescaling de l'image IRM"<<endl;
    RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
    //rescaler->SetInput(image_IRM);
    rescaler->SetInput(MRI_shrunk);
    rescaler->SetOutputMinimum(0);
    rescaler->SetOutputMaximum(255);
    try {
        rescaler->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error while rescaling image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    
    ImageType::Pointer rescaled_IRM = rescaler->GetOutput();
    
//    //ecrire res
//    WriterType::Pointer writer10 = WriterType::New();
//    writer10->SetImageIO(io);
//    writer10->SetInput(rescaled_IRM);
//    string out10 = outputPath+"/rescaled_MRI.nii.gz";
//    writer10->SetFileName(out10);
//    try {
//        writer10->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"Error while writing rescaled MRI"<<endl;
//        cerr<<e<<endl;
//    }

    //min max image IRM
    MinMaxCalculatorType::Pointer minMaxIRM2 = MinMaxCalculatorType::New();
    minMaxIRM2->SetImage(rescaled_IRM);
    minMaxIRM2->Compute();

    cout<<"rescaled intensity range IRM image : "<<"[ "<<minMaxIRM2->GetMinimum()<<","<<minMaxIRM2->GetMaximum()<<" ]"<<endl;
 //*****************************************************
	ImageType::Pointer mask_US = ImageType::New();
    
        
    ReaderType::Pointer reader_m = ReaderType::New();
    itk::NiftiImageIO::Pointer mask_io = itk::NiftiImageIO::New();
    reader_m->SetImageIO(mask_io);
    reader_m->SetFileName("C:/Users/Marc-Antoine/Documents/Imagecode/input/02-05_t1w.nii");
    try {
        reader_m->Update();
    } catch (itk::ExceptionObject &e) {
        cout<<"Error while reading US mask"<<endl;
        cout<<e<<endl;
        EXIT_FAILURE;
    }
    
    mask_US = reader_m->GetOutput();


	WriterType::Pointer   writer00 =  WriterType::New();
		itk::NiftiImageIO::Pointer ioimagenifti00=itk::NiftiImageIO::New();

		 string outputfile00="C:/Users/Marc-Antoine/Documents/Imagecode/output/US_avant.nii";
		

		 writer00->SetImageIO(ioimagenifti00);
		writer00->SetFileName( outputfile00);
		writer00->SetInput(image_US);
		 writer00->Update();
//*************************************************************************	
//SimilarityMeasure
//*************************************************************************
		//Similaritymeasure similaritymeasure;
		//similaritymeasure.computeoverlap();
		//float dicecoefficient=similaritymeasure.getindex();
		//cout<<"Similarity index (Dice): "<< float(dicecoefficient)<<endl;







//*************************************************************************	
//Align image center
//************************************************************************* 
	
	//trouve les centres
	ImageType::SizeType size_irm=rescaled_IRM->GetLargestPossibleRegion().GetSize();
	ImageType::SizeType size_us=US_shrunk->GetLargestPossibleRegion().GetSize();
	ImageType::IndexType centre_irm;
	centre_irm[0]=size_irm[0]/2;
	centre_irm[1]=size_irm[1]/2;
	centre_irm[2]=size_irm[2]/2;
	ImageType::IndexType centre_us;
	centre_us[0]=size_us[0]/2;
	centre_us[1]=size_us[1]/2;
	centre_us[2]=size_us[2]/2;

	ImageType::PointType centre_spatial_irm;
	rescaled_IRM->TransformIndexToPhysicalPoint(centre_irm,centre_spatial_irm);
	ImageType::PointType centre_spatial_us;
	US_shrunk->TransformIndexToPhysicalPoint(centre_us,centre_spatial_us);

	cout<<"us centre"<<centre_spatial_us<<endl;
	cout<<"irm centre"<<centre_spatial_irm<<endl;

	

	EulerTransformType::Pointer translation =  EulerTransformType::New();
	EulerTransformType::ParametersType parametre(6);
	
	parametre[0]=0;
	parametre[1]=0;
	parametre[2]=0;
	parametre[3] = centre_spatial_us[0]-centre_spatial_irm[0];
	parametre[4] = centre_spatial_us[1]-centre_spatial_irm[1];
	parametre[5] = centre_spatial_us[2]-centre_spatial_irm[2];
	parametre[3]=parametre[3]*-1;
	parametre[4]=parametre[4]*-1;
	parametre[5]=parametre[5]*-1;
	translation->SetParameters(parametre);
	
	cout<< "parametre"<<parametre<<endl;
	
  
	//transformation
	
	
	ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
	resampleFilter->SetTransform(translation);
	resampleFilter->SetInput(rescaled_IRM);
	ImageType::SizeType   size = rescaled_IRM->GetLargestPossibleRegion().GetSize();
	resampleFilter->SetSize( size );
	resampleFilter->SetOutputSpacing(rescaled_IRM->GetSpacing());
	resampleFilter->SetOutputDirection(translation->GetInverseMatrix()*rescaled_IRM->GetDirection());
	resampleFilter->SetOutputOrigin(translation->GetInverseTransform()->TransformPoint(rescaled_IRM->GetOrigin()));
	resampleFilter->Update();
	rescaled_IRM = resampleFilter->GetOutput();

	if(useLiverMask){
	ResamplerBinaryType::Pointer resampleFilter2 = ResamplerBinaryType::New();
	resampleFilter2->SetTransform(translation);
	resampleFilter2->SetInput(LiverMask);
	ImageType::SizeType   sizemask = LiverMask->GetLargestPossibleRegion().GetSize();
	resampleFilter2->SetSize( sizemask );
	resampleFilter2->SetOutputSpacing(LiverMask->GetSpacing());
	resampleFilter2->SetOutputDirection(translation->GetInverseMatrix()*LiverMask->GetDirection());
	resampleFilter2->SetOutputOrigin(translation->GetInverseTransform()->TransformPoint(LiverMask->GetOrigin()));
	resampleFilter2->Update();
	LiverMask = resampleFilter2->GetOutput();
	}
	//ImageType::Pointer sizeImage = ImageType::New();
    
	


//*************************************************************************	
//Rotation a US
//************************************************************************* 
	
	EulerTransformType::Pointer rotation2 =  EulerTransformType::New();
	EulerTransformType::ParametersType parametrerot2(6);
	EulerTransformType::FixedParametersType rotationcenter2(3);
	int typeorientation=2; //1 v3 de base--- 2 v3 a lenvers
	double angle90rad2=3.1416/2;
	double angle180rad=3.1416;

	if(typeorientation==1){
	parametrerot2[0]=0; //sagital vert reste
	parametrerot2[1]=angle90rad2; //coronal bleu reste
	parametrerot2[2]=-angle90rad2;//axial (rouge) reste
	} else if (typeorientation==2){
		parametrerot2[0]=angle180rad; //sagital vert reste
		parametrerot2[1]=angle90rad2; //coronal bleu reste
		parametrerot2[2]=-angle90rad2;


	}


	parametrerot2[3] = 0;
	parametrerot2[4] = 0;
	parametrerot2[5] = 0;
	rotationcenter2[0]=centre_spatial_us[0];
	rotationcenter2[1]=centre_spatial_us[1];
	rotationcenter2[2]=centre_spatial_us[2];
	rotation2->SetFixedParameters(rotationcenter2);
	rotation2->SetParameters(parametrerot2);
	
	

	cout<<"centre rot"<<centre_spatial_us<<endl;
	
	
	ResampleFilterType::Pointer resampler2 = ResampleFilterType::New();
    resampler2->SetInput(US_shrunk);
	resampler2->SetTransform(rotation2);

	resampler2->SetSize(US_shrunk->GetLargestPossibleRegion().GetSize());
    resampler2->SetOutputSpacing(US_shrunk->GetSpacing());
   resampler2->SetOutputDirection(US_shrunk->GetDirection());
   resampler2->SetOutputOrigin((US_shrunk->GetOrigin()));
   resampler2->Update();

   US_shrunk = resampler2->GetOutput();


  

	

	

    
		
		 
		 
		 
		 
		 //*************************************************************************
  // Selectionne la region of interest du crane avec le gradient
  //*************************************************************************
	
	Roi_crane crane(US_shrunk,1);
	
	
	
	crane.setdim();
	
	
	for (int z=60;z<180;z++){
		//crane.setlimitexy(z);
		//crane.limiteparcercle(z);
		//crane.regionexterieur(z);
	}
	crane.regioninferieur(62);
	crane.regionssuperieur();

	
	
	
	//*************************************************************************
  // technique avec le masque
  //*************************************************************************
	/*
	ImageType::Pointer mask_US = ImageType::New();
    
        
    ReaderType::Pointer reader_m = ReaderType::New();
    itk::NiftiImageIO::Pointer mask_io = itk::NiftiImageIO::New();
    reader_m->SetImageIO(mask_io);
    reader_m->SetFileName("C:/Users/Marc-Antoine/Documents/Imagecode/input/02-05_mask.nii");
    try {
        reader_m->Update();
    } catch (itk::ExceptionObject &e) {
        cout<<"Error while reading US mask"<<endl;
        cout<<e<<endl;
        EXIT_FAILURE;
    }
    
    mask_US = reader_m->GetOutput();
	crane.maskcrane(mask_US);
	
	*/

	crane.sauvegardeimage();
	
	
	US_shrunk=crane.getcraneROI();
	cout<<"region selectionné"<<endl;
	//***********************************************
	//rotation vers l'avant
	//*************************************************

	 EulerTransformType::Pointer rotation3 =  EulerTransformType::New();
	EulerTransformType::ParametersType parametrerot3(6);
	EulerTransformType::FixedParametersType rotationcenter3(3);
	double angle90rad3=3.1416/6;
	parametrerot3[0]=-angle90rad3; //sagital vert reste
	parametrerot3[1]=0; //coronal bleu reste
	parametrerot3[2]=0;//axial (rouge) reste
	parametrerot3[3] = 0;
	parametrerot3[4] = 0;
	parametrerot3[5] = 0;
	rotationcenter3[0]=centre_spatial_us[0];
	rotationcenter3[1]=centre_spatial_us[1];
	rotationcenter3[2]=centre_spatial_us[2];
	rotation3->SetFixedParameters(rotationcenter3);
	rotation3->SetParameters(parametrerot3);
	
	

	cout<<"centre rot"<<centre_spatial_us<<endl;
	
	
	ResampleFilterType::Pointer resampler3 = ResampleFilterType::New();
    resampler3->SetInput(US_shrunk);
	resampler3->SetTransform(rotation3);

	resampler3->SetSize(US_shrunk->GetLargestPossibleRegion().GetSize());
    resampler3->SetOutputSpacing(US_shrunk->GetSpacing());
   resampler3->SetOutputDirection(US_shrunk->GetDirection());
   resampler3->SetOutputOrigin((US_shrunk->GetOrigin()));
   resampler3->Update();

   US_shrunk = resampler3->GetOutput();

	Roi_crane crane2(rescaled_IRM,2);
	
	crane2.setdim();
	crane2.maskcrane(LiverMask);
	rescaled_IRM=crane2.getcraneROI();

	 WriterType::Pointer   writer19 =  WriterType::New();
		itk::NiftiImageIO::Pointer ioimagenifti19=itk::NiftiImageIO::New();

		 string outputfile29="C:/Users/Marc-Antoine/Documents/Imagecode/output/IRMav.nii";
		

		 writer19->SetImageIO(ioimagenifti19);
		writer19->SetFileName( outputfile29);
		writer19->SetInput(rescaled_IRM);
		 writer19->Update();

	cout<<"region selectionné"<<endl;
	
	
	
	//*************************************************************************
  // Calcul de volume
  //*************************************************************************
	
	Roi_crane volumeirm(rescaled_IRM,2);
	volumeirm.setdim();
	volumeirm.calculvolumecrane();
	int largueurirm=volumeirm.getlargueurcerveau();
	
	
	Roi_crane volumeus(US_shrunk,1);
	volumeus.setdim();
	volumeus.regioninferieur(65);
	volumeus.calculvolumecrane();
	int largueurus=volumeus.getlargueurcerveau();
	US_shrunk=volumeus.getcraneROI();

	
	//Trouve le scale factor pour la mise à l'échelle
	float scalefactor_irm=float(largueurirm)/float(largueurus);

	cout<<"scalefactor"<<scalefactor_irm<<endl;

	//correction
	scalefactor_irm=scalefactor_irm*1.04;

//*************************************************************************
  // Mise a l'échelle
  //*************************************************************************


	AffineTransformfixedcenter::Pointer scaletrans=AffineTransformfixedcenter::New();
	AffineTransformfixedcenter::OutputVectorType scalefactortransform1;

		scalefactortransform1[0]=scalefactor_irm;
		scalefactortransform1[1]=scalefactor_irm;
		scalefactortransform1[2]=scalefactor_irm;

		scaletrans->Scale(scalefactortransform1,true);

		// ImageType::PointType centerscale;
       
		//Patient 2
       // center[0] =98;
        //center[1] =68;
        //center[2] =98;

		//Patient 5
		//centerscale[0]=111;
		//centerscale[1]=80;
		//centerscale[2]=111;

		scaletrans->SetCenterOfRotationComponent(centre_spatial_us);

		ResamplerType::Pointer transscaleimage = ResamplerType::New();
		transscaleimage->SetInput(rescaled_IRM);
		 transscaleimage->SetSize(rescaled_IRM->GetLargestPossibleRegion().GetSize());
       transscaleimage->SetOutputSpacing(rescaled_IRM->GetSpacing());
        transscaleimage->SetOutputDirection(rescaled_IRM->GetDirection());
        transscaleimage->SetOutputOrigin(rescaled_IRM->GetOrigin());
        transscaleimage->SetTransform(scaletrans);
		transscaleimage->Update();
		rescaled_IRM=transscaleimage->GetOutput();

//sauvegarde


		WriterType::Pointer   writer0 =  WriterType::New();
		itk::NiftiImageIO::Pointer ioimagenifti0=itk::NiftiImageIO::New();

		 string outputfile="C:/Users/Marc-Antoine/Documents/Imagecode/output/US.nii";
		

		 writer0->SetImageIO(ioimagenifti0);
		writer0->SetFileName( outputfile);
		writer0->SetInput(US_shrunk);
		 writer0->Update();

		 WriterType::Pointer   writer1 =  WriterType::New();
		itk::NiftiImageIO::Pointer ioimagenifti1=itk::NiftiImageIO::New();

		 string outputfile2="C:/Users/Marc-Antoine/Documents/Imagecode/output/IRM.nii";
		

		 writer1->SetImageIO(ioimagenifti1);
		writer1->SetFileName( outputfile2);
		writer1->SetInput(rescaled_IRM);
		 writer1->Update();

//*************************************************************************	
//Filtre
//*************************************************************************
		// Filtre filtre(US_shrunk);
		 //filtre.FFT_transform();

		// Filtre filtrespatial(US_shrunk);
		//filtrespatial.filtrepassehaut();
		//filtrespatial.createkernel(3);
		//filtrespatial.filtregaussian();
		//filtrespatial.threshold();
		//filtrespatial.sauvegardeimage();
		//US_shrunk=filtrespatial.getresult();

//*************************************************************************
  // Translation
  //***********************************************************************

Roi_crane volumeirm2(rescaled_IRM,2);
	volumeirm2.setdim();
	volumeirm2.calculvolumecrane();
	
	ImageType::IndexType centre_irm_init=volumeirm2.getcentrecerveau();
	
	Roi_crane volumeus2(US_shrunk,1);
	volumeus2.setdim();
	volumeus2.regioninferieur(65);
	volumeus2.calculvolumecrane();
	
	ImageType::IndexType centre_us_init=volumeus2.getcentrecerveau();
	US_shrunk=volumeus2.getcraneROI();

	cout<< "centre US"<<centre_us_init<<endl;
	cout<<"centre IRM"<<centre_irm_init<<endl;

	//Trouve la translation pour aligner les bords
	ImageType::PointType pcentre;
	US_shrunk->TransformIndexToPhysicalPoint(centre_us_init,pcentre);
	ImageType::PointType icentre;
	rescaled_IRM->TransformIndexToPhysicalPoint(centre_irm_init,icentre);

	cout<< "centre US spatial"<<pcentre<<endl;
	cout<<"centre IRM spatial"<<icentre<<endl;

	VectorGeo vecteurtranslation;
	vecteurtranslation[0]=icentre[0]-pcentre[0];
	vecteurtranslation[1]=icentre[1]-pcentre[1];
	vecteurtranslation[2]=icentre[2]-pcentre[2];


	VectorGeo vecteurtranslationus;
	VectorGeo vecteurtranslationirm;

	vecteurtranslationirm[0]=icentre[0];
	vecteurtranslationirm[1]=icentre[1];
	vecteurtranslationirm[2]=icentre[2];

	vecteurtranslationus[0]=pcentre[0];
	vecteurtranslationus[1]=pcentre[1];
	vecteurtranslationus[2]=pcentre[2];


	cout<<"Vecteur translation: "<< vecteurtranslation<<endl;

	//en z dimension il manque le haut du crane
	//if (abs(vecteurtranslation[2])==vecteurtranslation[2])

	//correction
		vecteurtranslation[2]=vecteurtranslation[2]+2;
		vecteurtranslation[1]=vecteurtranslation[1]-5;
	//else
		//vecteurtranslation[2]=vecteurtranslation[2]-5;


	

		 cout<<"US et IRM initialise"<<endl;


	//volumeirm.sauvegardeimage();
	//volumeirm.maskcrane(mask_US);
	
	//volumeirm.calculvolumecrane();
		 
		 AffineTransformType::Pointer transla= AffineTransformType::New();

//translate
		 AffineTransformType::OutputVectorType vecteurtranslation_x_scale_ini;
		 AffineTransformType::OutputVectorType vecteurtranslation_y_scale_ini;
		 AffineTransformType::OutputVectorType vecteurtranslation_z_scale_ini;
		 vecteurtranslation_x_scale_ini[0]=1;
		 vecteurtranslation_x_scale_ini[1]=0;
		 vecteurtranslation_x_scale_ini[2]=0;
		 vecteurtranslation_y_scale_ini[0]=0;
		 vecteurtranslation_y_scale_ini[1]=1;
		 vecteurtranslation_y_scale_ini[2]=0;
		 vecteurtranslation_z_scale_ini[0]=0;
		 vecteurtranslation_z_scale_ini[1]=0;
		 vecteurtranslation_z_scale_ini[2]=1;


		 AffineTransformType::OutputVectorType vecteurtranslation_x_scale=vecteurtranslation_x_scale_ini*vecteurtranslation[0]; 
		 AffineTransformType::OutputVectorType vecteurtranslation_y_scale=vecteurtranslation_y_scale_ini*vecteurtranslation[1];
		 AffineTransformType::OutputVectorType vecteurtranslation_z_scale=vecteurtranslation_z_scale_ini*vecteurtranslation[2];

		transla->Translate(vecteurtranslation_x_scale);
		transla->Translate(vecteurtranslation_y_scale);
		transla->Translate(vecteurtranslation_z_scale);

		

		ResamplerType::Pointer transimage = ResamplerType::New();
		transimage->SetInput(rescaled_IRM);
		 transimage->SetSize(rescaled_IRM->GetLargestPossibleRegion().GetSize());
        transimage->SetOutputSpacing(rescaled_IRM->GetSpacing());
        transimage->SetOutputDirection(rescaled_IRM->GetDirection());
        transimage->SetOutputOrigin(rescaled_IRM->GetOrigin());
        transimage->SetTransform(transla);
		transimage->Update();
		rescaled_IRM=transimage->GetOutput();
		 /*
//translation
		 AffineTransformType::Pointer transla= AffineTransformType::New();

//translate
		 AffineTransformType::OutputVectorType vecteurtranslation_x_scale_ini;
		 AffineTransformType::OutputVectorType vecteurtranslation_y_scale_ini;
		 AffineTransformType::OutputVectorType vecteurtranslation_z_scale_ini;
		 vecteurtranslation_x_scale_ini[0]=1;
		 vecteurtranslation_x_scale_ini[1]=0;
		 vecteurtranslation_x_scale_ini[2]=0;
		 vecteurtranslation_y_scale_ini[0]=0;
		 vecteurtranslation_y_scale_ini[1]=1;
		 vecteurtranslation_y_scale_ini[2]=0;
		 vecteurtranslation_z_scale_ini[0]=0;
		 vecteurtranslation_z_scale_ini[1]=0;
		 vecteurtranslation_z_scale_ini[2]=1;


		 AffineTransformType::OutputVectorType vecteurtranslation_x_scale=vecteurtranslation_x_scale_ini*vecteurtranslationirm[0]; 
		 AffineTransformType::OutputVectorType vecteurtranslation_y_scale=vecteurtranslation_y_scale_ini*vecteurtranslationirm[1];
		 AffineTransformType::OutputVectorType vecteurtranslation_z_scale=vecteurtranslation_z_scale_ini*vecteurtranslationirm[2];

		transla->Translate(vecteurtranslation_x_scale);
		transla->Translate(vecteurtranslation_y_scale);
		transla->Translate(vecteurtranslation_z_scale);

		

		ResamplerType::Pointer transimage = ResamplerType::New();
		transimage->SetInput(rescaled_IRM);
		 transimage->SetSize(rescaled_IRM->GetLargestPossibleRegion().GetSize());
        transimage->SetOutputSpacing(rescaled_IRM->GetSpacing());
        transimage->SetOutputDirection(rescaled_IRM->GetDirection());
        transimage->SetOutputOrigin(rescaled_IRM->GetOrigin());
        transimage->SetTransform(transla);
		transimage->Update();
		rescaled_IRM=transimage->GetOutput();


		 AffineTransformType::Pointer transla2 = AffineTransformType::New();

		 AffineTransformType::OutputVectorType vecteurtranslation_x_scale2=vecteurtranslation_x_scale_ini*vecteurtranslationus[0]; 
		 AffineTransformType::OutputVectorType vecteurtranslation_y_scale2=vecteurtranslation_y_scale_ini*vecteurtranslationus[1];
		 AffineTransformType::OutputVectorType vecteurtranslation_z_scale2=vecteurtranslation_z_scale_ini*vecteurtranslationus[2];

		transla2->Translate(vecteurtranslation_x_scale2);
		transla2->Translate(vecteurtranslation_y_scale2);
		transla2->Translate(vecteurtranslation_z_scale2);

		

		ResamplerType::Pointer transimage2 = ResamplerType::New();
		transimage2->SetInput(US_shrunk);
		 transimage2->SetSize(US_shrunk->GetLargestPossibleRegion().GetSize());
        transimage2->SetOutputSpacing(US_shrunk->GetSpacing());
        transimage2->SetOutputDirection(US_shrunk->GetDirection());
        transimage2->SetOutputOrigin(US_shrunk->GetOrigin());
        transimage2->SetTransform(transla2);
		transimage2->Update();
		US_shrunk=transimage2->GetOutput();

		*/

		 WriterType::Pointer   writer11 =  WriterType::New();
		

		 string outputfile22="C:/Users/Marc-Antoine/Documents/Imagecode/output/IRM2.nii";
		

		 writer11->SetImageIO(ioimagenifti1);
		writer11->SetFileName( outputfile22);
		writer11->SetInput(rescaled_IRM);
		 writer11->Update();

		 WriterType::Pointer   writer12 =  WriterType::New();
		

		 string outputfile23="C:/Users/Marc-Antoine/Documents/Imagecode/output/US2.nii";
		

		 writer12->SetImageIO(ioimagenifti1);
		writer12->SetFileName( outputfile23);
		writer12->SetInput(US_shrunk);
		 writer12->Update();

		 cout<<"US et IRM initialise"<<endl;
	

		
	//*************************************************************************
  // Transform de mise à l'échelle (scalable) sur IRM
  //*************************************************************************
		 /*
		 AffineTransformfixedcenter::Pointer transformscalable = AffineTransformfixedcenter::New();
	AffineTransformType::OutputVectorType scalefactortransform;

	scalefactor_irm=1.14;

	scalefactortransform[0]=scalefactor_irm;
	scalefactortransform[1]=scalefactor_irm;
	scalefactortransform[2]=scalefactor_irm;

	transformscalable->Scale(scalefactortransform,true);

	//Définis l'axe de rotation avec un vecteur
		AffineTransformType::OutputVectorType axederotation_x0;
		axederotation_x0[0]=1;
		axederotation_x0[1]=0;
		axederotation_x0[2]=0;
		AffineTransformType::OutputVectorType axederotation_y0;
		axederotation_y0[0]=0;
		axederotation_y0[1]=1;
		axederotation_y0[2]=0;
		AffineTransformType::OutputVectorType axederotation_z0;
		axederotation_z0[0]=0;
		axederotation_z0[1]=0;
		axederotation_z0[2]=1;

		double anglerotation_x0=0.1;
		double anglerotation_y0=0.12;
		double anglerotation_z0=-0.05;

		 AffineTransformType::OutputVectorType  vecteurtranslation_x_scalef=vecteurtranslation_x_scale_ini-5; //p5 ish
		  AffineTransformType::OutputVectorType vecteurtranslation_y_scalef=vecteurtranslation_y_scale_ini*-6;
		  AffineTransformType::OutputVectorType vecteurtranslation_z_scalef=vecteurtranslation_z_scale_ini*5;

		 transformscalable->Translate(vecteurtranslation_x_scalef,true);
		 transformscalable->Translate(vecteurtranslation_y_scalef,true);
		 transformscalable->Translate(vecteurtranslation_z_scalef,true);

		// transformscalable->Rotate3D(axederotation_x0,anglerotation_x0);
		// transformscalable->Rotate3D(axederotation_y0,anglerotation_y0);
		// transformscalable->Rotate3D(axederotation_z0,anglerotation_z0);

		// ImageType::PointType center;
       
		//center[0]=111;
		//center[1]=80;
		//center[2]=111;

		transformscalable->SetCenterOfRotationComponent(centre_spatial_us);

		 //rotation pas le bon centre
		 AffineTransformType::OutputVectorType axederotation_y_scale;
		axederotation_y_scale[0]=1;
		axederotation_y_scale[1]=0;
		axederotation_y_scale[2]=0;

		double anglerotation_y_scale=0.3;

		//transformscalable->Rotate3D(axederotation_y_scale,anglerotation_y_scale,true);

	ResamplerType::Pointer scaleimage = ResamplerType::New();
	scaleimage->SetInput(rescaled_IRM);
      scaleimage->SetSize(rescaled_IRM->GetLargestPossibleRegion().GetSize());
        scaleimage->SetOutputSpacing(rescaled_IRM->GetSpacing());
        scaleimage->SetOutputDirection(rescaled_IRM->GetDirection());
        scaleimage->SetOutputOrigin(rescaled_IRM->GetOrigin());
        scaleimage->SetTransform(transformscalable);
        //resamplefilter->SetTransform(transform);

		//InterpolatorNearestNeighbor::Pointer interpolator = InterpolatorNearestNeighbor::New();
		//resamplefilteraff->SetInterpolator( interpolator );

        
        try {
           // scaleimage->Update();
        } catch (itk::ExceptionObject &e) {
            std::cerr<<"error while scaling irm image"<<std::endl;
            std::cerr<<e<<std::endl;
        }
        
      // rescaled_IRM= scaleimage->GetOutput();
	  */
/*******************
  * US sauvegarde avant recalage
 ******************/

	   WriterType::Pointer   writerus =  WriterType::New();
		itk::NiftiImageIO::Pointer ioimagenifti01=itk::NiftiImageIO::New();

		 string outputfilefinal="C:/Users/Marc-Antoine/Documents/Imagecode/output/US_avantregist.nii";
		

		 writerus->SetImageIO(ioimagenifti01);
		writerus->SetFileName( outputfilefinal);
		writerus->SetInput(US_shrunk);
		 writerus->Update();
	/*******************
     * DOWNSAMPLING US
     ******************/
    //RESOLUTION ADAPTATION US
	
    ShrinkFilterType::Pointer shrinkFilter2 = ShrinkFilterType::New();
    shrinkFilter2->SetInput(US_shrunk);
    
    //recuperation des spacing des deux images pour savoir quel facteur utiliser
    //ImageType::SpacingType spacingUS = image_US->GetSpacing();
    //ImageType::SpacingType spacingIRM = image_IRM->GetSpacing();
   
	//calcul du spacing factor
    //int shrinkX = int(spacingIRM[0]/spacingUS[0]);
    //int shrinkY = int(spacingIRM[1]/spacingUS[1]);
    //int shrinkZ = int(spacingIRM[2]/spacingUS[2]);

	int shrinkX2 = 2;
    int shrinkY2 = 2;
    int shrinkZ2 = 2;
    
    cout<<"shrinking factors : "<<shrinkX2<<", "<<shrinkY2<<", "<<shrinkZ2<<endl;

    //porc 6 : (3,3,2)/ porc 1 (5,4,3)
    shrinkFilter2->SetShrinkFactor(0, shrinkX2);
    shrinkFilter2->SetShrinkFactor(1, shrinkY2);
    shrinkFilter2->SetShrinkFactor(2, shrinkZ2);
    try {
        shrinkFilter2->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error while downsampling US image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
	
	
    
   US_shrunk = shrinkFilter2->GetOutput();
   



	//*************************************************************************	

//sauvegarde l'image modifier
//*************************************************************************
  WriterType::Pointer   writer =  WriterType::New();
  itk::NiftiImageIO::Pointer ioimagenifti=itk::NiftiImageIO::New();

  writer->SetImageIO(ioimagenifti);
  writer->SetFileName( "C:/Users/Marc-Antoine/Documents/Imagecode/output/US_shrunk.nii");
  writer->SetInput(US_shrunk);
  writer->Update();

  WriterType::Pointer   writer2 =  WriterType::New();
  itk::NiftiImageIO::Pointer ioimagenifti2=itk::NiftiImageIO::New();

  writer2->SetImageIO(ioimagenifti2);
  writer2->SetFileName( "C:/Users/Marc-Antoine/Documents/Imagecode/output/rescaled_IRM.nii");
  writer2->SetInput(rescaled_IRM);
  writer2->Update();



  // AffineRegistration(US_shrunk,rescaled_IRM);
  //*************************************************************************	

//ouverture des images déjà recalés
//*************************************************************************
  /*
  ReaderType::Pointer reader_us_reca = ReaderType::New();
  
    reader_us_reca->SetImageIO(ioimagenifti);
    reader_us_reca->SetFileName("C:/Users/Marc-Antoine/Documents/Imagecode/output/Affine_auto_centre/Patient 10/US_shrunk.nii");
    try {
        reader_m->Update();
    } catch (itk::ExceptionObject &e) {
        cout<<"Error while reading"<<endl;
        cout<<e<<endl;
        EXIT_FAILURE;
    }
    	
	
	US_shrunk=reader_us_reca->GetOutput();
	*/
  /*
	ReaderType::Pointer reader_irm_reca = ReaderType::New();
  
    reader_irm_reca->SetImageIO(ioimagenifti);
    reader_irm_reca->SetFileName("C:/Users/Marc-Antoine/Documents/Imagecode/output/finalregistreredUSBOBYQA.nii.gz");
    try {
        reader_irm_reca->Update();
    } catch (itk::ExceptionObject &e) {
        cout<<"Error while reading"<<endl;
        cout<<e<<endl;
        EXIT_FAILURE;
    }
    	
	
	rescaled_IRM=reader_irm_reca->GetOutput();
	*/
	

    /*************************
     * RECALAGE
    *************************/
    
    cout<<"registration routine"<<endl;
    
    //test cropping MRI
    
//    LC2_function LC2 = LC2_function(rescaled_IRM, US_shrunk);
//    LC2.limitMRI();
    
    //they must be column_vector !!!
    
    //BOBYQA initiaux FULL RES
    //Porc 1
    //best score : 0.183876
    //final parameters : [-0.03514169419515513, 0.04780031905516767, 0.08510374453376125, -3.16002229212777, -0.7603111683282886, -1.5771778577960676]
    
    //porc 2
    //MRI Shrunk : best score 1ere etape : 0.149085
    //final parameters : [-0.006769331840125886, 0.09820386001984022, -0.01826162744483251, 2.7224875212753967, 0.7292929664137704, 3.0486718768199808]
    
    //full resolution MRI
//    best score 1ere etape : 0.156658
//    final parameters : [0.0036297748801233246, 0.013143829260491974, 0.0010727709811333918, 0.08008875517288504, -0.04607054405119398, 1.0090378953697132]
    
    //MRI rescaled
    
    //porc4
//    best score : 0.185428
//    final parameters : [-0.025424753944333096, 0.39962575090364716, 0.013544113893079517, -0.047444368532859084, -0.739287807247723, -0.5011674690804339]
    
    //porc6
//    best score : 0.150082
//    final parameters : [-0.0150802635531933, -0.010100620088642457, 0.004725328308930741, 8.30720384677868, 0.3776652489827957, 0.9237480980253407]


    //MRI Shrunk + mask
//    best score 1ere etape : 0.195807
//    final parameters : [-0.004576003822187878, -0.018662112915280166, 0.021674120903660443, 1.335316780176604, 0.012001995175122276, -0.08281897027573802]
    
    /**********
     *IDEAS TO TEST
     ************/

    //using mask of liver to limit ROI
    //try with arterial phased MRA (p6)
    //try playing on parameters to understand BOBYQA better, what about CMA-ES ?
    //GPU programming ?
  
    
    matrix<double> initialParameters (6,1);
    initialParameters(0) = 0;
    initialParameters(1) = 0;
    initialParameters(2) = 0;
    initialParameters(3) = 0;
    initialParameters(4) = 0;
    initialParameters(5) = 0;
    
    
    //cout<<"test params intial : "<<initialParameters<<endl;
    
    matrix<double> x_lower (6,1); //-0.5 rot, -10 trans
    x_lower(0) = -1;
    x_lower(1) = -1;
    x_lower(2) = -1;
    x_lower(3) = -1;
    x_lower(4) = -1;
    x_lower(5) = -1;
    
    matrix<double> x_upper (6,1); //0.5 rot, 10 trans
    x_upper(0) = 1;
    x_upper(1) = 1;
    x_upper(2) = 1;
    x_upper(3) = 1;
    x_upper(4) = 1;
    x_upper(5) = 1;
    
    double m = 13;
    
    //rho begin
    double radius = 0.9; //0.9 //zero makes no sense !!!
    
    //rho end
    double precision = 0.1; //0.01 0.001
    
    //n iter
    double nombreIteration = 200;

	//euler
	int type=1;

	//hemisphere central
	int coderegion=1;

	for (int i=0;i<8;i++){


				

		double rot_init=0.12/(1+i*0.35);
		double trans_init=3/(1+i*0.35);
		
		if(i>4)
			coderegion=3;
		
    
	LC2_function LC2 = LC2_function(US_shrunk,rescaled_IRM,outputPath,type,centre_spatial_us,coderegion);//MRI_shrunk or rescaled
    //make sure that the mask is computed on US_Shrunk but that we use the rescaled image to compute LC2
    //LC2.setMovingImage(rescaled_US);
    LC2.setMaxRot(rot_init); //0.3 for best initialisation
    LC2.setMaxTrans(trans_init);//10
    LC2.setRadius(radius);
    
    if(useLiverMask) LC2.setLiverMask(LiverMask);
    
	double best_score=1;
  // double best_score = find_min_bobyqa(LC2, initialParameters, m, x_lower, x_upper, LC2.getRadius(), precision, nombreIteration);
    //double best_score = find_min_bobyqa(LC2_function(rescaled_IRM, US_shrunk), initialParameters, m, x_lower, x_upper, radius, precision, nombreIteration);
    
    cout<<"best score 1ere etape : "<<best_score<<endl;
    //cout<<"rough parameters : "<<initialParameters<<endl;

//    EulerTransformType::ParametersType Step1Parameters(6);
//    
//    Step1Parameters[0] = initialParameters(0)*LC2.getMaxRot()/(LC2.getRadius());
//    Step1Parameters[1]= initialParameters(1)*LC2.getMaxRot()/(LC2.getRadius());
//    Step1Parameters[2] = initialParameters(2)*LC2.getMaxRot()/(LC2.getRadius());
//    Step1Parameters[3] = initialParameters(3)*LC2.getMaxTrans()/(LC2.getRadius());
//    Step1Parameters[4] = initialParameters(4)*LC2.getMaxTrans()/(LC2.getRadius());
//    Step1Parameters[5] = initialParameters(5)*LC2.getMaxTrans()/(LC2.getRadius());
//    cout<<"test best parameters 1ere etape: "<<Step1Parameters<<endl;
    
    
//************************************************************
//Euler 3D transform (Rigid)
//**********************************************


    //enregistrement resultats
    EulerTransformType::ParametersType finalParameters(6);
    finalParameters[0] = initialParameters(0)*LC2.getMaxRot()/(LC2.getRadius());
    finalParameters[1] = initialParameters(1)*LC2.getMaxRot()/(LC2.getRadius());
    finalParameters[2] = initialParameters(2)*LC2.getMaxRot()/(LC2.getRadius());
    finalParameters[3] = initialParameters(3)*LC2.getMaxTrans()/(LC2.getRadius());
    finalParameters[4] = initialParameters(4)*LC2.getMaxTrans()/(LC2.getRadius());
    finalParameters[5] = initialParameters(5)*LC2.getMaxTrans()/(LC2.getRadius());

	//on bloque ces operations geometriques pour la region 1
		if (coderegion==1){
		finalParameters[0]=0;
		finalParameters[4]=0;
		finalParameters[5]=0;
		}
	    
    
    //FINAL TSF
    
    EulerTransformType::Pointer finalTsf = EulerTransformType::New();
    finalTsf->SetParameters(finalParameters);
    
    cout<<"final parameters : "<<finalTsf->GetParameters()<<endl;
    
    //ecritire des parametres dans un fichier txt
    
    string outP =outputPath+"parameters.txt";
    ofstream fichier(outP.c_str(),ios::out | ios::trunc );
    
    if(fichier)
    {
        fichier<<"Parameters for rigid transform : "<<endl;
        fichier<<finalTsf->GetParameters()<<endl;
		fichier<<"6  initial Parameters for rigid transform : "<<endl;
		fichier<<initialParameters<<endl;
        fichier<<" Score for this position : "<<endl;
        fichier<<best_score<<endl;
        fichier.close();
    }
    
    else
    {
        cerr<<"Error in opening txt file for parameters"<<endl;
    }
    
    cout<<"Writing final registered US image"<<endl;
    
    ImageType::SizeType sizeUS2 = US_shrunk->GetLargestPossibleRegion().GetSize();
    ImageType::PointType origin2 = US_shrunk->GetOrigin();
    ImageType::SpacingType spacing2 = US_shrunk->GetSpacing();
    ImageType::PointType center2;
    //center2[0] = origin2[0]+spacing2[0]*sizeUS2[0]/2;
    //center2[1] = origin2[1]+spacing2[1]*sizeUS2[1]/2;
    //center2[2] = origin2[2]+spacing2[2]*sizeUS2[2]/2;
	
	
	//center2[0] = 0;
    //center2[1] = 0;
    //center2[2] = 0;

	//Patient 5
	//center2[0]=111;
	//center2[1]=80;
	//center2[2]=111;
    
    
    EulerTransformType::ParametersType eulerFixedParameters2(3);
    eulerFixedParameters2[0] =centre_spatial_us[0];
    eulerFixedParameters2[1] =centre_spatial_us[1];
    eulerFixedParameters2[2] =centre_spatial_us[2];
    
    finalTsf->SetFixedParameters(eulerFixedParameters2);
    
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
	resampler->SetInput(rescaled_IRM);
    resampler->SetTransform(finalTsf);
    resampler->SetSize(image_US->GetLargestPossibleRegion().GetSize());
    resampler->SetOutputSpacing(image_US->GetSpacing());
    resampler->SetOutputDirection(finalTsf->GetInverseMatrix()*image_US->GetDirection());
    resampler->SetOutputOrigin(finalTsf->GetInverseTransform()->TransformPoint(image_US->GetOrigin()));
	resampler->Update();
    //est ce que c'est une bonne idee de garder l'image full res alors qu'on fait le recalage base sur l'us downsampled ?
    
    cout<<"verification origine : "<<endl;
    cout<<"avant tsf : "<<image_US->GetOrigin()<<endl;
    cout<<"apres : "<<finalTsf->GetInverseTransform()->TransformPoint(image_US->GetOrigin());
    
    ImageType::Pointer finalImage = ImageType::New();
    finalImage = resampler->GetOutput();
	rescaled_IRM=resampler->GetOutput();
    cout<<"writing final result"<<endl;

	
	//integrer le numero diteration dans le nom du fichier
	string Result;          
	ostringstream convert;  
	convert << i;      
	Result = convert.str();
	string nomiteration="iterationx_finalregistreredUSBOBYQA.nii.gz";
	nomiteration.replace(9,1,Result);

	
    
    WriterType::Pointer writer8 = WriterType::New();
    string out8 =outputPath+nomiteration;

	cout<<"string fichier"<<out8<<endl;

    writer8->SetImageIO(io);
    writer8->SetInput(finalImage);
    writer8->SetFileName(out8);
    try {
        writer8->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error whilte writing registered image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
	int aaa=8;
	cout <<"Iteration "<<i<<endl;
	}

	//**********************
	//Affine transform
	//****************************
	
	//nombre de parametres (ndimensions=3) (ndim+1)*(ndim)
	
	matrix<double> initialParameters_affineT (15,1);


	//rotation
    initialParameters_affineT(0) = 0;
    initialParameters_affineT(1) = 0;
    initialParameters_affineT(2) = 0;
	//scale
    initialParameters_affineT(3) = 0;
    initialParameters_affineT(4) = 0;
    initialParameters_affineT(5) = 0;
	//shear
	initialParameters_affineT(6) = 0;
    initialParameters_affineT(7) = 0;
    initialParameters_affineT(8) = 0;
    initialParameters_affineT(9) = 0;
    initialParameters_affineT(10) = 0;
    initialParameters_affineT(11) = 0;
	//trans
    initialParameters_affineT(12) = 0;
    initialParameters_affineT(13) = 0;
    initialParameters_affineT(14) = 0;
    
    //cout<<"test params intial : "<<initialParameters<<endl;
    
    matrix<double> x_lower_affineT (15,1); //-0.5 rot, -10 trans
   x_lower_affineT(0) = -1;
    x_lower_affineT(1) = -1;
    x_lower_affineT(2) = -1;
    x_lower_affineT(3) = -1;
    x_lower_affineT(4) = -1;
    x_lower_affineT(5) = -1;
	 x_lower_affineT(6) = -1;
    x_lower_affineT(7) = -1;
    x_lower_affineT(8) = -1;
    x_lower_affineT(9) = -1;
    x_lower_affineT(10) = -1;
    x_lower_affineT(11) = -1;
	x_lower_affineT(12) = -1;
    x_lower_affineT(13) = -1;
    x_lower_affineT(14) = -1;
    
    matrix<double> x_upper_affineT (15,1); //0.5 rot, 10 trans
    x_upper_affineT(0) = 1;
    x_upper_affineT(1) = 1;
    x_upper_affineT(2) = 1;
    x_upper_affineT(3) = 1;
    x_upper_affineT(4) = 1;
    x_upper_affineT(5) = 1;
	x_upper_affineT(6) = 1;
    x_upper_affineT(7) = 1;
    x_upper_affineT(8) = 1;
    x_upper_affineT(9) = 1;
    x_upper_affineT(10) = 1;
    x_upper_affineT(11) = 1;
	x_upper_affineT(12) = 1;
    x_upper_affineT(13) = 1;
    x_upper_affineT(14) = 1;

	double m_af = 31;
    
    //rho begin
    double radius_af = 0.9; //0.9 //zero makes no sense !!!
    
    //rho end
    double precision_af = 0.1; //0.01 0.001
    
    //niter
    double nombreIteration_af = 200;
    
	cout<<"initialise LC2"<<endl;

	//ImageType::PointType origine;
	//origine[0]=0;
	//origine[1]=0;
	//origine[2]=0;
	//rescaled_IRM->SetOrigin(origine);
	//US_shrunk->SetOrigin(origine);
	
	int typeaf=2;
	 coderegion=2;

    LC2_function LC2_af = LC2_function(US_shrunk,rescaled_IRM,outputPath,typeaf,centre_spatial_us,coderegion);//MRI_shrunk or rescaled
    //make sure that the mask is computed on US_Shrunk but that we use the rescaled image to compute LC2
    //LC2.setMovingImage(rescaled_US);
	double ini_rot_af=0.08;
	double ini_trans_af=3;
	double ini_scale_af=1.07;
	double ini_shear_af=0.15;

    LC2_af.setMaxRot(ini_rot_af); //0.3 for best initialisation
    LC2_af.setMaxTrans(ini_trans_af);//10
	LC2_af.setMaxScale(ini_scale_af); //valeur arbitraire a verifier
	//plus petit shear
	LC2_af.setMaxShear(ini_shear_af);
    LC2_af.setRadius(radius_af);
    
    if(useLiverMask) LC2_af.setLiverMask(LiverMask);

	cout<<"calcul LC2"<<endl;

	//Optimisation des paramètres avec Bobyqa
	double best_score_af = find_min_bobyqa(LC2_af, initialParameters_affineT, m_af, x_lower_affineT, x_upper_affineT, LC2_af.getRadius(), precision_af, nombreIteration_af);
    //double best_score = find_min_bobyqa(LC2_function(rescaled_IRM, US_shrunk), initialParameters, m, x_lower, x_upper, radius, precision, nombreIteration);
    
    cout<<"best score 1ere etape : "<<best_score_af<<endl;

	
	
	//FINAL TSF
    
    AffineTransformType::Pointer finalAffTsf = AffineTransformType::New();
    ///Rotation

		//Définis l'axe de rotation avec un vecteur
		AffineTransformType::OutputVectorType axederotation_x;
		axederotation_x[0]=1;
		axederotation_x[1]=0;
		axederotation_x[2]=0;
		AffineTransformType::OutputVectorType axederotation_y;
		axederotation_y[0]=0;
		axederotation_y[1]=1;
		axederotation_y[2]=0;
		AffineTransformType::OutputVectorType axederotation_z;
		axederotation_z[0]=0;
		axederotation_z[1]=0;
		axederotation_z[2]=1;
		
		
		///3er parametre sur les angles de rotation 1 pour chaque axe de rotation
		//double anglerotation=0.3;
		double anglerotation_x=initialParameters_affineT(0)*LC2_af.getMaxRot()/(LC2_af.getRadius());
		double anglerotation_y=initialParameters_affineT(1)*LC2_af.getMaxRot()/(LC2_af.getRadius());
		double anglerotation_z=initialParameters_affineT(2)*LC2_af.getMaxRot()/(LC2_af.getRadius());

	
	
		std::cout<<"Affine tsf parameters : "<<finalAffTsf->GetParameters()<<std::endl;

		finalAffTsf->Rotate3D(axederotation_x,anglerotation_x,true);
		finalAffTsf->Rotate3D(axederotation_y,anglerotation_y,true);
		finalAffTsf->Rotate3D(axederotation_z,anglerotation_z,true);

       // std::cout<<"Affine tsf parameters : "<<transformaff->GetParameters()<<std::endl;

		//essayer de scaler isotropique aussi!!!!!!

	//Scale
		AffineTransformType::OutputVectorType scalefactor;
		scalefactor[0]=1+initialParameters_affineT(3)*(LC2_af.getMaxScale()-1)/(LC2_af.getRadius());
		scalefactor[1]=1+initialParameters_affineT(4)*(LC2_af.getMaxScale()-1)/(LC2_af.getRadius());
		scalefactor[2]=1+initialParameters_affineT(5)*(LC2_af.getMaxScale()-1)/(LC2_af.getRadius());

		//finalAffTsf->Scale(scalefactor,true);

		
		
		 std::cout<<"Affine tsf parameters : "<<finalAffTsf->GetParameters()<<std::endl;
		 
	//Shear
		


		 double coefficient[6];
		coefficient[0]=initialParameters_affineT(6)*LC2_af.getMaxShear()/(LC2_af.getRadius());
		coefficient[1]=initialParameters_affineT(7)*LC2_af.getMaxShear()/(LC2_af.getRadius());
		coefficient[2]=initialParameters_affineT(8)*LC2_af.getMaxShear()/(LC2_af.getRadius());
		coefficient[3]=initialParameters_affineT(9)*LC2_af.getMaxShear()/(LC2_af.getRadius());
		coefficient[4]=initialParameters_affineT(10)*LC2_af.getMaxShear()/(LC2_af.getRadius());
		coefficient[5]=initialParameters_affineT(11)*LC2_af.getMaxShear()/(LC2_af.getRadius());

		//a verifier peut besoin de 6 parametre pcq inverser 0 et 1
		 finalAffTsf->Shear(0,1,coefficient[0]);
		 finalAffTsf->Shear(1,0,coefficient[1]);
		 finalAffTsf->Shear(2,0,coefficient[2]);
		finalAffTsf->Shear(0,2,coefficient[3]);
		finalAffTsf->Shear(1,2,coefficient[4]);
		finalAffTsf->Shear(2,1,coefficient[5]);
		 
		

		 std::cout<<"Affine tsf parameters : "<<finalAffTsf->GetParameters()<<std::endl;
		 
	//translate
		 AffineTransformType::OutputVectorType vecteurtranslation_x;
		 AffineTransformType::OutputVectorType vecteurtranslation_y;
		 AffineTransformType::OutputVectorType vecteurtranslation_z;
		 vecteurtranslation_x[0]=1;
		 vecteurtranslation_x[1]=0;
		 vecteurtranslation_x[2]=0;
		 vecteurtranslation_y[0]=0;
		 vecteurtranslation_y[1]=1;
		 vecteurtranslation_y[2]=0;
		 vecteurtranslation_z[0]=0;
		 vecteurtranslation_z[1]=0;
		 vecteurtranslation_z[2]=1;

		 vecteurtranslation_x=vecteurtranslation_x*initialParameters_affineT(12)*LC2_af.getMaxTrans()/(LC2_af.getRadius());
		 vecteurtranslation_y=vecteurtranslation_y*initialParameters_affineT(13)*LC2_af.getMaxTrans()/(LC2_af.getRadius());
		 vecteurtranslation_z=vecteurtranslation_z*initialParameters_affineT(14)*LC2_af.getMaxTrans()/(LC2_af.getRadius());

		 finalAffTsf->Translate(vecteurtranslation_x,true);
		 finalAffTsf->Translate(vecteurtranslation_y,true);
		 finalAffTsf->Translate(vecteurtranslation_z,true);


		  //Affine transform avec centre pré déterminé

		AffineTransformfixedcenter::Pointer AffineTransformfixefinal=AffineTransformfixedcenter::New();

		AffineTransformfixefinal->Rotate3D(axederotation_x,anglerotation_x,true);
		AffineTransformfixefinal->Rotate3D(axederotation_y,anglerotation_y,true);
		AffineTransformfixefinal->Rotate3D(axederotation_z,anglerotation_z,true);

		AffineTransformfixefinal->Scale(scalefactor,true);

	
		 AffineTransformfixefinal->Shear(0,1,coefficient[0]);
		 AffineTransformfixefinal->Shear(1,0,coefficient[1]);
		 AffineTransformfixefinal->Shear(2,0,coefficient[2]);
		 AffineTransformfixefinal->Shear(0,2,coefficient[3]);
		 AffineTransformfixefinal->Shear(1,2,coefficient[4]);
		 AffineTransformfixefinal->Shear(2,1,coefficient[5]);

		  AffineTransformfixefinal->Translate(vecteurtranslation_x,true);
		 AffineTransformfixefinal->Translate(vecteurtranslation_y,true);
		 AffineTransformfixefinal->Translate(vecteurtranslation_z,true);
     
       // ImageType::PointType center_af;
       
		//Patient 2
       // center[0] =98;
        //center[1] =68;
        //center[2] =98;

		//Patient 5
		//center_af[0]=111;
		//center_af[1]=80;
		//center_af[2]=111;

		cout<<"center of rotation"<<centre_spatial_us<<endl;

		AffineTransformfixefinal->SetCenterOfRotationComponent(centre_spatial_us);
    
    cout<<"final parameters : "<<initialParameters_affineT<<endl;
    
    //ecritire des parametres dans un fichier txt
    
    string outP2 =outputPath+"parametersaffine.txt";
    ofstream fichier2(outP2.c_str(),ios::out | ios::trunc );
    
    if(fichier2)
    {
        fichier2<<"Matrix parameters for Affine transform  : "<<endl;
        fichier2<<finalAffTsf->GetParameters()<<endl;
		fichier2<<"Transform optimize parameters for Affine transform  : "<<endl;
        fichier2<<initialParameters_affineT<<endl;
        fichier2<<" Score for this position : "<<endl;
        fichier2<<best_score_af<<endl;
		fichier2<<"Initialisation value: "<<endl;
		fichier2<<"Rotation: "<<ini_rot_af<<"Translation: "<<ini_trans_af<<"Scale: "<<ini_scale_af<<"Shear: "<<ini_shear_af<<endl;
        fichier2.close();
    }
    
    else
    {
        cerr<<"Error in opening txt file for parameters"<<endl;
    }
    
    cout<<"Writing final registered US image"<<endl;
    
    
    
    ResampleFilterType::Pointer resampler_af = ResampleFilterType::New();
	resampler_af->SetInput(rescaled_IRM);
	resampler_af->SetTransform(AffineTransformfixefinal);
    //resampler->SetTransform(finalAffTsf);
    //resampler->SetSize(image_US->GetLargestPossibleRegion().GetSize());
    //resampler->SetOutputSpacing(image_US->GetSpacing());
   // resampler->SetOutputDirection(finalAffTsf->GetInverseMatrix()*image_US->GetDirection());
   // resampler->SetOutputOrigin(finalAffTsf->GetInverseTransform()->TransformPoint(image_US->GetOrigin()));

	resampler_af->SetSize(image_US->GetLargestPossibleRegion().GetSize());
    resampler_af->SetOutputSpacing(image_US->GetSpacing());
    resampler_af->SetOutputDirection(image_US->GetDirection());
    resampler_af->SetOutputOrigin(image_US->GetOrigin());
    //est ce que c'est une bonne idee de garder l'image full res alors qu'on fait le recalage base sur l'us downsampled ?
    try {
            resampler_af->Update();
        } catch (itk::ExceptionObject &e) {
            cerr<<"error while transforming mesh"<<endl;
            cerr<<e<<endl;
    }

    cout<<"verification origine : "<<endl;
    cout<<"avant tsf : "<<image_IRM->GetOrigin()<<endl;
    cout<<"apres : "<<finalAffTsf->GetInverseTransform()->TransformPoint(image_IRM->GetOrigin());
	cout<<"apres IRM: "<<image_US->GetOrigin();
    
    ImageType::Pointer finalImage_af = ImageType::New();
    finalImage_af = resampler_af->GetOutput();
    
    cout<<"writing final result"<<endl;
    
    WriterType::Pointer writer9 = WriterType::New();
    string out9 =outputPath+"finalregistreredAffineUSBOBYQA.nii.gz";
    writer9->SetImageIO(io);
    writer9->SetInput(finalImage_af);
    writer9->SetFileName(out9);
    try {
        writer9->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error whilte writing registered image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
	
	
//    
//    //target points evaluation
//    ImageType::PointType target1;
//    target1[0]=-49.000;
//    target1[1]=-23.689;
//    target1[2]=-47.019;
//    
//    ImageType::PointType target2;
//    target2[0]=-22.620;
//    target2[1]=-23.424;
//    target2[2]=-11.070;
//    
//    ImageType::PointType target3;
//    target3[0]=-44.518;
//    target3[1]=-21.495;
//    target3[2]=-50.116;
//    
//    ImageType::PointType target4;
//    target4[0]=-31.493;
//    target4[1]=1.394;
//    target4[2]=-16.948;
//    
//    ImageType::PointType target5;
//    target5[0]=-41.194;
//    target5[1]=-1.120;
//    target5[2]=-12.452;
//    
//    //transformation des targets points
//    
//    ImageType::PointType resTarget1 = finalTsf->GetInverseTransform()->TransformPoint(target1);
//    ImageType::PointType resTarget2 = finalTsf->GetInverseTransform()->TransformPoint(target2);
//    ImageType::PointType resTarget3 = finalTsf->GetInverseTransform()->TransformPoint(target3);
//    ImageType::PointType resTarget4 = finalTsf->GetInverseTransform()->TransformPoint(target4);
//    ImageType::PointType resTarget5 = finalTsf->GetInverseTransform()->TransformPoint(target5);
//    
//    //enregistrement dans fichier txt
//    
//    //ecritire des parametres dans un fichier txt
//    
//    string outT =outputPath+"/results_targets.txt";
//    ofstream fichier2(outT.c_str(),ios::out | ios::trunc );
//    
//    if(fichier)
//    {
//        fichier2<<"Transformed target points : "<<endl;
//        fichier2<<"resTarget 1 : "<<resTarget1<<endl;
//        fichier2<<"resTarget 2 : "<<resTarget2<<endl;
//        fichier2<<"resTarget 3 : "<<resTarget3<<endl;
//        fichier2<<"resTarget 4 : "<<resTarget4<<endl;
//        fichier2<<"resTarget 5 : "<<resTarget5<<endl;
//        fichier2.close();
//    }
//    
//    else
//    {
//        cerr<<"Error in opening txt file for parameters"<<endl;
//    }
    
    /*
    ////////////////
    //1. TRANSFORM /
    ////////////////
    
    cout<<"creation tsf"<<endl;
    
    //euler tsf
    EulerTransformType::Pointer transform = EulerTransformType::New();
    transform->SetIdentity();
    
    
    /////////////////
    // 2. OPTIMIZER /
    /////////////////
    

    //Amoeba optimizer
    
    cout<<"creation optimizer"<<endl;
    
    AmoebaOptimizerType::Pointer optimizer = AmoebaOptimizerType::New();
    
    
    //interpolateur et metrique
    
    /////////////////////
    // 3. INTERPOLATOR //
    /////////////////////
    
    cout<<"creation interpolateur"<<endl;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    
    
    /////////////////
    //  4. REGISTOR /
    ////////////////
    //pipeline de recalage
    
    RegistrationType::Pointer registration = RegistrationType::New();
    
    //settings des blocs constitutifs
   
    cout<<"settings des blocs constitutifs"<<endl;
  
    registration->SetOptimizer(optimizer);
    registration->SetTransform(transform);
    registration->SetInterpolator(interpolator);

//
//    ///////////////
//    // 4. METRIC //
//    ///////////////
    cout<<"setting metrique"<<endl;
    LC2MetricType::Pointer metric = LC2MetricType::New();

    
  //Amoeba optimizer
    
    registration->SetMetric(metric);
    
    
    //setting des images
    cout<<"setting images"<<endl;
    registration->SetFixedImage(rescaled_IRM);
    registration->SetMovingImage(US_shrunk);
    //just to be sure for the metric
    metric->SetFixed(rescaled_IRM);
    metric->ComputeGradImage();
    //metric->ComputeVesselnessImage();
    metric->SetMoving(US_shrunk);
    metric->ComputeMask();
    if(useLiverMask) metric->setLiverMask(LiverMask);
    
    registration->SetFixedImageRegion(rescaled_IRM->GetBufferedRegion());
    
    
    RegistrationType::ParametersType initialParameters = transform->GetParameters();
    //setting des parametres optimizable
    initialParameters[0] = 0;
    initialParameters[1] = 0;
    initialParameters[2] = 0;
    initialParameters[3] = 0; //euler tsf = 6 parameters
    initialParameters[4] = 0;
    initialParameters[5] = 0;
    
    registration->SetInitialTransformParameters(initialParameters);
    
    cout<<"Initial transform param = "<<registration->GetTransform()->GetParameters()<<endl;
    
    /************
     * OPTIMIZER
     **************/
    /*
    //setting des params de l'optimizer
    const unsigned int numberOfParameters = transform->GetNumberOfParameters();
    AmoebaOptimizerType::ParametersType simplexDelta(numberOfParameters);
    simplexDelta[0] =0.3;
    simplexDelta[1] =0.3;
    simplexDelta[2] =0.3;
    simplexDelta[3] = 10;
    simplexDelta[4] = 10;
    simplexDelta[5] = 10;
    
    cout<<"verification simplex delta structure : "<<endl;
    cout<<simplexDelta<<endl;
    
    optimizer->AutomaticInitialSimplexOff();
    optimizer->SetInitialSimplexDelta(simplexDelta);
    
    //tolerance setting on param units and LC2 units
    optimizer->SetParametersConvergenceTolerance(0.1);
    optimizer->SetFunctionConvergenceTolerance(0.001);
    
    //max number of iterations to stop if no convergence
    optimizer->SetMaximumNumberOfIterations(3);
    
    //on le met en maximisation
    optimizer->SetMaximize(true);
    
    //Command observer
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    optimizer->AddObserver(itk::IterationEvent(), observer);
    

    cout<<"registration"<<endl;
    
        try {
            registration->Update();
            cout<<"Optimizer stop condition : "<<registration->GetOptimizer()->GetStopConditionDescription()<<endl;
        } catch (itk::ExceptionObject &e) {
            cerr<<"erreur dans le recalage"<<endl;
            cerr<<e<<endl;
            return EXIT_FAILURE;
        }
    

    RegistrationType::ParametersType finalParameters = registration->GetLastTransformParameters();
    
    double bestValue = optimizer->GetValue();
    
    //Print out results
    
    cout<<"Results : "<<endl;
    //cout<<"Translation vector : "<<"[ "<<TX<<", "<<TY<<", "<<TZ<<" ]"<<endl;
    cout<<"parametres : "<<finalParameters<<endl;
    cout<<"metric value : "<<bestValue<<endl;
    
    EulerTransformType::Pointer finalTsf = EulerTransformType::New();
    finalTsf->SetParameters(finalParameters);
    
    cout<<"final parameters : "<<finalTsf->GetParameters()<<endl;
    
    //ecritire des parametres dans un fichier txt
    
    string outP =outputPath+"\parameters.txt";
    ofstream fichier(outP.c_str(),ios::out | ios::trunc );
    
    if(fichier)
    {
        fichier<<"Parameters for rigid transform : "<<endl;
        fichier<<finalTsf->GetParameters()<<endl;
        fichier<<" Score for this position : "<<endl;
        fichier<<bestValue<<endl;
        fichier.close();
    }
    
    else
    {
        cerr<<"Error in opening txt file for parameters"<<endl;
    }


    cout<<"Writing final registered US image"<<endl;
    
    ImageType::SizeType sizeUS2 = US_shrunk->GetLargestPossibleRegion().GetSize();
    ImageType::PointType origin2 = US_shrunk->GetOrigin();
    ImageType::SpacingType spacing2 = US_shrunk->GetSpacing();
    ImageType::PointType center2;
    center2[0] = origin2[0]+spacing2[0]*sizeUS2[0]/2;
    center2[1] = origin2[1]+spacing2[1]*sizeUS2[1]/2;
    center2[2] = origin2[2]+spacing2[2]*sizeUS2[2]/2;
    
    
    EulerTransformType::ParametersType eulerFixedParameters2(3);
    eulerFixedParameters2[0] =center2[0];
    eulerFixedParameters2[1] =center2[1];
    eulerFixedParameters2[2] =center2[2];
    
    finalTsf->SetFixedParameters(eulerFixedParameters2);
    
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetInput(image_US);
    resampler->SetTransform(finalTsf);
    resampler->SetSize(image_US->GetLargestPossibleRegion().GetSize());
    resampler->SetOutputSpacing(image_US->GetSpacing());
    resampler->SetOutputDirection(finalTsf->GetInverseMatrix()*image_US->GetDirection());
    resampler->SetOutputOrigin(finalTsf->GetInverseTransform()->TransformPoint(image_US->GetOrigin()));
    //est ce que c'est une bonne idee de garder l'image full res alors qu'on fait le recalage base sur l'us downsampled ?
    
    ImageType::Pointer finalImage = ImageType::New();
    finalImage = resampler->GetOutput();
    
    cout<<"writing final result"<<endl;
    
    WriterType::Pointer writer8 = WriterType::New();
	string out8 =outputPath+"\finalregistreredUSBOBYQA.nii";
    writer8->SetImageIO(io);
    writer8->SetInput(finalImage);
    writer8->SetFileName(out8);
    try {
        writer8->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error whilte writing registered image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    
    /***********************
     * TEST TSF DEFORMABLE
     **********************/
    /*
    BSplineTransformType::Pointer Btransform = BSplineTransformType::New();
    //space dimension est a 3
	//spline order est a 3 (peut le changer les 2)
    int SpaceDimension=3;
	int SplineOrder=3;

	unsigned int numberOfGridNodesInOneDimension = 7;
	BSplineTransformType::PhysicalDimensionsType   fixedPhysicalDimensions;
  BSplineTransformType::MeshSizeType             meshSize;
  BSplineTransformType::OriginType               fixedOrigin;

  for( unsigned int i=0; i< SpaceDimension; i++ )
    {
    fixedOrigin[i] = rescaled_IRM->GetOrigin()[i];
    fixedPhysicalDimensions[i] = rescaled_IRM->GetSpacing()[i] * static_cast<double>(rescaled_IRM->GetLargestPossibleRegion().GetSize()[i] - 1 );
    }
  meshSize.Fill( numberOfGridNodesInOneDimension - SplineOrder );
	
	
  Btransform->SetTransformDomainOrigin( fixedOrigin );
  Btransform->SetTransformDomainPhysicalDimensions(
    fixedPhysicalDimensions );
  Btransform->SetTransformDomainMeshSize( meshSize );
  Btransform->SetTransformDomainDirection(rescaled_IRM->GetDirection() );
  
  const unsigned int numberOfParameters = Btransform->GetNumberOfParameters();
  cout<<"number of parametre B SPline"<<numberOfParameters<<endl;
  BSParametersType parameters( numberOfParameters );
  parameters.Fill( 0.0 );
  */
  
  //definition de la region de deformation
    



    
    
    
    
    std::cout << "Done running routine!\n";
  
    tend = std::time(NULL);
    texec = std::difftime(tend,tbegin);
    
    cout<<"temps en s : "<<texec<<endl;
    return 0;
}
