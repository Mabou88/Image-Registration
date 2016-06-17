#include "AffineRegistration.h"

//typedef itk::Image<double,3> ImageType;

void AffineRegistration(ImageType::Pointer itkimageirm ,ImageType::Pointer itkimageus){


	/*
//Convertit l'image IRM de mitk en image itk
	mitk::ImageToItk<ImageType>::Pointer toItkFilter =mitk::ImageToItk<ImageType>::New();
    toItkFilter->SetInput(imagecloneirm);
    toItkFilter->Update();
	 ImageType::Pointer itkimageirm = toItkFilter->GetOutput();

	 mitk::ImageToItk<ImageType>::Pointer toItkFilter2 =mitk::ImageToItk<ImageType>::New();
    toItkFilter2->SetInput(imagecloneus);
    toItkFilter2->Update();
	 ImageType::Pointer itkimageus = toItkFilter2->GetOutput();
	 */


  //transformation
	 typedef itk::AffineTransform< double,3> TransformType;
  //optimisateur
  typedef itk::RegularStepGradientDescentOptimizer  OptimizerType;

  typedef itk::QuaternionRigidTransformGradientDescentOptimizer  OptimizerType2;
  //metric
   typedef itk::MeanSquaresImageToImageMetric< ImageType, ImageType >   MetricType;
   //interpolateur
    typedef itk:: LinearInterpolateImageFunction< ImageType,  double>    InterpolatorType;
	//registrationmethod
	typedef itk::ImageRegistrationMethod< ImageType, ImageType >    RegistrationType;
	
	// Create components
  MetricType::Pointer         metric        = MetricType::New();
  TransformType::Pointer      transform     =TransformType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();
  
  // Each component is now connected to the instance of the registration method.
  registration->SetMetric ( metric);
  registration->SetOptimizer(optimizer);
  registration->SetTransform( transform );
  registration->SetInterpolator(interpolator);
  
  // Write the two synthetic inputs
  typedef itk::ImageFileWriter< ImageType >  WriterType;
 
  WriterType::Pointer      fixedWriter =  WriterType::New();
  fixedWriter->SetFileName("fixed.nii");
  fixedWriter->SetInput( itkimageirm);
  fixedWriter->Update();
 //verifier  les setinput
  WriterType::Pointer      movingWriter =  WriterType::New();
  movingWriter->SetFileName("moving.nii");
  movingWriter->SetInput(itkimageus);
  movingWriter->Update();
  
  
 
  // Set the registration inputs
  registration->SetFixedImage(itkimageirm);
  registration->SetMovingImage(itkimageus);
  
  //creation d'une region pour appliquer la metric
  typedef itk::ImageRegion<3> RegionType;
  itk::Index<3> index={{ 60,80,60}};
  itk::Size<3> size = {50, 50, 40};
  RegionType region(index,size);

  //itkimageirm->SetRequestedRegion(region);
  //RegionType roi;
  //roi=itkimageirm->GetRequestedRegion();
  metric->SetFixedImageRegion(region);
  //registration->SetFixedImageRegion(region);
  //registration->SetFixedImageRegionDefined(true);
  // a verifier
  //registration->SetFixedImageRegion(itkimageirm->GetLargestPossibleRegion());

  //initialise le centre de la rotation pour la transformation
  typedef itk::CenteredTransformInitializer< TransformType, ImageType,ImageType >  TransformInitializerType;
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();
   initializer->SetTransform(   transform );
  initializer->SetFixedImage(itkimageirm );
  initializer->SetMovingImage(itkimageus);
  initializer->MomentsOn();
  initializer->InitializeTransform();

  //  Now we pass the parameters of the current transform as the initial
  //  parameters to be used when the registration process starts.
  registration->SetInitialTransformParameters( transform->GetParameters());





	  //  Initialize the transform
  double translationScale = 1.0 / 1000.0;
  typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
  optimizerScales[0] =  1.0;
  optimizerScales[1] =  1.0;
  optimizerScales[2] =  1.0;
  optimizerScales[3] =  1.0;
  optimizerScales[4] =  1.0;
  optimizerScales[5] =  1.0;
  optimizerScales[6] =  1.0;
  optimizerScales[7] =  1.0;
  optimizerScales[8] =  1.0;
  optimizerScales[9]  =  translationScale;
  optimizerScales[10] =  translationScale;
  optimizerScales[11] =  translationScale;
  optimizer->SetScales( optimizerScales );

  int nombre=transform->GetNumberOfParameters();
  cout<<"parametre"<<endl;
  cout<<nombre<<endl;
  
  

  //registration->SetInitialTransformParameters( initialParameters );

  optimizer->SetMaximumStepLength( 0.1 );
  optimizer->SetMinimumStepLength( 0.001 );
 
  // Set a stopping criterion
  optimizer->SetNumberOfIterations(2);

  try
  {
    registration->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
   // return EXIT_FAILURE;
  }
 
  //  The result of the registration process is an array of parameters that
  //  defines the spatial transformation in an unique way. This final result is
  //  obtained using the \code{GetLastTransformParameters()} method.
 
  OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();

   
  //  The optimizer can be queried for the actual number of iterations
  //  performed to reach convergence.  The \code{GetCurrentIteration()}
  //  method returns this value. A large number of iterations may be an
  //  indication that the maximum step length has been set too small, which
  //  is undesirable since it results in long computational times.
 
  const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
 
  //  The value of the image metric corresponding to the last set of parameters
  //  can be obtained with the \code{GetValue()} method of the optimizer.
 
  const double bestValue = optimizer->GetValue();
 
  // Print out results
  //
  cout << "Result = " << std::endl;
  cout << " Iterations    = " << numberOfIterations << endl;
  cout << " Metric value  = " << bestValue          << endl;
  
  //Crée la transformation avec les bon paramètres pour le recalage
  TransformType::Pointer finalTransform = TransformType::New();
   finalTransform->SetParameters( finalParameters );
  finalTransform->SetFixedParameters( transform->GetFixedParameters() );

  //Crée un resample filter pour faire la transformation
  typedef itk::ResampleImageFilter<ImageType,ImageType >    ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetTransform( finalTransform );
  resampler->SetInput( itkimageus);
  ImageType::Pointer fixedImage = itkimageirm;
  resampler->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
  resampler->SetOutputSpacing( fixedImage->GetSpacing() );
  resampler->SetOutputDirection( fixedImage->GetDirection() );
  resampler->SetDefaultPixelValue( 0 );

  //ecris le resultats dans un fichier results.nii
   typedef  double  OutputPixelType;
  typedef itk::Image< OutputPixelType,3> OutputImageType;
  typedef itk::CastImageFilter<ImageType, OutputImageType > CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();
  itk::NiftiImageIO::Pointer ioimagenifti=itk::NiftiImageIO::New();

  writer->SetImageIO(ioimagenifti);
  writer->SetFileName( "C:/im/results_roi5.nii");
  caster->SetInput( resampler->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();

}