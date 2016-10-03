#include "Filtre.h"

Filtre::Filtre(ImageType::Pointer imagetofiltre){
	
	//copy image
	imagefiltre=imagetofiltre;
	kernel=ImageType::New();

}
void Filtre::fft_transform(){
	//version 1
	
	//FFTFilterType::Pointer fftinput = FFTFilterType::New();
	//fftinput->SetInput(imagefiltre);
	//fftinput->Update();

	//imagefiltreFFT=fftinput->GetOutput();
	

	//version 2
	//VnlFFTFilterType::Pointer fftfiltre=VnlFFTFilterType::New();
	//fftfiltre->SetInput(imagefiltre);
	//fftfiltre->Update();

}

void Filtre::selectROI(){

	 int indMaxX =0;
     int indMinX = 1000;
     int indMaxY =0;
     int indMinY = 1000;
     int indMaxZ =0;
     int indMinZ = 1000;

	 ImageConstIteratorType us_it(imagefiltre,imagefiltre->GetLargestPossibleRegion());
      us_it.GoToBegin();
	int iterateurparcours=0;
      while(!us_it.IsAtEnd())
        {
			iterateurparcours+=1;
            if(us_it.Get()>0)
            {
                ImageType::IndexType ind = us_it.GetIndex();
                
                //X
                if(ind[0]>indMaxX)
                {
                    indMaxX=ind[0];
                }
                
                else if(ind[0]<indMinX)
                {
                    indMinX=ind[0];
                }
                
                //Y
                if(ind[1]>indMaxY)
                {
                    indMaxY=ind[1];
                }
                
                else if(ind[1]<indMinY)
                {
                    indMinY=ind[1];
                }
                
                //Z
                if(ind[2]>indMaxZ)
                {
                    indMaxZ=ind[2];
                }
                
                else if(ind[2]<indMinZ)
                {
                    indMinZ=ind[2];
                }
            }
            
            ++us_it;
            
        }

	  ImageType::IndexType startCropped;
        startCropped[0] = indMinX ;
        startCropped[1] = indMinY;
        startCropped[2] = indMinZ;
        
        
        
        ImageType::IndexType endCropped;
        endCropped[0] = indMaxX;
        endCropped[1] = indMaxY;
        endCropped[2] = indMaxZ;
        
        std::cout<<"indices minimum X,Y,Z "<<std::endl;
        std::cout<<"X : "<<indMinX<<" "<<indMaxX<<std::endl;
        std::cout<<"Y : "<<indMinY<<" "<<indMaxY<<std::endl;
        std::cout<<"Z : "<<indMinZ<<" "<<indMaxZ<<std::endl;
        
        ImageType::SizeType sizeCropped;
        sizeCropped[0] = endCropped[0]-startCropped[0]+1;
        sizeCropped[1] = endCropped[1]-startCropped[1]+1;
        sizeCropped[2] = endCropped[2]-startCropped[2]+1;
        
        ImageType::RegionType regionCropped;
        regionCropped.SetIndex(startCropped);
        regionCropped.SetSize(sizeCropped);

		//Cropping de l'image originale
        
        CropFilter::Pointer CroppingFilter = CropFilter::New();
        CroppingFilter->SetExtractionRegion(regionCropped);
        CroppingFilter->SetInput(imagefiltre);
        CroppingFilter->SetDirectionCollapseToIdentity();
        CroppingFilter->Update();
        ImageType::Pointer imagefiltre2 = CroppingFilter->GetOutput();
        std::cout<<"verification image size : "<<imagefiltre2->GetLargestPossibleRegion().GetSize()<<std::endl;


}

void Filtre::filtregaussian(){

	FilterType::Pointer convolutionFilter = FilterType::New();
	selectROI();
	convolutionFilter->SetInput(imagefiltre);
	convolutionFilter->SetKernelImage(kernel);
	convolutionFilter->NormalizeOn();
	convolutionFilter->Update();
	ImageType::Pointer filtre0=convolutionFilter->GetOutput();

	

	//AddFilter::Pointer additionfilter=AddFilter::New();
	//additionfilter->SetInput1(imagefiltre);
	//additionfilter->SetInput2(filtre);
	//additionfilter->Update();
	//imageresult=additionfilter->GetOutput();

	//GaussianType::Pointer gausfiltre=GaussianType::New();
	//gausfiltre->SetInput(imagefiltre);
	//gausfiltre->SetVariance(1);
	//gausfiltre->SetMaximumKernelWidth(1);
	//gausfiltre->Update();
	//imageresult=gausfiltre->GetOutput();

	//SubtractFilter::Pointer subtractionfilter=SubtractFilter::New();
	//subtractionfilter->SetInput1(imagefiltre);
	//subtractionfilter->SetInput2(imageresult);
	//subtractionfilter->Update();
	//filtre=subtractionfilter->GetOutput();

	//AddFilter::Pointer additionfilter=AddFilter::New();
	//additionfilter->SetInput1(imagefiltre);
	//additionfilter->SetInput2(filtre);
	//additionfilter->Update();
	//imageresult=additionfilter->GetOutput();

	//RescaleFilterType::Pointer rescaleimage=RescaleFilterType::New();
	//rescaleimage->SetInput(imageresult);
	//rescaleimage->SetOutputMinimum(0);
	//rescaleimage->SetOutputMaximum(255);
	//rescaleimage->Update();
	//imageresult=rescaleimage->GetOutput();

}

void Filtre::filtrepassehaut(){

	//lissage de l'image
	GaussianType2::Pointer gausfiltre=GaussianType2::New();
	gausfiltre->SetInput(imagefiltre);
	gausfiltre->SetVariance(2);
	gausfiltre->SetMaximumKernelWidth(7);
	gausfiltre->Update();
	ImageType::Pointer imagegaus=gausfiltre->GetOutput();

	GaussianType2::Pointer gausfiltre2=GaussianType2::New();
	gausfiltre2->SetInput(imagefiltre);
	gausfiltre2->SetVariance(1);
	gausfiltre2->SetMaximumKernelWidth(1);
	gausfiltre2->Update();
	ImageType::Pointer imagegaus2=gausfiltre->GetOutput();

	 

  //gradient de l'image
	Gradient imagegrad;
	imagegrad.compute_gradient(imagegaus);
	GradientMFilterType::Pointer m_gradientMagnitudeFilter;
    VectorGradientFilterType::Pointer m_gradientMapFilter;
	ImageType::Pointer imagegradient=imagegrad.getImageGradient();

	WriterType::Pointer   writer0 =  WriterType::New();
  itk::NiftiImageIO::Pointer ioimagenifti0=itk::NiftiImageIO::New();

  writer0->SetImageIO(ioimagenifti0);
  writer0->SetFileName( "C:/Users/Marc-Antoine/Documents/Imagecode/output/gradient.nii");
  writer0->SetInput(imagegradient);
  writer0->Update();
	
	//accentue les contours
	MultiplyFilter::Pointer multiplyfilter=MultiplyFilter::New();
	multiplyfilter->SetInput(imagegradient);
	multiplyfilter->SetConstant(1.5);
	multiplyfilter->Update();
	imagegradient=multiplyfilter->GetOutput();
	
	//additionne le gradient sur l'image
	AddFilter::Pointer additionfilter=AddFilter::New();
	additionfilter->SetInput1(imagefiltre);
	additionfilter->SetInput2(imagegradient);
	additionfilter->Update();
	imageresult=additionfilter->GetOutput();
	
	//soustrait constante
	SubtractFilter::Pointer subtractionfilter=SubtractFilter::New();
	subtractionfilter->SetInput1(imageresult);
	subtractionfilter->SetConstant2(20);
	subtractionfilter->Update();
	imageresult=subtractionfilter->GetOutput();
	
	/*
	//ImageType::Pointer vectorgradient=imagegrad.getImageVectorGradient();
	 WriterType::Pointer   writer =  WriterType::New();
  itk::NiftiImageIO::Pointer ioimagenifti=itk::NiftiImageIO::New();

  writer->SetImageIO(ioimagenifti);
  writer->SetFileName( "C:/Users/Marc-Antoine/Documents/Imagecode/output/gradient.nii");
  writer->SetInput(imagegradient);
  writer->Update();

  GaussianType::Pointer gausfiltrelap=GaussianType::New();
	gausfiltrelap->SetInput(imagefiltre);
	gausfiltrelap->SetVariance(1);
	gausfiltrelap->SetMaximumKernelWidth(5);
	gausfiltrelap->Update();
	imagelaplacian=gausfiltrelap->GetOutput();
 

  LaplacianFilter::Pointer lapFilter = LaplacianFilter::New();
  lapFilter->SetInput(imagelaplacian);
  lapFilter->Update();
  ImageTypeLap::Pointer laplacian=lapFilter->GetOutput();

  ZeroCrossingFilter::Pointer zeroFilter = ZeroCrossingFilter::New();
  zeroFilter->SetInput(laplacian);
  ImageTypeLap::Pointer EdgeLap=zeroFilter->GetOutput();

  WriterLapType::Pointer   writer2 =  WriterLapType::New();
  itk::NiftiImageIO::Pointer ioimagenifti2=itk::NiftiImageIO::New();

  writer2->SetImageIO(ioimagenifti2);
  writer2->SetFileName( "C:/Users/Marc-Antoine/Documents/Imagecode/output/laplacianedge.nii");
  writer2->SetInput(EdgeLap);
  writer2->Update();

  */

}

void Filtre::createkernel(unsigned int largueur){
	ImageType::RegionType masque;
	ImageType::RegionType::IndexType start;
	start.Fill(0);

	ImageType::SizeType size;
	size.Fill(largueur);

	
	 masque.SetSize(size);
	 masque.SetIndex(start);

	
	 kernel->SetRegions(masque);
	 kernel->Allocate();

	 ImageIteratorType imageIterator(kernel, masque);
	
	 int compteur=0;
	 int grandeur=largueur*largueur*largueur;
	 int milieu=(grandeur+1)/2;
	 while(!imageIterator.IsAtEnd())
    {
		compteur++;

		imageIterator.Set(1);

		if(compteur==milieu)
			imageIterator.Set(grandeur*0.25);
    
		
 
	    ++imageIterator;
		
    }
	 
	 cout<<"grandeur patch"<<compteur<<endl;
	 cout<<"milieu"<<milieu<<endl;
}

void Filtre::threshold(){
	ThresholdFilter::Pointer threshimage=ThresholdFilter::New();
	threshimage->SetInput(imagefiltre);
	threshimage->SetOutsideValue(0);
	//threshimage->SetUpper(120);
	//threshimage->SetLower(70);
	threshimage->ThresholdOutside(60,140);
	threshimage->Update();
	//imageresult=threshimage->GetOutput();

	BinaryThresholdFilter::Pointer binthreshimage=BinaryThresholdFilter::New();
	binthreshimage->SetInput(imagefiltre);
	binthreshimage->SetUpperThreshold(160);
	binthreshimage->SetLowerThreshold(73);
	binthreshimage->SetInsideValue(0);
	binthreshimage->SetOutsideValue(250);
	binthreshimage->Update();
	imageresult=binthreshimage->GetOutput();
	

}

void Filtre::sauvegardeimage(){
	//sauvegarde l'image modifier
	cout<<"sauvegarde de limage"<<endl;
  WriterType::Pointer   writer =  WriterType::New();
  itk::NiftiImageIO::Pointer ioimagenifti=itk::NiftiImageIO::New();

  writer->SetImageIO(ioimagenifti);
  writer->SetFileName( "C:/Users/Marc-Antoine/Documents/Imagecode/output/resultthresh.nii");
  writer->SetInput(imageresult);
  writer->Update();

  WriterType::Pointer   writer2 =  WriterType::New();
 

 // writer2->SetImageIO(ioimagenifti);
 // writer2->SetFileName( "C:/Users/Marc-Antoine/Documents/Imagecode/output/filtreconvolution.nii");
 // writer2->SetInput(filtre);
 // writer2->Update();

}

ImageType::Pointer Filtre::getresult(){
	return imageresult;

}