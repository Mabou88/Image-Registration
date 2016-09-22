#include "Similaritymeasure.h"

Similaritymeasure::Similaritymeasure(){

	BinaryReaderType::Pointer reader_m = BinaryReaderType::New();
    itk::NiftiImageIO::Pointer mask_io = itk::NiftiImageIO::New();
    reader_m->SetImageIO(mask_io);
    reader_m->SetFileName("C:/Users/Marc-Antoine/Documents/Imagecode/test/2e/Us_ventricule.nii");
    try {
        reader_m->Update();
    } catch (itk::ExceptionObject &e) {
        cout<<"Error while reading US mask"<<endl;
        cout<<e<<endl;
        EXIT_FAILURE;
    }
    	
	image_us_seg=BinaryImageType::New();
	image_us_seg=reader_m->GetOutput();
 
  
	BinaryReaderType::Pointer reader_irm = BinaryReaderType::New();
    reader_irm->SetImageIO(mask_io);
    reader_irm->SetFileName("C:/Users/Marc-Antoine/Documents/Imagecode/output/Affine_auto_centre/Patient 2/PH Gaus+grad/IRM_ventricule.nii");
    try {
        reader_irm->Update();
    } catch (itk::ExceptionObject &e) {
        cout<<"Error while reading US mask"<<endl;
        cout<<e<<endl;
        EXIT_FAILURE;
    }
    	
	image_irm_seg=BinaryImageType::New();
	image_irm_seg=reader_irm->GetOutput();
	
	
	
	
}

void Similaritymeasure::computeoverlap(){
	BinaryImageType::RegionType ROI_us=image_us_seg->GetLargestPossibleRegion();
	ImageBinaryConstIteratorType Us_it(image_us_seg,ROI_us);

	

	Us_it.GoToBegin();

	int comptetotalus=0;
	int compteoverlap=0;
	int compteoverlap1=0;
	int compteoverlap2=0;
	
	int compteiterateur=0;

	//parcours de la ROI avec l'iterateur
	while (!Us_it.IsAtEnd()){
		BinaryImageType::IndexType index=Us_it.GetIndex();


		//compte la grandeur des ventricules latéraux
		if (image_us_seg->GetPixel(index)==1){
			comptetotalus++;

			ImageType::PointType p;
			 //get the equivalent in physical point to evaluate whether it's within the MRI image
			image_us_seg->TransformIndexToPhysicalPoint(index,p);
			ImageType::IndexType i;
			image_irm_seg->TransformPhysicalPointToIndex(p, i);
		
			if (image_irm_seg->GetPixel(i)==1)
				compteoverlap++;
		}
		
		ImageType::PointType p;
		//get the equivalent in physical point to evaluate whether it's within the MRI image
		image_us_seg->TransformIndexToPhysicalPoint(index,p);
		ImageType::IndexType i;
		image_irm_seg->TransformPhysicalPointToIndex(p, i);
		
		if (image_irm_seg->GetPixel(i)==1)
				compteoverlap1++;
		//incrementation de l'iterateur
		Us_it++;
		compteiterateur++;
	}

	cout<<"comptage :"<<comptetotalus<<" 2e overlap: "<<compteoverlap<<" 3e overlap: "<<compteoverlap1<<endl;

	BinaryImageType::RegionType ROI_irm=image_irm_seg->GetLargestPossibleRegion();
	ImageBinaryConstIteratorType Irm_it(image_irm_seg,ROI_irm);

	Irm_it.GoToBegin();

	int comptetotalirm=0;
	
	int compteiterateur2=0;

	//parcours de la ROI avec l'iterateur
	while (!Irm_it.IsAtEnd()){
		BinaryImageType::IndexType index=Irm_it.GetIndex();


		//compte la grandeur des ventricules latéraux
		if (image_irm_seg->GetPixel(index)==1){
			comptetotalirm++;

			ImageType::PointType p;
			 //get the equivalent in physical point to evaluate whether it's within the MRI image
			image_irm_seg->TransformIndexToPhysicalPoint(index,p);
			ImageType::IndexType i;
			image_us_seg->TransformPhysicalPointToIndex(p, i);
		
			if (image_us_seg->GetPixel(i)==1)
				compteoverlap2++;
		}
			

				
		//incrementation de l'iterateur
		Irm_it++;
		compteiterateur2++;
	}

	//mise à l'échelle car us est shrink de 8 en volume
	//Les 2e devraient etre égale avec le meme shrink, à essayer
	comptetotalirm=comptetotalirm;
	compteoverlap2=compteoverlap2;

	cout<<"comptage2 :"<<comptetotalirm<<" 2e: "<<compteoverlap2<<endl;
	cout<<"Les 2e devraient etre égale avec le meme shrink, à essayer"<<endl;

	similarityindex=float(compteoverlap+compteoverlap2)/float(comptetotalus+comptetotalirm);


}

float Similaritymeasure::getindex(){
	/*
	SimilarityIndex::Pointer computeindex=SimilarityIndex::New(); 

	const ImageType::RegionType input1=image_us_seg->GetLargestPossibleRegion();
	ImageType::RegionType input2=image_irm_seg->GetLargestPossibleRegion();

	computeindex->SetInput1( image_us_seg);
	computeindex->SetInput2(image_irm_seg);

	SimilarityIndex::RealType valeurindex=computeindex->GetSimilarityIndex();
	similarityindex=computeindex->GetSimilarityIndex();
	cout<<"valeur index: "<<valeurindex<<endl;
	cout<<"valeur similarite index: "<<similarityindex<<endl;
	
	//ImageType::RegionType imageregion=image_us_seg->GetLargestPossibleRegion();

	LabelMapFilter::Pointer filterlabel=LabelMapFilter::New();
	filterlabel->SetInput(image_us_seg);
	filterlabel->Update();
	
	LabelMapType::Pointer Us_label;

	Us_label=filterlabel->GetOutput();
	cout<<"nombre object"<<filterlabel->GetOutput()->GetNumberOfLabelObjects()<<endl;
	LabelMapFilter::Pointer filterlabel2=LabelMapFilter::New();
	filterlabel2->SetInput(image_irm_seg);
	filterlabel2->Update();
	
	LabelMapType::Pointer Irm_label;

	Irm_label=filterlabel2->GetOutput();
	*/
	
	
	//OverlapMeasurement::Pointer computeoverlap;
	
	//OverlapMeasurement::Pointer computeoverlap=OverlapMeasurement::New();
	/*
	computeoverlap->SetSourceImage( image_us_seg);
	computeoverlap->SetTargetImage( image_irm_seg);
	//computeoverlap->SetInput(Us_label);
	//computeoverlap->SetInput(Irm_label);
	
	RealType dice=computeoverlap->GetDiceCoefficient();
	float dice2=computeoverlap->GetDiceCoefficient();

	RealType jac=computeoverlap->GetJaccardCoefficient();

	RealType overlap=computeoverlap->GetTotalOverlap();

	RealType error=computeoverlap->GetFalseNegativeError();

	cout<<"valeur dice"<<dice<<endl;
	cout<<"valeur dice"<<dice2<<endl;
	cout<<"valeur jaccard"<<jac<<endl;
	cout<<"valeur overlap"<<overlap<<endl;
	cout<<"valeur erreur"<<error<<endl;
	
	*/

	return similarityindex;
  }