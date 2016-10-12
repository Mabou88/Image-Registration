//
//  LC2_function.cpp
//  LC2_M
//
//  Created by Maxime Gérard on 10/02/16.
//  Copyright © 2016 Maxime Gérard. All rights reserved.
//

#include "LC2_function.hpp"

LC2_function::LC2_function(ImageType::Pointer im_Fixed, ImageType::Pointer im_Moving,string out, int type,ImageType::PointType centre,int coderegion)
{
    //std::cout<<"creating measure object"<<std::endl;
    m_outputPath=out;
    m_FixedImage = im_Fixed;
    m_MovingImage = im_Moving;
    computeGradient();
    //limitMRI();
    computeMask();
    
    m_maxTrans =0;
    m_maxRot = 0;
	m_maxScale=0;
	m_maxShear=0;
    m_radius = 0;
    
    m_useLiverMask=false;
    nombreloop=0;


	//type 1:Euler/rigide 2: Affine 
	typetransform=type;
	centrerot=centre;
	//coderegion 1:hemisphere avant central 2: ventricule complet
	coderegion_=coderegion;

	//coord
	//hemisphere centrale
	if (coderegion_==1){
		coordregion_[0]=122;
		coordregion_[1]=138;
		coordregion_[2]=65;
		coordregion_[3]=90;
		coordregion_[4]=82;
		coordregion_[5]=87;

	}
	//ventricule complet
	else if (coderegion_==2){
		coordregion_[0]=95;
		coordregion_[1]=170;
		coordregion_[2]=50;
		coordregion_[3]=145;
		coordregion_[4]=64;
		coordregion_[5]=84;

	//ventricule avant
	}else if(coderegion_==3){
		coordregion_[0]=118;
		coordregion_[1]=150;
		coordregion_[2]=55;
		coordregion_[3]=85;
		coordregion_[4]=66;
		coordregion_[5]=72;
	//haut des ventricules
	}else if (coderegion_==4){
		coordregion_[0]=105;
		coordregion_[1]=160;
		coordregion_[2]=60;
		coordregion_[3]=140;
		coordregion_[4]=78;
		coordregion_[5]=83;


	}
	//bout posterieur ventricules
	else if(coderegion_==5){
		coordregion_[0]=96;
		coordregion_[1]=160;
		coordregion_[2]=110;
		coordregion_[3]=140;
		coordregion_[4]=76;
		coordregion_[5]=70;


	}
    
}

void LC2_function::setLiverMask(MaskType::Pointer liver)
{
    m_useLiverMask=true;
    m_LiverMask=liver;
    
    //we than transform the fixed MRI with that mask to keep only the liver info
    //in the computation of LC2 we are going to have to exclude all points outside that mask
}

void LC2_function::limitMRI()
{
    std::cout<<"test limitating MRI image to interesting area of the liver"<<endl;
    ImageType::IndexType startPoint;
    startPoint[0] = 119;
    startPoint[1] = 180;
    startPoint[2] = 0;
    
    ImageType::SizeType sizeR;
    sizeR[0] = 80;
    sizeR[1] = 86;
    sizeR[2] = 79;
    
    ImageType::RegionType regionR;
    regionR.SetIndex(startPoint);
    regionR.SetSize(sizeR);
    
    ExtractorType::Pointer extractFilter = ExtractorType::New();
    extractFilter->SetExtractionRegion(regionR);
    extractFilter->SetInput(m_FixedImage);
    extractFilter->SetDirectionCollapseToIdentity();
    try {
        extractFilter->Update();
    } catch (itk::ExceptionObject &e) {
        std::cerr<<"error while cropping MRI image"<<std::endl;
        std::cerr<<e<<std::endl;
    }
    
    m_accessibleMRI = extractFilter->GetOutput();
    
    //enregistrement pour verification
    WriterType::Pointer writer10 = WriterType::New();
    itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
    std::string out10 = m_outputPath+"/testCroppingMRI.nii.gz";
    writer10->SetImageIO(io);
    writer10->SetInput(m_accessibleMRI);
    writer10->SetFileName(out10);
    try {
        writer10->Update();
    } catch (itk::ExceptionObject &e) {
        std::cerr<<"error while writing cropped MRI"<<std::endl;
        std::cerr<<e<<std::endl;
    }
}

void LC2_function::computeGradient()
{
    std::cout<<"compute gradient image of fixed MRI image"<<std::endl;
    
    GradientFilterType::Pointer filterG = GradientFilterType::New();
    filterG->SetInput(m_MovingImage);
    try {
        filterG->Update();
    } catch (itk::ExceptionObject &e) {
        std::cerr<<"error while computing gradient image"<<std::endl;
        std::cerr<<e<<std::endl;
        EXIT_FAILURE;
    }
    
    m_grad = filterG->GetOutput();
    
    //write image for test
//    
//            typename WriterType::Pointer writer1 = WriterType::New();
//            string out1 = m_outputPath+"/testgrad.nii.gz";
//            itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
//            writer1->SetInput(m_grad);
//            writer1->SetImageIO(io);
//            writer1->SetFileName(out1);
//            try {
//                writer1->Update();
//            } catch (itk::ExceptionObject &e) {
//                cerr<<"error while writing image file"<<endl;
//                cerr<<e<<endl;
//                EXIT_FAILURE;
//            }
//    
//            cout<<"done writing gradient image"<<endl;

    
}

void LC2_function::computeMask()
{
    std::cout<<"computing cropping mask"<<std::endl;
    //binarisation
    
    BinaryThresholdFilterType::Pointer thresholder = BinaryThresholdFilterType::New();
    thresholder->SetInput(m_MovingImage);
    thresholder->SetOutsideValue(255);
    thresholder->SetInsideValue(0);
    thresholder->SetLowerThreshold(0);
    thresholder->SetUpperThreshold(1);
    
    try {
        thresholder->Update();
    } catch (itk::ExceptionObject &e) {
        std::cerr<<"error while binarizing US image"<<std::endl;
        std::cerr<<e<<std::endl;
    }
    
    MaskType::Pointer mask1 =thresholder->GetOutput() ;
    
    
    std::cout<<"done writing initial mask image"<<std::endl;
    
    std::cout<<"test closing"<<std::endl;
    //operation morphologiques
    
    /*
    kernelType::RadiusType radius;
    radius.Fill(3);
    
    kernelType kernel = kernelType::Ball(radius);
    std::cout<<"radius kernel : "<<kernel.GetRadius()<<std::endl;
    std::cout<<"kernel size : "<<kernel.GetSize()<<std::endl;
    
    CloserType::Pointer closer = CloserType::New();
    closer->SetInput(mask1);
    closer->SetKernel(kernel);
    closer->Update();
    
    m_mask = closer->GetOutput();
   */
	m_mask= mask1;
    //writing mask images
    
    BinaryWriterType::Pointer writer3 = BinaryWriterType::New();
	std::string out3 = m_outputPath+"testmask3.nii.gz";
    itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
    writer3->SetInput(m_mask);
    writer3->SetImageIO(io);
    writer3->SetFileName(out3);
    try {
        writer3->Update();
    } catch (itk::ExceptionObject &e) {
        std::cerr<<"error while writing image file"<<std::endl;
        std::cerr<<e<<std::endl;
        EXIT_FAILURE;
    }
    
    std::cout<<"done writing final mask image"<<std::endl;

    
}

ImageType::Pointer LC2_function::TransformImage(const dlib::matrix<double> &params, int ind) const
{
    ImageType::Pointer imageTransformed;
    //nombreloop+=1;
    //1 represents euler tsf
    if(ind==1)
    {
        EulerTransformType::Pointer transform = EulerTransformType::New();
        EulerTransformType::ParametersType parameters(6);
        //mise a l'echelle des parametres
        parameters[0] = params(0)*(m_maxRot/m_radius);
        parameters[1] = params(1)*(m_maxRot/m_radius);
        parameters[2] = params(2)*(m_maxRot/m_radius);
        parameters[3] = params(3)*(m_maxTrans/m_radius);
        parameters[4] = params(4)*(m_maxTrans/m_radius);
        parameters[5] = params(5)*(m_maxTrans/m_radius);


		//on bloque ces operations geometriques pour la region
		if (coderegion_==1){
		parameters[0]=0;
		parameters[1]=0;
		parameters[4]=0;
		parameters[5]=0;
		}
		if(coderegion_==4){
			parameters[0]=0;
			parameters[2]=0;
			parameters[3]=0;
			parameters[4]=0;
			//parameters[5]=0;
		}
		if(coderegion_==3){
			parameters[0]=0;
			parameters[1]=0;
			parameters[2]=0;

		}
        transform->SetParameters(parameters);

		cout<<"parametre bobyqa_"<<params<<endl;
		cout<<"rot_"<<m_maxRot<<"trans_"<<m_maxTrans<<"radius_"<<m_radius<<endl;


        std::cout<<"euler tsf parameters : "<<transform->GetParameters()<<std::endl;
        
        ImageType::SizeType sizeUS = m_MovingImage->GetLargestPossibleRegion().GetSize();
        ImageType::PointType origin = m_MovingImage->GetOrigin();
        ImageType::SpacingType spacing = m_MovingImage->GetSpacing();
        ImageType::PointType center;
       // center[0] = origin[0]+spacing[0]*sizeUS[0]/2;
        //center[1] = origin[1]+spacing[1]*sizeUS[1]/2;
        //center[2] = origin[2]+spacing[2]*sizeUS[2]/2;
       
		//center[0] =0;
        //center[1] =0;
        //center[2] =0;


		//Patient 5
		//center[0]=111;
		//center[1]=80;
		//center[2]=111;

		cout<<"center of rotation"<<center<<endl;
        
        EulerTransformType::ParametersType eulerFixedParameters(3);
        eulerFixedParameters[0] =centrerot[0];
        eulerFixedParameters[1] =centrerot[1];
        eulerFixedParameters[2] =centrerot[2];
        
        transform->SetFixedParameters(eulerFixedParameters);
        //std::cout<<"tsf fixed param : "<<transform->GetFixedParameters()<<std::endl;
        
        
        
        ResamplerType::Pointer resamplefilter = ResamplerType::New();
        resamplefilter->SetInput(m_MovingImage);
        resamplefilter->SetSize(m_FixedImage->GetLargestPossibleRegion().GetSize());
        resamplefilter->SetOutputSpacing(m_FixedImage->GetSpacing());
        resamplefilter->SetOutputDirection(m_FixedImage->GetDirection());
        resamplefilter->SetOutputOrigin(m_FixedImage->GetOrigin());
        resamplefilter->SetTransform(transform);
        //resamplefilter->SetTransform(transform);
        
        try {
            resamplefilter->Update();
        } catch (itk::ExceptionObject &e) {
            std::cerr<<"error while transforming moving image"<<std::endl;
            std::cerr<<e<<std::endl;
        }
        
       imageTransformed= resamplefilter->GetOutput();

        
    }
	else if (ind=2){
		AffineTransformType::Pointer transformaff = AffineTransformType::New();
		AffineTransformType::Pointer transformaff2 = AffineTransformType::New();
		AffineTransformType::Pointer transformaff3 = AffineTransformType::New();
		AffineTransformType::Pointer transformaff4 = AffineTransformType::New();
		
		

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
		double anglerotation_x=params(0)*(m_maxRot/m_radius);
		double anglerotation_y=params(1)*(m_maxRot/m_radius);
		double anglerotation_z=params(2)*(m_maxRot/m_radius);

	
	
		std::cout<<"Affine tsf parameters : "<<transformaff->GetParameters()<<std::endl;

		transformaff->Rotate3D(axederotation_x,anglerotation_x,true);
		transformaff->Rotate3D(axederotation_y,anglerotation_y,true);
		transformaff->Rotate3D(axederotation_z,anglerotation_z,true);

       // std::cout<<"Affine tsf parameters : "<<transformaff->GetParameters()<<std::endl;

		//essayer de scaler isotropique aussi!!!!!!
		
		
		//scale isotropique
	//Scale
		AffineTransformType::OutputVectorType scalefactor;
		scalefactor[0]=1+params(3)*((m_maxScale-1)/m_radius);
		scalefactor[1]=1+params(4)*((m_maxScale-1)/m_radius);
		scalefactor[2]=1+params(5)*((m_maxScale-1)/m_radius);

		//transformaff->Scale(scalefactor,true);

		
		
		 std::cout<<"Affine tsf parameters : "<<transformaff->GetParameters()<<std::endl;
		 
	//Shear
		


		 double coefficient[6];
		coefficient[0]=params(6)*((m_maxShear)/m_radius);
		coefficient[1]=params(7)*((m_maxShear)/m_radius);
		coefficient[2]=params(8)*((m_maxShear)/m_radius);
		coefficient[3]=params(9)*((m_maxShear)/m_radius);
		coefficient[4]=params(10)*((m_maxShear)/m_radius);
		coefficient[5]=params(11)*((m_maxShear)/m_radius);

		//a verifier peut besoin de 6 parametre pcq inverser 0 et 1
		 transformaff->Shear(0,1,coefficient[0]);
		 transformaff->Shear(1,0,coefficient[1]);
		 transformaff->Shear(2,0,coefficient[2]);
		 transformaff->Shear(0,2,coefficient[3]);
		 transformaff->Shear(1,2,coefficient[4]);
		 transformaff->Shear(2,1,coefficient[5]);
		 
		

		 std::cout<<"Affine tsf parameters : "<<transformaff->GetParameters()<<std::endl;
		 
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

		 vecteurtranslation_x=vecteurtranslation_x*params(12)*((m_maxTrans)/m_radius);
		 vecteurtranslation_y=vecteurtranslation_y*params(13)*((m_maxTrans)/m_radius);
		 vecteurtranslation_z=vecteurtranslation_z*params(14)*((m_maxTrans)/m_radius);

		 transformaff->Translate(vecteurtranslation_x,true);
		 transformaff->Translate(vecteurtranslation_y,true);
		 transformaff->Translate(vecteurtranslation_z,true);

		 std::cout<<"Affine tsf parameters : "<<transformaff->GetParameters()<<std::endl;


		 //Affine transform avec centre pré déterminé

		AffineTransformfixedcenter::Pointer AffineTransformfixe=AffineTransformfixedcenter::New();

		AffineTransformfixe->Rotate3D(axederotation_x,anglerotation_x,true);
		AffineTransformfixe->Rotate3D(axederotation_y,anglerotation_y,true);
		AffineTransformfixe->Rotate3D(axederotation_z,anglerotation_z,true);

		AffineTransformfixe->Scale(scalefactor,true);

	
		 AffineTransformfixe->Shear(0,1,coefficient[0]);
		 AffineTransformfixe->Shear(1,0,coefficient[1]);
		 AffineTransformfixe->Shear(2,0,coefficient[2]);
		 AffineTransformfixe->Shear(0,2,coefficient[3]);
		 AffineTransformfixe->Shear(1,2,coefficient[4]);
		 AffineTransformfixe->Shear(2,1,coefficient[5]);

		  AffineTransformfixe->Translate(vecteurtranslation_x,true);
		 AffineTransformfixe->Translate(vecteurtranslation_y,true);
		 AffineTransformfixe->Translate(vecteurtranslation_z,true);
     
       // ImageType::PointType center;
       
		//Patient 2
        //center[0] =98;
        //center[1] =68;
        //center[2] =98;

		//Patient 5
		//center[0]=111;
		//center[1]=80;
		//center[2]=111;


		cout<<"center of rotation"<<centrerot<<endl;

		AffineTransformfixe->SetCenterOfRotationComponent(centrerot);
        
        
                     
		ImageType::SizeType sizefixed=m_FixedImage->GetLargestPossibleRegion().GetSize();
		//sizefixed[1]+=100;
        
        ResamplerType::Pointer resamplefilteraff = ResamplerType::New();
        resamplefilteraff->SetInput(m_MovingImage);
        resamplefilteraff->SetSize(sizefixed);
        resamplefilteraff->SetOutputSpacing(m_FixedImage->GetSpacing());
        resamplefilteraff->SetOutputDirection(m_FixedImage->GetDirection());
        resamplefilteraff->SetOutputOrigin(m_FixedImage->GetOrigin());
        //resamplefilteraff->SetTransform(transformaff);
		resamplefilteraff->SetTransform(AffineTransformfixe);
        

		//InterpolatorNearestNeighbor::Pointer interpolator = InterpolatorNearestNeighbor::New();
		//resamplefilteraff->SetInterpolator( interpolator );

        
        try {
            resamplefilteraff->Update();
        } catch (itk::ExceptionObject &e) {
            std::cerr<<"error while transforming moving image"<<std::endl;
            std::cerr<<e<<std::endl;
        }
        
       imageTransformed= resamplefilteraff->GetOutput();

	}
    
    return imageTransformed;
}

MaskType::Pointer LC2_function::TransformMask(const dlib::matrix<double> &params, int ind) const
{
    MaskType::Pointer maskTsf;
    
    //1 represents euler tsf
    if(ind==1)
    {
        EulerTransformType::Pointer transform = EulerTransformType::New();
        EulerTransformType::ParametersType parameters(6);
        //mise a l'echelle des parametres
        parameters[0] = params(0)*(m_maxRot/m_radius);
		parameters[1] = params(1)*(m_maxRot/m_radius);
        parameters[2] = params(2)*(m_maxRot/m_radius);
        parameters[3] = params(3)*(m_maxTrans/m_radius);
        parameters[4] = params(4)*(m_maxTrans/m_radius);
		parameters[5] = params(5)*(m_maxTrans/m_radius);
		


		//on bloque ces operations geometriques pour la region
		if (coderegion_==1){
		parameters[0]=0;
		parameters[1]=1;
		parameters[4]=0;
		parameters[5]=0;
		}
		if(coderegion_==4){
			parameters[0]=0;
			parameters[2]=0;
			parameters[3]=0;
			parameters[4]=0;
			//parameters[5]=0;
		}
		if(coderegion_==3){
			parameters[0]=0;
			parameters[1]=0;
			parameters[2]=0;

		}

        transform->SetParameters(parameters);
        //std::cout<<"euler tsf parameters : "<<transform->GetParameters()<<std::endl;
        
        ImageType::SizeType sizeUS = m_MovingImage->GetLargestPossibleRegion().GetSize();
        ImageType::PointType origin = m_MovingImage->GetOrigin();
        ImageType::SpacingType spacing = m_MovingImage->GetSpacing();
        ImageType::PointType center;
        //center[0] = origin[0]+spacing[0]*sizeUS[0]/2;
        //center[1] = origin[1]+spacing[1]*sizeUS[1]/2;
        //center[2] = origin[2]+spacing[2]*sizeUS[2]/2;
		
		//center[0] =0;
        //center[1] =0;
        //center[2] =0;

		//Patient 5
		//center[0]=111;
		//center[1]=80;
		//center[2]=111;
        
        
        EulerTransformType::ParametersType eulerFixedParameters(3);
        eulerFixedParameters[0] =centrerot[0];
        eulerFixedParameters[1] =centrerot[1];
        eulerFixedParameters[2] =centrerot[2];
        
        transform->SetFixedParameters(eulerFixedParameters);
        //std::cout<<"tsf fixed param : "<<transform->GetFixedParameters()<<std::endl;
        
        ResamplerBinaryType::Pointer maskResampler = ResamplerBinaryType::New();
        maskResampler->SetInput(m_mask);
        maskResampler->SetOutputDirection(m_FixedImage->GetDirection());
        maskResampler->SetOutputOrigin(m_FixedImage->GetOrigin());
        maskResampler->SetOutputSpacing(m_FixedImage->GetSpacing());
        maskResampler->SetSize(m_FixedImage->GetLargestPossibleRegion().GetSize());
        maskResampler->SetTransform(transform);
        
        try {
            maskResampler->Update();
        } catch (itk::ExceptionObject &e) {
            cerr<<"error while transforming mesh"<<endl;
            cerr<<e<<endl;
        }
        
        maskTsf = maskResampler->GetOutput();
        
    }
	else if (ind=2){
		AffineTransformType::Pointer transformaff = AffineTransformType::New();
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
		double anglerotation_x=params(0)*(m_maxRot/m_radius);
		double anglerotation_y=params(1)*(m_maxRot/m_radius);
		double anglerotation_z=params(2)*(m_maxRot/m_radius);

	
	
	//	std::cout<<"Affine tsf parameters : "<<transformaff->GetParameters()<<std::endl;

		transformaff->Rotate3D(axederotation_x,anglerotation_x,true);
		transformaff->Rotate3D(axederotation_y,anglerotation_y,true);
		transformaff->Rotate3D(axederotation_z,anglerotation_z,true);

       // std::cout<<"Affine tsf parameters : "<<transformaff->GetParameters()<<std::endl;

		//essayer de scaler isotropique aussi!!!!!!

	//Scale
		AffineTransformType::OutputVectorType scalefactor;
		scalefactor[0]=1+params(3)*((m_maxScale-1)/m_radius);
		scalefactor[1]=1+params(4)*((m_maxScale-1)/m_radius);
		scalefactor[2]=1+params(5)*((m_maxScale-1)/m_radius);

		//transformaff->Scale(scalefactor,true);

		
		
		 //std::cout<<"Affine tsf parameters : "<<transformaff->GetParameters()<<std::endl;
		 
	//Shear
		


		 double coefficient[6];
		coefficient[0]=params(6)*((m_maxShear)/m_radius);
		coefficient[1]=params(7)*((m_maxShear)/m_radius);
		coefficient[2]=params(8)*((m_maxShear)/m_radius);
		coefficient[3]=params(9)*((m_maxShear)/m_radius);
		coefficient[4]=params(10)*((m_maxShear)/m_radius);
		coefficient[5]=params(11)*((m_maxShear)/m_radius);

		//a verifier peut besoin de 6 parametre pcq inverser 0 et 1
		 transformaff->Shear(0,1,coefficient[0]);
		 transformaff->Shear(1,0,coefficient[1]);
		 transformaff->Shear(2,0,coefficient[2]);
		 transformaff->Shear(0,2,coefficient[3]);
		 transformaff->Shear(1,2,coefficient[4]);
		 transformaff->Shear(2,1,coefficient[5]);
		 
		

		// std::cout<<"Affine tsf parameters : "<<transformaff->GetParameters()<<std::endl;
		 
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

		 vecteurtranslation_x=vecteurtranslation_x*params(12)*((m_maxTrans)/m_radius);
		 vecteurtranslation_y=vecteurtranslation_y*params(13)*((m_maxTrans)/m_radius);
		 vecteurtranslation_z=vecteurtranslation_z*params(14)*((m_maxTrans)/m_radius);

		 transformaff->Translate(vecteurtranslation_x,true);
		 transformaff->Translate(vecteurtranslation_y,true);
		 transformaff->Translate(vecteurtranslation_z,true);
        
		 //Affine transform en fixant le centre de rotation
		 AffineTransformfixedcenter::Pointer AffineTransformfixe=AffineTransformfixedcenter::New();

		AffineTransformfixe->Rotate3D(axederotation_x,anglerotation_x,true);
		AffineTransformfixe->Rotate3D(axederotation_y,anglerotation_y,true);
		AffineTransformfixe->Rotate3D(axederotation_z,anglerotation_z,true);

		AffineTransformfixe->Scale(scalefactor,true);

	
		 AffineTransformfixe->Shear(0,1,coefficient[0]);
		 AffineTransformfixe->Shear(1,0,coefficient[1]);
		 AffineTransformfixe->Shear(2,0,coefficient[2]);
		 AffineTransformfixe->Shear(0,2,coefficient[3]);
		 AffineTransformfixe->Shear(1,2,coefficient[4]);
		 AffineTransformfixe->Shear(2,1,coefficient[5]);

		  AffineTransformfixe->Translate(vecteurtranslation_x,true);
		 AffineTransformfixe->Translate(vecteurtranslation_y,true);
		 AffineTransformfixe->Translate(vecteurtranslation_z,true);
     
        //ImageType::PointType center;
       
		//Patient 2
       // center[0] =98;
        //center[1] =68;
        //center[2] =98;

		//Patient 5
		//center[0]=111;
		//center[1]=80;
		//center[2]=111;

		cout<<"center of rotation"<<centrerot<<endl;

		AffineTransformfixe->SetCenterOfRotationComponent(centrerot);


       
        ImageType::SizeType sizefixed=m_FixedImage->GetLargestPossibleRegion().GetSize();
		//sizefixed[1]+=100;
        
        
        ResamplerBinaryType::Pointer maskResampleraff = ResamplerBinaryType::New();
        maskResampleraff->SetInput(m_mask);
        maskResampleraff->SetOutputDirection(m_FixedImage->GetDirection());
        maskResampleraff->SetOutputOrigin(m_FixedImage->GetOrigin());
        maskResampleraff->SetOutputSpacing(m_FixedImage->GetSpacing());
        maskResampleraff->SetSize(sizefixed);
       // maskResampleraff->SetTransform(transformaff);
		maskResampleraff->SetTransform(AffineTransformfixe);
        
        try {
            maskResampleraff->Update();
        } catch (itk::ExceptionObject &e) {
            cerr<<"error while transforming mesh"<<endl;
            cerr<<e<<endl;
        }
        
        maskTsf = maskResampleraff->GetOutput();

	}
    
    return maskTsf;
}

void LC2_function::choosecoordregion(){
	
	//hemisphere centrale
	if (coderegion_==1){
		coordregion_[0]=122;
		coordregion_[1]=137;
		coordregion_[2]=77;
		coordregion_[3]=107;
		coordregion_[4]=82;
		coordregion_[5]=93;

	}
	//ventricule complet
	else if (coderegion_==2){
		coordregion_[0]=95;
		coordregion_[1]=170;
		coordregion_[2]=48;
		coordregion_[3]=145;
		coordregion_[4]=64;
		coordregion_[5]=84;


	}



}

