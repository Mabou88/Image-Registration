#include "Normalisationimage.h"
#include <iostream>

Normalisationimage::Normalisationimage(){
	maxirm_=0;
	maxus_=0;
	facteurconversion_=0;
}

void Normalisationimage::calculefacteur(){
	facteurconversion_=float(maxirm_)/float(maxus_);
}

float Normalisationimage::getfacteur(){
	return facteurconversion_;
}

void Normalisationimage::calculemax(mitk::Image::Pointer image,int modalite){
	


	//Prend les grandeurs des dimensions
	unsigned int dimx=image->GetDimension(0);
	unsigned int dimy=image->GetDimension(1);
	unsigned int dimz=image->GetDimension(2);

	//Initialise le maximum
	double maximum=0;
	
	mitk::ImageDataItem::Pointer volume=image->GetVolumeData(0);
	//crée un accès à l'image pour lire les pixels
	mitk::ImagePixelReadAccessor<short,3> readAccess(image, volume);
	for (int i=0;i<dimx;i++){
		for (int j=0;j<dimy;j++){
			for (int k=0;k<dimz;k++){
				itk::Index<3> idx = {{ i, j,k}};
				short value = readAccess.GetPixelByIndex(idx);
				if (double(value)>maximum){
					maximum=double(value);
				}
			}
		}
	}
	cout<<"maximum"<<endl;
	cout<<maximum<<endl;
	//dependamment de l'image de US ou IRM 
	switch(modalite){
	case 0:
		 maxus_=maximum;
	case 1:
		maxirm_=maximum;
	}
	
}