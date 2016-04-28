#include "Matrice.h"
#include "gradient.h"
#include <iostream>
#include <mitkITKImageImport.h>
//fichier itk
#include "itkimage.h"
#include "itkImageFileReader.h"
#include "itkImageIOBase.h"
#include "itkNiftiImageIO.h"

//typedef itk::Image<short,3> ImageType;
//typedef itk::Index<3> Index;


Matricelc2::Matricelc2(int nbrvoisins, Index index, mitk::BaseGeometry::Pointer geo_irm, mitk::BaseGeometry::Pointer geo_us, mitk::Image::Pointer image1,mitk::Image::Pointer image2,float facteurconversion,mitk::Image::Pointer mitkgradient)
{
	facteurconversion_=facteurconversion;
	//calcule la grandeur de la patch donc des matrices
	grandeurpatch=pow((2*nbrvoisins)+1,3);

	//incrementation pour remplir les matrices
	int position=0;
	
	//accède à au volume de l'image
	mitk::ImageDataItem::Pointer volume=image1->GetVolumeData(0);
	mitk::ImageDataItem::Pointer volumeus=image2->GetVolumeData(0);
	mitk::ImageDataItem::Pointer volumegrad=mitkgradient->GetVolumeData(0);

	//crée un accès à l'image pour lire les pixels
	mitk::ImagePixelReadAccessor<short,3> readAccess(image1, volume);
	mitk::ImagePixelReadAccessor<short,3> readAccessus(image2, volumeus);
	mitk::ImagePixelReadAccessor<short,3> readAccessgrad(mitkgradient, volumegrad);
	
	for (int i=-nbrvoisins;i<nbrvoisins+1;i++){
		for (int j=-nbrvoisins;j<nbrvoisins+1;j++){
			for (int k=-nbrvoisins;k<nbrvoisins+1;k++){

				//les index pour chacun des pixels de la patch
				itk::Index<3> index_irm_temp;
				itk::Index<3> index_us_temp;

				//index calculé à partir du pixel central
				index_us_temp[0]=index[0]+k;
				index_us_temp[1]=index[1]+j;
				index_us_temp[2]=index[2]+i;

				//convertit en index de l'US pour la même coordonné spatial (mm)
				mitk::Point3D positionspatial;
				geo_us->IndexToWorld(index_us_temp,positionspatial);
				geo_irm->WorldToIndex(positionspatial,index_irm_temp);
				
				//accède à la valeur pointé par l'index
				short value = readAccess.GetPixelByIndex(index_irm_temp);
				short valuegrad=readAccessgrad.GetPixelByIndex(index_irm_temp);
				short valueus=readAccessus.GetPixelByIndex(index_us_temp);
				
				//construit la matrice et le vecteur
				matrice[position][0]=double(value);
				matrice[position][1]=double(valuegrad);
				matrice[position][2]=double(1);
				vecteur[position]=double(valueus);
				position++;
				//*double(facteurconversion_)
			}

		}


	}
	
};


int Matricelc2::getgrandeurpatch(){
	return grandeurpatch;
}

double Matricelc2::getmat(){
return matrice[27][3];
}
double Matricelc2::getmattrans(){
return matricetrans[3][27];
};
double Matricelc2::getmatavantinv(){
return matriceavantinv[3][3];
};
double Matricelc2::getmatinv(){
return matriceinv[3][3];
};
double Matricelc2::getvecteur(){
return vecteur[27];
}
double Matricelc2::getalpha(){
return coefficient[0];
}
double Matricelc2::getbeta(){
return coefficient[1];
}
double Matricelc2::getgamma(){
return coefficient[2];
}
double Matricelc2::getlc2local(){
	return lc2local;
}
double Matricelc2::getvariance(){
	return variance;
}

void Matricelc2::affichemat(){
	cout<<"matrice"<<endl;
	for (int i=0;i<grandeurpatch;i++){
		for (int j=0;j<3;j++){
			cout<<matrice[i][j]<<endl;
		}
	}
}
void Matricelc2::affichevect(){
	cout<<"vecteur"<<endl;
	for (int i=0;i<grandeurpatch;i++){
		cout<<vecteur[i]<<endl;
	}
}
void Matricelc2::affichemattrans(){
	cout<<"transpose"<<endl;
	for (int i=0;i<3;i++){
		for (int j=0;j<grandeurpatch;j++){
			cout<<matricetrans[i][j]<<endl;
		}
	}
}
void Matricelc2::affichematavantinv(){
	cout<<"matriceavantinverse"<<endl;
	for (int i=0;i<3;i++){
		for (int j=0;j<3;j++){
			cout<<matriceavantinv[i][j]<<endl;
		}
	}
}
void Matricelc2::affichematinv(){
	cout<<"matinverse"<<endl;
	for (int i=0;i<3;i++){
		for (int j=0;j<3;j++){
			cout<<matriceinv[i][j]<<endl;
		}
	}
}
void Matricelc2::affichematfinal(){
	cout<<"matfinal"<<endl;
	for (int i=0;i<3;i++){
		for (int j=0;j<grandeurpatch;j++){
			cout<<matricefinal[i][j]<<endl;
		}
	}
}

void Matricelc2::affichecoef(){
	cout<<"coefficient"<<endl;
	for (int j=0;j<3;j++){
			cout<<coefficient[j]<<endl;
	}
}

void Matricelc2::domatricetransverse()
{
	for (int i=0;i<grandeurpatch;i++){
		for (int j=0;j<3;j++){
			matricetrans[j][i]=matrice[i][j];
		}
	}
}

void Matricelc2::matricemultipli()
{
		for (int i=0;i<3;i++){
			for (int k=0;k<3;k++){
				matriceavantinv[i][k]=0;
				for (int j=0;j<grandeurpatch;j++){
					matriceavantinv[i][k]+=matricetrans[i][j]*matrice[j][k];
				}
			}
		}
		
}
void Matricelc2::matricemultipli2()
{
		for (int i=0;i<3;i++){
			for (int k=0;k<grandeurpatch;k++){
				matricefinal[i][k]=0;
				for (int j=0;j<3;j++){
					matricefinal[i][k]+=matriceinv[i][j]*matricetrans[j][k];
				}
			}
		}
		
}
void Matricelc2::calculecoefficient()
{
		for (int i=0;i<3;i++){
			coefficient[i]=0;
			for (int k=0;k<grandeurpatch;k++){
				coefficient[i]+=matricefinal[i][k]*vecteur[k];
			}
		}
		
}
void Matricelc2::domatriceinverse()
{
	//calcul du determinant
	double det1=(matriceavantinv[0][0]*matriceavantinv[1][1]*matriceavantinv[2][2]+matriceavantinv[1][0]*matriceavantinv[2][1]*matriceavantinv[0][2]+matriceavantinv[0][1]*matriceavantinv[1][2]*matriceavantinv[2][0]);
    double det2=(matriceavantinv[0][2]*matriceavantinv[1][1]*matriceavantinv[2][0] + matriceavantinv[0][1]*matriceavantinv[1][0]*matriceavantinv[2][2] +matriceavantinv[1][2]*matriceavantinv[2][1]*matriceavantinv[0][0]);
	double det=det1-det2;
	//méthode alternative
	double detv1=matriceavantinv[0][0]*((matriceavantinv[1][1]*matriceavantinv[2][2])-(matriceavantinv[1][2]*matriceavantinv[2][1]));
	double detv2=-1*matriceavantinv[0][1]*((matriceavantinv[1][0]*matriceavantinv[2][2])-(matriceavantinv[1][2]*matriceavantinv[2][0]));
	double detv3=matriceavantinv[0][2]*((matriceavantinv[1][0]*matriceavantinv[2][1])-(matriceavantinv[2][0]*matriceavantinv[1][1]));
	double detv0=detv1+detv2+detv3;
	
	matriceinv[0][0]=1/detv0*(matriceavantinv[1][1]*matriceavantinv[2][2]-matriceavantinv[1][2]*matriceavantinv[2][1]);
	matriceinv[1][0]=1/detv0*(matriceavantinv[1][2]*matriceavantinv[2][0]-matriceavantinv[1][0]*matriceavantinv[2][2]);
	matriceinv[2][0]=1/detv0*(matriceavantinv[1][0]*matriceavantinv[2][1]-matriceavantinv[2][0]*matriceavantinv[1][1]);
	matriceinv[0][1]=1/detv0*(matriceavantinv[0][2]*matriceavantinv[2][1]-matriceavantinv[0][1]*matriceavantinv[2][2]);
	matriceinv[1][1]=1/detv0*(matriceavantinv[0][0]*matriceavantinv[2][2]-matriceavantinv[2][0]*matriceavantinv[0][2]);
	matriceinv[2][1]=1/detv0*(matriceavantinv[0][1]*matriceavantinv[2][0]-matriceavantinv[0][0]*matriceavantinv[2][1]);
	matriceinv[0][2]=1/detv0*(matriceavantinv[0][1]*matriceavantinv[1][2]-matriceavantinv[0][2]*matriceavantinv[1][1]);
	matriceinv[1][2]=1/detv0*(matriceavantinv[0][2]*matriceavantinv[1][0]-matriceavantinv[0][0]*matriceavantinv[1][2]);
	matriceinv[2][2]=1/detv0*(matriceavantinv[0][0]*matriceavantinv[1][1]-matriceavantinv[0][1]*matriceavantinv[1][0]);
		
		
}

void Matricelc2::calculelc2local(){
	lc2local=0;
	cout<<"simulation lc2"<<endl;
	for (int i=0;i<grandeurpatch;i++){
	double transIRM=coefficient[0]*matrice[i][0]+coefficient[1]*matrice[i][1]+coefficient[2]*matrice[i][2];
	double diff=vecteur[i]-transIRM;
	lc2local+=pow(diff,2);
	}
	
	//lc2local=1-lc2local/variance;
	lc2local=lc2local/variance;
	cout<<lc2local<<endl;
}

void Matricelc2::calculevariance(){
	double moyenne=0;
	variance=0;
	cout<<"variance"<<endl;
	for (int i=0;i<grandeurpatch;i++){
	moyenne+=vecteur[i];
	}
	moyenne=moyenne/grandeurpatch;
	for (int j=0;j<grandeurpatch;j++){
		double diff=vecteur[j]-moyenne;
		variance+=pow(diff,2);
	}
	variance=variance/grandeurpatch;
	cout<<variance<<endl;
}