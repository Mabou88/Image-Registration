#include "Lc2.h"

Lc2::Lc2(mitk::BaseGeometry::Pointer geo_irm,mitk::BaseGeometry::Pointer geo_us,mitk::Image::Pointer image_irm,mitk::Image::Pointer image_us)
{
	//Clone l'image avant de la convertit en itk
	mitk::Image::Pointer imageclone = image_irm->Clone();

	//Convertit l'image IRM de mitk en image itk
	mitk::ImageToItk<ImageType>::Pointer toItkFilter =mitk::ImageToItk<ImageType>::New();
    toItkFilter->SetInput(imageclone);
    toItkFilter->Update();
	 ImageType::Pointer itkImage = toItkFilter->GetOutput();
	
	 //transfert des geometries pour construire les matrices
	geo_irm_=geo_irm;
	geo_us_=geo_us;
	image_irm_=image_irm;
	image_us_=image_us;
	 
	 //Normalisation de matrice en fonction de l'image
	Normalisationimage maximum;
	maximum.calculemax(image_us,0);
	maximum.calculemax(image_irm,1);
	maximum.calculefacteur();
	float facteurconversion=maximum.getfacteur();
	cout<<"facteur intensité"<<endl;
	cout<<facteurconversion<<endl;

	//calcule le gradient l'image
	Gradient gradient;
	gradient.compute_gradient(itkImage);
	ImageType::Pointer gradientimage=gradient.getImageGradient();
	//Reconvertit en image mitk
	 mitk::Image::Pointer mitkgradient = mitk::ImportItkImage(gradientimage);

	 //variable pour le calcul final de lc2
	 double variancesum=0;
	 double variancelc2poids=0;

	 //inscrire l'index du pixel central de la patch de départ sur l'US
	 itk::Index<3> idx = {{ 291,127,104}};


	 //section de l'image qui sera analysé à partir de l'index de départ (volume=longueur^3)
	 int longueur_section=3;

	 for (int x=0;x<longueur_section;x++){
		for (int y=0;y<longueur_section;y++){
			for (int z=0;z<longueur_section;z++){

				//defini les pixels centrals de chaquep patch
				itk::Index<3> idx_central;
				idx_central[2]=idx[2]+x;
				idx_central[1]=idx[1]+y;
				idx_central[0]=idx[0]+z;
				
				//Le nombre de voisins du pixel central détermine la grandeur de la patch dans Matrice
				  int nombrevoisins=1;

				  //construire matrice
				  Matricelc2 matricedelc2(nombrevoisins,idx_central,geo_irm_,geo_us_,image_irm_,image_us_,facteurconversion,mitkgradient);
				  //matricedelc2.affichemat();
				  //matricedelc2.affichevect();

				   //construire matrice transverse
				  matricedelc2.domatricetransverse();
				  //matricedelc2.affichemattrans();

				  //multiplication des matrices
				  matricedelc2.matricemultipli();
				   //matricedelc2.affichematavantinv();

				  //inverse du résultats
				  matricedelc2.domatriceinverse();
				  //matricedelc2.affichematinv();

				  matricedelc2.matricemultipli2();
				  //matricedelc2.affichematfinal();

				  matricedelc2.calculecoefficient();
   
  
				  //affiche alpha beta gamma
				  matricedelc2.affichecoef();

				  //calculer transformée de l'IRM
				  matricedelc2.calculevariance();
				  double variance=matricedelc2.getvariance();
				  matricedelc2.calculelc2local();
				  double lc2=matricedelc2.getlc2local();
				  //cout<<lc2<<endl;
					

				  //moyenne pondéré
				  variancesum+=variance;
				  variancelc2poids+=variance*lc2;
				  
  			}
		}
	}
	 //moyenne pondéré
	 lc2final=variancelc2poids/variancesum;
	cout<<"lc2final"<<endl;
	//cout<<lc2final<<endl;
	
}

double Lc2::getlc2(){
	return lc2final;
}