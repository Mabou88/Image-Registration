#include "mitkImagePixelReadAccessor.h"
#include "mitkImagePixelWriteAccessor.h"
#include "mitkImage.h"
#include "mitkBaseData.h"
#include "mitkBaseGeometry.h"
#include "mitkGeometry3D.h"

typedef itk::Image<double,3> ImageType;
typedef itk::Index<3> Index;

#ifndef Matrice_H
#define Matrice_H

class Matricelc2 {
	public:
		Matricelc2(int nbrvoisins, Index index_irm, mitk::BaseGeometry::Pointer geo_irm, mitk::BaseGeometry::Pointer geo_us, mitk::Image::Pointer image1,mitk::Image::Pointer image2,float facteurconversion,mitk::Image::Pointer mitkgradient);
		void domatricetransverse();
		void matricemultipli();
		void matricemultipli2();
		void calculecoefficient();
		void domatriceinverse();
		void calculelc2local();
		void calculevariance();

				
		//accédeur
		double getmat();
		double getmattrans();
		double getmatavantinv();
		double getmatinv();
		double getvecteur();
		double getalpha();
		double getbeta();
		double getgamma();
		int getgrandeurpatch();
		double getlc2local();
		double getvariance();

		//affiche
		void affichemat();
		void affichevect();
		void affichemattrans();
		void affichematavantinv();
		void affichematinv();
		void affichematfinal();
		void affichecoef();

	private:
		//les matrices sont assez grosses pour contenir au maximum 3 voisins du pixel central
		double matrice[350][3];
		double matricetrans[3][350];
		double matriceavantinv[3][3];
		double matriceinv[3][3];
		double matricefinal[3][350];
		double vecteur[350];
		double coefficient[3];
		int grandeurpatch;
		double lc2local;
		double variance;
		float facteurconversion_;

};

#endif