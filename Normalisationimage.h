#include "mitkImagePixelReadAccessor.h"
#include "mitkImagePixelWriteAccessor.h"
#include "mitkImage.h"
#include "mitkBaseData.h"
#include "mitkBaseGeometry.h"
#include "mitkGeometry3D.h"

#ifndef NORMALISATIONIMAGE_H
#define NORMALISATIONIMAGE_H

class Normalisationimage {
	public:
		Normalisationimage();
		void calculemax(mitk::Image::Pointer image, int modalite);
		void calculefacteur();
		float getfacteur();
		
		

	private:
		int maxirm_;
		int maxus_;
		float facteurconversion_;

};

#endif