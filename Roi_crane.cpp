#include "Roi_crane.h"
typedef itk::ImageDuplicator< ImageType > DuplicatorType;

Roi_crane::Roi_crane(ImageType::Pointer itkimageus ){
		cout<<"1"<<endl;
 // DuplicatorType::Pointer duplicator = DuplicatorType::New();
  cout<<"2"<<endl;
  //duplicator->SetInputImage(itkimageus);
  cout<<"3"<<endl;
  //duplicator->Update();
  cout<<"4"<<endl;
  roi_crane=ImageType::New();
  cout<<"5"<<endl;
  roi_crane = itkimageus;
  cout<<"6"<<endl;
  
  
}
ImageType::Pointer Roi_crane::getcraneROI(){
	return roi_crane;
}

void Roi_crane::sauvegardeimage(){
	//sauvegarde l'image modifier
	cout<<"sauvegarde de limage"<<endl;
  WriterType::Pointer   writer =  WriterType::New();
  itk::NiftiImageIO::Pointer ioimagenifti=itk::NiftiImageIO::New();

  writer->SetImageIO(ioimagenifti);
  writer->SetFileName( "C:/Users/Marc-Antoine/Documents/Imagecode/output/roi_crane_P2V2_cercle.nii");
  writer->SetInput(roi_crane);
  writer->Update();

}

void Roi_crane::setdim(){
	ImageType::SizeType dimension=roi_crane->GetLargestPossibleRegion().GetSize();
	dim_x=dimension[0];
	dim_y=dimension[1];
	dim_z=dimension[2];
	cout<< dim_x<<endl;
	cout<< dim_y<<endl;
	cout<< dim_z<<endl;
}



void Roi_crane::setlimitexy(int z){
	lim_min_x=dim_x;
	lim_max_x=0;
	lim_min_y=dim_y;
	lim_max_y=0;
	

	
	for (int i=0;i<dim_x;i++){
		for (int j=0;j<dim_y;j++){
			
				ImageType::IndexType index ={{i,j,z}};
				ImageType::PixelType pixelvalue=roi_crane->GetPixel(index);
				if (pixelvalue!=0){
					if (i>lim_max_x)
						lim_max_x=i;	
					if (i<lim_min_x)
						lim_min_x=i;
					if (j>lim_max_y)
						lim_max_y=j;			
					if (j<lim_min_y)
						lim_min_y=j;
					
										
				}

				
			
		}
	}


	cout<<"x"<<lim_min_x<<lim_max_x<<endl;
	cout<<"y"<<lim_min_y<<lim_max_y<<endl;
	
}

int Roi_crane::max(int a,int b){
	int max=0;
	if (a<b)
		max=b;
	if (a>b)
		max=a;
	if (a==b)
		max=a;
		
	return max;
}
void Roi_crane::setlimitepost_value(int a){

	limitepost_value=a;
}
//elimine les zones grises dans les coins de l'image
void Roi_crane::zonegrise(){
	//au dessus du crane
	for (int j=0;j<217;j++){
		for (int i=0;i<180;i++){
			//au dessus du crane
			for (int k=180;k>0;k--){
				itk::Index<3> curseur={{i,j,k}};
				ImageType::PixelType pixelvalue;
				pixelvalue=roi_crane->GetPixel(curseur);
				if (pixelvalue==0){
					break;
				}
				if (pixelvalue==100){
					roi_crane->SetPixel(curseur,0);
				}
				

			}
			//en-dessous
			for (int k=0;k<180;k++){
				itk::Index<3> curseur={{i,j,k}};
				ImageType::PixelType pixelvalue;
				pixelvalue=roi_crane->GetPixel(curseur);
				if (pixelvalue==0){
					break;
				}
				if (pixelvalue==100){
					roi_crane->SetPixel(curseur,0);
				}
				

			}
		}
	}
	
	//en-dessous
	for (int j=0;j<217;j++){
		for (int i=0;i<180;i++){
			for (int k=0;k<45;k++){
				itk::Index<3> curseur2={{i,j,k}};		
				roi_crane->SetPixel(curseur2,0);
			}
		}
	}
	//posterieur
	for (int j=0;j<40;j++){
		for (int i=0;i<181;i++){
			for (int k=0;k<181;k++){
				itk::Index<3> curseur2={{i,j,k}};		
				roi_crane->SetPixel(curseur2,0);
			}
		}
	}
	//anterieur
	for (int j=200;j<217;j++){
		for (int i=0;i<181;i++){
			for (int k=0;k<181;k++){
				itk::Index<3> curseur2={{i,j,k}};		
				roi_crane->SetPixel(curseur2,0);
			}
		}
	}
	//gauche droite
	for (int j=0;j<217;j++){
		for (int i=0;i<30;i++){
			for (int k=0;k<181;k++){
				itk::Index<3> curseur2={{i,j,k}};		
				roi_crane->SetPixel(curseur2,0);
			}
		}
	}
	for (int j=0;j<217;j++){
		for (int i=160;i<181;i++){
			for (int k=0;k<181;k++){
				itk::Index<3> curseur2={{i,j,k}};		
				roi_crane->SetPixel(curseur2,0);
			}
		}
	}

}


void Roi_crane::maskcrane(BinaryImageType::Pointer mask_us){
	
/*
 typedef itk::ImageConstIterator<ImageType> iterator;
 iterator it_m(mask_us,mask_us->GetBufferedRegion);

 it_m.GoToBegin();

 while(!it_m.IsAtEnd()){
    //on recupere le point dans l'espace correspondant pour aller le chercher dans l'IRM
     ImageType::IndexType index_m = it_m.GetIndex();
	 ImageType::PixelType pixelvalue_m;
	 pixelvalue_m=mask_us->GetPixel(index_m);
	 ImageType::PointType pt;
                     
     mask_us->TransformIndexToPhysicalPoint(index_m, pt);
	
	 ImageType::IndexType index_us;
	 roi_crane->TransformPhysicalPointToIndex(pt,index_us);
	  ImageType::PixelType pixelvalue_us;
	  pixelvalue_us=roi_crane->GetPixel(index_us);

	if (pixelvalue_m==1){
		roi_crane->SetPixel(index_us,0);
	}
	

	
	
 }
 */
	//si le masque superpose la partie on la garde l'image US sinon on l'élimine
	typedef itk::LinearInterpolateImageFunction<BinaryImageType> LinearInterpolatorFilterType;
	LinearInterpolatorFilterType::Pointer interpolator = LinearInterpolatorFilterType::New();
     interpolator->SetInputImage(mask_us);
	 
	ImageType::SizeType dimension=roi_crane->GetLargestPossibleRegion().GetSize();
	for (int i=0;i<dimension[0]-1;i++){
		for (int j=0;j<dimension[1]-1;j++){
			for (int k=0;k<dimension[2]-1;k++){
				ImageType::IndexType index_us ={{i,j,k}};
				 ImageType::PointType pt;
                     
				 roi_crane->TransformIndexToPhysicalPoint(index_us, pt);
	
				 
				  ImageType::PixelType pixelvalue_m;
				  ImageType::IndexType index_m;
				  if(mask_us->TransformPhysicalPointToIndex(pt,index_m)){
				  
				  
						  pixelvalue_m=mask_us->GetPixel(index_m);
						  //pixelvalue_m=interpolator->Evaluate(pt);

						if (pixelvalue_m!=1){
							roi_crane->SetPixel(index_us,0);
						}
				  }else{
					 // roi_crane->SetPixel(index_us,0);
				  }
			}
		}
	}

}

int Roi_crane::getlargueurcerveau(){
	return largueurcerveau;

}
void Roi_crane::calculvolumecrane(){

	//*******************************************
	//calcul du volume ellipsoide principal
	//*******************************************
	int min_x=dim_x;
	int max_x=0;
	int min_y=dim_y;
	int max_y=0;
	int min_z=dim_z;
	int max_z=0;
	int y_max_z;

	//pour le volume de l'ellipse principal
	for (int i=0;i<dim_x;i++){
		for (int j=0;j<dim_y;j++){
			for (int k=0;k<dim_z;k++){
				ImageType::IndexType index ={{i,j,k}};
				ImageType::PixelType pixelvalue=roi_crane->GetPixel(index);
				if (pixelvalue>10){
					


					//compteur dans chaque direction
					int c1=0;
					int c2=0;
					int c3=0;
					int c4=0;
					int c5=0;
					int c6=0;
					
					//compteur pour verifier que la surface est significative et non un point égaré
					for (int g=1;g<10;g++){
						ImageType::PixelType pixelvalue1=0;
						if(i<dim_x-10){
						ImageType::IndexType index1 ={{i+g,j,k}};
						pixelvalue1=roi_crane->GetPixel(index1);
							}
						ImageType::PixelType pixelvalue2=0;
						if(i>10){
						ImageType::IndexType index2 ={{i-g,j,k}};
						pixelvalue2=roi_crane->GetPixel(index2);
						}
						ImageType::PixelType pixelvalue3=0;
						if(j<dim_y-10){
						ImageType::IndexType index3 ={{i,j+g,k}};
						pixelvalue3=roi_crane->GetPixel(index3);
						}
						ImageType::PixelType pixelvalue4=0;
						if(j>10){
						ImageType::IndexType index4 ={{i,j-g,k}};
						pixelvalue4=roi_crane->GetPixel(index4);
						}
						ImageType::PixelType pixelvalue5=0;
						if(k<dim_z-10){
						ImageType::IndexType index5 ={{i,j,k+g}};
						pixelvalue5=roi_crane->GetPixel(index5);
						}
						ImageType::PixelType pixelvalue6=0;
						if(k>10){
						ImageType::IndexType index6 ={{i,j,k-g}};
						pixelvalue6=roi_crane->GetPixel(index6);
						}
					
							if(pixelvalue1>10)
							c1++;
							if(pixelvalue2>10)
							c2++;
							if(pixelvalue3>10)
							c3++;
							if(pixelvalue4>10)
							c4++;
							//tronc cérébrale vaut moins
							if(pixelvalue5>10)
							c5+=0.25;
							if(pixelvalue6>10)
							c6++;
						
						
					}
					int ct=c1+c2+c3+c4+c5+c6;
					if (ct>30){
						if (i>max_x)
							max_x=i;	
						if (i<min_x)
							min_x=i;
						if (j>max_y)
							max_y=j;			
						if (j<min_y)
							min_y=j;
						if (k>max_z){
							max_z=k;
							y_max_z=j;
						}
						if (k<min_z)
							min_z=k;
					}
										
				}

				
			}
		}
	}
	double centrage_haut_crane=double(y_max_z-min_y)/double(max_y-min_y);
	cout<<"y max z"<<y_max_z<<endl;
	cout<<"centrage"<<centrage_haut_crane<<endl;
	//*******************************************
	//calcul du volume ellipsoide principal
	//*******************************************
	ImageType::SpacingType spacing_us=roi_crane->GetSpacing();
	
	cout<<"limite x"<<min_x<<max_x<<endl;
	cout<<"limite y"<<min_y<<max_y<<endl;
	cout<<"limite z"<<min_z<<max_z<<endl;

	largueurcerveau=((max_x-min_x))*spacing_us[0];
	cout<<"largueur"<<largueurcerveau<<endl;

	int a=((max_x-min_x)/2)*spacing_us[0];
	int b=((max_y-min_y)/2)*spacing_us[1];
	int c=((max_z-min_z)/2)*spacing_us[2];

	cout<<"a"<<a<<"b"<<b<<"c"<<c<<endl;
	double fraction=1.333333;
	double pi=3.14159265359;
	double volume=fraction*pi*a*b*c;
	
	cout<<"volume"<<volume<<endl;
	
	//**********************************************
	//calcul des droites de la secondes ellipses
	//**********************************************
	
	ImageType::IndexType indexmin_z;
	ImageType::IndexType indexmax_y;

	int x_central=(min_x)+(max_x-min_x)/2;
	cout<<"x central"<<x_central<<endl;

	min_x=dim_x;
	max_x=0;
	min_y=dim_y;
	max_y=0;
	min_z=dim_z;
	max_z=0;

	//pour les valeurs sur la tranche sagittal central
	for (int i=92;i<101;i++){
		for (int j=0;j<dim_y;j++){
			for (int k=0;k<dim_z;k++){
				ImageType::IndexType index ={{x_central,j,k}};
				ImageType::PixelType pixelvalue=roi_crane->GetPixel(index);
				if (pixelvalue!=0){
					
					if (j>max_y){
						max_y=j;
						indexmax_y[0]=x_central;
						indexmax_y[1]=j;
						indexmax_y[2]=k;
					}
					
						
					
					if (k<min_z){
						min_z=k;
						indexmin_z[0]=x_central;
						indexmin_z[1]=j;
						indexmin_z[2]=k;
					}
					
				}

				
			}
		}
	}

	ImageType::IndexType point1=indexmax_y;
	ImageType::IndexType point2=indexmin_z;
	//Point 1 sera à y-1 au plus bas z
	point1[1]=indexmax_y[1]-1;

		
	for (int k=0;k<dim_z;k++){
				ImageType::IndexType index ={{point1[0],point1[1],k}};
				ImageType::PixelType pixelvalue=roi_crane->GetPixel(index);
				if (pixelvalue!=0){
									
					//prend le point le plus bas à y-1
					if (k<point1[2]){
						point1[2]=k;
					}
					
				}

				
	}
		
	
	cout<<"point1"<<point1<<endl;
	cout<<"point2"<<point2<<endl;

	int milieuz=(point1[2]-point2[2])/2+point2[2];
	int milieuy=(point1[1]-point2[1])/2+point2[1];
	
	cout<<"milieuz"<<milieuz<<endl;
	//int max_x_milieuz=dim_x;
	//int min_x_milieuz=0;

	//pour la largueur de la seconde ellipse
	for (int j=0;j<dim_y;j++){
			for (int i=0;i<dim_x;i++){
				ImageType::IndexType index ={{i,j,milieuz}};
				ImageType::PixelType pixelvalue=roi_crane->GetPixel(index);
				if (pixelvalue!=0){
					
					if (i>max_x){
						max_x=i;
						
					}
					
						
					
					if (i<min_x){
						min_x=i;
						
					}
					
				}

				
			}
	}
	int point3x=max_x;
	int point4x=min_x;
	cout<<"min x"<<min_x<<endl;
	cout<<"max x"<<max_x<<endl;
	//min max en x y et z


	

	//centre de l'ellipse
	ImageType::IndexType centre;
	centre[0]=min_x+a;
	centre[1]=min_y+b;
	centre[2]=min_z+c;
	//centre de la droite principal
	ImageType::IndexType centredroite;
	centredroite[0]=103;
	centredroite[1]=milieuy;
	centredroite[2]=milieuz;

	//creation des points
	ImageType::IndexType point5=centredroite;
	ImageType::IndexType point6;

	//2e ellipse
	//point 1 et 2
	int droiteprincipal=0;
	//point 3 et 4
	int droitelargueur=0;
	//point 5 et 6
	int petitedroite=0;
	
	

	//pente
	float pente=abs(double(point1[2]-point2[2]))/abs(double(point1[1]-point2[1]));
	//intersection entre ellipse et droite
	//ellipse centré à (122,77)

	//translation de -122,-77 pour que le centre soit a lorigine
	point5[2]=point5[2]-77;
	point5[1]=point5[1]-122;
	double m=-1/pente;
	double ord=point5[2]-m*point5[1];
	cout<<"ordonne"<<ord<<endl;
	cout<<"pente"<<m<<endl;
	cout<<"centre"<<centre<<endl;
	cout<<"centredroit"<<point5<<endl;

	//calcul des termes quadratiques
	//a
	double premier=b*b*m*m+c*c;
	//b
	double deuxieme=2*b*b*m*ord;
	//c
	double troisieme=b*b*(ord*ord-c*c);

	//equation quadratique
	point6[0]=x_central;
	point6[1]=(-deuxieme+sqrt(deuxieme*deuxieme-4*premier*troisieme))/(2*premier);
	point6[2]=m*point6[1]+ord;
	point6[1]=point6[1]+122;
	point6[2]=point6[2]+77;
	cout<<"point6"<<point6<<endl;
	//point 5 avant translation
	point5[1]=point5[1]+122;
	point5[2]=point5[2]+77;
	cout<<"point5"<<point5<<endl;
	
	//equation ellipse y,z
	//int point_y=sqrt((b*b)*(1-double((point_z-centre[2])*(point_z-centre[2]))/double(c*c)))+centre[1];

	

	cout<<point1<<point2<<endl;
	cout<<point3x<<point4x<<endl;
	cout<<milieuz<<endl;
	
	//calcul des droites de l'ellipsoide
	petitedroite=sqrt((point5[2]-point6[2])*(point5[2]-point6[2])+(point5[1]-point6[1])*(point5[1]-point6[1]));
	droiteprincipal=sqrt(abs(point1[1]-point2[1])*abs(point1[1]-point2[1])+abs(point1[2]-point2[2])*abs(point1[2]-point2[2]));
	droitelargueur=point3x-point4x;
	cout<<"droiteprincipal"<<droiteprincipal<<endl;
	cout<<"droitelargueur"<<droitelargueur<<endl;
	cout<<"petitedroite"<<petitedroite<<endl;
	//calcul du volume de l'ellipsoide à enlever
	double f2=0.125;
	double V2=fraction*pi*droitelargueur*droiteprincipal*petitedroite*f2;
	
	//volume du cerveau final
	double volumefinal=volume-V2;
	cout<<"volumefinal"<<volumefinal<<endl;
}
//elimine les reflets de l'ultrason
void Roi_crane::bordureexterieur(){
	//du niveau inferieur
	for (int i=0;i<181;i++){

		for (int j=0;j<211;j++){
			int noirceur=0;
			int estdansimage=0;

			for (int z=0;z<83;z++){
				itk::Index<3> curseurread={{i,j,z}};
				ImageType::PixelType pixelvalue;
				pixelvalue=roi_crane->GetPixel(curseurread);

				if (pixelvalue>10){
					estdansimage+=1;
					noirceur=0;
				}
				if (estdansimage>5){
					if (pixelvalue==0){
						noirceur+=1;
					}
				}
				//plus loin en y que 195, il n'y a que des reflets
				if (noirceur>3 || j>194){
		
		
						for (int k=z;k>0;k--){
									itk::Index<3> curseurwrite2={{i,j,k}};
									roi_crane->SetPixel(curseurwrite2,0);
							}
		
		
				}
			}
		}

	}
		
}
// appeler la fonction dans un for pour chaque z de limage
void Roi_crane::limiteparcercle(int z){
	
	//setlimitexy(z);
	
	// coord x et y du crane
	int cranecercle[10000][2];
	int cranecercleav[10000][2];
	int cranecerclepost[10000][2];

	//centre du crane
	int centre[2];
	double moitie=2;
	centre[0]=int(double(dim_x)/moitie);
	centre[1]=int(double(dim_y)/moitie);

	if(z<140){
		centre[0]=centre[0]+(150-z)/5;
		centre[1]=centre[1]+(150-z)/2;
	}

	cout<<"centre"<<centre[0]<<centre[1]<<endl;

	//pour le cercle
	
	int limitesuperieur=200-(z-60)/1.5;
	int limiteinferieur=200-(z-60)/1.5;
	cout<<"limite"<<limitesuperieur<<endl;

	//determination des longueurs dellipse du crane pour le calcul du rayon
	int largueur=130;
	int longueur=180;

	if (z<150){
		longueur=longueur-(150-z);
	}
	if (z<120){
		largueur=largueur-(120-z)/0.9;
	}
	
	

	
	//****************************************************
	//limite au niveau anterieur
	//****************************************************
	int position_a=0;
	for (int i=lim_min_x;i<lim_max_x;i++){
		

		//redefinition des limites superieurs et inferieurs
		int limitesuperieurant=limitesuperieur;
		int limiteinferieurant=limiteinferieur;
		if(z<85){
			limitesuperieurant=limitesuperieurant/1.25;
			limiteinferieurant=limiteinferieurant/1.25;
			//cout<<"limite"<<limitesuperieurant<<endl;
		}
		if(110<z){
			limitesuperieurant=limitesuperieurant*1.25;
			limiteinferieurant=limiteinferieurant*1.25;
			//cout<<"limite"<<limitesuperieurant<<endl;
		}
		
		bool limite0=false;
		int epaisseur0=0;
		for (int j=0;j<centre[1];j++){

			itk::Index<3> curseurread={{i,j,z}};
			ImageType::PixelType pixelvalue;
			pixelvalue=roi_crane->GetPixel(curseurread);

			if (pixelvalue>limitesuperieurant){
				
				epaisseur0+=1;
			}
			if (pixelvalue<limiteinferieurant){
				//limite=false;
				epaisseur0=0;
			}
			if (epaisseur0>3){
				
				limite0=true;
			}
			if (limite0==true && pixelvalue<limiteinferieurant){
				cranecercleav[position_a][0]=i;
				cranecercleav[position_a][1]=j;
				position_a+=1;
				epaisseur0=0;
				limite0=false;
				
			}
		}

	}

	//les coordonnees doivent rentrer dans un intervalle de rayon
	int cranecercle_av2[10000][2];
	int position_a2=0;
	
	cout<<"position a"<<position_a<<endl;
	for (int p=0;p<position_a;p++){
		//calcul du rayon
		double rayon=(double((cranecercleav[p][0]-centre[0])*(cranecercleav[p][0]-centre[0]))/double(largueur*largueur)+double((cranecercleav[p][1]-centre[1])*(cranecercleav[p][1]-centre[1]))/double((longueur)*(longueur)));
		if(z<148 && z>126){
			int longueur_ant=longueur/((148-z)*0.025+1);
			rayon=(double((cranecercleav[p][0]-centre[0])*(cranecercleav[p][0]-centre[0]))/double(largueur*largueur)+double((cranecercleav[p][1]-centre[1])*(cranecercleav[p][1]-centre[1]))/double((longueur_ant)*(longueur_ant)));
			cout<<"rayon ant"<<rayon<<endl;
		}
		if(z<127 && z>120){
			int longueur_ant=longueur/((127-z)*0.03+1.5);
			rayon=(double((cranecercleav[p][0]-centre[0])*(cranecercleav[p][0]-centre[0]))/double(largueur*largueur)+double((cranecercleav[p][1]-centre[1])*(cranecercleav[p][1]-centre[1]))/double((longueur_ant)*(longueur_ant)));
			cout<<"rayon ant"<<rayon<<endl;
		}
		if(z<121 && z>110){
			int longueur_ant=longueur/((121-z)*0.03+1.65);
			rayon=(double((cranecercleav[p][0]-centre[0])*(cranecercleav[p][0]-centre[0]))/double(largueur*largueur)+double((cranecercleav[p][1]-centre[1])*(cranecercleav[p][1]-centre[1]))/double((longueur_ant)*(longueur_ant)));
			cout<<"rayon ant"<<rayon<<endl;
		}
		if(z<111 && z>89){
			int longueur_ant=longueur/((110-z)*0.025+1.95);
			rayon=(double((cranecercleav[p][0]-centre[0])*(cranecercleav[p][0]-centre[0]))/double(largueur*largueur)+double((cranecercleav[p][1]-centre[1])*(cranecercleav[p][1]-centre[1]))/double((longueur_ant)*(longueur_ant)));
			cout<<"rayon ant"<<rayon<<endl;
		}
		if(z<90){
			rayon=(double((cranecercleav[p][0]-centre[0])*(cranecercleav[p][0]-centre[0]))/double(largueur/1.25*largueur/1.25)+double((cranecercleav[p][1]-centre[1])*(cranecercleav[p][1]-centre[1]))/double((longueur/3.5)*(longueur/3.5)));
			cout<<"rayon ant"<<rayon<<endl;
		}
		if(z<72){
			rayon=(double((cranecercleav[p][0]-centre[0])*(cranecercleav[p][0]-centre[0]))/double(largueur/1.25*largueur/1.25)+double((cranecercleav[p][1]-centre[1])*(cranecercleav[p][1]-centre[1]))/double((longueur/3.5)*(longueur/3.5)));
			cout<<"rayon ant"<<rayon<<endl;
		}
		//cout<<"rayon"<<rayon<<endl;
		if (0.75<rayon && rayon<20){
			cranecercle_av2[position_a2][0]=cranecercleav[p][0];
			cranecercle_av2[position_a2][1]=cranecercleav[p][1];
			position_a2+=1;
			
		}
	}
	//ajuste ceux qui rentre dans les limites du rayon
	cout<<"position 2a"<<position_a2<<endl;

	for (int a=0;a<position_a2;a++){
		int x=cranecercle_av2[a][0];
		int y=cranecercle_av2[a][1];

		
		itk::Index<3> curseur={{x,y,z}};
		roi_crane->SetPixel(curseur,0);
		
		
		
			for (int n=0;n<y;n++){
				itk::Index<3> curseurwrite={{x,n,z}};
				//cout<<"x"<<n<<endl;
				roi_crane->SetPixel(curseurwrite,0);
				
			}
			
			
	}
	//****************************************************
	//limite au niveau postérieur
	//****************************************************
	int positionpost=0;
	for (int i=lim_min_x;i<lim_max_x;i++){
		bool limitepost=false;
		int epaisseurpost=0;
		int limitepost_y=centre[1];

		int limitesuperieur_p=limitesuperieur;
		int limiteinferieur_p=limiteinferieur;


		if(z>=150){
			limitesuperieur_p=limitesuperieur/1.25;
			limiteinferieur_p=limiteinferieur/1.25;

		}

		for (int j=lim_max_y;j>limitepost_y;j--){

			itk::Index<3> curseurread={{i,j,z}};
			ImageType::PixelType pixelvalue;
			pixelvalue=roi_crane->GetPixel(curseurread);

			if (pixelvalue>limitesuperieur_p){
				
				epaisseurpost+=1;
			}
			if (pixelvalue<limiteinferieur_p){
				//limite=false;
				epaisseurpost=0;
			}
			if (epaisseurpost>3){
				
				limitepost=true;
			}
			if (limitepost==true && pixelvalue<limiteinferieur_p){
				cranecerclepost[positionpost][0]=i;
				cranecerclepost[positionpost][1]=j;
				positionpost+=1;
				epaisseurpost=0;
				limitepost=false;
				
			}
		}

	}

	//les coordonnees doivent rentrer dans un intervalle de rayon
	int cranecerclepost2[10000][2];
	int positionpost2=0;
	
	
	for (int p=0;p<positionpost;p++){
		//calcul du rayon
		double rayon=(double((cranecerclepost[p][0]-centre[0])*(cranecerclepost[p][0]-centre[0]))/double(largueur*largueur)+double((cranecerclepost[p][1]-centre[1])*(cranecerclepost[p][1]-centre[1]))/double((longueur)*(longueur)));
		if(z<80){
			rayon=(double((cranecerclepost[p][0]-centre[0])*(cranecerclepost[p][0]-centre[0]))/double((largueur)*(largueur))+double((cranecerclepost[p][1]-centre[1])*(cranecerclepost[p][1]-centre[1]))/double((longueur/2.6)*(longueur/2.6)));
			//cout<<"rayon"<<rayon<<endl;		
		}
		//cout<<"rayon"<<rayon<<endl;
		if (0.85<rayon && rayon<20){
			cranecerclepost2[positionpost2][0]=cranecerclepost[p][0];
			cranecerclepost2[positionpost2][1]=cranecerclepost[p][1];
			positionpost2+=1;
			
		}
	}
	//ajuste ceux qui rentre dans les limites du rayon
	cout<<"longueur"<<longueur<<endl;

	for (int a=0;a<positionpost2;a++){
		int x=cranecerclepost2[a][0];
		int y=cranecerclepost2[a][1];

		
		itk::Index<3> curseur={{x,y,z}};
		roi_crane->SetPixel(curseur,0);
		
		
		
		for (int n=lim_max_y;n>y;n--){
				itk::Index<3> curseurwrite={{x,n,z}};
				//cout<<"x"<<n<<endl;
				roi_crane->SetPixel(curseurwrite,0);
				
			}
			
			
	}
	//***********************************************
	// par la gauche et la droite
	//***********************************************
	int position=0;
	for (int j=lim_min_y;j<lim_max_y;j++){
		//1ere moitie
		bool limite=false;
		int epaisseur=0;
		
		int limitesuperieur_gd=limitesuperieur;
		int limiteinferieur_gd=limiteinferieur;
		
		if(z>140){
			int limitesuperieur_gd=limitesuperieur/1.75;
			int limiteinferieur_gd=limiteinferieur/1.75;

		}
		
		for (int i=lim_max_x;i>centre[0];i--){
			
			
			
			itk::Index<3> curseurread={{i,j,z}};
			ImageType::PixelType pixelvalue;
			pixelvalue=roi_crane->GetPixel(curseurread);

			if (pixelvalue>limitesuperieur_gd){
				
				epaisseur+=1;
			}
			if (pixelvalue<limiteinferieur_gd){
				//limite=false;
				epaisseur=0;
			}
			if (epaisseur>3){
				
				limite=true;
			}
			if (limite==true && pixelvalue<limiteinferieur_gd){
				cranecercle[position][0]=i;
				cranecercle[position][1]=j;
				position+=1;
				epaisseur=0;
				limite=false;
				
			}


		
		}

		bool limite2=false;
		int epaisseur2=0;
		//2e moitie

		for (int i2=lim_min_x;i2<centre[0];i2++){
			
			
			
			itk::Index<3> curseurread={{i2,j,z}};
			ImageType::PixelType pixelvalue;
			pixelvalue=roi_crane->GetPixel(curseurread);

			if (pixelvalue>limitesuperieur)
				epaisseur2+=1;
			if (pixelvalue<limiteinferieur){
				//limite=false;
				epaisseur2=0;
			}
			if (epaisseur2>3)
				limite2=true;
			if (limite2==true && pixelvalue<limiteinferieur){
				cranecercle[position][0]=i2;
				cranecercle[position][1]=j;
				position+=1;
				epaisseur2=0;
				limite2=false;
			}
		
		}
	}
	cout<<"position"<<position<<endl;
	//cout<<"calcul du rayon"<<endl;

	//les coordonnees doivent rentrer dans un intervalle de rayon
	int cranecercle_2[10000][2];
	int position_2=0;
	
	
	for (int p=0;p<position;p++){
		//calcul du rayon
		
		double rayon=(double((cranecercle[p][0]-centre[0])*(cranecercle[p][0]-centre[0]))/double(largueur*largueur)+double((cranecercle[p][1]-centre[1])*(cranecercle[p][1]-centre[1]))/double(longueur*longueur));
		if (z<100){
			rayon=(double((cranecercle[p][0]-centre[0])*(cranecercle[p][0]-centre[0]))/double((largueur/1.1)*(largueur/1.1))+double((cranecercle[p][1]-centre[1])*(cranecercle[p][1]-centre[1]))/double((longueur/1.2)*(longueur/1.2)));
		}
		if (z<70){
			rayon=(double((cranecercle[p][0]-centre[0])*(cranecercle[p][0]-centre[0]))/double((largueur/1.5)*(largueur/1.5))+double((cranecercle[p][1]-centre[1])*(cranecercle[p][1]-centre[1]))/double((longueur/1.25)*(longueur/1.25)));
		}
		//cout<<"rayon"<<rayon<<endl;
		if (0.85<rayon && rayon<20){
			cranecercle_2[position_2][0]=cranecercle[p][0];
			cranecercle_2[position_2][1]=cranecercle[p][1];
			position_2+=1;
			
		}
	}
	
	//cout<<"position 2"<<position_2<<endl;
	//cout<<"ajustement des limites"<<endl;


	// les coordonnees restantes sont utilisées pour faire le masque
	for (int a=0;a<position_2;a++){
		int x=cranecercle_2[a][0];
		int y=cranecercle_2[a][1];

		
		itk::Index<3> curseur={{x,y,z}};
		roi_crane->SetPixel(curseur,0);
		

		if (x<centre[0]){
			for (int n=x;n>lim_min_x;n--){
				itk::Index<3> curseurwrite={{n,y,z}};
				//cout<<"x"<<n<<endl;
				roi_crane->SetPixel(curseurwrite,0);
				
			}
		}
		if (x>centre[0]){
			for (int n=x;n<lim_max_x+1;n++){
				itk::Index<3> curseurwrite={{n,y,z}};
				
				roi_crane->SetPixel(curseurwrite,0);
				
			}
		}
		
	}
	cout<<"z"<<z<<"L"<<longueur<<"y_centre"<<centre[1]<<endl;

}

void Roi_crane::regionexterieur(int z){
	int centre[2];
	double moitie=2;
	centre[0]=int(double(dim_x)/moitie);
	centre[1]=int(double(dim_y)/moitie);

	int a=180-(180-z)*1.75;
	int b=130-(180-z)/2;
	//pour les tranches sur le côtés de l'image
	
	//par les côtés en x
	for (int j=0;j<dim_y;j++){
		for (int i=0;i<centre[0]-b/1.3;i++){
		itk::Index<3> curseurread={{i,j,z}};
		ImageType::PixelType pixelvalue;
		pixelvalue=roi_crane->GetPixel(curseurread);
			if (pixelvalue==0){
				for (int ii=i;ii>=0;ii--){
					itk::Index<3> curseurwrite={{ii,j,z}};
					roi_crane->SetPixel(curseurwrite,0);
				}
			}

		}
		for (int i=lim_max_x;i>centre[0]+b/1.3;i--){
		itk::Index<3> curseurread={{i,j,z}};
		ImageType::PixelType pixelvalue;
		pixelvalue=roi_crane->GetPixel(curseurread);
			if (pixelvalue==0){
				for (int ii=i;ii<lim_max_x;ii++){
					itk::Index<3> curseurwrite={{ii,j,z}};
					roi_crane->SetPixel(curseurwrite,0);
				}
			}

		}
		itk::Index<3> curseurread={{lim_max_x-1,j,z}};
		ImageType::PixelType pixelvalue;
		pixelvalue=roi_crane->GetPixel(curseurread);
		if(pixelvalue==0){
			//itk::Index<3> curseurwrite={{lim_max_x,j,z}};
			//roi_crane->SetPixel(curseurwrite,0);
		}
	}
	
	
	//par les cotés en y
	
	for (int i=lim_min_x;i<lim_max_x;i++){
		for (int j=0;j<centre[1]/1.25;j++){
		itk::Index<3> curseurread0={{i,j,z}};
		ImageType::PixelType pixelvalue;
		pixelvalue=roi_crane->GetPixel(curseurread0);
			if (pixelvalue==0){
				for (int jj=j;jj>=0;jj--){
					itk::Index<3> curseurwrite0={{i,jj,z}};
					roi_crane->SetPixel(curseurwrite0,0);
				}
			}

		}
	}

	for (int i=0;i<lim_max_x;i++){
		for (int j=lim_max_y;j>centre[1]*1.25;j--){
		itk::Index<3> curseurread00={{i,j,z}};
		ImageType::PixelType pixelvalue00;
		pixelvalue00=roi_crane->GetPixel(curseurread00);
			if (pixelvalue00==0){
				for (int jj=j;jj<lim_max_y;jj++){
					itk::Index<3> curseurwrite00={{i,jj,z}};
					roi_crane->SetPixel(curseurwrite00,0);
				}
			}

		}
	}
	
	if (z<90){
		//par les côtés en x
		for (int j=0;j<dim_y;j++){
			for (int i=0;i<centre[0]/1.1;i++){
			itk::Index<3> curseurread={{i,j,z}};
			ImageType::PixelType pixelvalue;
			pixelvalue=roi_crane->GetPixel(curseurread);
				if (pixelvalue==0){
					for (int ii=i;ii>=0;ii--){
						itk::Index<3> curseurwrite={{ii,j,z}};
						roi_crane->SetPixel(curseurwrite,0);
					}
				}

			}
			for (int i=lim_max_x;i>centre[0]*1.1;i--){
			itk::Index<3> curseurread={{i,j,z}};
			ImageType::PixelType pixelvalue;
			pixelvalue=roi_crane->GetPixel(curseurread);
				if (pixelvalue==0){
					for (int ii=i;ii<lim_max_x;ii++){
						itk::Index<3> curseurwrite={{ii,j,z}};
						roi_crane->SetPixel(curseurwrite,0);
					}
				}

			}
		}
		//par les cotés en y
	
		for (int i=lim_min_x;i<lim_max_x+1;i++){
			for (int j=0;j<centre[1]/1.1;j++){
			itk::Index<3> curseurread0={{i,j,z}};
			ImageType::PixelType pixelvalue;
			pixelvalue=roi_crane->GetPixel(curseurread0);
				if (pixelvalue==0){
					for (int jj=j;jj>=0;jj--){
						itk::Index<3> curseurwrite0={{i,jj,z}};
						roi_crane->SetPixel(curseurwrite0,0);
					}
				}

			}
		}

		for (int i=0;i<lim_max_x+1;i++){
			for (int j=lim_max_y;j>centre[1]*1.1;j--){
			itk::Index<3> curseurread00={{i,j,z}};
			ImageType::PixelType pixelvalue00;
			pixelvalue00=roi_crane->GetPixel(curseurread00);
				if (pixelvalue00==0){
					for (int jj=j;jj<lim_max_y;jj++){
						itk::Index<3> curseurwrite00={{i,jj,z}};
						roi_crane->SetPixel(curseurwrite00,0);
					}
				}

			}
		}
	}

	//Pour les tranches au milieu qui ont encore leur intensite encadré par du noir
	if(z<150){

		//limite anterieur
		for (int i=lim_min_x;i<lim_max_x;i++){
			//bool trou=false;
			//bool bord=false;
			//int bordy=centre[1];
			int longueurtrou=0;
			int largueurtrou=60+(180-z)/4;
			for (int j=0;j<lim_max_y;j++){

				if(j>lim_min_y+40 && longueurtrou<10){
				break;
				}
				itk::Index<3> curseurread0={{i,j,z}};
				ImageType::PixelType pixelvalue;
				pixelvalue=roi_crane->GetPixel(curseurread0);
				if(pixelvalue!=0){
					longueurtrou++;
					bool trouaccepte=false;
					for (int ii=0;ii<largueurtrou/2;ii++){
						itk::Index<3> curseurread1={{i+ii,j,z}};
						ImageType::PixelType pixelvalue1;
						pixelvalue1=roi_crane->GetPixel(curseurread1);
						itk::Index<3> curseurread2={{i-ii,j,z}};
						ImageType::PixelType pixelvalue2;
						pixelvalue2=roi_crane->GetPixel(curseurread2);
						if (pixelvalue1==0 && pixelvalue2==0){
							trouaccepte=true;
						}
					}
					if(trouaccepte){
						for (int jj=0;jj<largueurtrou/2;jj++){
						itk::Index<3> curseurwrite0={{i+jj,j,z}};
						roi_crane->SetPixel(curseurwrite0,0);
						itk::Index<3> curseurwrite1={{i-jj,j,z}};
						roi_crane->SetPixel(curseurwrite1,0);
						}


					}
				}
			}
		}

		//limite posterieur
		for (int i=lim_min_x;i<lim_max_x;i++){
			//bool trou=false;
			//bool bord=false;
			//int bordy=centre[1];
			int longueurtrou=0;
			for (int j=lim_max_y;j>0;j--){

				if(j<lim_max_y-30 && longueurtrou<10){
				break;
				}
				itk::Index<3> curseurread0={{i,j,z}};
				ImageType::PixelType pixelvalue;
				pixelvalue=roi_crane->GetPixel(curseurread0);
				if(pixelvalue!=0){
					longueurtrou++;
					bool trouaccepte=false;
					for (int ii=0;ii<30;ii++){
						itk::Index<3> curseurread1={{i+ii,j,z}};
						ImageType::PixelType pixelvalue1;
						pixelvalue1=roi_crane->GetPixel(curseurread1);
						itk::Index<3> curseurread2={{i-ii,j,z}};
						ImageType::PixelType pixelvalue2;
						pixelvalue2=roi_crane->GetPixel(curseurread2);
						if (pixelvalue1==0 && pixelvalue2==0){
							trouaccepte=true;
						}
					}
					if(trouaccepte){
						for (int jj=0;jj<30;jj++){
						itk::Index<3> curseurwrite0={{i+jj,j,z}};
						roi_crane->SetPixel(curseurwrite0,0);
						itk::Index<3> curseurwrite1={{i-jj,j,z}};
						roi_crane->SetPixel(curseurwrite1,0);
						}


					}
				}
			}
		}
		
		//limite droite
		for (int j=lim_min_y;j<lim_max_y;j++){
			//bool trou=false;
			//bool bord=false;
			//int bordy=centre[1];
			int longueurtrou=0;
			for (int i=lim_max_x;i>lim_min_x;i--){

				if(i<lim_max_x-30 && longueurtrou<10){
				break;
				}
				itk::Index<3> curseurread0={{i,j,z}};
				ImageType::PixelType pixelvalue;
				pixelvalue=roi_crane->GetPixel(curseurread0);
				if(pixelvalue!=0){
					longueurtrou++;
					bool trouaccepte=false;
					for (int ii=0;ii<40;ii++){
						itk::Index<3> curseurread1={{i,j+ii,z}};
						ImageType::PixelType pixelvalue1;
						pixelvalue1=roi_crane->GetPixel(curseurread1);
						itk::Index<3> curseurread2={{i,j-ii,z}};
						ImageType::PixelType pixelvalue2;
						pixelvalue2=roi_crane->GetPixel(curseurread2);
						if (pixelvalue1==0 && pixelvalue2==0){
							trouaccepte=true;
						}
					}
					if(trouaccepte){
						for (int jj=0;jj<40;jj++){
						itk::Index<3> curseurwrite0={{i,j+jj,z}};
						roi_crane->SetPixel(curseurwrite0,0);
						itk::Index<3> curseurwrite1={{i,j-jj,z}};
						roi_crane->SetPixel(curseurwrite1,0);
						}
					}


				}
			}
		}
		

		//limite gauche
		for (int j=lim_min_y;j<lim_max_y;j++){
			//bool trou=false;
			//bool bord=false;
			//int bordy=centre[1];
			int longueurtrou=0;
			for (int i=0;i<lim_max_x;i++){

				if(i>lim_min_x+30 && longueurtrou<10){
				break;
				}
				itk::Index<3> curseurread0={{i,j,z}};
				ImageType::PixelType pixelvalue;
				pixelvalue=roi_crane->GetPixel(curseurread0);
				if(pixelvalue!=0){
					longueurtrou++;
					bool trouaccepte=false;
					for (int ii=0;ii<40;ii++){
						itk::Index<3> curseurread1={{i,j+ii,z}};
						ImageType::PixelType pixelvalue1;
						pixelvalue1=roi_crane->GetPixel(curseurread1);
						itk::Index<3> curseurread2={{i,j-ii,z}};
						ImageType::PixelType pixelvalue2;
						pixelvalue2=roi_crane->GetPixel(curseurread2);
						if (pixelvalue1==0 && pixelvalue2==0){
							trouaccepte=true;
						}
					}
					if(trouaccepte){
						for (int jj=0;jj<40;jj++){
						itk::Index<3> curseurwrite0={{i,j+jj,z}};
						roi_crane->SetPixel(curseurwrite0,0);
						itk::Index<3> curseurwrite1={{i,j-jj,z}};
						roi_crane->SetPixel(curseurwrite1,0);
						}
					}


				}
			}
		}
	}
		
		
	
		
		
	
	

	
	

}

void Roi_crane::regioninferieur(int k){
	for (int i=0;i<dim_x;i++){
		for (int j=0;j<dim_y;j++){
			for (int z=0;z<k;z++){
				itk::Index<3> curseurwrite0={{i,j,z}};
				roi_crane->SetPixel(curseurwrite0,0);
			}

		}

	}
	for (int i=400;i<dim_x;i++){
		for (int j=0;j<dim_y;j++){
			for (int z=0;z<dim_z;z++){
				itk::Index<3> curseurwrite0={{i,j,z}};
				roi_crane->SetPixel(curseurwrite0,0);
			}

		}

	}

}

void Roi_crane::regionssuperieur(){

	
	int min_z=dim_z;
	int max_z=0;
	

	//pour la hauteur de la limite superieur du crane
	for (int i=0;i<dim_x;i++){
		for (int j=0;j<dim_y;j++){
			for (int k=0;k<dim_z;k++){
				ImageType::IndexType index ={{i,j,k}};
				ImageType::PixelType pixelvalue=roi_crane->GetPixel(index);
				if (pixelvalue>10){
														
						
						if (k>max_z)
							max_z=k;
							
						if (k<min_z)
							min_z=k;
					
										
				}

				
			}
		}
	}
	cout<<"max z hauteur"<<max_z<<endl;
	ImageType::SpacingType spacingus=roi_crane->GetSpacing();

	int hauteurindex=double(8.5)/(double(spacingus[2]));

	cout<<"hauteur index"<<hauteurindex<<endl;
	int nouvellehauteur=max_z-hauteurindex;

	
	setlimitexy(nouvellehauteur);
	int centre_y=(lim_min_y+lim_max_y)/2;

	//en ligne droite
	for (int z=max_z;z>nouvellehauteur;z--){
		for (int i=lim_min_x;i<lim_max_x;i++){
			for (int j=lim_min_y;j<lim_max_y;j++){
			
				
				ImageType::IndexType index ={{i,j,z}};
				roi_crane->SetPixel(index,0);
				

				
			
			}
		}

	}
	//pour arrondir les coins
	
	for (int i=lim_min_x-20;i<lim_max_x+20;i++){
		for (int j=lim_min_y-20;j<lim_max_y+20;j++){
			int différence=abs(j-centre_y)/10;
			int nouvellehauteurcourbe=nouvellehauteur-différence;
			for (int z=nouvellehauteur;z>nouvellehauteurcourbe;z--){
				
				ImageType::IndexType index ={{i,j,z}};
				roi_crane->SetPixel(index,0);
				

				
			
			}
		}

	}
	
}

void Roi_crane::limiteposterieur(){
	
	//variable pour une trache axiale
	int limite_y[200];
	int limite_z[200];
	int position=0;
	
	for (int z=30;z<135;z++){

	
		
	//variable pour définir la limite du crane
	bool limite=false;
	//int limite_y[200];
	int epaisseur=0;
	int coord=0;
	//int position;


		//Analyse la partie posterieur
		
		for (int a=0;a<90;a++){
			itk::Index<3> curseurread={{90,a,z}};
			ImageType::PixelType pixelvalue;
			pixelvalue=roi_crane->GetPixel(curseurread);

			//on rencontre le crâne
				if (pixelvalue>160){
					//compteur de l'epaisseur du crane
					epaisseur+=1;
					//verifier que les pixels du crane sont collé
					if (epaisseur>1){
						if (a-coord>1){
							epaisseur=0;
						}
					}
					coord=a;
					
					
					
				}
				
				//si la distance est de 5 au plus, c'est le crane
				if (epaisseur>3){
					limite=true;
				}
				
				//lorsqu'on sort de l'épaisseur du crane
				if (limite==true && pixelvalue<141){
					
					//mettre à 0 les pixels extérieurs et ceux du crâne
						for (int y=a;y>0;y--){
							for (int x=0;x<181;x++){
								itk::Index<3> curseurwrite={{x,y,z}};
								roi_crane->SetPixel(curseurwrite,0);
							}
						}
						limite_y[position]=a;
						limite_z[position]=z;
						position+=1;
					
					break;
				}
		
		}

	//bouche les trous
	for (int i=0;i<199;i++){
		int difference=limite_z[i+1]-limite_z[i];
		//s'il y a des trous entre les lignes noires
		if (difference>1){
			
			for (int k=limite_z[i]+1;k<limite_z[i+1];k++){
				//y est le milieu de deux lignes a 0 a ajuster!!!!!
				int y=(limite_y[i]+limite_y[i+1])/2;
				//on va les remplir de 0 (noire) aussi
				for (int j=y;j>0;j--){
					for (int x=0;x<181;x++){
							itk::Index<3> curseurwrite2={{x,j,k}};
							roi_crane->SetPixel(curseurwrite2,0);
					}
				}

			}
		}


	}
		
							
	}
}



void Roi_crane::limitedroit(int z){
	
		
	//variable pour une trache axiale
	int limite_y[200];
	int limite_x_droite[200];
	int position=0;

	
	for (int y=limitepost_value;y<200;y++){
		//Analyse la partie droite
		
		//variable pour une tranche coronale
		
		//int estdansimage=0;
		//int noirceur=0;
		bool limite=false;
		int epaisseur_droite=0;
		
	    int trou_coor=0;
		int trou=0;
		int coord=0;

		bool coin=false;

		int patchdelete=0;
		int compteurnoir=0;

		for (int x=0;x<90;x++){
			itk::Index<3> curseurread={{x,y,z}};
			ImageType::PixelType pixelvalue;
			pixelvalue=roi_crane->GetPixel(curseurread);

			//on rencontre le crâne
				if (pixelvalue>150){
					//compteur de l'epaisseur du crane			
					epaisseur_droite+=1;
					//verifier que les pixels du crane sont collé
					if (epaisseur_droite>1){
						if (x-coord>1){
							epaisseur_droite=0;
						}
					}
					coord=x;

					
														
				}
				if (z<90 && x<55)
					coin=true;
				if(z<100 && y>145)
					coin=true;
				//eliminer les rectangles d'intensité dans les coins (zone d'intensité entre les zones noires)
				if (coin){
					if (pixelvalue>10){
						patchdelete+=1;
					}
					if (patchdelete>4 && pixelvalue==0){
						compteurnoir+=1;
					}
					if (compteurnoir>4){
								//mettre à 0 les pixels extérieurs et ceux du crâne
								for (int i=x;i>0;i--){
									itk::Index<3> curseurwrited1={{i,y,z}};
									roi_crane->SetPixel(curseurwrited1,0);
								}
								//reset
								compteurnoir=0;
								patchdelete=0;
								break;
					}
				}
				
				//si la distance est de 5 au plus, c'est le crane
				if (epaisseur_droite>3){
					limite=true;
				}
				
				//lorsqu'on sort de l'épaisseur du crane
				if (limite==true && pixelvalue<141){
					
					
					//dependamment si cest la première limite de la slice
					if (position!=0){

						int différence_y=y-limite_y[position-1];
						int différence_x=x-limite_x_droite[position-1];
						//verifie que les limites sont collées
						if (différence_x<10 && différence_y<10){
							limite_x_droite[position]=x;
							limite_y[position]=y;

							//mettre à 0 les pixels extérieurs et ceux du crâne
							for (int i=x;i>0;i--){
								itk::Index<3> curseurwrite1={{i,y,z}};
								roi_crane->SetPixel(curseurwrite1,0);
							}
							position+=1;
							
						}
						break;
					//si cest la première limite de la slice
					}else{
						limite_x_droite[position]=x;
						limite_y[position]=y;
						
						//mettre à 0 les pixels extérieurs et ceux du crâne
						for (int i=x;i>0;i--){
							for (int i=x;i>0;i--){
								itk::Index<3> curseurwrite2={{i,y,z}};
								roi_crane->SetPixel(curseurwrite2,0);
							}
						}
						position+=1;
						break;

					}
					
									
											
									     				 
				
				}
			
		}
		
	}
	//reevalue toutes les lignes mises a zero
	for (int i=0;i<199;i++){
		int difference=limite_y[i+1]-limite_y[i];
		//s'il y a des trous entre les lignes noires
		if (difference>1){
			for (int j=limite_y[i]+1;j<limite_y[i+1];j++){
				//x est le milieu de deux lignes a 0 a ajuster!!!!!
				int x=(limite_x_droite[i]+limite_x_droite[i+1])/2;
				//on va les remplir de 0 (noire) aussi
				for (int i=x;i>0;i--){
							itk::Index<3> curseurwrite2={{i,j,z}};
							roi_crane->SetPixel(curseurwrite2,0);
				}

			}
		}


	}

}
/*

void Roi_crane::limitedroit(int z){
	
		
	//variable pour une trache axiale
	int limite_y[200];
	int limite_x_droite[200];
	int position=0;

	
	for (int y=limitepost_value;y<200;y++){
		//Analyse la partie droite
		
		//variable pour une tranche coronale
		
		//int estdansimage=0;
		//int noirceur=0;
		bool limite=false;
		int epaisseur_droite=0;
		
	    int trou_coor=0;
		int trou=0;
		int coord=0;

		int patchdelete=0;
		int compteurnoir=0;

		for (int x=0;x<100;x++){
			itk::Index<3> curseurread={{x,y,z}};
			ImageType::PixelType pixelvalue;
			pixelvalue=roi_crane->GetPixel(curseurread);

			//on rencontre le crâne
				if (pixelvalue>150){
										
					epaisseur_droite+=1;
					//verifier que les pixels du crane sont collé
					if (epaisseur_droite>1){
						if (x-coord>1){
							epaisseur_droite=0;
						}
					}
					coord=x;

					
														
				}

				//eliminer les rectangles d'intensité dans les coins
				if (z<80){
					if (pixelvalue>30){
						patchdelete+=1;
					}
					if (patchdelete>10){
						compteurnoir+=1;
					}
					if (compteurnoir>15){
								//mettre à 0 les pixels extérieurs et ceux du crâne
								for (int i=x;i>0;i--){
									itk::Index<3> curseurwrite11={{i,y,z}};
									roi_crane->SetPixel(curseurwrite11,0);
								}
					}
				}
				
				//si la distance est de 5 au plus, c'est le crane
				if (epaisseur_droite>3){
					limite=true;
				}
				
				//lorsqu'on sort de l'épaisseur du crane
				if (limite==true && pixelvalue<141){
					
					
					//dependamment si cest la première limite de la slice
					if (position!=0){

						int différence_y=y-limite_y[position-1];
						int différence_x=x-limite_x_droite[position-1];
						//verifie que les limites sont collées
						if (différence_x<10 && différence_y<10){
							limite_x_droite[position]=x;
							limite_y[position]=y;

							//mettre à 0 les pixels extérieurs et ceux du crâne
							for (int i=x;i>0;i--){
								itk::Index<3> curseurwrite1={{i,y,z}};
								roi_crane->SetPixel(curseurwrite1,0);
							}
							position+=1;
							
						}
						break;
					//si cest la première limite de la slice
					}else{
						limite_x_droite[position]=x;
						limite_y[position]=y;
						
						//mettre à 0 les pixels extérieurs et ceux du crâne
						for (int j=y;j>0;j--){
							for (int i=x;i>0;i--){
								itk::Index<3> curseurwrite2={{i,j,z}};
								roi_crane->SetPixel(curseurwrite2,0);
							}
						}
						position+=1;
						break;

					}
					
									
											
									     				 
				
				}
			
		}
		
	}
	//reevalue toutes les lignes mises a zero
	for (int i=0;i<199;i++){
		int difference=limite_y[i+1]-limite_y[i];
		//s'il y a des trous entre les lignes noires
		if (difference>1){
			for (int j=limite_y[i]+1;j<limite_y[i+1];j++){
				//x est le milieu de deux lignes a 0 a ajuster!!!!!
				int x=(limite_x_droite[i]+limite_x_droite[i+1])/2;
				//on va les remplir de 0 (noire) aussi
				for (int i=x;i>0;i--){
							itk::Index<3> curseurwrite2={{i,j,z}};
							roi_crane->SetPixel(curseurwrite2,0);
				}

			}
		}


	}

}
*/
void Roi_crane::limitegauche(int z){
	
		
	//variable pour une trache axiale
	int limite_y[200];
	int limite_x_gauche[200];
	int position=0;

	for (int p=0;p<200;p++){
			//limite_x_droite[p]=0;
			//limite_y[p]=0;
		

	}
	for (int y=limitepost_value;y<200;y++){
		//Analyse la partie droite
		
		//variable pour une tranche coronale
		bool limite=false;
		int epaisseur=0;
	    int trou_coor=0;
		int trou=0;
		int coord=0;
		int patchdelete=0;
		int compteurnoir=0;

		bool coin=false;

		for (int x=180;x>90;x--){
			itk::Index<3> curseurread={{x,y,z}};
			ImageType::PixelType pixelvalue;
			pixelvalue=roi_crane->GetPixel(curseurread);

			//on rencontre le crâne
				if (pixelvalue>141){
					//compteur de l'epaisseur du crane					
					epaisseur+=1;
					//verifier que les pixels du crane sont collé
					if (epaisseur>1){
						if (x-coord>1){
							epaisseur=0;
						}
					}
					coord=x;

																
				}
				if(z<90 && x>120)
					coin=true;
				if(z<100 && y<145)
					coin=true;

				//eliminer les rectangles d'intensité dans les coins (zone d'intensité entre les zones noires)
				if (coin){
					if (pixelvalue>10){
						patchdelete+=1;
					}
					if (patchdelete>4 && pixelvalue==0){
						compteurnoir+=1;
					}
					if (compteurnoir>4){
								//mettre à 0 les pixels extérieurs et ceux du crâne
								for (int i=x;i<181;i++){
									itk::Index<3> curseurwriteg1={{i,y,z}};
									roi_crane->SetPixel(curseurwriteg1,0);
								}
								//reset
								compteurnoir=0;
								patchdelete=0;
								break;
					}
				}
				
				//si la distance est de 5 au plus, c'est le crane
				if (epaisseur>2){
					limite=true;
				}
				if (x==107 && y==161 && z==136){
					
					cout<< epaisseur<<endl;
				}
				if (x==108 && y==121 && z==136){
					
					cout<< epaisseur<<endl;
				}
				//lorsqu'on sort de l'épaisseur du crane
				if (limite==true && pixelvalue<135){
					
					
					//dependamment si cest la première limite de la slice
					if (position!=0){

						int différence_y=y-limite_y[position-1];
						int différence_x=limite_x_gauche[position-1]-x;
						//verifie que les limites sont collées
						if (différence_x<10 && différence_y<10){
							limite_x_gauche[position]=x;
							limite_y[position]=y;

							//mettre à 0 les pixels extérieurs et ceux du crâne
							for (int i=x;i<181;i++){
								itk::Index<3> curseurwrite1={{i,y,z}};
								roi_crane->SetPixel(curseurwrite1,0);
							}
							position+=1;
							
						}
						break;
					//si cest la première limite de la slice
					}else{
						limite_x_gauche[position]=x;
						limite_y[position]=y;
						
						//mettre à 0 les pixels extérieurs et ceux du crâne
						for (int i=x;i<181;i++){
							for (int i=x;i<181;i++){
								itk::Index<3> curseurwrite2={{i,y,z}};
								roi_crane->SetPixel(curseurwrite2,0);
							}
						}
						position+=1;
						break;

					}
					
									
											
									     				 
				
				}
			
		}
		
	}	

	//bouche les trous
	for (int i=0;i<199;i++){
		int difference=limite_y[i+1]-limite_y[i];
		//s'il y a des trous entre les lignes noires
		if (difference>1){
			for (int j=limite_y[i]+1;j<limite_y[i+1];j++){
				//x est le milieu de deux lignes a 0 a ajuster!!!!!
				int x=(limite_x_gauche[i]+limite_x_gauche[i+1])/2;
				//on va les remplir de 0 (noire) aussi
				for (int i=x;i<181;i++){
							itk::Index<3> curseurwrite2={{i,j,z}};
							roi_crane->SetPixel(curseurwrite2,0);
				}

			}
		}


	}

}
/*
void Roi_crane::limitegauche(int z){
	
		
	//variable pour une trache axiale
	int limite_y[200];
	int limite_x_gauche[200];
	int position=0;

	for (int p=0;p<200;p++){
			//limite_x_droite[p]=0;
			//limite_y[p]=0;
		

	}
	for (int y=limitepost_value;y<150;y++){
		//Analyse la partie droite
		
		//variable pour une tranche coronale
		bool limite=false;
		int epaisseur=0;
	    int trou_coor=0;
		int trou=0;
		int coord=0;
		int patchdelete=0;
		int compteurnoir=0;

		for (int x=180;x>90;x--){
			itk::Index<3> curseurread={{x,y,z}};
			ImageType::PixelType pixelvalue;
			pixelvalue=roi_crane->GetPixel(curseurread);

			//on rencontre le crâne
				if (pixelvalue>150){
										
					epaisseur+=1;
					//verifier que les pixels du crane sont collé
					if (epaisseur>1){
						if (x-coord>1){
							epaisseur=0;
						}
					}
					coord=x;

																
				}
				//eliminer les rectangles d'intensité dans les coins
				if (z<80){
					if (pixelvalue>30){
						patchdelete+=1;
					}
					if (patchdelete>10){
						compteurnoir+=1;
					}
					if (compteurnoir>15){
								//mettre à 0 les pixels extérieurs et ceux du crâne
								for (int i=x;i<181;i++){
									itk::Index<3> curseurwrite11={{i,y,z}};
									roi_crane->SetPixel(curseurwrite11,0);
								}
					}
				}
				
				//si la distance est de 5 au plus, c'est le crane
				if (epaisseur>3){
					limite=true;
				}
				
				//lorsqu'on sort de l'épaisseur du crane
				if (limite==true && pixelvalue<141){
					
					
					//dependamment si cest la première limite de la slice
					if (position!=0){

						int différence_y=y-limite_y[position-1];
						int différence_x=limite_x_gauche[position-1]-x;
						//verifie que les limites sont collées
						if (différence_x<10 && différence_y<10){
							limite_x_gauche[position]=x;
							limite_y[position]=y;

							//mettre à 0 les pixels extérieurs et ceux du crâne
							for (int i=x;i<181;i++){
								itk::Index<3> curseurwrite1={{i,y,z}};
								roi_crane->SetPixel(curseurwrite1,0);
							}
							position+=1;
							
						}
						break;
					//si cest la première limite de la slice
					}else{
						limite_x_gauche[position]=x;
						limite_y[position]=y;
						
						//mettre à 0 les pixels extérieurs et ceux du crâne
						for (int j=y;j>0;j--){
							for (int i=x;i<181;i++){
								itk::Index<3> curseurwrite2={{i,j,z}};
								roi_crane->SetPixel(curseurwrite2,0);
							}
						}
						position+=1;
						break;

					}
					
									
											
									     				 
				
				}
			
		}
		
	}	

	//bouche les trous
	for (int i=0;i<199;i++){
		int difference=limite_y[i+1]-limite_y[i];
		//s'il y a des trous entre les lignes noires
		if (difference>1){
			for (int j=limite_y[i]+1;j<limite_y[i+1];j++){
				//x est le milieu de deux lignes a 0 a ajuster!!!!!
				int x=(limite_x_gauche[i]+limite_x_gauche[i+1])/2;
				//on va les remplir de 0 (noire) aussi
				for (int i=x;i<181;i++){
							itk::Index<3> curseurwrite2={{i,j,z}};
							roi_crane->SetPixel(curseurwrite2,0);
				}

			}
		}


	}

}

*/
void Roi_crane::limiteinferieuravant(){
	
	//on analyse avec une tranche coronal à la fois		
	for (int j=130;j<200;j++){

		//variable pour une trache coronal
		int limite_z[200];
		int limite_x[200];
		int position=0;

		for (int p=0;p<200;p++){
				limite_x[p]=0;
				limite_z[p]=0;
		

		}
		//premiere partie	
		for (int i=40;i<150;i++){
			//Analyse la partie inferieur
		
			//variable pour une tranche coronale
			bool limite=false;
			int epaisseur=0;
			int coord=0;

			//analyse du haut vers le bas la limite du crane
				for (int z=95;z>0;z--){
					itk::Index<3> curseurread={{i,j,z}};
					ImageType::PixelType pixelvalue;
					pixelvalue=roi_crane->GetPixel(curseurread);

					//on rencontre le crâne (130 aussi....?)
						if (pixelvalue>180){
						
							epaisseur+=1;
							//verifier que les pixels du crane sont collé
							if (epaisseur>1){
								if (z-coord>1){
									epaisseur=0;
									
								}
							}
							coord=z;
									
													
						}
					
						//verifier qu'il n'y a pas de baisse d'intensité dans le crane
							if (epaisseur>0 && limite==false){
								if (pixelvalue<150){
									epaisseur=0;
									
								}
							}
						
				
						//si la distance est de 5 au plus, c'est le crane
						if (epaisseur>2){
							limite=true;
						}
						if (i==125 && j==165 && z==89){
							cout<<"epaisseur"<<endl;
							cout<<epaisseur<<endl;
						}
						//lorsqu'on sort de l'épaisseur du crane
						if (limite==true){
					
						
									//dependamment si cest la première limite de la slice
							if (position!=0){

								int différence_z=z-limite_z[position-1];
						
								//verifie que le prochain pixel du crane est proche de l'autre
								if (différence_z<7){
							
									limite_z[position]=z;
									limite_x[position]=i;

									//mettre à 0 les pixels extérieurs et ceux du crâne
									for (int k=z+epaisseur;k>0;k--){
												itk::Index<3> curseurwrite1={{i,j,k}};
												roi_crane->SetPixel(curseurwrite1,0);
											}
									position+=1;
									break;
							//si le prochain pixel est plus loin, on fait une 2e condition
								}else if( 6 < différence_z && différence_z < 20){
									
									
									for (int p=-1;p<2;p++){
										int longueur=0;
										for (int l=0;l<10;l++){
											
												itk::Index<3> curseurread2={{i+l,j,z+p}};
												ImageType::PixelType pixelvalue2;
												pixelvalue2=roi_crane->GetPixel(curseurread2);
										
												if (pixelvalue2>160){
													longueur+=1;
												}
									
										}
										//alors on l'accepte comme étant la limite du crane
										if (longueur>8){
											limite_z[position]=z;
											limite_x[position]=i;

											//mettre à 0 les pixels extérieurs et ceux du crâne
											for (int k=z+epaisseur;k>0;k--){
												itk::Index<3> curseurwrite1={{i,j,k}};
												roi_crane->SetPixel(curseurwrite1,0);
											}
											position+=1;
											break;
										}
									}
								}
								//break;
							//si cest la première limite de la slice, on l'accepte
							
							}else{
								limite_z[position]=z;
								limite_x[position]=i;
						
								//mettre à 0 les pixels extérieurs et ceux du crâne
								for (int k=z+epaisseur;k>0;k--){
												itk::Index<3> curseurwrite2={{i,j,k}};
												roi_crane->SetPixel(curseurwrite2,0);
											}
								position+=1;
								break;

							}
							
								     				 
				}
			
				
			
				
			}
		
		}
		
		limite_z[position]=limite_z[position-1];
		limite_x[position]=150;
		//reevalue toutes les lignes mises a zero
		
		for (int a=0;a<199;a++){
			int difference=limite_x[a+1]-limite_x[a];
			//s'il y a des trous entre les lignes noires
			if (difference>1){
				for (int i=limite_x[a]+1;i<limite_x[a+1];i++){
					//x est le milieu de deux lignes a 0 a ajuster!!!!!
					
					int z=(limite_z[a]+limite_z[a+1])/2;
					//on va les remplir de 0 (noire) aussi
					int maxi=max(limite_z[a],limite_z[a+1]);
					for (int k=maxi;k>0;k--){
								itk::Index<3> curseurwrite2={{i,j,k}};
								roi_crane->SetPixel(curseurwrite2,0);
					}

				}
			}
		}
		
		

	}

}
void Roi_crane::limiteinferieurarriere(){

	//on analyse avec une tranche coronal à la fois	
	for (int j=75;j<130;j++){

		//variable pour une trache coronal
		int limite_z[200];
		int limite_x[200];
		int position=0;
		int limite_z2[200];
		int limite_x2[200];
		int position2=0;

		for (int p=0;p<200;p++){
				limite_x[p]=0;
				limite_z[p]=0;
				limite_x2[p]=0;
				limite_z2[p]=0;
		

		}
		//moitie de droite	
		for (int i=38;i<91;i++){
			//Analyse la partie inferieur
		
			//variable pour une tranche coronale
			bool limite=false;
			int epaisseur=0;
			int coord=0;
			int* pointernoirceur;
			int* pointerestdansimage;
			int noirceur=0;
			int estdansimage=0;
			pointernoirceur=&noirceur;
			pointerestdansimage=&estdansimage;

			ImageType::PixelType max=0;

			//analyse du bas vers le haut la limite de crane
				for (int z=0;z<90;z++){

					
					itk::Index<3> curseurread={{i,j,z}};
					ImageType::PixelType pixelvalue;
					pixelvalue=roi_crane->GetPixel(curseurread);

					//on rencontre le crâne
						if (pixelvalue>160){
						
							epaisseur+=1;
							//sauvegarde le maximum
							if (max<pixelvalue)
								max=pixelvalue;

							//verifier que les pixels du crane sont collé
							if (epaisseur>1){
								if (z-coord>1){
									epaisseur=0;
									
								}
							}
							coord=z;
									
													
						}
					
						//verifier qu'il n'y a pas de baisse d'intensité dans le crane
							if (epaisseur>0 && limite==false){
								if (pixelvalue<140){
									epaisseur=0;
									
								}
							}
						
				
						//si la distance est de 5 au plus, c'est le crane
						if (epaisseur>3){
							limite=true;
						}
				
						//lorsqu'on sort de l'épaisseur du crane
						if (limite==true && pixelvalue<155){
					
						
									//dependamment si cest la première limite de la slice
							if (position!=0){

								//difference de hauteur avec le pixel limite precedent
								int différence_z=z-limite_z[position-1];
						
								//verifie que les limites sont collées
								if (différence_z<7 ){
									//sauvegarde
									limite_z[position]=z;
									limite_x[position]=i;

									//mettre à 0 les pixels extérieurs et ceux du crâne
									for (int k=z;k>0;k--){
												itk::Index<3> curseurwrite1={{i,j,k}};
												roi_crane->SetPixel(curseurwrite1,0);
											}
									position+=1;
								break;
								}else if( 6 < différence_z && différence_z < 20){
									
									int longueur1=0;
									for (int p=-1;p<2;p++){
										int longueur2=0;
										for (int l=-3;l<4;l++){
											
												itk::Index<3> curseurread2={{i+l,j,z+p}};
												ImageType::PixelType pixelvalue2;
												pixelvalue2=roi_crane->GetPixel(curseurread2);
										
												if (pixelvalue2>160){
													longueur2+=1;
												}
									
										}
										//alors on l'accepte comme étant la limite du crane
										if (longueur2>8){
											limite_z[position]=z;
											limite_x[position]=i;

											//mettre à 0 les pixels extérieurs et ceux du crâne
											for (int k=z;k>0;k--){
												itk::Index<3> curseurwrite1={{i,j,k}};
												roi_crane->SetPixel(curseurwrite1,0);
											}
											position+=1;
											break;
										}
									}
								}
								//break;
							//si cest la première limite de la slice
							}else{
								limite_z[position]=z;
								limite_x[position]=i;
						
								//mettre à 0 les pixels extérieurs et ceux du crâne
								for (int k=z;k>0;k--){
												itk::Index<3> curseurwrite2={{i,j,k}};
												roi_crane->SetPixel(curseurwrite2,0);
											}
								position+=1;
								break;

							}
							
								     				 
						}
			
				
			
				
				}
		
		}	
		
		//moitie de gauche	
		for (int i=144;i>89;i--){

			//variable pour une trache coronal
			

			
			//Analyse la partie inferieur
		
			//variable pour une tranche coronale
			bool limite=false;
			int epaisseur=0;
			int coord=0;
			int* pointernoirceur;
			int* pointerestdansimage;
			int noirceur=0;
			int estdansimage=0;
			pointernoirceur=&noirceur;
			pointerestdansimage=&estdansimage;

			//analyse du bas vers le haut la limite de crane
				for (int z=0;z<90;z++){

					
					itk::Index<3> curseurread={{i,j,z}};
					ImageType::PixelType pixelvalue;
					pixelvalue=roi_crane->GetPixel(curseurread);

					//on rencontre le crâne
						if (pixelvalue>160){
						
							
							epaisseur+=1;
							//verifier que les pixels du crane sont collé
							if (epaisseur>1){
								if (z-coord>1){
									epaisseur=0;
									
								}
							}
							coord=z;
									
													
						}
					
						//verifier qu'il n'y a pas de baisse d'intensité dans le crane
							if (epaisseur>0 && limite==false){
								if (pixelvalue<140){
									epaisseur=0;
									
								}
							}
						
				
						//si la distance est de 5 au plus, c'est le crane
						if (epaisseur>3){
							limite=true;
						}
				
						//lorsqu'on sort de l'épaisseur du crane
						if (limite==true && pixelvalue<155){
					
						
									//dependamment si cest la première limite de la slice
							if (position2!=0){

								//difference de hauteur avec le pixel limite precedent
								int différence_z=z-limite_z2[position2-1];
						
								//verifie que les limites sont collées
								if (différence_z<7 ){
									//sauvegarde
									limite_z2[position2]=z;
									limite_x2[position2]=i;

									//mettre à 0 les pixels extérieurs et ceux du crâne
									for (int k=z;k>0;k--){
												itk::Index<3> curseurwrite1={{i,j,k}};
												roi_crane->SetPixel(curseurwrite1,0);
											}
									position2+=1;
								break;
								}else if( 6 < différence_z && différence_z < 20){
									
									int longueur1=0;
									for (int p=-1;p<2;p++){
										int longueur2=0;
										for (int l=-3;l<4;l++){
											
												itk::Index<3> curseurread2={{i+l,j,z+p}};
												ImageType::PixelType pixelvalue2;
												pixelvalue2=roi_crane->GetPixel(curseurread2);
										
												if (pixelvalue2>160){
													longueur2+=1;
												}
									
										}
										//alors on l'accepte comme étant la limite du crane
										if (longueur2>8){
											limite_z2[position2]=z;
											limite_x2[position2]=i;

											//mettre à 0 les pixels extérieurs et ceux du crâne
											for (int k=z;k>0;k--){
												itk::Index<3> curseurwrite1={{i,j,k}};
												roi_crane->SetPixel(curseurwrite1,0);
											}
											position2+=1;
											break;
										}
									}
								}
								//break;
							//si cest la première limite de la slice
							}else{
								limite_z2[position2]=z;
								limite_x2[position2]=i;
						
								//mettre à 0 les pixels extérieurs et ceux du crâne
								for (int k=z;k>0;k--){
												itk::Index<3> curseurwrite2={{i,j,k}};
												roi_crane->SetPixel(curseurwrite2,0);
											}
								position2+=1;
								break;

							}
							
								     				 
						}
			
				
			
				
				}
		
		}	
		

		limite_z[position]=limite_z[position-1];
		limite_x[position]=90;
		//reevalue toutes les lignes mises a zero
		for (int a=0;a<199;a++){
			int difference=limite_x[a+1]-limite_x[a];
			//s'il y a des trous entre les lignes noires
			if (difference>1){
				for (int i=limite_x[a]+1;i<limite_x[a+1];i++){
					//x est le milieu de deux lignes a 0 a ajuster!!!!!
					int z=(limite_z[a]+limite_z[a+1])/2;
					//on va les remplir de 0 (noire) aussi
					for (int k=z;k>0;k--){
								itk::Index<3> curseurwrite2={{i,j,k}};
								roi_crane->SetPixel(curseurwrite2,0);
					}

				}
			}
		}
		//2e moitié
		limite_z2[position2]=limite_z2[position2-1];
		limite_x2[position2]=90;
		//reevalue toutes les lignes mises a zero
		for (int a=0;a<199;a++){
			int difference=limite_x2[a+1]-limite_x2[a];
			//s'il y a des trous entre les lignes noires
			if (difference>1){
				for (int i=limite_x2[a]+1;i<limite_x2[a+1];i++){
					//x est le milieu de deux lignes a 0 a ajuster!!!!!
					int z=(limite_z2[a]+limite_z2[a+1])/2;
					//on va les remplir de 0 (noire) aussi
					for (int k=z;k>0;k--){
								itk::Index<3> curseurwrite2={{i,j,k}};
								roi_crane->SetPixel(curseurwrite2,0);
					}

				}
			}
		}
		

	}

}


