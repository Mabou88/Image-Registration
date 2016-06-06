#include "Roi_crane.h"
typedef itk::ImageDuplicator< ImageType > DuplicatorType;

Roi_crane::Roi_crane(ImageType::Pointer itkimageus ){
		
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(itkimageus);
  duplicator->Update();
  roi_crane=ImageType::New();
  roi_crane = duplicator->GetOutput();
  
  
}
ImageType::Pointer Roi_crane::getcraneROI(){
	return roi_crane;
}

void Roi_crane::sauvegardeimage(){
	//sauvegarde l'image modifier
  WriterType::Pointer   writer =  WriterType::New();
  itk::NiftiImageIO::Pointer ioimagenifti=itk::NiftiImageIO::New();

  writer->SetImageIO(ioimagenifti);
  writer->SetFileName( "C:/Users/Marc-Antoine/Documents/Imagecode/output/roi_crane_P2V2_3.nii");
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


void Roi_crane::maskcrane(ImageType::Pointer mask_us){
	
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
	typedef itk::LinearInterpolateImageFunction<ImageType> LinearInterpolatorFilterType;
	LinearInterpolatorFilterType::Pointer interpolator = LinearInterpolatorFilterType::New();
     interpolator->SetInputImage(mask_us);
	 
	ImageType::SizeType dimension=roi_crane->GetLargestPossibleRegion().GetSize();
	for (int i=0;i<dimension[0];i++){
		for (int j=0;j<dimension[1];j++){
			for (int k=0;k<dimension[2];k++){
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
					  roi_crane->SetPixel(index_us,0);
				  }
			}
		}
	}

}


void Roi_crane::calculvolumecrane(){
	int min_x=dim_x;
	int max_x=0;
	int min_y=dim_y;
	int max_y=0;
	int min_z=dim_z;
	int max_z=0;
	
	//min max en x y et z
	
	for (int i=0;i<dim_x;i++){
		for (int j=0;j<dim_y;j++){
			for (int k=0;k<dim_z;k++){
				ImageType::IndexType index ={{i,j,k}};
				ImageType::PixelType pixelvalue=roi_crane->GetPixel(index);
				if (pixelvalue!=0){
					if (i>max_x)
						max_x=i;
					if (i<min_x)
						min_x=i;
					if (j>max_y)
						max_y=j;
					if (j<min_y)
						min_y=j;
					if (k>max_z)
						max_z=k;
					if (k<min_z)
						min_z=k;
					
				}

				
			}
		}
	}
	
	cout<<"limite x"<<min_x<<max_x<<endl;
	cout<<"limite y"<<min_y<<max_y<<endl;
	cout<<"limite x"<<min_z<<max_z<<endl;

	int a=(max_x-min_x)/2;
	int b=(max_y-min_y)/2;
	int c=(max_z-min_z)/2;

	double fraction=1.333333;
	double pi=3.14159265359;
	double volume=fraction*pi*a*b*c;
	
	cout<<"volume"<<volume<<endl;


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


