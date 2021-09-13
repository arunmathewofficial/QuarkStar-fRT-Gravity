//[=============================================================================================================]
//[ Source library for Quark star in f(R) = R + alpha R^2  gravity. 											]
//[																												]																
//[	The main purpose of the code is to identify the central value of Ricci scalar R_c correspoding to a given	]	
//[ central density. The equation of state of the matter is given by the bag model, with pressure expressed as  ]
//[ p = k(epsilon - 4B), where B representing the bag.															]
//[ 																											]																																						]				
//[ The code solves the interior and exterior field equations iterratively until the correct value of R_c		]
//[ satifing the boundary condition is obtained. 																]
//[																												]
//[ Author: Arun Mathew, New Numerical Lab, Dept. of Physics, IIT Guwahati.                      				]
//[																												]
//[=============================================================================================================]


//[Included library] 
	#include<stdio.h>
	#include<math.h>
	#include<string.h>
	#include<stdlib.h>
	#include<time.h>
	#include"starobinsky.h"   //[ User defined library]

//[Definition for colors]                                                                          				
	#define white   "\x1B[37m"
	#define red     "\x1b[31m"
	#define green   "\x1b[32m"
	#define yellow  "\x1b[33m"
	#define blue    "\x1b[34m"
	#define violet  "\x1b[35m"
	#define cyan    "\x1b[36m"
	#define reset   "\x1b[0m"

//[Physical Constants in SI units]  																			]
	#define c    	2.99792458*pow(10,8)		   	 //[  Speed of light in meter/sec      						]	
	#define pi   	3.141592654				 	 	 //[  Pi	         										]
	#define pl   	6.626070040*pow(10,-34)		 	 //[  Planck constant       								]
	#define G    	6.6743015*pow(10,-11)	 		 //[  Gravitational Constant           	 				    ]
	#define S    	1.98847*pow(10,30)  		     //[  One Solar Mass in the units of kilograms				]
	#define B       9.613059804*pow(10,33.0)		 //[  The Bag constant in the units of Joule/meter^3		]	
 
 //[Constants used in the code] 
	#define r_g		G*S*pow(c,-2.0)
	#define alpha   10.0*pow(r_g,2.0)
	#define beta	0.01
	#define kappa	8*pi*G*pow(c,-4.0)	
	#define q1 		4.0*kappa*B*pow(r_g,2.0)
	#define q2		pow(r_g,2.0)*pow(6.0*alpha,-1.0)
	#define q3		beta
	#define k 		0.28

//[ Scaling constants 																							]
	#define rho_scale   4.0*B*pow(c,-2.0)*pow(10,-3.0)		//[ Density scale in the units of gm/cm^3 			]
	#define R_scale     G*S*pow(c,-2.0)*pow(10,-3.0)		//[ Radial scale in the units of km 				]


//[Interior Functions ==========================================================================================]

/*  Interior field equations for starobinsky gravity. "0" stands for starobinsky quatities and 
	"I" stands for interior.                                                              
*/

long double gamma0_I (long double eta, long double chi0, long double xi0)
{
	//[ Gamma function]
   long double RV = pow( 6.0*q2 , -1.0  )*eta*pow( phi(chi0) , -1.0)*xi0;
   return RV;
}

long double Dnu0_I (long double eta, long double the0, long double lam0, long double chi0, long double xi0)
	{
		long double Eq_1 =   q1*eta*exp(lam0)*pressure(the0)*pow( phi(chi0) , -1.0 );
		long double Eq_2 =   ( exp(lam0) -1.0 )/eta;
		long double Eq_3 = - eta*exp(lam0)*pow( chi0, 2.0 )*pow( 12.0*q2*phi(chi0) , -1.0 );
		long double Eq_4 = - 4.0*gamma0_I(eta, chi0, xi0)/eta;

		long double RV   =   (Eq_1 + Eq_2 + Eq_3 + Eq_4)*pow( 1.0 + gamma0_I(eta, chi0, xi0) , -1.0 );

	return RV;
	}

long double Dlam0_I (long double eta, long double the0, long double lam0, long double chi0, long double xi0)
	{
		long double Eq_1 =   ( 1.0 - exp(lam0) )/eta ; 
		long double Eq_2 =   eta*exp(lam0)*pow( 6.0*phi(chi0), -1.0 )*( 2.0 + 0.5*chi0/q2 )*chi0;
		long double Eq_3 = - gamma0_I(eta, chi0, xi0)*Dnu0_I(eta, the0, lam0, chi0, xi0);
		long double Eq_4 =   q1*eta*exp(lam0)*( 3.0*epsilon(the0) + trace(the0) )*pow( 3.0*phi(chi0), -1.0 );

		long double RV   =   Eq_1 + Eq_2 + Eq_3 + Eq_4;

	return RV;
	}

long double Dchi0_I (long double xi0)
	{
		long double RV = xi0;      

		return RV;
	}

long double Dxi0_I (long double eta, long double the0, long double lam0, long double chi0, long double xi0)
	{
		long double Eq_1 = - 2.0*xi0/eta + q1*q2*trace(the0)*exp(lam0) + q2*exp(lam0)*chi0;
		long double Eq_2 =   0.5*( Dlam0_I(eta, the0, lam0, chi0, xi0) - Dnu0_I(eta, the0, lam0, chi0, xi0) )*xi0;

		long double RV   =   Eq_1 + Eq_2 ;

		return RV;
	}

long double Dthe0_I (long double eta, long double the0, long double lam0, long double chi0, long double xi0)
	{
		long double RV = - 0.5*( the0 + pressure(the0) )*pow(k,-1.0)*Dnu0_I(eta, the0, lam0, chi0, xi0); 

		return RV;
	}


/*  Interior field equations for the actual f(R,T) gravity. "A" stands for the actual quatities and 
	"I" stands for interior. "F" represent the terms without omega and "S" represent terms with omega.                                                            
*/


long double gammaA_I (long double eta, long double chi, long double xi)
	{
   		long double RV = pow( 6.0*q2 , -1.0  )*eta*pow( phi(chi) , -1.0)*xi;
   		return RV;
	}


long double DnuA_I (long double eta, long double the, long double lam, long double chi, long double xi, long double the0, long double lam0, long double chi0, long double xi0)
	{
		long double FEq_1 =   q1*eta*exp(lam)*pressure(the)*pow( phi(chi) , -1.0 );
		long double FEq_2 =   ( exp(lam) -1.0 )/eta;
		long double FEq_3 = - eta*exp(lam)*pow( chi, 2.0 )*pow( 12.0*q2*phi(chi) , -1.0 );
		long double FEq_4 = - 4.0*gammaA_I(eta, chi, xi)/eta;

		long double SEq_1 =   ( trace(the0) + 0.5*eta*Dtrace(eta, the0, lam0, chi0, xi0) )*Dnu0_I(eta, the0, lam0, chi0, xi0)*pow( phi(chi0) , -1.0 );
		long double SEq_2 =   ( 1.0 - exp(lam0) )*trace(the0)*pow( phi(chi0) , -1.0 );
		long double SEq_3 =   2.0*Dtrace(eta, the0, lam0, chi0, xi0)*pow( phi(chi0) , -1.0 );
		

		long double RV   =   (FEq_1 + FEq_2 + FEq_3 + FEq_4)*pow( 1.0 + gammaA_I(eta, chi, xi) , -1.0 ) - q3*( SEq_1 + SEq_2 + SEq_3 )*pow( 1.0 + gamma0_I(eta, chi0, xi0) , -1.0 );

	return RV;
	}

long double DlamA_I (long double eta, long double the, long double lam, long double chi, long double xi, long double the0, long double lam0, long double chi0, long double xi0)
	{

		long double FEq_1 =   ( 1.0 - exp(lam) )/eta ; 
		long double FEq_2 =   eta*exp(lam)*pow( 6.0*phi(chi), -1.0 )*( 2.0 + 0.5*chi/q2 )*chi;
		long double FEq_3 = - gammaA_I(eta, chi, xi)*DnuA_I(eta, the, lam, chi, xi, the0, lam0, chi0, xi0);
		long double FEq_4 =   q1*eta*exp(lam)*( 3.0*epsilon(the) + trace(the) )*pow( 3.0*phi(chi), -1.0 );

		long double SEq_1 =   eta*exp(lam0)*( 3.0*epsilon(the0) - pressure(the0) + 2.0*trace(the0)  )*chi0*pow( 3.0*phi(chi0), -1.0 );
		long double SEq_2 = - ( 0.5*eta*Dtrace( eta, the0, lam0, chi0, xi0)*Dnu0_I( eta, the0, lam0, chi0, xi0) + trace(the0)*Dlam0_I( eta, the0, lam0, chi0, xi0) )*pow( phi(chi0) , -1.0);
		long double SEq_3 =   (1.0 - exp(lam0) )*trace(the0)*pow( eta*phi(chi0) , -1.0 );

		long double RV   =   FEq_1 + FEq_2 + FEq_3 + FEq_4 + q3*( SEq_1 + SEq_2 + SEq_3 );

	return RV;
	}

long double DchiA_I (long double xi)
	{
		long double RV = xi;      

		return RV;
	}

long double DxiA_I (long double eta, long double the, long double lam, long double chi, long double xi, long double the0, long double lam0, long double chi0, long double xi0, long double h_I)
	{
		long double FEq_1 = - 2.0*xi/eta + q1*q2*trace(the)*exp(lam) + q2*exp(lam)*chi;
		long double FEq_2 =   0.5*( DlamA_I(eta, the, lam, chi, xi, the0, lam0, chi0, xi0) - DnuA_I(eta, the, lam, chi, xi, the0, lam0, chi0, xi0) )*xi;

		long double SEq_1 =   2.0*( trace(the0) - 2.0*pressure(the0) )*chi0*exp(lam0);
		long double SEq_2 = - 3.0*( 0.5*( Dnu0_I(eta, the0, lam0, chi0, xi0) - Dlam0_I(eta, the0, lam0, chi0, xi0) ) + 2.0/eta)*Dtrace( eta, the0, lam0, chi0, xi0);
		long double SEq_3 = - 3.0*DDtrace( eta, the0, lam0, chi0, xi0, h_I );						

		long double RV   =   FEq_1 + FEq_2 + q2*q3*( SEq_1 + SEq_2 +SEq_3 );

		return RV;
	}

long double DtheA_I (long double eta, long double the, long double lam, long double chi, long double xi, long double the0, long double lam0, long double chi0, long double xi0)
	{
		long double FEq_1 =  0.5*( the + pressure(the) )*DnuA_I(eta, the, lam, chi, xi, the0, lam0, chi0, xi0); 
		long double SEq_1 =  0.5*( the0 + pressure(the0) )*chi0*Dnu0_I( eta, the0, lam0, chi0, xi0);

		long double RV     = - ( FEq_1 + q3*SEq_1/q1 )/k;

		return RV;
	}

//[Exterior functions ==========================================================================================]

/*  Exterior field equations for the actual f(R,T) gravity. "A" stands for the actual quatities and 
	"E" stands for exterior.                                                              
*/
	 											

long double gammaA_E (long double eta, long double chi, long double xi)
{
	//[ Gamma function]
   long double RV = pow( 6.0*q2 , -1.0  )*eta*pow( phi(chi) , -1.0)*xi;
   return RV;
}


long double DnuA_E (long double eta, long double lam, long double chi, long double xi)
	{
		long double Eq_1 =   ( exp(lam) - 1.0 )/eta;
		long double Eq_2 = - eta*exp(lam)*pow( chi, 2.0 )*pow( 12.0*q2*phi(chi) , -1.0 );
		long double Eq_3 = - 4.0*gammaA_E(eta, chi, xi )/eta;

		long double RV   =   (Eq_1 + Eq_2 + Eq_3)*pow( 1.0 + gammaA_E(eta, chi, xi ) , -1.0 );

	return RV;
	}

long double DlamA_E (long double eta, long double lam, long double chi, long double xi)
	{
		long double Eq_1 =   ( 1.0 - exp(lam) )/eta ; 
		long double Eq_2 =   eta*exp(lam)*pow( 6.0*phi(chi), -1.0 )*( 2.0 + 0.5*chi/q2 )*chi;
		long double Eq_3 = - gammaA_E(eta, chi, xi)*DnuA_E(eta, lam, chi, xi);
		
		long double RV   =   Eq_1 + Eq_2 + Eq_3;

	return RV;
	}

long double DchiA_E (long double xi)
	{
		long double RV = xi;      

		return RV;
	}

long double DxiA_E (long double eta, long double lam, long double chi, long double xi)
	{
		long double Eq_1 = - 2.0*xi/eta + q2*exp(lam)*chi;
		long double Eq_2 =   0.5*( DlamA_E(eta, lam, chi, xi) - DnuA_E(eta, lam, chi, xi) )*xi;

		long double RV   =   Eq_1 + Eq_2 ;

		return RV;
	}



//[Main program ================================================================================================]

int main () 
{	
	int n;
	long double the_c = 5.94;
	long double the_step;

  	FILE* main_datafile  = fopen("main_data.dat","w");

	//for(n=1; n<=103; n++)
	{
		//if(the_c<2.0){the_step = 0.05;}
		//else the_step = 0.2; 

		//the_c = the_c + the_step;
  	

  		printf("\n ["white"%d"reset"]\t[theta_c: %0.3Lf]\n",n,the_c);
		

  		//[ Declaration of the central variables 
  		long double lam0_c, nu0_c, chi0_c, xi0_c;															
  		long double lam_c, nu_c, chi_c, xi_c;
  	     		        

  		//[ Integration variables of interior loop																	]
  		long double h_I			 = 0.0001;           //[ Interior Step Size   								
  		long double eta;    


  		long double I10, I20, I30, I40;
		long double J10, J20, J30, J40;
		long double K10, K20, K30, K40;
		long double M10, M20, M30, M40;
		long double N10, N20, N30, N40;
		long double the0, lam0, nu0, chi0, xi0;


		long double I1A, I2A, I3A, I4A;
		long double J1A, J2A, J3A, J4A;
		long double K1A, K2A, K3A, K4A;
		long double M1A, M2A, M3A, M4A;
		long double N1A, N2A, N3A, N4A;
		long double the, lam, nu, chi, xi;

		//[ Declaration of the surface variables 																	]
		long double eta_s;															
  		long double lam_s;
  		long double nu_s; 
  		long double chi_s; 			
  		long double xi_s;      

  		long double eta_H;          //[Halo radius]              		        

  		//[ Integration variables of exterior loop																	]
  		long double h_E			 = 0.001;           //[ Exterior Step Size   									]         						
  		long double R1, R2, R3, R4;
		long double S1, S2, S3, S4;
		long double T1, T2, T3, T4;
		long double U1, U2, U3, U4;
	

		//[ Temporary Varibale 																						]
	 	long double temp1_I, temp2_I, temp3_I, temp4_I, temp5_I;
		long double temp1_E, temp2_E, temp3_E, temp4_E, temp5_E;

		//[ To reduce the file size 																			]
 		int loop_counter;

		//[Boundary Conditions at the center ===================================================================]
	
  		int flag;
  		int count;
  		long double chi_c_Eienstein;
		long double chi_c_1, chi_c_2;
		long double chi_c_intermediate;

  		chi0_c   	=   0.0058032369228479;	

  		printf("\n");	
  		printf(yellow" \t[:: Calculating chi_c value in f_RT Gravity ::]\n");
  		printf(" \tParameters: [alpha: %0.3f r_g^2] [omega: %0.3f (4B)^-1]\n",alpha*pow(r_g,-2.0),beta);
		printf(" \tCentral Density: [%Le  gcm^-3][theta_c: %0.3Lf]\n",rho_scale*the_c, the_c);


		chi_c   	=   0.006391010003;	

		printf(" \tCenter Starobinsky chi_c: %0.16Lf\n", chi_c);
  		
  		/*	flag = 0;
  			count= 0;
			chi_c_Eienstein = - q1*trace(the_c); 

			chi_c_1 = 0.0; 
			chi_c_2 = chi_c_Eienstein;

			chi_c_intermediate = (chi_c_1+chi_c_2)/2.0;


			jump:

				if( flag == 1)
					{
					chi_c_1 = chi_c_intermediate;
					chi_c_intermediate = (chi_c_1+chi_c_2)/2.0;
					}

			 
				if( flag == 2)
					{
					chi_c_2 = chi_c_intermediate;
   					chi_c_intermediate = (chi_c_1+chi_c_2)/2.0;
					} 

  			chi_c   	= chi_c_intermediate;

  		   //[ Progress Bar ]
  			printf("\r \tTuning chi_c value: [%0.12Lf] ["red"%d"yellow"]", chi_c, count);
			fflush(stdout);

		*/	


  		lam0_c	 	=   0.0;
  		nu0_c     	=   1.0; 		
  		xi0_c       =   0.0;

  		lam_c	 	=   0.0;
  		nu_c     	=   1.0; 		
  		xi_c     	=   0.0;

  		eta_H 		= 	0.0;

		//[End of Boundary conditions ==========================================================================]

  		//[initializing the loop counter																		]
  		loop_counter = 0;


  		//[ Declaration of data files 

  		char title1[50];																				
  		sprintf(title1,"./f_RT_DATA/f_RT[%0.3Lf].dat",the_c);
  		FILE* fRT_datafile  = fopen(title1,"w");

  		char title2[50];																				
  		sprintf(title2,"./f_RT_EC/f_RT_EC[%0.3Lf].dat",the_c);
		FILE* fRT_EC_data  = fopen(title2,"w");

  	
  		//[ Initializing the central values of the variables 														]
  			eta 		= h_I;

			the0  		= the_c;
  			lam0 		= lam0_c; 	
  			nu0   	  	= nu0_c;
  			chi0 	   	= chi0_c;
  			xi0      	= xi0_c;

  			the  		= the_c;
  			lam 		= lam_c; 	
  			nu   	  	= nu_c;
  			chi 	   	= chi_c;
  			xi      	= xi_c;



  		fprintf(fRT_datafile,"%Lf\t%Lf\t%Lf\t%Lf\t%0.12Lf\t%0.12Lf\t%Lf\n",eta*R_scale,the,lam,nu,chi,xi,mass(eta,lam));
 		fprintf(fRT_EC_data,"%Lf\t%Lf\t%Lf\t%Lf\t%Lf\n",eta*R_scale,the,the+pressure(the),the-pressure(the),the+3.0*pressure(the));	
					
  		while( the > 1.0 ) //[ INTERIOR LOOP ===================================================================]
  																	
  		{ 
         
        	if(loop_counter % 100 == 0)
        	//{  	   
 			fprintf(fRT_datafile,"%Lf\t%Lf\t%Lf\t%Lf\t%0.12Lf\t%0.12Lf\t%Lf\n",eta*R_scale,the,lam,nu,chi,xi,mass(eta,lam));
 			fprintf(fRT_EC_data,"%Lf\t%Lf\t%Lf\t%Lf\t%Lf\n",eta*R_scale,the,the+pressure(the),the-pressure(the),the+3.0*pressure(the));
 			//}
 		printf("%Lf\t%Lf\t%Lf\t%Lf\t%0.12Lf\t%0.12Lf\n",eta*R_scale,the,lam,nu,chi,xi);

 			temp1_I = eta; temp2_I = lam; temp3_I = nu, temp4_I = chi; temp5_I = xi;

    		//[ Evaluation of the first slope in RK4 method for the Starobinsky case	 						]
    		I10 =  h_I*Dlam0_I(eta, the0, lam0, chi0, xi0);
    		J10 =  h_I*Dnu0_I(eta, the0, lam0, chi0, xi0);
    		K10 =  h_I*Dxi0_I(eta, the0, lam0, chi0, xi0);
    		M10 =  h_I*Dchi0_I(xi0);
    		N10 =  h_I*Dthe0_I(eta, the0, lam0, chi0, xi0);


    		//[ Evaluation of the first slope in RK4 method for the Actual case									]
    		I1A =  h_I*DlamA_I(eta, the, lam, chi, xi, the0, lam0, chi0, xi0);
    		J1A =  h_I*DnuA_I(eta, the, lam, chi, xi, the0, lam0, chi0, xi0);
    		K1A =  h_I*DxiA_I(eta, the, lam, chi, xi, the0, lam0, chi0, xi0, h_I);
    		M1A =  h_I*DchiA_I(xi);
    		N1A =  h_I*DtheA_I(eta, the, lam, chi, xi, the0, lam0, chi0, xi0);


    		//[ Evaluation of the second slope in RK4 method for the Starobinsky case							]
			I20 =  h_I*Dlam0_I(eta+0.5*h_I, the0+0.5*N10, lam0+0.5*I10, chi0+0.5*M10, xi0+0.5*K10); 
			J20 =  h_I*Dnu0_I(eta+0.5*h_I, the0+0.5*N10, lam0+0.5*I10, chi0+0.5*M10, xi0+0.5*K10); 
			K20 =  h_I*Dxi0_I(eta+0.5*h_I, the0+0.5*N10, lam0+0.5*I10, chi0+0.5*M10, xi0+0.5*K10);
			M20 =  h_I*Dchi0_I(xi0+0.5*K10);
			N20 =  h_I*Dthe0_I(eta+0.5*h_I, the0+0.5*N10, lam0+0.5*I10, chi0+0.5*M10, xi0+0.5*K10);

			//[ Evaluation of the second slope in RK4 method for the Actual case								]							]
			I2A =  h_I*DlamA_I(eta+0.5*h_I, the+0.5*N1A, lam+0.5*I1A, chi+0.5*M1A, xi+0.5*K1A, the0+0.5*N10, lam0+0.5*I10, chi0+0.5*M10, xi0+0.5*K10); 
			J2A =  h_I*DnuA_I(eta+0.5*h_I, the+0.5*N1A, lam+0.5*I1A, chi+0.5*M1A, xi+0.5*K1A, the0+0.5*N10, lam0+0.5*I10, chi0+0.5*M10, xi0+0.5*K10); 
			K2A =  h_I*DxiA_I(eta+0.5*h_I, the+0.5*N1A, lam+0.5*I1A, chi+0.5*M1A, xi+0.5*K1A, the0+0.5*N10, lam0+0.5*I10, chi0+0.5*M10, xi0+0.5*K10, h_I);
			M2A =  h_I*DchiA_I(xi+0.5*K1A);
			N2A =  h_I*DtheA_I(eta+0.5*h_I, the+0.5*N1A, lam+0.5*I1A, chi+0.5*M1A, xi+0.5*K1A, the0+0.5*N10, lam0+0.5*I10, chi0+0.5*M10, xi0+0.5*K10);

    		//[ Evaluation of the third slope in RK4 method for the Starobinsky case							]
    		I30 =  h_I*Dlam0_I(eta+0.5*h_I, the0+0.5*N20, lam0+0.5*I20, chi0+0.5*M20, xi0+0.5*K20);
    		J30 =  h_I*Dnu0_I(eta+0.5*h_I, the0+0.5*N20, lam0+0.5*I20, chi0+0.5*M20, xi0+0.5*K20);
    		K30 =  h_I*Dxi0_I(eta+0.5*h_I, the0+0.5*N20, lam0+0.5*I20, chi0+0.5*M20, xi0+0.5*K20);
    		M30 =  h_I*Dchi0_I(xi0+0.5*K20);
    		N30 =  h_I*Dthe0_I(eta+0.5*h_I, the0+0.5*N20, lam0+0.5*I20, chi0+0.5*M20, xi0+0.5*K20);



    		//[ Evaluation of the third slope in RK4 method for the Actual case									]							]
    		I3A =  h_I*DlamA_I(eta+0.5*h_I, the+0.5*N2A, lam+0.5*I2A, chi+0.5*M2A, xi+0.5*K2A, the0+0.5*N20, lam0+0.5*I20, chi0+0.5*M20, xi0+0.5*K20);
    		J3A =  h_I*DnuA_I(eta+0.5*h_I, the+0.5*N2A, lam+0.5*I2A, chi+0.5*M2A, xi+0.5*K2A, the0+0.5*N20, lam0+0.5*I20, chi0+0.5*M20, xi0+0.5*K20);
    		K3A =  h_I*DxiA_I(eta+0.5*h_I, the+0.5*N2A, lam+0.5*I2A, chi+0.5*M2A, xi+0.5*K2A, the0+0.5*N20, lam0+0.5*I20, chi0+0.5*M20, xi0+0.5*K20, h_I);
    		M3A =  h_I*DchiA_I(xi+0.5*K2A);
    		N3A =  h_I*DtheA_I(eta+0.5*h_I, the+0.5*N2A, lam+0.5*I2A, chi+0.5*M2A, xi+0.5*K2A, the0+0.5*N20, lam0+0.5*I20, chi0+0.5*M20, xi0+0.5*K20);

    		//[ Evaluation of the forth slope in RK4 method for the Starobinsky case							]	
    		I40 =  h_I*Dlam0_I(eta+h_I, the0+N30, lam0+I30, chi0+M30, xi0+K30);
    		J40 =  h_I*Dnu0_I(eta+h_I, the0+N30, lam0+I30, chi0+M30, xi0+K30);
    		K40 =  h_I*Dxi0_I(eta+h_I, the0+N30, lam0+I30, chi0+M30, xi0+K30);
    		M40 =  h_I*Dchi0_I(xi0+K30);
    		N40 =  h_I*Dthe0_I(eta+h_I, the0+N30, lam0+I30, chi0+M30, xi0+K30);


    		//[ Evaluation of the forth slope in RK4 method for the Actual case								     ]	
            I4A =  h_I*DlamA_I(eta+h_I, the+N3A, lam+I3A, chi+M3A, xi+K3A, the0+N30, lam0+I30, chi0+M30, xi0+K30);
    		J4A =  h_I*DnuA_I(eta+h_I, the+N3A, lam+I3A, chi+M3A, xi+K3A, the0+N30, lam0+I30, chi0+M30, xi0+K30);
    		K4A =  h_I*DxiA_I(eta+h_I, the+N3A, lam+I3A, chi+M3A, xi+K3A, the0+N30, lam0+I30,  chi0+M30, xi0+K30, h_I);
    		M4A =  h_I*DchiA_I(xi+K3A);
    		N4A =  h_I*DtheA_I(eta+h_I, the+N3A, lam+I3A, chi+M3A, xi+K3A, the0+N30, lam0+I30, chi0+M30, xi0+K30);

   
		    lam0    = lam0     +  (I10 + 2.0*I20 + 2.0*I30 + I40)/6.0;
			nu0     = nu0      +  (J10 + 2.0*J20 + 2.0*J30 + J40)/6.0; 
			xi0     = xi0      +  (K10 + 2.0*K20 + 2.0*K30 + K40)/6.0; 
			chi0    = chi0     +  (M10 + 2.0*M20 + 2.0*M30 + M40)/6.0; 
		    the0    = the0     +  (N10 + 2.0*N20 + 2.0*N30 + N40)/6.0;



		    lam   	 = lam     +  (I1A + 2.0*I2A + 2.0*I3A + I4A)/6.0;
			nu    	 = nu      +  (J1A + 2.0*J2A + 2.0*J3A + J4A)/6.0; 
			xi    	 = xi      +  (K1A + 2.0*K2A + 2.0*K3A + K4A)/6.0; 
			chi   	 = chi     +  (M1A + 2.0*M2A + 2.0*M3A + M4A)/6.0; 
		    the   	 = the     +  (N1A + 2.0*N2A + 2.0*N3A + N4A)/6.0;

		    eta    = eta     +   h_I;

		    loop_counter++;

  			//flag = flag_check(temp4_I,chi); if( flag != 0 ) { count++; fclose(fRT_datafile); fclose(fRT_EC_data); goto jump;}
		  
		    	
		}	

   		//[ Boundary Conditions at the surface	
   		eta_s		= temp1_I;	
   		lam_s	 	= temp2_I;
  		nu_s     	= temp3_I; 														
		chi_s   	= temp4_I;        																							
  		xi_s        = temp5_I;
  			
  		//[ Initializing the surface values of the variables 														]
  		eta 		= eta_s;
  		lam 		= lam_s; 	
  		nu     		= nu_s;
  		chi    		= chi_s;
  		xi     		= xi_s;

  		while( eta < 12.0*eta_s ) //[ EXTERIOR LOOP ============================================================]	
																													
  		{ 
          //  if(loop_counter % 100 == 0)
            {             
 		    fprintf(fRT_datafile,"%Lf\t%f\t%Lf\t%Lf\t%0.12Lf\t%0.12Lf\t%Lf\n",eta*R_scale,0.000000,lam,nu,chi,xi,mass(eta,lam));
 			}
 		printf("%Lf\t%f\t%Lf\t%Lf\t%0.12Lf\t%0.12Lf\n",eta*R_scale,0.000000,lam,nu,chi,xi);
 		
			temp1_E = eta; temp2_E = lam; temp3_E = nu, temp4_E = chi; temp5_E = xi;

			if(chi == pow(10.0,-7.0) ) {eta_H = eta;}

    		//[ Evaluation of the first slope in RK4 method in the given Starobinsky 	 						]
    		R1 =  h_E*DlamA_E(eta, lam, chi, xi);
    		S1 =  h_E*DnuA_E(eta, lam, chi, xi);
    		T1 =  h_E*DxiA_E(eta, lam, chi, xi);
    		U1 =  h_E*DchiA_E(xi);
    			

    		//[ Evaluation of the second slope in RK4 method in the given Starobinsky							]
			R2 =  h_E*DlamA_E(eta+0.5*h_E, lam+0.5*R1, chi+0.5*U1, xi+0.5*T1); 
			S2 =  h_E*DnuA_E(eta+0.5*h_E, lam+0.5*R1, chi+0.5*U1, xi+0.5*T1); 
			T2 =  h_E*DxiA_E(eta+0.5*h_E, lam+0.5*R1, chi+0.5*U1, xi+0.5*T1);
			U2 =  h_E*DchiA_E(xi+0.5*T1);
				

    		//[ Evaluation of the third slope in RK4 method in the given Starobinsky							]
    		R3 =  h_E*DlamA_E(eta+0.5*h_E, lam+0.5*R2, chi+0.5*U2, xi+0.5*T2);
    		S3 =  h_E*DnuA_E(eta+0.5*h_E, lam+0.5*R2, chi+0.5*U2, xi+0.5*T2);
    		T3 =  h_E*DxiA_E(eta+0.5*h_E, lam+0.5*R2, chi+0.5*U2, xi+0.5*T2);
    		U3 =  h_E*DchiA_E(xi+0.5*T2);
    			

    		//[ Evaluation of the forth slope in RK4 method in the given Starobinsky							]	
    		R4 =  h_E*DlamA_E(eta+h_E, lam+R3, chi+U3, xi+T3);
    		S4 =  h_E*DnuA_E(eta+h_E, lam+R3, chi+U3, xi+T3);
    		T4 =  h_E*DxiA_E(eta+h_E, lam+R3, chi+U3, xi+T3);
    		U4 =  h_E*DchiA_E(xi+T3);
    			

   
		    lam    = lam     +  (R1 + 2.0*R2 + 2.0*R3 + R4)/6.0;
			nu     = nu      +  (S1 + 2.0*S2 + 2.0*S3 + S4)/6.0; 
			xi     = xi      +  (T1 + 2.0*T2 + 2.0*T3 + T4)/6.0; 
			chi    = chi     +  (U1 + 2.0*U2 + 2.0*U3 + U4)/6.0;
		    eta    = eta     +   h_E;

		    loop_counter++;
		    
			//if( temp4_E < pow( 10, - 12.0 ) ) { break; }

		   //flag = flag_check(temp4_E,chi); if( flag != 0 ){ count++; fclose(fRT_datafile); fclose(fRT_EC_data); goto jump; }

		}

		fprintf(fRT_datafile,"%Lf\t%f\t%Lf\t%Lf\t%0.12Lf\t%0.12Lf\t%Lf\n",temp1_E*R_scale,0.000000,temp2_E,temp2_E,temp4_E,temp5_E,mass(temp1_E,temp2_E));

		fprintf(fRT_datafile,"\n");	
		fprintf(fRT_datafile," R + alpha R^2 +omega R T Gravity: [alpha: %0.3f r_g^2] [omega: %f 4B]\n",alpha*pow(r_g,-2.0), beta);
		fprintf(fRT_datafile," Central Density: [%Le  gcm^-3] [theta_c = %Lf]\n",rho_scale*the_c, the_c);
		fprintf(fRT_datafile," Stellar Radius R_s: [%Lf  km] \n", R_scale*eta_s);	
		fprintf(fRT_datafile," Stellar Mass   M_s: [%Lf  SM] \n", mass( temp1_I, temp2_I) );
		fprintf(fRT_datafile," Halo Radius    R_H: [%Lf  km] \n", R_scale*eta_H);	
		fprintf(fRT_datafile," Total Mass     M  : [%Lf  SM] ", mass( temp1_E, temp2_E) );
		fclose(fRT_datafile);
		fclose(fRT_EC_data);


		printf("\n");	
		printf(" \tStellar Radius R_s: [%Lf  km] \n", R_scale*eta_s);	
		printf(" \tStellar Mass   M_s: [%Lf  SM] \n", mass( temp1_I, temp2_I) );
		printf(" \tHalo Radius    R_H: [%Lf  km] \n", R_scale*eta_H);	
		printf(" \tTotal Mass     M  : [%Lf  SM] \n", mass( temp1_E, temp2_E) );
		printf(" \t[:: Done ::]\n\n"reset);
		printf(" [:: Done ::]");
		printf("\n");

	fprintf(main_datafile,"%0.4Lf\t%Le\t%0.12Lf\t%0.12Lf\t%0.12Lf\t%Lf\t%Lf\t%Lf\t%Lf\n",the_c,rho_scale*the_c,chi_c_Eienstein,chi0_c,chi_c,R_scale*eta_s,mass( temp1_I, temp2_I),R_scale*eta_H,mass( temp1_E, temp2_E));

	}

fclose(main_datafile);
return (0);
}
