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

	#include<stdio.h>
	#include<math.h>
	#include<string.h>
	#include<stdlib.h>
	#include<time.h>
	#include"common.h"

//[=============================================================================================================]
//[Definite the colors                                                                          				]
	#define white   "\x1B[37m"
	#define red     "\x1b[31m"
	#define green   "\x1b[32m"
	#define yellow  "\x1b[33m"
	#define blue    "\x1b[34m"
	#define violet  "\x1b[35m"
	#define cyan    "\x1b[36m"
	#define reset   "\x1b[0m"

//[=============================================================================================================]
//[Physical Constants in SI units  																				]


	#define c    	2.99792458*pow(10,8)		   	 //[  Speed of light in meter/sec      						]	
	#define pi   	3.141592654				 	 	 //[  Pi	         										]
	#define pl   	6.626070040*pow(10,-34)		 	 //[  Planck constant       								]
	#define G    	6.6743015*pow(10,-11)	 		     //[  Gravitational Constant           	 				    ]
	#define S    	1.98847*pow(10,30)  		     //[  One Solar Mass in the units of kilograms				]
	#define B       9.613059804*pow(10,33.0)		 //[  The Bag constant in the units of Joule/meter^3		]	


//[=============================================================================================================]
//[ Constants used in the code ]

																							
	#define r_g		G*S*pow(c,-2.0)
	#define alpha   10.0*pow(r_g,2.0)
	#define kappa	8*pi*G*pow(c,-4.0)	
	#define q1 		4.0*kappa*B*pow(r_g,2.0)
	#define q2		pow(r_g,2.0)*pow(6.0*alpha,-1.0)
	#define k 		0.28

//[ Scaling constants 																							]
	#define rho_scale   4.0*B*pow(c,-2.0)*pow(10,-3.0)		//[ Density scale in the units of gm/cm^3 			]
	#define R_scale     G*S*pow(c,-2.0)*pow(10,-3.0)		//[ Radial scale in the units of km 				]



//[=============================================================================================================]

//[INTERIOR FIELD EQUATIONS FOR STAROBINSKY GRAVITY		 														]

long double gamma_I (long double eta, long double chi, long double xi)
{
	//[ Gamma function]
   long double RV = pow( 6.0*q2 , -1.0  )*eta*pow( phi(chi) , -1.0)*xi;
   return RV;
}

long double Dnu_I (long double eta, long double the, long double lam, long double chi, long double xi)
	{
		long double Eq_1 =   q1*eta*exp(lam)*pressure(the)*pow( phi(chi) , -1.0 );
		long double Eq_2 =   ( exp(lam) -1.0 )/eta;
		long double Eq_3 = - eta*exp(lam)*pow( chi, 2.0 )*pow( 12.0*q2*phi(chi) , -1.0 );
		long double Eq_4 = - 4.0*gamma_I(eta, chi, xi)/eta;

		long double RV   =   (Eq_1 + Eq_2 + Eq_3 + Eq_4)*pow( 1.0 + gamma_I(eta, chi, xi) , -1.0 );

	return RV;
	}

long double Dlam_I (long double eta, long double the, long double lam, long double chi, long double xi)
	{
		long double Eq_1 =   ( 1.0 - exp(lam) )/eta ; 
		long double Eq_2 =   eta*exp(lam)*pow( 6.0*phi(chi), -1.0 )*( 2.0 + 0.5*chi/q2 )*chi;
		long double Eq_3 = - gamma_I(eta, chi, xi)*Dnu_I(eta, the, lam, chi, xi);
		long double Eq_4 =   q1*eta*exp(lam)*( 3.0*epsilon(the) + trace(the) )*pow( 3.0*phi(chi), -1.0 );

		long double RV   =   Eq_1 + Eq_2 + Eq_3 + Eq_4;

	return RV;
	}

long double Dchi_I (long double xi)
	{
		long double RV = xi;      

		return RV;
	}

long double Dxi_I (long double eta, long double the, long double lam, long double chi, long double xi)
	{
		long double Eq_1 = - 2.0*xi/eta + q1*q2*trace(the)*exp(lam) + q2*exp(lam)*chi;
		long double Eq_2 =   0.5*( Dlam_I(eta, the, lam, chi, xi) - Dnu_I(eta, the, lam, chi, xi) )*xi;

		long double RV   =   Eq_1 + Eq_2 ;

		return RV;
	}

long double Dthe_I (long double eta, long double the, long double lam, long double chi, long double xi)
	{
		long double RV = - 0.5*( the + pressure(the) )*pow(k,-1.0)*Dnu_I(eta, the, lam, chi, xi); 

		return RV;
	}

//[=============================================================================================================]


//[EXTERIOR FIELD EQUATIONS FOR STAROBINSKY GRAVITY		 														]

long double gamma_E (long double eta, long double chi, long double xi)
{
	//[ Gamma function]
   long double RV = pow( 6.0*q2 , -1.0  )*eta*pow( phi(chi) , -1.0)*xi;
   return RV;
}


long double Dnu_E (long double eta, long double lam, long double chi, long double xi)
	{
		long double Eq_1 =   ( exp(lam) - 1.0 )/eta;
		long double Eq_2 = - eta*exp(lam)*pow( chi, 2.0 )*pow( 12.0*q2*phi(chi) , -1.0 );
		long double Eq_3 = - 4.0*gamma_E(eta, chi, xi )/eta;

		long double RV   =   (Eq_1 + Eq_2 + Eq_3)*pow( 1.0 + gamma_E(eta, chi, xi ) , -1.0 );

	return RV;
	}

long double Dlam_E (long double eta, long double lam, long double chi, long double xi)
	{
		long double Eq_1 =   ( 1.0 - exp(lam) )/eta ; 
		long double Eq_2 =   eta*exp(lam)*pow( 6.0*phi(chi), -1.0 )*( 2.0 + 0.5*chi/q2 )*chi;
		long double Eq_3 = - gamma_E(eta, chi, xi)*Dnu_E(eta, lam, chi, xi);
		
		long double RV   =   Eq_1 + Eq_2 + Eq_3;

	return RV;
	}

long double Dchi_E (long double xi)
	{
		long double RV = xi;      

		return RV;
	}

long double Dxi_E (long double eta, long double lam, long double chi, long double xi)
	{
		long double Eq_1 = - 2.0*xi/eta + q2*exp(lam)*chi;
		long double Eq_2 =   0.5*( Dlam_E(eta, lam, chi, xi) - Dnu_E(eta, lam, chi, xi) )*xi;

		long double RV   =   Eq_1 + Eq_2 ;

		return RV;
	}

//[=============================================================================================================]
//[ Main Function: 																								]

long double Starobinsky_chi_c ( long double the_c )
//[ This fcuntion evaluate the correct valye of R_c in the ure Starobinsky for a given input density  			]
{

  		printf("\n");	
  		printf(cyan" \t[:: Calculating chi_c value in Starobinsky Gravity ::]\n");
  		printf(" \tParameter: [alpha: %0.3f r_g^2]\n",alpha*pow(r_g,-2.0));
		printf(" \tCentral Density: [%Le  gcm^-3][theta_c: %0.3Lf]\n",rho_scale*the_c, the_c);
		

  	//[ Declaration of the central variables 																	]
  		long double lam_c;
  		long double nu_c; 
  		long double chi_c; 			
  		long double xi_c;                    		        

  	//[ Integration variables of INTIRIOR LOOP																	]
  		long double h_I			 = 0.0001;           //[ Interior Step Size   									]
  		long double eta;            						
  		long double I1, I2, I3, I4;
		long double J1, J2, J3, J4;
		long double K1, K2, K3, K4;
		long double M1, M2, M3, M4;
		long double N1, N2, N3, N4;
		long double the, lam, nu, chi, xi;

	//[ Declaration of the surface variables 																	]
		long double eta_s;															
  		long double lam_s;
  		long double nu_s; 
  		long double chi_s; 			
  		long double xi_s;     

  		long double eta_H;          //[Halo radius]                 		        

  	//[ Integration variables of exterior loop																	]
  		long double h_E			 = 0.0001;           //[ Exterior Step Size   									]         						
  		long double R1, R2, R3, R4;
		long double S1, S2, S3, S4;
		long double T1, T2, T3, T4;
		long double U1, U2, U3, U4;
	

	//[ Temporary Varibale 																						]
	 	long double temp_1_I, temp_2_I, temp_3_I, temp_4_I, temp_5_I;
		long double temp_1_E, temp_2_E, temp_3_E, temp_4_E, temp_5_E;


	//[ To reduce the file size 																			]
 		int loop_counter;
 

//[=============================================================================================================]
	//[ Boundary Conditions at the center	

  		int flag;
  		int count;
  		long double chi_c_Eienstein;
		long double chi_c_1, chi_c_2;
		long double chi_c_intermediate;

  			flag = 0;
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
  			printf("\r \tTuning chi_c value: [%0.12Lf] ["red"%d"cyan"]", chi_c, count);
			fflush(stdout);
			
																       														
  		lam_c	 	=   0.0;
  		nu_c     	=   1.0; 			//[ The central value of nu is arbitrary, we shall choose nu_c = 1.0. 	]
  		xi_c        =   0.0;

  		eta_H 		= 	0.0;

//[End of boundary condition ===================================================================================]

  	//[initializing the loop counter																			]
  		loop_counter = 0;

  	//[ Declaration of data files 

 		char title1[50];																				
  		sprintf(title1,"./STAROB_DATA/STAROBINSKY[%0.3Lf].dat",the_c);
		FILE* datafile  = fopen(title1,"w");

		char title2[50];																				
  		sprintf(title2,"./STAROB_EC/STAROB_EC[%0.3Lf].dat",the_c);
		FILE* EC_data  = fopen(title2,"w");

  	
  	//[ Initializing the central values of the variables 														]
  		eta 		= h_I;
		the  		= the_c;
  		lam 		= lam_c; 	
  		nu   	  	= nu_c;
  		chi 	   	= chi_c;
  		xi      	= xi_c;


  		fprintf(datafile,"%Lf\t%Lf\t%Lf\t%Lf\t%0.12Lf\t%0.12Lf\t%Lf\n",eta*R_scale,the,lam,nu,chi,xi,mass(eta,lam));
 		fprintf(EC_data,"%Lf\t%Lf\t%Lf\t%Lf\t%Lf\n",eta*R_scale,the,the+pressure(the),the-pressure(the),the+3.0*pressure(the));
					
  	while( the > 1.0 ) //[ INTERIOR LOOP =======================================================================]
  																	
  		{ 

  			if(loop_counter % 100 == 0)
        	{    	   
 			fprintf(datafile,"%Lf\t%Lf\t%Lf\t%Lf\t%0.12Lf\t%0.12Lf\t%Lf\n",eta*R_scale,the,lam,nu,chi,xi,mass(eta,lam));
 			fprintf(EC_data,"%Lf\t%Lf\t%Lf\t%Lf\t%Lf\n",eta*R_scale,the,the+pressure(the),the-pressure(the),the+3.0*pressure(the));
 			}
 		//printf("%Lf\t%Lf\t%Lf\t%Lf\t%0.12Lf\t%0.12Lf\n",eta*R_scale,the,lam,nu,chi,xi);

 			temp_1_I = eta; temp_2_I = lam; temp_3_I = nu, temp_4_I = chi; temp_5_I = xi;

    		//[ Evaluation of the first slope in RK4 method in the given Starobinsky 	 						]
    		I1 =  h_I*Dlam_I(eta, the, lam, chi, xi);
    		J1 =  h_I*Dnu_I(eta, the, lam, chi, xi);
    		K1 =  h_I*Dxi_I(eta, the, lam, chi, xi);
    		M1 =  h_I*Dchi_I(xi);
    		N1 =  h_I*Dthe_I(eta, the, lam, chi, xi);

    		//[ Evaluation of the second slope in RK4 method in the given Starobinsky							]
			I2 =  h_I*Dlam_I(eta+0.5*h_I, the+0.5*N1, lam+0.5*I1, chi+0.5*M1, xi+0.5*K1); 
			J2 =  h_I*Dnu_I(eta+0.5*h_I, the+0.5*N1, lam+0.5*I1, chi+0.5*M1, xi+0.5*K1); 
			K2 =  h_I*Dxi_I(eta+0.5*h_I, the+0.5*N1, lam+0.5*I1, chi+0.5*M1, xi+0.5*K1);
			M2 =  h_I*Dchi_I(xi+0.5*K1);
			N2 =  h_I*Dthe_I(eta+0.5*h_I, the+0.5*N1, lam+0.5*I1, chi+0.5*M1, xi+0.5*K1);

    		//[ Evaluation of the third slope in RK4 method in the given Starobinsky							]
    		I3 =  h_I*Dlam_I(eta+0.5*h_I, the+0.5*N2, lam+0.5*I2, chi+0.5*M2, xi+0.5*K2);
    		J3 =  h_I*Dnu_I(eta+0.5*h_I, the+0.5*N2, lam+0.5*I2, chi+0.5*M2, xi+0.5*K2);
    		K3 =  h_I*Dxi_I(eta+0.5*h_I, the+0.5*N2, lam+0.5*I2, chi+0.5*M2, xi+0.5*K2);
    		M3 =  h_I*Dchi_I(xi+0.5*K2);
    		N3 =  h_I*Dthe_I(eta+0.5*h_I, the+0.5*N2, lam+0.5*I2, chi+0.5*M2, xi+0.5*K2);

    		//[ Evaluation of the forth slope in RK4 method in the given Starobinsky							]	
    		I4 =  h_I*Dlam_I(eta+h_I, the+N3, lam+I3, chi+M3, xi+K3);
    		J4 =  h_I*Dnu_I(eta+h_I, the+N3, lam+I3, chi+M3, xi+K3);
    		K4 =  h_I*Dxi_I(eta+h_I, the+N3, lam+I3, chi+M3, xi+K3);
    		M4 =  h_I*Dchi_I(xi+K3);
    		N4 =  h_I*Dthe_I(eta+h_I, the+N3, lam+I3, chi+M3, xi+K3);

   
		    lam    = lam     +  (I1 + 2.0*I2 + 2.0*I3 + I4)/6.0;
			nu     = nu      +  (J1 + 2.0*J2 + 2.0*J3 + J4)/6.0; 
			xi     = xi      +  (K1 + 2.0*K2 + 2.0*K3 + K4)/6.0; 
			chi    = chi     +  (M1 + 2.0*M2 + 2.0*M3 + M4)/6.0; 
		    the    = the     +  (N1 + 2.0*N2 + 2.0*N3 + N4)/6.0;
		    eta    = eta     +   h_I;

		    loop_counter++;

		    flag = flag_check(temp_4_I,chi); if( flag != 0 ) { count++; fclose(datafile); fclose(EC_data); goto jump;} //[check]
		    	
		}	

   	//[ Boundary Conditions at the surface	
   		eta_s		= temp_1_I;	
   		lam_s	 	= temp_2_I;
  		nu_s     	= temp_3_I; 														
		chi_s   	= temp_4_I;        																							
  		xi_s        = temp_5_I;
  			
  	//[ Initializing the surface values of the variables 														]
  		eta 		= eta_s;
  		lam 		= lam_s; 	
  		nu     		= nu_s;
  		chi    		= chi_s;
  		xi     		= xi_s;

  	
  	while( eta < 15.0*eta_s ) //[ EXTERIOR LOOP ================================================================]	
																													
  		{ 

        if(loop_counter % 100 == 0)           
 		fprintf(datafile,"%Lf\t%f\t%Lf\t%Lf\t%0.12Lf\t%0.12Lf\t%Lf\n",eta*R_scale,0.000000,lam,nu,chi,xi,mass(eta,lam));	
 		//printf("%Lf\t%f\t%Lf\t%Lf\t%0.12Lf\t%0.12Lf\n",eta*R_scale,0.000000,lam,nu,chi,xi);
 		
			temp_1_E = eta; temp_2_E = lam; temp_3_E = nu, temp_4_E = chi; temp_5_E = xi;

			if(chi == pow(10.0,-7.0) ) {eta_H = eta;}

    		//[ Evaluation of the first slope in RK4 method in the given Starobinsky 	 						]
    		R1 =  h_E*Dlam_E(eta, lam, chi, xi);
    		S1 =  h_E*Dnu_E(eta, lam, chi, xi);
    		T1 =  h_E*Dxi_E(eta, lam, chi, xi);
    		U1 =  h_E*Dchi_E(xi);
    			

    		//[ Evaluation of the second slope in RK4 method in the given Starobinsky							]
			R2 =  h_E*Dlam_E(eta+0.5*h_E, lam+0.5*R1, chi+0.5*U1, xi+0.5*T1); 
			S2 =  h_E*Dnu_E(eta+0.5*h_E, lam+0.5*R1, chi+0.5*U1, xi+0.5*T1); 
			T2 =  h_E*Dxi_E(eta+0.5*h_E, lam+0.5*R1, chi+0.5*U1, xi+0.5*T1);
			U2 =  h_E*Dchi_E(xi+0.5*T1);
				

    		//[ Evaluation of the third slope in RK4 method in the given Starobinsky							]
    		R3 =  h_E*Dlam_E(eta+0.5*h_E, lam+0.5*R2, chi+0.5*U2, xi+0.5*T2);
    		S3 =  h_E*Dnu_E(eta+0.5*h_E, lam+0.5*R2, chi+0.5*U2, xi+0.5*T2);
    		T3 =  h_E*Dxi_E(eta+0.5*h_E, lam+0.5*R2, chi+0.5*U2, xi+0.5*T2);
    		U3 =  h_E*Dchi_E(xi+0.5*T2);
    			

    		//[ Evaluation of the forth slope in RK4 method in the given Starobinsky							]	
    		R4 =  h_E*Dlam_E(eta+h_E, lam+R3, chi+U3, xi+T3);
    		S4 =  h_E*Dnu_E(eta+h_E, lam+R3, chi+U3, xi+T3);
    		T4 =  h_E*Dxi_E(eta+h_E, lam+R3, chi+U3, xi+T3);
    		U4 =  h_E*Dchi_E(xi+T3);
    			

   
		    lam    = lam     +  (R1 + 2.0*R2 + 2.0*R3 + R4)/6.0;
			nu     = nu      +  (S1 + 2.0*S2 + 2.0*S3 + S4)/6.0; 
			xi     = xi      +  (T1 + 2.0*T2 + 2.0*T3 + T4)/6.0; 
			chi    = chi     +  (U1 + 2.0*U2 + 2.0*U3 + U4)/6.0;
		    eta    = eta     +   h_E;

		    loop_counter++;
		    
			if( temp_4_E < pow( 10, - 12.0 ) ) { break; }

		    flag = flag_check(temp_4_E,chi); if( flag != 0 ){ count++; fclose(datafile); fclose(EC_data); goto jump; }

		}


		fprintf(datafile,"%Lf\t%f\t%Lf\t%Lf\t%0.12Lf\t%0.12Lf\t%Lf\n",temp_1_E*R_scale,0.000000,temp_2_E,temp_3_E,temp_4_E,temp_5_E,mass(temp_1_E,temp_2_E));

		fprintf(datafile,"\n");	
		fprintf(datafile," R + alpha R^2 Gravity: [alpha: %0.3f r_g^2]\n",alpha*pow(r_g,-2.0));
		fprintf(datafile," Central Density: [%Le  gcm^-3] [theta_c = %Lf]\n",rho_scale*the_c, the_c);
		fprintf(datafile," Stellar Radius R_s: [%Lf  km] \n", R_scale*eta_s);	
		fprintf(datafile," Stellar Mass   M_s: [%Lf  SM] \n", mass( temp_1_I, temp_2_I) );
		fprintf(datafile," Halo Radius    R_H: [%Lf  km] \n", R_scale*eta_H);	
		fprintf(datafile," Total Mass     M  : [%Lf  SM] ", mass( temp_1_E, temp_2_E) );
		fclose(datafile);
		fclose(EC_data);


		printf("\n");	
		printf(" \tStellar Radius R_s: [%Lf  km] \n", R_scale*eta_s);	
		printf(" \tStellar Mass   M_s: [%Lf  SM] \n", mass( temp_1_I, temp_2_I) );
		printf(" \tHalo Radius    R_H: [%Lf  km] \n", R_scale*temp_1_E);	
		printf(" \tTotal Mass     M  : [%Lf  SM] \n", mass( temp_1_E, temp_2_E) );
		printf(" \t[:: Done ::]"reset);
		printf("\n");


return chi_c;
}
