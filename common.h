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
//[ Common functions ]  

    /* The arguments in follwing functions define whether they are related to perturbed or 	
	   unperturbed fucntions 
	*/   										 					


long double trace (long double the )         
	{
		//[ Returns the Trace of Energy Momemtum Tensor for the given equation of state 						] 
   		long double RV = - 0.16*the - 0.84;
   		return RV;										//[ RV represents the return value 						]
	}

long double pressure (long double the )		
	{
		//[ Returns the value of dimensionless pressure of the given matter     								]
   		long double RV = 0.28*( the - 1.0 );
   		return RV;
	}

long double epsilon (long double the )    
	{
		//[ Returns the diemntionless energy density for the  matter    										]
   		long double RV = the;
   		return RV;
	}
	
long double phi (long double chi)			 
	{
		//[ Returns the value ofphi(R) = 1 + 2*alpha*R 															]
   		long double RV = 1.0 + pow( 3.0*q2, -1.0 )*chi;
   		return RV;
	}

long double mass (long double eta, long double lam)
	{
		//[Returns the mass at a radial distance r in the units of Solar masses.								]
		long double RV = 0.5*eta*( 1.0 - exp(-lam) );
   		return RV;													
	}


int flag_check (long double chi_1, long double chi_2)
	{
		int flag = 0;

		if( chi_2 < 0.0 )         	{flag = 1;}
		if( chi_2 - chi_1 > 0.0 )   {flag = 2;}

		return flag;
	}



//[=============================================================================================================]

    /* The following fucntions are defined only for the actual field equations 
       of the f(R,T) gravity.
	*/

	/* To determine these fcuntions we need to redefine the interior field 
	   equations for the starobinsky
	*/   

    //[  Interior field equations for starobinsky gravity, RD stands for redefinition	 							     					]

 
long double gamma_I_RD (long double eta, long double chi, long double xi)
{
	//[ Gamma function]
   long double RV = pow( 6.0*q2 , -1.0  )*eta*pow( phi(chi) , -1.0)*xi;
   return RV;
}

long double Dnu_I_RD (long double eta, long double the, long double lam, long double chi, long double xi)
	{
		long double Eq_1 =   q1*eta*exp(lam)*pressure(the)*pow( phi(chi) , -1.0 );
		long double Eq_2 =   ( exp(lam) -1.0 )/eta;
		long double Eq_3 = - eta*exp(lam)*pow( chi, 2.0 )*pow( 12.0*q2*phi(chi) , -1.0 );
		long double Eq_4 = - 4.0*gamma_I_RD(eta, chi, xi)/eta;

		long double RV   =   (Eq_1 + Eq_2 + Eq_3 + Eq_4)*pow( 1.0 + gamma_I_RD(eta, chi, xi) , -1.0 );

	return RV;
	}

long double Dlam_I_RD (long double eta, long double the, long double lam, long double chi, long double xi)
	{
		long double Eq_1 =   ( 1.0 - exp(lam) )/eta ; 
		long double Eq_2 =   eta*exp(lam)*pow( 6.0*phi(chi), -1.0 )*( 2.0 + 0.5*chi/q2 )*chi;
		long double Eq_3 = - gamma_I_RD(eta, chi, xi)*Dnu_I_RD(eta, the, lam, chi, xi);
		long double Eq_4 =   q1*eta*exp(lam)*( 3.0*epsilon(the) + trace(the) )*pow( 3.0*phi(chi), -1.0 );

		long double RV   =   Eq_1 + Eq_2 + Eq_3 + Eq_4;

	return RV;
	}

long double Dchi_I_RD (long double xi)
	{
		long double RV = xi;      

		return RV;
	}

long double Dxi_I_RD (long double eta, long double the, long double lam, long double chi, long double xi)
	{
		long double Eq_1 = - 2.0*xi/eta + q1*q2*trace(the)*exp(lam) + q2*exp(lam)*chi;
		long double Eq_2 =   0.5*( Dlam_I_RD(eta, the, lam, chi, xi) - Dnu_I_RD(eta, the, lam, chi, xi) )*xi;

		long double RV   =   Eq_1 + Eq_2 ;

		return RV;
	}

long double Dthe_I_RD (long double eta, long double the, long double lam, long double chi, long double xi)
	{
		long double RV = - 0.5*( the + pressure(the) )*pow(k,-1.0)*Dnu_I_RD(eta, the, lam, chi, xi); 

		return RV;
	}



    /* The first derivative and second derivative of trace of the given matter is 
        evaluted as follows.
	*/


long double Dtrace (long double eta, long double the, long double lam, long double chi, long double xi)         
	{

		long double Dthe;
		long double RV;

		Dthe 	= Dthe_I_RD(eta,the,lam,chi,xi);
		RV 		= - 0.16*Dthe;

	return RV;			
	}

long double DDtrace (long double eta, long double the, long double lam, long double chi, long double xi,
					 long double h)         
	{

		long double DDthe;
		long double Dthe_1, Dthe_2;
		long RV;
          						
  		long double I1, I2, I3, I4;
		long double J1, J2, J3, J4;
		long double K1, K2, K3, K4;
		long double M1, M2, M3, M4;
		long double N1, N2, N3, N4;

		Dthe_1 = Dthe_I_RD(eta, the, lam, chi, xi);

		//[ Evaluation of the first slope in RK4 method in the given Starobinsky 	 						    ]
    		I1 =  h*Dlam_I_RD(eta, the, lam, chi, xi);
    		K1 =  h*Dxi_I_RD(eta, the, lam, chi, xi);
    		M1 =  h*Dchi_I_RD(xi);
    		N1 =  h*Dthe_I_RD(eta, the, lam, chi, xi);

    		//[ Evaluation of the second slope in RK4 method in the given Starobinsky							]
			I2 =  h*Dlam_I_RD(eta+0.5*h, the+0.5*N1, lam+0.5*I1, chi+0.5*M1, xi+0.5*K1); 
			K2 =  h*Dxi_I_RD(eta+0.5*h, the+0.5*N1, lam+0.5*I1, chi+0.5*M1, xi+0.5*K1);
			M2 =  h*Dchi_I_RD(xi+0.5*K1);
			N2 =  h*Dthe_I_RD(eta+0.5*h, the+0.5*N1, lam+0.5*I1, chi+0.5*M1, xi+0.5*K1);

    		//[ Evaluation of the third slope in RK4 method in the given Starobinsky							]
    		I3 =  h*Dlam_I_RD(eta+0.5*h, the+0.5*N2, lam+0.5*I2, chi+0.5*M2, xi+0.5*K2);
    		K3 =  h*Dxi_I_RD(eta+0.5*h, the+0.5*N2, lam+0.5*I2, chi+0.5*M2, xi+0.5*K2);
    		M3 =  h*Dchi_I_RD(xi+0.5*K2);
    		N3 =  h*Dthe_I_RD(eta+0.5*h, the+0.5*N2, lam+0.5*I2, chi+0.5*M2, xi+0.5*K2);

    		//[ Evaluation of the forth slope in RK4 method in the given Starobinsky							]	
    		I4 =  h*Dlam_I_RD(eta+h, the+N3, lam+I3, chi+M3, xi+K3);
    		K4 =  h*Dxi_I_RD(eta+h, the+N3, lam+I3, chi+M3, xi+K3);
    		M4 =  h*Dchi_I_RD(xi+K3);
    		N4 =  h*Dthe_I_RD(eta+h, the+N3, lam+I3, chi+M3, xi+K3);

   
		    lam    = lam     +  (I1 + 2.0*I2 + 2.0*I3 + I4)/6.0;
			xi     = xi      +  (K1 + 2.0*K2 + 2.0*K3 + K4)/6.0; 
			chi    = chi     +  (M1 + 2.0*M2 + 2.0*M3 + M4)/6.0; 
		    the    = the     +  (N1 + 2.0*N2 + 2.0*N3 + N4)/6.0;
		   

		   Dthe_2 = Dthe_I_RD(eta, the, lam, chi, xi);
		   DDthe  = ( Dthe_2 - Dthe_1 )/h;
		   RV 	  = - 0.16*DDthe;

		return RV;
	}




