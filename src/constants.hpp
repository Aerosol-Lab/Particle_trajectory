#pragma once

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <random>
#include <stdlib.h>
#include <time.h>

using namespace std;

//	Constant
const double kb=1.38e-23;
const double e=1.6e-19;

//	Constant (Air)
const double myu0_Air=1.81e-5;  // Sutherland's equation parameter
const double ST0_Air=288.15;    // Sutherland's equation parameter
const double SC0_Air=110.4;     // Sutherland's equation parameter
const double lamda0_Air=67.0e-9;// Mean free path at normal condition
const double mgas_Air=(28*0.8+32*0.2)*0.001/6.02e23;  // Mass of a molecule
const double Cp_Air=1.00;       // Specific heat at constant pressure
const double Cv_Air=0.718;      // Specific heat at constant volume

//	Constant (He)
const double myu0_He=3.93267e-5;  // Sutherland's equation parameter
const double ST0_He=809;          // Sutherland's equation parameter
const double SC0_He=147.2;        // Sutherland's equation parameter
const double lamda0_He=285.0e-9;  // Sutherland's equation parameter
const double mgas_He=(4.0)*0.001/6.02e23;  // Mass of a molecule
const double Cp_He=5.1926;       // Specific heat at constant pressure
const double Cv_He=3.1156;       // Specific heat at constant volume

const double myu0=myu0_Air;
const double ST0=ST0_Air;
const double SC0=SC0_Air;
const double lamda0=lamda0_Air;
const double mgas=mgas_Air;
const double Cp=Cp_Air;
const double Cv=Cv_Air;
const double gam=Cp/Cv;

/*const double myu0=myu0_He;
const double ST0=ST0_He;
const double SC0=SC0_He;
const double lamda0=lamda0_He;
const double mgas=mgas_He;
const double Cp=Cp_He;
const double Cv=Cv_He;
const double gam=Cp/Cv;*/

const double lamda_coeff=sqrt(M_PI*mgas*0.5/kb);
const double Kn_coeff=sqrt(M_PI*gam*0.5);
const double gamkb_m=gam*kb/mgas;
const double meshScale=1e-7;

// Parameters for slip correction
const double A1=2.514;
const double A2=0.8;
const double A3=0.55;
