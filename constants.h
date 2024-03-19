#ifndef constants_h
#define constants_h 1
#include <Rtypes.h>
#include <TString.h>
#include <TMath.h>

namespace cs
{
	//*************************DETS DISPLACEMENT*******************************************
	static const float		widthStripX	=	1.8125;
	static const float		widthStripY	=	3.625;

	static const float		sqlXzero	=	widthStripX * 15.5;//28.09375;
	static const float		sqlYzero	=	-widthStripY * 7.5;//-27.1875;

	static const float		sqlXstart	=	sqlXzero;
	static const float		sqlYstart	=	sqlYzero;

	static const float		sqrXzero	=	widthStripX * 15.5;//28.09375;
	static const float		sqrYzero	=	-widthStripY * 7.5;//-27.1875;

	static const float		sqrXstart	=	sqrXzero;
	static const float		sqrYstart	=	sqrYzero;

	//************************* CUTS FOR FOR 0-2*******************************************
	static const float		tc_SQX_L	=	200.0;
	static const float		tc_SQX_R	=	200.0;
	static const float		tc_SQY_L	=	200.0;
	static const float		tc_SQY_R	=	200.0;
	static const float		tc_CsI_L	=	1500.0;
	static const float		tc_CsI_R	=	1500.0;

	static const float		ec_SQX_L_5	=	1.5;	//a lot of dirt below 1.5
	static const float		ec_SQY_L_5	=	1.5;	//a lot of dirt below 1.5
	static const float		ec_CsI_L_5	=	0.5;	//it's actually SSD_L 

	static const float		ec_SQX_R_5	=	5.0;	//since it's for He6, I can use 5 MeV
	static const float		ec_SQY_R_5	=	5.0;	//since it's for He6, I can use 5 MeV
	static const float		ec_CsI_R_5	=	150.0;

	static const float		sqRDeadLayer = 3.45;
	static const float		sqLDeadLayer = 3.5;

	//*************************MWPC******************************************************

	static const float	MWPC1_X_displacement = 0.0;
	static const float	MWPC1_Y_displacement = 0.0;
	static const float	MWPC2_X_displacement = 0.0;
	static const float	MWPC2_Y_displacement = 0.0;

	static const float	leftDetShift = 0.0;
	static const float	rightDetShift = 0.0;
	static const float	leftDetShiftX = 0.0;
	static const float	rightDetShiftX = 0.0;

	static const float	leftDetDist = 0.0;
	static const float	rightDetDist = 0.0;

	static const float	leftAngShift = 0.0;
	static const float	rightAngShift = 0.0;

	static const float	tarThicknessShift = 36.8251;
	static const float	tarPos = 0.0;
	static const float	tarAngle = 0.0;


	static const float		MWPC1_X_zero_position	=	15.5*1.25;
	static const float		MWPC1_Y_zero_position	=	-15.5*1.25;
	static const float		MWPC2_X_zero_position	=	15.5*1.25;
	static const float		MWPC2_Y_zero_position	=	-15.5*1.25;

	static const short		MWPC_1_X_id	=	0;
	static const short		MWPC_1_Y_id	=	1;
	static const short		MWPC_2_X_id	=	2;
	static const short		MWPC_2_Y_id	=	3;

	//*************************PHYSICS**************************************************
	static const float		alpha_from_Ra226[4]{4.751, 5.459, 5.972, 7.661};
	static const float		u_to_MeV	=	931.494028;

	static const float		massN	= 939.565;
	static const float 		mass1H 	= 938.782;
	static const float 		mass2H 	= 1876.12;
	static const float 		mass3H 	= 2809.43;
	static const float 		mass3He = 2809.4;
	static const float 		mass4He = 3728.4;
	static const float 		mass5He = 4668.7;
	static const float 		mass6He = 5606.556;
	static const float 		mass7He = 6546.53;
	static const float		Qpt		= (mass1H + mass6He) - (mass3H + mass4He);
	static const float		Qdt		= (mass2H + mass6He) - (mass3H + mass5He);

	static const float		c			=	299.792;	// mm/ns
	
	static const float		tofBase			=	12320.0;
	static const float		dist_Tar_to_F5	=	-953.0;
	static const float		dist_Tar_to_F6	=	478.0;
	static const float		tof_const		=	89.165;
	static const float		tof_const_5		=	68.475;


	static const TString 	dir_runs("/home/zalewski/dataTmp/");
	static const TString 	dir_CsI("/home/zalewski/data/he6_d/miscroot/CsI/parts");
	static const TString 	dir_params("/home/zalewski/dataTmp/calibrationParameters/geo123/");	
	static const TString 	s_inFname("run00_12");
	static const TString	dir_gcut("/home/zalewski/aku/wrk/GCuts/");

	static const int		runNo = 2;
	static const TString	inDir = "raw";
	static const bool 		fixedMWPC = false;
	static const int		tarMass = 2;
}
#endif