R__LOAD_LIBRARY(libgsl.so);
R__LOAD_LIBRARY(/home/zalewski/aku/ELC/build/libEloss.so);

#include "/home/zalewski/aku/ELC/AELC.h"
#include "/home/zalewski/aku/ELC/ELC.h"

void mini()
{
    AELC *li7_Ne, *he6_Si;
    int nEl=1;
    double elA[nEl], elZ[nEl], elW[nEl];
    double SiA[nEl], SiZ[nEl], SiW[nEl];

    SiA[0] = 14;
	SiZ[0] = 7;
	SiW[0] = 1;
    he6_Si = new ELC(6, 2, nEl, 2.35, SiA, SiZ, SiW, 200,1500);

	elA[0] = 20;
	elZ[0] = 10;
	elW[0] = 1;
    li7_Ne = new ELC(7, 3, nEl, 2.74e-9, elA, elZ, elW, 1.0,500);



    double startEnergy = 7.0*0.0145;
    double endEnergy = li7_Ne->GetE(startEnergy, 100.0);

    //printf("%f\n", he6_Si->GetE(150.0, 500.0));

    printf("End energy: %f[keV]\tEnergy difference: [keV]%f\n", endEnergy*1000.0, (startEnergy - endEnergy)*1000.0);
}