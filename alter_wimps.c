#include "src/include.h"

/*-----------------------------------------------------------*/
/* Calculation of the observables in the light WIMP scenario */
/*-----------------------------------------------------------*/

int main(int argc,char** argv)
{
        struct relicparam paramrelic;
        double ratioH[NNUC+1],sigma_ratioH[NNUC+1];
        double H2_H,He3_H,Yp,Li7_H,Li6_H,Be7_H;
        double sigma_H2_H,sigma_He3_H,sigma_Yp,sigma_Li7_H,sigma_Li6_H,sigma_Be7_H;
        double m_chi;
        int type_wimp, SMC_wimp, beta, selfConjugate;

        if(argc<5)
        {
                printf(" This program needs 3 parameters:\n"
                       "   type_wimp   type of WIMP\n"
                       "                   1 - real scalar\n"
                       "                   2 - complex scalar\n"
                       "                   3 - Majorana fermion\n"
                       "                   4 - Dirac fermion\n"
                       "   SMC_wimp    SM coupling of light WIMP\n"
                       "                   1 - coupled to SM neutrinos\n"
                       "                   2 - electromagnetically coupled\n"
                       "                   3 - coupled to SM and equivalent neutrinos\n"
                       "   m_chi       mass of light WIMP (in MeV)\n");
                exit(1);
        }
        else
        {
                sscanf(argv[1],"%d",&type_wimp);
                sscanf(argv[2],"%d",&SMC_chi);
                sscanf(argv[3],"%lf",&m_chi);
        }

        if (type_wimp == 1) // Real scalar
        {
            beta = 0;           // Boson
            selfConjugate = 1;  // True
        }
        else if (config.type_wimp == 2) // Complex scalar
        {
            beta = 0;
            selfConjugate = 0;
        }
        else if (config.type_wimp == 3) // Majorana fermion
        {
            beta = 1;           // Fermion
            selfConjugate = 1;
        }
        else if (config.type_wimp == 4) // Dirac fermion
        {
            beta = 1;
            selfConjugate = 0;
        }
        else
        {
            printf("\t [ERROR] Incorrect input value for parameter 'type_wimp'.");
            exit(1);
        }

        Init_cosmomodel(&paramrelic);

        Init_wimp(m_chi, SMC_wimp, beta, selfConjugate, ,&paramrelic);

        printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
        nucl(2,paramrelic,ratioH);
        H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
        printf("  low:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);

        nucl(0,paramrelic,ratioH);
        H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
        printf(" cent:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);

        nucl(1,paramrelic,ratioH);
        H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
        printf(" high:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);

        if(nucl_witherrors(3,paramrelic,ratioH,sigma_ratioH))
        {
                printf("With uncertainties:\n");
        H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
                sigma_H2_H=sigma_ratioH[3];sigma_Yp=sigma_ratioH[6];sigma_Li7_H=sigma_ratioH[8];sigma_Be7_H=sigma_ratioH[9];sigma_He3_H=sigma_ratioH[5];sigma_Li6_H=sigma_ratioH[7];
                printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");

                printf("value:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);
                printf(" +/- :\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n\n",sigma_Yp,sigma_H2_H,sigma_He3_H,sigma_Li7_H,sigma_Li6_H,sigma_Be7_H);
        }

        int compat=bbn_excluded(0,paramrelic);

        if(compat==1) printf("Excluded by BBN constraints\n");
        else if(compat==0) printf("Compatible with BBN constraints\n");
        else printf("Computation failed\n");

        return 1;
}
