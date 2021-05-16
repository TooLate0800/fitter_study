#include<iostream>
#include<fstream>
#include <iomanip>
#include <string>
#include <sstream>
using namespace std;

#include "stdlib.h"
#include "TROOT.h"
#include "TApplication.h"
#include "Rtypes.h"
#include "math.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TObject.h"
#include "TKey.h"
#include "TTree.h"

#include "TChain.h"
#include "TString.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TEventList.h"

//#include "data_one.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"

int roofit_loop(){
    gStyle->SetOptFit(111111);
    for(int i = 0; i <  100; i++){
        //define fitting function
        TF1 *f1 = new TF1 ("f1","[1]/(1+[0]*x/12)/(1+[0]*x/12)",0.,25);//dipole
        //TF1 *f1 = new TF1 ("f1","[2]*(1+[0]*x)/(1+[1]*x)",0.,25);//rational_1_1
        //auto c=new TCanvas();
        //TF1 *f1 = new TF1 ("f1","(1+[0]*x)/(1+[1]*x+[2]*x)",0.,1.5);
        //TF1 *f1 = new TF1 ("f1","[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x+[7]*x*x*x*x*x*x*x+[8]*x*x*x*x*x*x*x*x+[9]*x*x*x*x*x*x*x*x*x+[10]*x*x*x*x*x*x*x*x*x*x)",0.,26);
        //TF1 *f1 = new TF1 ("f1","[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x+[7]*x*x*x*x*x*x*x+[8]*x*x*x*x*x*x*x*x+[9]*x*x*x*x*x*x*x*x*x+[10]*x*x*x*x*x*x*x*x*x*x+[11]*x*x*x*x*x*x*x*x*x*x*x)",0.,26);
        //TF1 *f1 = new TF1 ("f1","[2]*(1+[0]*x)/(1+[1]*x)",0.,26);
        //TF1 *f1 = new TF1 ("f1","[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x+[7]*x*x*x*x*x*x*x+[8]*x*x*x*x*x*x*x*x)",0.,26);
        //TF1 *f1 = new TF1 ("f1","[3]*(1+[0]*x+[1]*x*x)/(1+[2]*x)",0.,45);//rational_2_1
        //TF1 *f1 = new TF1 ("f1","[4]*(1+[0]*x)/(1+[1]*x+[2]*x*x+[3]*x*x*x)",0.,45);//rational_1_3
        //TF1 *f1 = new TF1 ("f1","[2]*(1+[0]*x)/(1+[1]*x)",0.,45);//modified_A_B
        //TF1 *f1 = new TF1 ("f1","real([4]*pow((1+[0]*x),[2])/pow((1+[1]*x),[3]))",0.,45);//modified_A_B
        //TF1 *f1 = new TF1 ("f1","real([3]*pow((1+[0]*x),[2])/(1+[1]*x))",0.,45);//modified_A
        //TF1 *f1 = new TF1 ("f1","(1+[0]*x)/(1+[1]*x+0.0227*x*x+0.005*x*x*x)",0.,1.5);
        //f1-> SetParameter(0,-0.047);
        //f1-> SetParameter(1,0.6626);
        //f1-> SetParameter(2,1.0);
        //f1-> SetParameter(3,0.80);
        //f1-> SetParameter(4,1.00);
        //f1-> SetParameter(0,-0.065);
        //f1-> SetParameter(1,0.6626);
        //f1-> SetParameter(2,0.0227);
        //f1-> SetParameter(3,0.005);
        TGraphErrors graph(Form("../generator/robust_table/1GeV_table_%d.txt",i+1),"%lg %lg %lg");//input file
        graph.Fit("f1","B");
        //cout<<f1->GetChisquare()<<endl;
        //cout<<f1->GetParError(1)<<endl;
        double p0 = f1->GetParameter(0);
        //double p1 = f1->GetParameter(1);
        //double R = 6*sqrt(p1*p1-p0*p0);
        //cout<<R<<endl;
        double Rfit = TMath::Sqrt(p0);//dipole 
        ofstream outFile1;
        outFile1.open("roofit_result/dipole_dipole_1e4.txt", std::ios_base::app);//output file
        outFile1<<Rfit<<endl;
        outFile1.close();


    }

    return 0;
}
