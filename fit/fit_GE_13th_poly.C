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

#include "data1.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"


double GE_fun(const double &Q2, const double &R, const double &pd, const double &pd1, const double &pd2 , const double &pd3, const double &pd4, const double &pd5, const double &pd6, const double &pd7, const double &pd8, const double &pd9, const double &pd10, const double &pd11, const double &fp)
{
  double GE;
  double R2=R*R;
  GE=fp*(1.-R2*Q2/6.+pd*Q2*Q2+pd1*pow(Q2,3)+pd2*pow(Q2,4)+pd3*pow(Q2,5)+pd4*pow(Q2,6)+pd5*pow(Q2,7)+pd6*pow(Q2,8)+pd7*pow(Q2,9)+pd8*pow(Q2,10)+pd9*pow(Q2,11)+pd10*pow(Q2,12)+pd11*pow(Q2,13));
  return GE;
}


double chi2(const double * fitpara){
  double R=fitpara[0];//radius, one of the fitting parameters
  double pd=fitpara[1];//the other parameter
  double fp1=fitpara[2];//the floating parameter
  double pd1=fitpara[3];//the other parameter
  double pd2=fitpara[4];//the other parameter
  double pd3=fitpara[5];//the other parameter
  double pd4=fitpara[6];//the other parameter
  double pd5=fitpara[7];//the other parameter
  double pd6=fitpara[8];//the other parameter
  double pd7=fitpara[9];//the other parameter
  double pd8=fitpara[10];//the other parameter
  double pd9=fitpara[11];//the other parameter
  double pd10=fitpara[12];//the other parameter
  double pd11=fitpara[13];//the other parameter
  
  double Q2v,GE_cal,GE0,GEd;
  double result=0.;

  for(int i=0;i<ndata;i++){
  //for(int i=1;i<ndata;i++){
    Q2v=Q2[i];GE0=GE[i];GEd=dGE[i];
    GE_cal=GE_fun(Q2v,R,pd,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,fp1);

    result=result+(GE_cal-GE0)*(GE_cal-GE0)/(GEd*GEd);
  }

  return result;
}


int xyminimizer(const char * minName = "Minuit", const char * algoName = ""){
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(1000000); // for Minuit
  min->SetMaxIterations(10000000);  // for GSL 
  min->SetTolerance(1);
  //min->SetTolerance(0.1);
  min->SetPrintLevel(1);   //output fit info
  // min->SetPrintLevel(0);   //not output fit info

  ROOT::Math::Functor f(&chi2,3);
  min->SetFunction(f);


  double step[14],variable[14];
  step[0]=0.001;step[1]=0.0001;step[2]=0.0001;step[3]=0.0001;step[4]=0.0001;step[5]=0.0001;step[6]=0.0001;step[7]=0.0001;step[8]=0.0001;step[9]=0.0001;step[10]=0.0001;step[11]=0.0001;step[12]=0.0001;step[13]=0.0001;
  variable[0]=0.87;variable[1]=7.19364e-01;variable[2]=1.0;variable[3]=7.19364e-01;variable[4]=7.19364e-01;variable[5]=7.19364e-01;variable[6]=7.19364e-01;variable[7]=7.19364e-01;variable[8]=7.19364e-01;variable[9]=7.19364e-01;variable[10]=7.19364e-01;variable[11]=7.19364e-01;variable[12]=7.19364e-01;variable[13]=7.19364e-01;

  //define the fitting range of each parameter
  min->SetLimitedVariable(0, "Rp", variable[0], step[0],  0.0, 2.5);//constrain of Rp
  min->SetLimitedVariable(1, "pd", variable[1], step[1], -500.0, 500.0);
  min->SetLimitedVariable(2, "fp1", variable[2], step[2], 0.9, 1.1);//constrain of floating parameter
  min->SetLimitedVariable(3, "pd1", variable[3], step[3], -500.0, 500.0);//constrain of fitting parameters, can be varied
  min->SetLimitedVariable(4, "pd2", variable[4], step[4], -500.0, 500.0);//constrain of fitting parameters, can be varied
  min->SetLimitedVariable(5, "pd3", variable[5], step[5], -500.0, 500.0);//constrain of fitting parameters, can be varied
  min->SetLimitedVariable(6, "pd4", variable[6], step[6], -500.0, 500.0);//constrain of fitting parameters, can be varied
  min->SetLimitedVariable(7, "pd5", variable[7], step[7], -500.0, 500.0);//constrain of fitting parameters, can be varied
  min->SetLimitedVariable(8, "pd6", variable[8], step[8], -500.0, 500.0);//constrain of fitting parameters, can be varied
  min->SetLimitedVariable(9, "pd7", variable[9], step[9], -500.0, 500.0);//constrain of fitting parameters, can be varied
  min->SetLimitedVariable(10, "pd8", variable[10], step[10], -500.0, 500.0);//constrain of fitting parameters, can be varied
  min->SetLimitedVariable(11, "pd9", variable[11], step[11], -500.0, 500.0);//constrain of fitting parameters, can be varied
  min->SetLimitedVariable(12, "pd10", variable[12], step[12], -500.0, 500.0);//constrain of fitting parameters, can be varied
  min->SetLimitedVariable(13, "pd11", variable[13], step[13], -500.0, 500.0);//constrain of fitting parameters, can be varied

  min->Minimize();

  const double * xsfit = min->X();
  Rfit=xsfit[0];
  pd_fit=xsfit[1];
  fp=xsfit[2];

  chi2fit=f(xsfit);


  const double * xsfiterr=min->Errors();
  Rfiterr=xsfiterr[0];
  pd_fiterr=xsfiterr[1];
  fperr=xsfiterr[2];


  return 0;
}


int main(Int_t argc, char *argv[]){
 
 for(int i = 0; i <  1; i++){//repetition times
  //string fname1=Form("/var/phy/project/mepg/jz271/runPRad/model_generator/robust_table/700_table_%d.txt",i+1);
  //string fname1=Form("/var/phy/project/mepg/jz271/runPRad/model_generator/robust_table/1400_table_%d.txt",i+1);
  //string fname1="../data/Carl-norm.dat";//input files
  string fname1="../generator/test_input/dipole_1000.txt";//input files
  //string fname1=Form("/var/phy/project/mepg/jz271/fitter_study/generator/robust_table/1GeV_table_%d.txt",i+1);//input files
  

  if (argc == 2){
      fname1 = argv[1];
  }
  
  double modiFactor = 1.;


  readdata(fname1, modiFactor);
  //readdata(fname1, fname2, modiFactor);

          double para[3];
          para[0]=0.87962e-00; para[1]=7.03484e-01; para[2]=9.99736e-01;para[3]=7.03484e-01;para[4]=7.03484e-01;para[5]=7.03484e-01;para[6]=7.03484e-01;para[7]=7.03484e-01;para[8]=7.03484e-01;para[9]=7.03484e-01;para[10]=7.03484e-01;para[11]=7.03484e-01;para[12]=7.03484e-01;para[13]=7.03484e-01;//define the initial value of the parameters
          double chi2_test=chi2(para);
          cout << "chi2 test= " << chi2_test << endl;


          cout << "minimize" << endl;
          xyminimizer();

          cout << "chi2fit= " << chi2fit << " , Rfit= " << Rfit << " , Rfiterr= " << Rfiterr << endl;
          ofstream outFile1, outFile2, outFile3, outFile4;
          outFile1.open("fit_result/test.txt", std::ios_base::app);//output file
          //outFile1<<Rfit<<endl;
          outFile1.close();
//      }
  }
  

  return 0;
}
