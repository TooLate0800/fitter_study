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


double GE_fun(const double &Q2, const double &R, const double &pd, const double &fp)
{
  double GE;
  double R2=R*R;
  //GE=fp*(1.+R*Q2)/(1.+pd*Q2);//rational_1_1
  GE=fp*(1.-R2*Q2/6.+pd*Q2)/(1.+pd*Q2);//definition of rational_1_1
  //GE=fp*(1.-R2*Q2/6.+pd*Q2*Q2);//poly_2nd
  //GE=fp*(1/(1+(R2*Q2/6)/(1+pd*Q2)));//CF2
  



  return GE;
}


double chi2(const double * fitpara){
  double R=fitpara[0];//radius, one of the fitting parameters
  double pd=fitpara[1];//the other parameter
  double fp1=fitpara[2];//the floating parameter
  
  double Q2v,GE_cal,GE0,GEd;
  double result=0.;

  for(int i=0;i<ndata;i++){
  //for(int i=1;i<ndata;i++){
    Q2v=Q2[i];GE0=GE[i];GEd=dGE[i];
    GE_cal=GE_fun(Q2v,R,pd,fp1);

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


  double step[3],variable[3];
  step[0]=0.001;step[1]=0.01;step[2]=0.001;
  variable[0]=0.85;variable[1]=7.19364e-03;variable[2]=1.0;

  //define the fitting range of each parameter
  min->SetLimitedVariable(0, "Rp", variable[0], step[0],  0.7, 1.0);
  min->SetLimitedVariable(1, "pd", variable[1], step[1], -500.0, 500.0);
  min->SetLimitedVariable(2, "fp1", variable[2], step[2], 0.9, 1.1);

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
  string fname1="../data/Carl-norm_PRad_range.dat";//input files
  //string fname1=Form("/var/phy/project/mepg/jz271/fitter_study/bootstrap/Mainz_bootstrap_PRad_60bins/1GeV_table_%d.txt",i+1);//input files
  //string fname1=Form("/var/phy/project/mepg/jz271/fitter_study/bootstrap/Mainz_bootstrap_weighted_6points/1GeV_table_%d.txt",i+1);//input files
  //string fname1=Form("/var/phy/project/mepg/jz271/fitter_study/generator/Mainz_PRad_range/Model_1/1GeV_table_%d.txt",i+1);//input files
  

  if (argc == 2){
      fname1 = argv[1];
  }
  
  double modiFactor = 1.;
  if (argc == 2){
      fname1 = argv[1];
      
      istringstream ss( argv[3] );
      double n;
      ss >> n;
      modiFactor += n/33;
      cout << modiFactor << endl;
  }



  readdata(fname1, modiFactor);
  //readdata(fname1, fname2, modiFactor);

          double para[3];
          para[0]=0.87962e-00; para[1]=7.03484e-01; para[2]=9.99736e-01;//define the initial value of the parameters
          double chi2_test=chi2(para);
          cout << "chi2 test= " << chi2_test << endl;


          cout << "minimize" << endl;
          xyminimizer();

          cout << "chi2fit= " << chi2fit << " , Rfit= " << Rfit << " , Rfiterr= " << Rfiterr << endl;
          ofstream outFile1, outFile2, outFile3, outFile4;
          outFile1.open("fit_result/Rfit_R11_bootstrap_MainzinPRad_60bins.txt", std::ios_base::app);//output file
          //outFile1<<Rfit<<endl;
          outFile1.close();
          outFile2.open("fit_result/Rfiterr_R11_bootstrap_MainzinPRad_60bins.txt", std::ios_base::app);//output file
          //outFile2<<Rfiterr<<endl;
          outFile2.close();
          outFile3.open("fit_result/chi2_R11_bootstrap_MainzinPRad_60bins.txt", std::ios_base::app);//output file
          //outFile3<<chi2fit/(120-3)<<endl;
          outFile3.close();
          //cout<<chi2fit/(1422-3)<<endl;
//      }
  }
  

  return 0;
}
