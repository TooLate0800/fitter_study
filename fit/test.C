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
  GE=fp*(1.-R2*Q2/6.+pd*Q2)/(1.+pd*Q2);//rational_1_1

  return GE;
}


double chi2(const double * fitpara){
  double R=fitpara[0];
  double pd=fitpara[1];
  double fp1=fitpara[2];
  //double fp2=fitpara[3];
  // double fp2=fitpara[2];
  
  double Q2v,GE_cal,GE0,GEd;
  double result=0.;

  for(int i=0;i<ndata;i++){
  //for(int i=1;i<ndata;i++){
    Q2v=Q2[i];GE0=GE[i];GEd=dGE[i];
    GE_cal=GE_fun(Q2v,R,pd,fp1);

    result=result+(GE_cal-GE0)*(GE_cal-GE0)/(GEd*GEd);
  }

  //for(int i=0;i<ndata_2;i++){
  //  Q2v=Q2_2[i];GE0=GE_2[i];GEd=dGE_2[i];
  //  GE_cal=GE_fun(Q2v,R,pd,fp2);

  //  result=result+(GE_cal-GE0)*(GE_cal-GE0)/(GEd*GEd);
  //}

  return result;
}


int xyminimizer(const char * minName = "Minuit", const char * algoName = ""){
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000); // for Minuit
  min->SetMaxIterations(10000);  // for GSL 
  //min->SetTolerance(1);
  min->SetTolerance(0.1);
  min->SetPrintLevel(1);   //output fit info
  // min->SetPrintLevel(0);   //not output fit info

  ROOT::Math::Functor f(&chi2,3);
  //ROOT::Math::Functor f(&chi2,4);
  min->SetFunction(f);


  double step[2],variable[2];
  step[0]=0.001;step[1]=0.0001;step[2]=0.0001;
  variable[0]=0.87;variable[1]=7.19364e-01;variable[2]=1.0;

  min->SetLimitedVariable(0, "Rp", variable[0], step[0],  0.0, 10.5);
  min->SetLimitedVariable(1, "pd", variable[1], step[1], -1.0, 1.0);
  min->SetLimitedVariable(2, "fp1", variable[2], step[2], 0.9, 1.1);
  //min->SetLimitedVariable(3, "fp2", variable[3], step[3], 0.9, 1.1);

  min->Minimize();

  const double * xsfit = min->X();
  Rfit=xsfit[0];
  pd_fit=xsfit[1];
  fp=xsfit[2];
  //fp_2=xsfit[3];

  chi2fit=f(xsfit);


  const double * xsfiterr=min->Errors();
  Rfiterr=xsfiterr[0];
  pd_fiterr=xsfiterr[1];
  fperr=xsfiterr[2];
  //fperr_2=xsfiterr[3];


  return 0;
}


int main(Int_t argc, char *argv[]){
 
 //for(int i = 0; i <  1; i++){
 for(int i = 0; i <  10000; i++){
  /*string fname1="no_syst_1GeV_full.txt";
  string fname2="no_syst_2GeV_full.txt";
  */
  //string fname1=Form("/var/phy/project/mepg/jz271/runPRad/model_generator/robust_table/700_table_%d.txt",i+1);
  //string fname1=Form("/var/phy/project/mepg/jz271/runPRad/model_generator/robust_table/1400_table_%d.txt",i+1);
  //string fname1=Form("/var/phy/project/mepg/jz271/runPRad/model_generator/robust_table/2100_table_%d.txt",i+1);
  string fname1=Form("/var/phy/project/mepg/jz271/fitter_study/generator/robust_table/1GeV_table_%d.txt",i+1);//input files
 
  //string fname1="/home/jz271/PRad/model_generator/0.7GeV_table.txt";
  //string fname1=Form("/home/jz271/PRad/model_generator/robust_table/700_table_%d.txt",i+1);

  //string fname1="/home/jz271/DRad/generator/PRad_II_10times_Kelly/1GeV_table_nosmear.txt";
  //string fname1=Form("/home/jz271/PRad/model_generator/robust_table_1per/0.7GeV_table_%d.txt",i+1);

  if (argc == 3){
      fname1 = argv[1];
      //fname2 = argv[2];
  }
  
  double modiFactor = 1.;
  if (argc == 4){
      fname1 = argv[1];
      //fname2 = argv[2];
      
      istringstream ss( argv[3] );
      double n;
      ss >> n;
      modiFactor += n/33;
      cout << modiFactor << endl;
  }



  readdata(fname1, modiFactor);
  //readdata(fname1, fname2, modiFactor);

          double para[3];
          para[0]=0.87962e-00; para[1]=7.03484e-01; para[2]=9.99736e-01;
          //para[0]=2.087962e-00; para[1]=7.03484e-01; para[2]=9.99736e-01; para[3]=9.98915e-01;
          double chi2_test=chi2(para);
          cout << "chi2 test= " << chi2_test << endl;


          cout << "minimize" << endl;
          xyminimizer();

          cout << "chi2fit= " << chi2fit << " , Rfit= " << Rfit << " , Rfiterr= " << Rfiterr << endl;
          ofstream outFile1, outFile2, outFile3, outFile4;
          outFile1.open("fit_result/PRad_z2_R11_2100_test20201230.txt", std::ios_base::app);
          //outFile1<<Rfit<<endl;
          outFile1.close();
//      }
  }
  

  return 0;
}
