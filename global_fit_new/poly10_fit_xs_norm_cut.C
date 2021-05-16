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

//common physical constants
const double mkb = 389.379404e-3;
const double alpha = 1. / 137.036;
const double Mp = 938.272046;
const double Me = 0.510998928;
const double mu = 2.79284736; // Magnetic moment of the proton
const double GeVtofm = 25.68189504;


double sig_fun(const double &theta, const double &Q2, const double &R_E, const double &pb_E1, const double &pb_E2, const double &pb_E3, const double &pb_E4, const double &pb_E5, const double &pb_E6, const double &pb_E7, const double &pb_E8, const double &pb_E9, const double &R_M, const double &pb_M1, const double &pb_M2, const double &pb_M3, const double &pb_M4, const double &pb_M5, const double &pb_M6, const double &pb_M7, const double &pb_M8, const double &pb_M9, const double &fp)
{
  double R2_E=R_E*R_E;
  double R2_M=R_M*R_M;
  double tau = Q2/4./pow(Mp*1e-3, 2);
  double t = theta / 180. * TMath::Pi();
  double epsilon = 1./(1. + 2.*(1.+tau)*pow(tan(t/2.), 2));
  double GE_dip = 1./pow(1. + Q2/0.71, 2); 
  double GM_dip = mu*1./pow(1. + Q2/0.71, 2); 

  // GE=1.-Q2*R2/6.;
  //GE=fp*(1.+R*Q2)/(1.+pb*Q2);
  double Q2_fm = Q2*GeVtofm;
  double GE=(1.-R2_E*Q2_fm/6.+pb_E1*Q2*Q2+pb_E2*pow(Q2,3)+pb_E3*pow(Q2,4)+pb_E4*pow(Q2,5)+pb_E5*pow(Q2,6)+pb_E6*pow(Q2,7)+pb_E7*pow(Q2,8)+pb_E8*pow(Q2,9)+pb_E9*pow(Q2,10));
  double GM=mu*(1.-R2_M*Q2_fm/6.+pb_M1*Q2*Q2+pb_M2*pow(Q2,3)+pb_M3*pow(Q2,4)+pb_M4*pow(Q2,5)+pb_M5*pow(Q2,6)+pb_M6*pow(Q2,7)+pb_M7*pow(Q2,8)+pb_M8*pow(Q2,9)+pb_M9*pow(Q2,10));
  double sig = fp*(epsilon*GE*GE+tau*GM*GM)/(epsilon*GE_dip*GE_dip+tau*GM_dip*GM_dip);

  return sig;
}


double chi2(const double * fitpara){
  double R_E=fitpara[0];
  double pb_E[9], pb_M[9];
  for (int j=1; j<10;j++){
      pb_E[j] = fitpara[j];
  }
  double R_M=fitpara[10];
  for (int j=1; j<10;j++){
      pb_M[j] = fitpara[j+10];
  }
  double fp[32];
  //for (int j=1; j<32;j++){
  //    fp[j] = fitpara[j+19];
  //}


  double thetav, Q2v, cs_cal, cs0, csd;
  double result=0.;
  double float_para;

//  double param_list[34]={fp[3], fp[1]*fp[3], fp[1]*fp[4], fp[1]*fp[5], fp[2]*fp[4], fp[2]*fp[5], fp[9], fp[7]*fp[9], fp[6]*fp[9], fp[8]*fp[9], fp[13], fp[14], fp[11]*fp[13], fp[10]*fp[13], fp[10]*fp[14], fp[10]*fp[15], fp[12]*fp[15], fp[18], fp[19], fp[16]*fp[18], fp[16]*fp[19], fp[16]*fp[20], fp[17]*fp[20], fp[25], fp[21]*fp[25], fp[21]*fp[26], fp[23]*fp[26], fp[22]*fp[26], fp[24]*fp[26], fp[29]*fp[30], fp[29], fp[27]*fp[29], fp[27]*fp[31], fp[28]*fp[31]};
// 
//  for (int l = 0; l<34; l++){
//      for(int i=0;i<ndata[l];i++){
//        thetav=theta[l][i];Q2v=Q2[l][i];cs0=cs[l][i];csd=dcs[l][i];
//        cs_cal=sig_fun(thetav,Q2v,R_E,pb_E[1], pb_E[2],pb_E[3],pb_E[4], pb_E[5], pb_E[6], pb_E[7], pb_E[8], pb_E[9], R_M, pb_M[1], pb_M[2], pb_M[3], pb_M[4], pb_M[5], pb_M[6], pb_M[7], pb_M[8], pb_M[9], param_list[l]);
//        result=result+(cs_cal-cs0)*(cs_cal-cs0)/(csd*csd);
//      }
//
//  //for(int i=0;i<ndata[1];i++){
//  //  thetav=theta[1][i];Q2v=Q2[1][i];cs0=cs[1][i];csd=dcs[1][i];
//  //  cs_cal=sig_fun(thetav,Q2v,R,pb_E, pa_M, pb_M,fp[1]*fp[3]);
//  //  result=result+(cs_cal-cs0)*(cs_cal-cs0)/(csd*csd);
//  //}
//  }
  if (number_of_norms<number_of_sets){
      for (int j=1; j<number_of_norms+1;j++){
          fp[j] = fitpara[j+19];
      }

      for (int l = 0; l<number_of_sets; l++){
          if (norm2[l] == -1){float_para = fp[norm1[l]];}
          else {float_para = fp[norm1[l]]*fp[norm2[l]];}
          for(int i=0;i<ndata[l];i++){
            thetav=theta[l][i];Q2v=Q2[l][i];cs0=cs[l][i];csd=dcs[l][i];
            cs_cal=sig_fun(thetav,Q2v,R_E,pb_E[1], pb_E[2],pb_E[3],pb_E[4], pb_E[5], pb_E[6], pb_E[7], pb_E[8], pb_E[9], R_M, pb_M[1], pb_M[2], pb_M[3], pb_M[4], pb_M[5], pb_M[6], pb_M[7], pb_M[8], pb_M[9], float_para);
            result=result+(cs_cal-cs0)*(cs_cal-cs0)/(csd*csd);
          }
      }
  }else{
      for (int j=1; j<number_of_sets+1;j++){
          fp[j] = fitpara[j+3];
      }
      for (int l = 0; l<number_of_sets; l++){
          for(int i=0;i<ndata[l];i++){
            thetav=theta[l][i];Q2v=Q2[l][i];cs0=cs[l][i];csd=dcs[l][i];
            cs_cal=sig_fun(thetav,Q2v,R_E,pb_E[1], pb_E[2],pb_E[3],pb_E[4], pb_E[5], pb_E[6], pb_E[7], pb_E[8], pb_E[9], R_M, pb_M[1], pb_M[2], pb_M[3], pb_M[4], pb_M[5], pb_M[6], pb_M[7], pb_M[8], pb_M[9], fp[l+1]);
            result=result+(cs_cal-cs0)*(cs_cal-cs0)/(csd*csd);
          }
      }

  }


  return result;
}


int xyminimizer(const char * minName = "Minuit", const char * algoName = ""){
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000); // for Minuit
  min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.01);
  min->SetPrintLevel(1);   //output fit info
  // min->SetPrintLevel(0);   //not output fit info

  int number_of_fp;
  if (number_of_norms<number_of_sets){number_of_fp=number_of_norms;
  }else{number_of_fp=number_of_sets;}

  ROOT::Math::Functor f(&chi2,20+number_of_fp);

  min->SetFunction(f);


  double step[6],variable[6];
  step[0]=0.001;step[1]=0.0001;step[2]=0.0001;step[3]=0.0001;step[4]=0.0001;step[5] = 0.001;
  variable[0]=0.87;variable[1]=9.19364e-02;variable[2]=0.70;variable[3]=9.19364e-02;variable[4]=1.0;variable[5] = 1.0;

  min->SetLimitedVariable(0, "R_E", variable[0], step[0], 0.6, 1.4);
  min->SetLimitedVariable(1, "pb_E1", 1.45487683e+01,  0.01, -1e5, 1e5);
  min->SetLimitedVariable(2, "pb_E2", -8.87959239e+01, 0.01, -1e5, 1e5);
  min->SetLimitedVariable(3, "pb_E3", 4.61097705e+02,  0.01, -1e5, 1e5);
  min->SetLimitedVariable(4, "pb_E4", -1.67562381e+03, 0.01, -1e5, 1e5);
  min->SetLimitedVariable(5, "pb_E5", 4.07646487e+03,  0.01, -1e5, 1e5);
  min->SetLimitedVariable(6, "pb_E6", -6.45411460e+03, 0.01, -1e5, 1e5);
  min->SetLimitedVariable(7, "pb_E7", 6.34035079e+03,  0.01, -1e5, 1e5);
  min->SetLimitedVariable(8, "pb_E8", -3.49373923e+03, 0.01, -1e5, 1e5);
  min->SetLimitedVariable(9, "pb_E9", 8.22601568e+02,  0.01, -1e5, 1e5);
  //for (int j=1; j<10; j++){
  //    min->SetLimitedVariable(j, Form("pb_E%d",j), 9.19364e-02, 0.1, -1e4, 1e4);
  //}
  min->SetLimitedVariable(10, "R_M", variable[2], step[2], 0.6, 1.4);
  min->SetLimitedVariable(11, "pb_M1", 1.0222000, 0.01, -1e5, 1e5);
  min->SetLimitedVariable(12, "pb_M2", 23.494500, 0.01, -1e5, 1e5);
  min->SetLimitedVariable(13, "pb_M3", -93.03720, 0.01, -1e5, 1e5);
  min->SetLimitedVariable(14, "pb_M4", 140.79840, 0.01, -1e5, 1e5);
  min->SetLimitedVariable(15, "pb_M5", -0.365600, 0.01, -1e5, 1e5);
  min->SetLimitedVariable(16, "pb_M6", -305.6759, 0.01, -1e5, 1e5);
  min->SetLimitedVariable(17, "pb_M7", 444.62510, 0.01, -1e5, 1e5);
  min->SetLimitedVariable(18, "pb_M8", -273.6688, 0.01, -1e5, 1e5);
  min->SetLimitedVariable(19, "pb_M9", 64.581100, 0.01, -1e5, 1e5);
  //for (int j=1; j<10; j++){
  //    min->SetLimitedVariable(j+10, Form("pb_M%d",j), 9.19364e-02, 0.1, -1e4, 1e4);
  //}

  for (int j=1; j<number_of_fp+1; j++){
      min->SetLimitedVariable(j+19, Form("fp_%d",j), 1.0, 0.001, 0.8, 1.2);
  }
  min->Minimize();

  const double * xsfit = min->X();
  REfit=xsfit[0];
  //pbE_fit=xsfit[1];
  RMfit=xsfit[10];
  //pbM_fit=xsfit[3];
  //for (int j=1; j<32; j++){
  //    fp_[j]=xsfit[j+3];
  //}
  //fp_1=xsfit[4];
  //fp_2=xsfit[5];

  chi2fit=f(xsfit);


  const double * xsfiterr=min->Errors();
  REfiterr=xsfiterr[0];
  //pbE_fiterr=xsfiterr[1];
  RMfiterr=xsfiterr[10];
  //pbM_fiterr=xsfiterr[3];
  //for (int j=1; j<32; j++){
  //    fperr[j]=xsfiterr[j+3];
  //}
  //fperr_1=xsfiterr[4];
  //fperr_2=xsfiterr[5];


  return 0;
}


int main(Int_t argc, char *argv[]){

  string fname[34];
  for (int j = 0; j <  1; j++){//models
       for(int i = 0; i <  1; i++){//table

           string tablename = "../data/cs_data/Cut_1GeV/norm_table.dat";
           readtable(tablename);
         
         
           for (int l=0; l<number_of_sets; l++){
               fname[l] = Form("../generator_31norms/xs_pseudo_GMfromRatio/Model8_norm%d_table_1.txt",l);
               //fname[l] = Form("../data/cs_data/Cut_009GeV/%d.txt",l);
               readdata(fname[l],l);
           }


           //for (int l=0; l<34; l++){
           //    //fname[l] = Form("../generator_31norms/xs_pseudo_GMfromRatio/Model%d_norm%d_table_%d.txt",j+1, l+1, i+1); 
           //    fname[l] = Form("../data/cs_data/Cut_08GeV/%d.txt",l); 
           //    readdata(fname[l],l);
           //    //readdata(fname[i],i);
           //}



           cout << "minimize" << endl;
           xyminimizer();

           cout << "chi2fit= " << chi2fit << " , Rfit= " << REfit << " , REfiterr= " << REfiterr << endl;

           
           ofstream outFile1,outFile2;
           outFile1.open(Form("result/poly10_xsfit_RE_GMfromRatio_check_Model%d.txt",j+1), std::ios_base::app);
           //outFile1.open(Form("result/R11_xsfit_RE_GMRandom_Model%d.txt",j+1), std::ios_base::app);
           //outFile1<<REfit<<endl;
           outFile1.close();
           outFile2.open(Form("result/poly_xsfit_all_GMfromRatio_check_Model%d.txt",j+1), std::ios_base::app);
           //outFile2.open(Form("result/R11_xsfit_all_GMRandom_Model%d.txt",j+1), std::ios_base::app);
           //outFile2<<chi2fit<<" "<<REfit<<" "<<REfiterr<<" "<<RMfit<<" "<<RMfiterr<<endl;
           outFile2.close();
      }
  }
  //cout<<Rfit<<" "<<Rfiterr<<" "<<fp<<" "<<fperr<<" "<<fp_2<<" "<<fperr_2<<" "<<chi2fit/(ndata+ndata_2-4)<<endl;
  return 0;
}
