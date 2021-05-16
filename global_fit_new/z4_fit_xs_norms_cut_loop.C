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

void z_trans(const int &ndata, const int &l)
{
  double num,den,Tc;
  Tc=2.;
  for(int i=0;i<ndata;i++){
    num=sqrt(Tc+Q2[l][i])-sqrt(Tc);
    den=sqrt(Tc+Q2[l][i])+sqrt(Tc);
    z_arr[l][i]=num/den;
  }
}


double sig_fun(const double &theta, const double &Q2, const double &z, const double &RE, const double &pb_E, const double &RM, const double &pb_M, const double &pd1_E, const double &pd1_M, const double &pd2_E, const double &pd2_M, const double &fp)

{
  double RE2=RE*RE;
  double RM2=RM*RM;
  double tau = Q2/4./pow(Mp*1e-3, 2);
  double t = theta / 180. * TMath::Pi();
  double epsilon = 1./(1. + 2.*(1.+tau)*pow(tan(t/2.), 2));
  double GE_dip = 1./pow(1. + Q2/0.71, 2); 
  double GM_dip = mu*1./pow(1. + Q2/0.71, 2); 

  // GE=1.-Q2*R2/6.;
  //GE=fp*(1.+R*Q2)/(1.+pd*Q2);
  double Q2_fm = Q2*GeVtofm;
  double GE=(1. - RE2*z*2.0*2./3.+pb_E*z*z + pd1_E*z*z*z + pd2_E*z*z*z*z);
  double GM=mu*(1. - RM2*z*2.0*2./3.+pb_M*z*z + pd1_M*z*z*z + pd2_M*z*z*z*z);
  //double GM=mu*(1. +RM*z+pb_M*z*z + pd1_M*z*z*z + pd2_M*z*z*z*z);
  double sig = fp*(epsilon*GE*GE+tau*GM*GM)/(epsilon*GE_dip*GE_dip+tau*GM_dip*GM_dip);

  return sig;
}


double chi2(const double * fitpara){
  double RE=fitpara[0];
  double pb_E=fitpara[1];
  double RM=fitpara[2];
  double pb_M=fitpara[3];
  double pd1_E=fitpara[4];
  double pd1_M=fitpara[5];
  double pd2_E=fitpara[6];
  double pd2_M=fitpara[7];
  //
  //double fp1=fitpara[4];
  //double fp2=fitpara[5];
  double fp[40];
  //double thetav, Q2v,cs_cal, cs0, csd;
  double thetav, zv, Q2v,cs_cal, cs0, csd;
  double result=0.;
  double float_para;
//cout<<number_of_norms<<" "<<number_of_sets<<endl;


  if (number_of_norms<number_of_sets){
      for (int j=1; j<number_of_norms+1;j++){
          fp[j] = fitpara[j+7];
      }
      //double param_list[34]={fp[3], fp[1]*fp[3], fp[1]*fp[4], fp[1]*fp[5], fp[2]*fp[4], fp[2]*fp[5], fp[9], fp[7]*fp[9], fp[6]*fp[9], fp[8]*fp[9], fp[13], fp[14], fp[11]*fp[13], fp[10]*fp[13], fp[10]*fp[14], fp[10]*fp[15], fp[12]*fp[15], fp[18], fp[19], fp[16]*fp[18], fp[16]*fp[19], fp[16]*fp[20], fp[17]*fp[20], fp[25], fp[21]*fp[25], fp[21]*fp[26], fp[23]*fp[26], fp[22]*fp[26], fp[24]*fp[26], fp[29]*fp[30], fp[29], fp[27]*fp[29], fp[27]*fp[31], fp[28]*fp[31]};

      for (int l = 0; l<number_of_sets; l++){
          if (norm2[l] == -1){float_para = fp[norm1[l]];}
          else {float_para = fp[norm1[l]]*fp[norm2[l]];}
          for(int i=0;i<ndata[l];i++){
            //zv = z_arr[l][i]; thetav=theta[l][i];Q2v=Q2[l][i];cs0=cs[l][i];csd=dcs[l][i];
            cs_cal=sig_fun(thetav, Q2v, zv, RE, pb_E, RM, pb_M, pd1_E, pd1_M, pd2_E, pd2_M, float_para);
            //cs_cal=sig_fun(thetav,Q2v,R,pb_E, RM, pb_M, param_list[l]);
            result=result+(cs_cal-cs0)*(cs_cal-cs0)/(csd*csd);
          }
      }
  }else{
      for (int j=1; j<number_of_sets+1;j++){
          fp[j] = fitpara[j+7];
      }
      for (int l = 0; l<number_of_sets; l++){
          for(int i=0;i<ndata[l];i++){
            //zv = z_arr[l][i]; thetav=theta[l][i];Q2v=Q2[l][i];cs0=cs[l][i];csd=dcs[l][i];
            cs_cal=sig_fun(thetav, Q2v, zv, RE, pb_E, RM, pb_M, pd1_E, pd1_M, pd2_E, pd2_M, fp[l+1]);
            //cs_cal=sig_fun(thetav,Q2v,R,pb_E, RM, pb_M, param_list[l]);
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
  min->SetTolerance(0.1);
  min->SetPrintLevel(1);   //output fit info
  // min->SetPrintLevel(0);   //not output fit info

  int number_of_fp;
  if (number_of_norms<number_of_sets){number_of_fp=number_of_norms;
  }else{number_of_fp=number_of_sets;}

  ROOT::Math::Functor f(&chi2,8+number_of_fp);
  min->SetFunction(f);


  double step[6],variable[6];
  step[0]=0.001;step[1]=0.01;step[2]=0.01;step[3]=0.01;step[4]=0.01;step[5] = 0.001;
  variable[0]=0.87;variable[1]=-7.40479e-01;variable[2]=0.80;variable[3]=-7.40479e-01;variable[4]=1.0;variable[5] = 1.0;

  min->SetLimitedVariable(0, "RE", variable[0], step[0], 0.7, 1.0);
  min->SetLimitedVariable(1, "pb_E", variable[1], step[1], -500., 500.);
  min->SetLimitedVariable(2, "RM", variable[2], step[2], 0.7, 1.0);
  min->SetLimitedVariable(3, "pb_M", variable[3], step[3], -500., 500.);
  min->SetLimitedVariable(4, "pd1_E", -1.85497e+00, 0.1,  -500., 500.);
  min->SetLimitedVariable(5, "pd2_E", 4.65666e+00, 0.1,  -500., 500.);
  min->SetLimitedVariable(6, "pd1_M", -1.85497e+00, 0.1,  -500., 500.);
  min->SetLimitedVariable(7, "pd2_M", 4.65666e+00, 0.1,  -500., 500.);

  //min->SetLimitedVariable(4, "fp_1", variable[4], step[4], 0.9, 1.1);
  //min->SetLimitedVariable(5, "fp_2", variable[5], step[5], 0.9, 1.1);


  for (int j=1; j<number_of_fp+1; j++){
      min->SetLimitedVariable(j+7, Form("fp_%d",j), 1.0, 0.001, 0.7, 1.3);
  }
  min->Minimize();

  const double * xsfit = min->X();
  Rfit=xsfit[0];
  pbE_fit=xsfit[1];
  paM_fit=xsfit[2];
  pbM_fit=xsfit[3];
  //for (int j=1; j<32; j++){
  //    fp_[j]=xsfit[j+3];
  //}
  //fp_1=xsfit[4];
  //fp_2=xsfit[5];

  chi2fit=f(xsfit);


  const double * xsfiterr=min->Errors();
  Rfiterr=xsfiterr[0];
  pbE_fiterr=xsfiterr[1];
  paM_fiterr=xsfiterr[2];
  pbM_fiterr=xsfiterr[3];
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

               string tablename = "../data/cs_data/Cut_05GeV/norm_table.dat";
               readtable(tablename); 


               for (int l=0; l<number_of_sets; l++){
                   fname[l] = Form("../generator_31norms/xs_pseudo_GMfromRatio_CutQ2_05GeV_1e3/Model%d_norm%d_table_%d.txt",j+1, l+1, i+1);
                   //fname[l] = Form("../data/cs_data/Cut_1GeV/%d.txt",l);
                   readdata(fname[l],l);
                   z_trans(ndata[l],l);
               }


               cout << "minimize" << endl;
               xyminimizer();

               cout << "chi2fit= " << chi2fit << " , Rfit= " << Rfit << " , Rfiterr= " << Rfiterr << endl;

               
               ofstream outFile1,outFile2;
               outFile1.open(Form("result/z4_xsfit_RE_GMfromRatio_Cut1GeV_1e3_Model%d.txt",j+1), std::ios_base::app);
               //outFile1<<Rfit<<endl;
               outFile1.close();
               outFile2.open(Form("result/z4_xsfit_all_GMfromRatio_Cut1GeV_1e3_Model%d.txt",j+1), std::ios_base::app);
               //outFile2<<chi2fit<<" "<<Rfit<<" "<<Rfiterr<<endl;
               outFile2.close();

               //cout<<Rfit<<" "<<Rfiterr<<" "<<fp<<" "<<fperr<<" "<<fp_2<<" "<<fperr_2<<" "<<chi2fit/(ndata+ndata_2-4)<<endl;
               //cout<<number_of_norms<<endl;
      }
  }
  return 0;
}
