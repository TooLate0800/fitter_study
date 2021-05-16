#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>

using namespace std;

//data
Int_t ndata[2000],ndata_2;
Double_t theta[40][2000], Q2[40][2000],cs[40][2000],dcs[40][2000];
Double_t GE[40][2000],dGE[40][2000];
Double_t Q2_2[2000],GE_2[2000],dGE_2[2000];
Double_t z_arr[40][2000];
//Double_t pd_fit;
//Double_t pd_fiterr;
Double_t chi2fit, Rfit, REfit, RMfit, pbE_fit, paM_fit, pbM_fit;
Double_t Rfiterr, REfiterr, RMfiterr, pbE_fiterr, paM_fiterr, pbM_fiterr;
Double_t fp_[40];
Double_t fperr[40];
Double_t fp_1, fp_2;
Double_t fperr_1, fperr_2;
Int_t fit_type;
Int_t npower;
Int_t dist_type;
Int_t number_of_sets;
Int_t number_of_norms = 0;
Int_t number_of_points[40];
Int_t norm1[40];
Int_t norm2[40];
Int_t norm[80];

int read_GEdata(string fname, int n)
{
  TString filename=fname;

  ifstream infile(filename);
  if(!infile.is_open()){
    cout << "File not exist." << endl;
    return 0;
  }
  infile >> ndata[n];

  for(int i=0;i<ndata[n];i++){
    infile >> Q2[n][i] >> GE[n][i] >> dGE[n][i];
  }
  infile.close();

  return 0;
}

int readdata(string fname, int n)
{
  TString filename=fname;

  ifstream infile(filename);
  if(!infile.is_open()){
    cout << "File not exist." << endl;
    return 0;
  }
  infile >> ndata[n];

  for(int i=0;i<ndata[n];i++){
    infile >> theta[n][i] >> Q2[n][i] >> cs[n][i] >> dcs[n][i];
  }
  infile.close();

  return 0;
}

int readtable(string fname)
{
  TString filename=fname;

  ifstream infile(filename);
  if(!infile.is_open()){
    cout << "File not exist." << endl;
    return 0;
  }
  infile >> number_of_sets;
 
  double num;
  for(int i=0;i<number_of_sets;i++){
    infile >> num >> number_of_points[i] >> norm1[i] >> norm2[i];
  }
  infile.close();
number_of_norms = number_of_sets*2;
for(int i=0;i<number_of_sets;i++){
    for(int j=0;j<i;j++){
        if (norm1[i] == norm1[j] || norm1[i] == norm2[j]){number_of_norms--;break;}//neglect repeated numbers
    }
    for(int j=0;j<i;j++){
        if (norm2[i] == norm2[j] || norm2[i] == norm1[j]){number_of_norms--;break;}
    }
    norm[i] = norm1[i];
    norm[i+number_of_sets] = norm2[i];//combine into 1 array
}
number_of_norms--;//neglect -1;
//sort numbers and renumber the norm1 norm2 arrays
sort(norm,norm+2*number_of_sets);
int new_norm = 0;
int start = 0;
for(int i=1;i<2*number_of_sets;i++){
    if (norm[i-1] != -1 && norm[i] > (norm[i-1]+1)){
        new_norm = norm[i-1] + 1;
        start = i+1;
        for (int j=0; j<number_of_sets; j++){
            if (norm1[j] == norm[i]){norm1[j]=new_norm;}
            if (norm2[j] == norm[i]){norm2[j]=new_norm;}
        }
        break;
    }
}
if (start != 0){
    for(int i=start;i<2*number_of_sets;i++){
        if (norm[i-1] != -1 && norm[i] > (norm[i-1])){
            new_norm++;//wrong
            for (int j=0; j<number_of_sets; j++){
                if (norm1[j] == norm[i]){norm1[j]=new_norm;}
                if (norm2[j] == norm[i]){norm2[j]=new_norm;}
            }
        }
    }
}





  return 0;
}



