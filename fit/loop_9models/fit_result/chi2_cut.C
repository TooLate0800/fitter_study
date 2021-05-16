#include <string>
#include <iostream>

void chi2_cut(){

    for (int j = 0; j <  9; j++){
    ifstream inf(Form("R11_xsfit_all_GMfromRatio_Cut_01GeV_Model%d.txt",j+1));
    ofstream outf(Form("R11_xsfit_RE_GMfromRatio_Cut_01GeV_chi2cut_Model%d.txt",j+1));
    double chi2, R, Rerr; 
    for (int i=0; i<10000; i++){
        inf>>chi2>>R>>Rerr;
        //cout<<chi2<<" "<<R<<endl;
        if (chi2<4000){
            outf<<R<<endl;
        }
    }
    inf.close();
    outf.close();
    }


}

