{
    double Q2[2000], R[2000],dR[2000], chi2[2000];
    double Q2_1[2000], Q2_2[2000], Q2_3[2000], Q2_4[2000];
    double R_1[2000], R_2[2000], R_3[2000], R_4[2000];
    double dR_1[2000], dR_2[2000], dR_3[2000], dR_4[2000];
    double chi2_1[2000], chi2_2[2000], chi2_3[2000], chi2_4[2000];
    int count[4];
    double weighted_center[4];
    double weighted_bar[4];


    ifstream inFile[5];
    TGraphErrors *gr[5];
    inFile[0].open("R_fit_Mainz_Carl_norm.dat");
    inFile[1].open("R_fit_Mainz_RS.dat");
    inFile[2].open("R11_fit_Mainz_global.dat");
    inFile[3].open("z2_fit_Mainz_global.dat");
    inFile[4].open("poly10_fit_Mainz_global.dat");
    //inFile[3].open("R_fit_PRad.dat");
    //inFile[4].open("R_fit_PRad_0.0038.dat");
    for (int i=0; i<5; i++){
        inFile[i]>>count[i];
        double sum_center = 0.0;
        double sum_bar = 0.0;
        for (int j=0; j<count[i]; j++){
            inFile[i]>>Q2[j]>>R[j]>>dR[j]>>chi2[j];
            sum_bar = sum_bar + 1/pow(dR[j],2);
            sum_center = sum_center + R[j]/pow(dR[j],2);
        }
        inFile[i].close();
        gr[i] = new TGraphErrors(count[i], Q2, R, 0, dR);
        weighted_bar[i] = TMath::Sqrt(1/sum_center);
        weighted_center[i] = sum_center/sum_bar;
    }


    TCanvas* c1 = new TCanvas("c1", "c1", 1000, 800);
    c1->SetBottomMargin(0.15);
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.05);


    //TGraphErrors *gr1 = new TGraphErrors(count[1], Q2, R, 0, dR);

    gr[2]->SetTitle("");
    gr[2]->GetXaxis()->SetLabelSize(0.045);
    gr[2]->GetYaxis()->SetLabelSize(0.045);
    gr[2]->GetXaxis()->SetTitleSize(0.055);
    gr[2]->GetYaxis()->SetTitleSize(0.055);
    gr[2]->GetYaxis()->SetTitle("r_{p} [fm]");
    gr[2]->GetXaxis()->SetTitle("Q^{2}_{max} [(GeV/c)^{2}]");
    gr[2]->GetYaxis()->CenterTitle();
    gr[2]->GetXaxis()->CenterTitle();
    gr[2]->GetXaxis()->SetTitleFont(62);
    gr[2]->GetXaxis()->SetLabelFont(62);
    gr[2]->GetYaxis()->SetTitleFont(62);
    gr[2]->GetYaxis()->SetLabelFont(62);
    gr[2]->GetYaxis()->SetTitleOffset(1.1);
    gr[2]->GetXaxis()->SetLimits(0.011, 0.11);
    gr[2]->GetYaxis()->SetRangeUser(0.78, 0.91);
    gr[0]->SetMarkerStyle(25);
    gr[0]->SetMarkerSize(2.0);
    gr[0]->SetMarkerColor(4);
    gr[0]->SetLineColor(4);
    gr[0]->SetLineWidth(3.0);
    gr[1]->SetMarkerStyle(23);
    gr[1]->SetMarkerSize(2.0);
    gr[1]->SetMarkerColor(2);
    gr[1]->SetLineColor(2);
    gr[1]->SetLineWidth(3.0);
    gr[2]->SetMarkerStyle(20);
    gr[2]->SetMarkerSize(2.0);
    gr[2]->SetMarkerColor(1);
    gr[2]->SetLineColor(1);
    gr[2]->SetLineWidth(3.0);
    gr[3]->SetMarkerStyle(22);
    gr[3]->SetMarkerSize(2.0);
    gr[3]->SetMarkerColor(6);
    gr[3]->SetLineColor(6);
    gr[3]->SetLineWidth(3.0);
    gr[4]->SetMarkerStyle(28);
    gr[4]->SetMarkerSize(2.0);
    gr[4]->SetMarkerColor(1);
    gr[4]->SetLineColor(1);
    gr[4]->SetLineWidth(3.0);
    gr[2]->Draw("AP");
    //gr[1]->Draw("Psame");
    gr[3]->Draw("Psame");
    //gr[4]->Draw("Psame");
    gr[0]->Draw("Psame");
    //gr[3]->SetTitle("");
    //gr[3]->GetXaxis()->SetLabelSize(0.045);
    //gr[3]->GetYaxis()->SetLabelSize(0.045);
    //gr[3]->GetXaxis()->SetTitleSize(0.055);
    //gr[3]->GetYaxis()->SetTitleSize(0.055);
    //gr[3]->GetYaxis()->SetTitle("r_{fit} [fm]");
    //gr[3]->GetXaxis()->SetTitle("Q^{2}_{max} [(GeV/c)^{2}]");
    //gr[3]->GetYaxis()->CenterTitle();
    //gr[3]->GetXaxis()->CenterTitle();
    //gr[3]->GetXaxis()->SetTitleFont(62);
    //gr[3]->GetXaxis()->SetLabelFont(62);
    //gr[3]->GetYaxis()->SetTitleFont(62);
    //gr[3]->GetYaxis()->SetLabelFont(62);
    //gr[3]->GetYaxis()->SetTitleOffset(1.1);
    //gr[3]->GetXaxis()->SetLimits(0.0, 0.062);
    //gr[3]->GetYaxis()->SetRangeUser(0.73, 0.90);
    //gr[3]->Draw("AP");
    ////gr[3]->Draw("Psame");
    //gr[3]->SetMarkerStyle(22);
    //gr[3]->SetMarkerSize(1.2);
    //gr[3]->SetMarkerColor(4);
    //gr[3]->SetLineColor(4);
    //gr[3]->SetLineWidth(2.5);
    //gr[4]->Draw("Psame");
    //gr[4]->SetMarkerStyle(25);
    //gr[4]->SetMarkerSize(1.2);
    //gr[4]->SetMarkerColor(6);
    //gr[4]->SetLineColor(6);
    //gr[4]->SetLineWidth(2.5);
    //TLegend* leg = new TLegend(0.20, 0.18, 0.85, 0.43);
    TLegend* leg = new TLegend(0.25, 0.18, 0.93, 0.43);
    leg->AddEntry(gr[0], "Rational (1,1) fit with normalized G_{E}^{p} data", "ep");
    //leg->AddEntry(gr[1], "Rational (1,1) fit with RS-separated G_{E}^{p} data", "ep");
    leg->AddEntry(gr[2], "Rational (1,1) fit with #sigma_{exp}/#sigma_{dipole} data", "ep");
    leg->AddEntry(gr[3], "Polynomial-Z (2) fit with normalized G_{E}^{p} data", "ep");
    //leg->AddEntry(gr[4], "polynomial(10) fit with 34 sets of #sigma_{exp}/#sigma_{dipole} data", "ep");
    //TLegend* leg = new TLegend(0.44, 0.23, 0.93, 0.43);
    //leg->AddEntry(gr[3], "PRad", "ep");
    //leg->AddEntry(gr[4], "PRad start from 0.0038 GeV^{2}", "ep");
    leg->SetMargin(0.1);
    leg->SetTextSize(0.04);
    leg->Draw("same");

   // TLine* line[2];
   // line[0] = new TLine(0.0, weighted_center[0], 0.6, weighted_center[0]);
   // line[1] = new TLine(0.0, weighted_center[1], 0.6, weighted_center[1]);
   // line[0]->SetLineWidth(1);
   // line[0]->SetLineColor(1);
   // //line[0]->Draw("same");
   // line[1]->SetLineWidth(1);
   // line[1]->SetLineColor(2);
   // //line[1]->Draw("same");
   // //line[0] = new TLine(0.0, weighted_center[3], 0.062, weighted_center[3]);
   // //line[1] = new TLine(0.0, weighted_center[4], 0.062, weighted_center[4]);
   // //line[0]->SetLineWidth(1);
   // //line[0]->SetLineColor(4);
   // //line[0]->Draw("same");
   // //line[1]->SetLineWidth(1);
   // //line[1]->SetLineColor(6);
   // //line[1]->Draw("same");
   // TBox* MainzBox1 = new TBox(0.0, weighted_center[0]-weighted_bar[0], 0.6, weighted_center[0]+weighted_bar[0]);
   // MainzBox1->SetFillStyle(3001);
   // MainzBox1->SetFillColorAlpha(14, 0.3);
   // TBox* MainzBox2 = new TBox(0.0, weighted_center[1]-weighted_bar[1], 0.6, weighted_center[1]+weighted_bar[1]);
   // MainzBox2->SetFillStyle(3001);
   // MainzBox2->SetFillColorAlpha(45, 0.3);
   // //MainzBox1->Draw("same");
   // //MainzBox2->Draw("same");

   // TBox* PRadBox1 = new TBox(0.0, weighted_center[2]-weighted_bar[2], 0.062, weighted_center[2]+weighted_bar[2]);
   // PRadBox1->SetFillStyle(3001);
   // PRadBox1->SetFillColorAlpha(4, 0.3);
   // TBox* PRadBox2 = new TBox(0.0, weighted_center[3]-weighted_bar[3], 0.062, weighted_center[3]+weighted_bar[3]);
   // PRadBox2->SetFillStyle(3001);
   // PRadBox2->SetFillColorAlpha(6, 0.3);
   // //PRadBox1->Draw("same");
   // //PRadBox2->Draw("same");

}
