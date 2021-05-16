double carlX[2000] = {0};
double carlEX[2000] = {0};
double carlY[2000] = {0};
double carlEY[2000] = {0};
double carlCount = 0;

double fitPara[6] = {8.30608e-01, 1.12207e-01, 1.00022e+00,  9.98316e-01, 0.0074757, 0.0116};

Double_t DipoleGE(Double_t* x, Double_t* par)
{
    Double_t xx = x[0];
    double dipolePara = 0.66;
    return 1./pow(1. + xx/dipolePara, 2);
}


void LoadTable()
{
    string tmp[3];
    ifstream inFile;
    inFile.open("Carl-norm.dat");
    inFile>>tmp[0]>>tmp[1]>>tmp[2];
    
    for (int i=0; i<1422; i++){
        inFile>>tmp[0]>>tmp[1]>>tmp[2];
        
        carlX[i] = stof(tmp[0]) / 25.68189504;
        carlY[i] = stof(tmp[1]);
        carlEY[i] = stof(tmp[2]);
        
        carlCount ++;
        
        cout<<carlX[i]<<" "<<carlY[i]<<" "<<carlEY[i]<<endl;
    }
}


void DrawCarlGE()
{
    gStyle->SetEndErrorSize(2);
    gStyle->SetFrameLineWidth(2);
    Double_t base = 0.83;
    Double_t top = 1.01;
    TFile* f1 = new TFile("tmp_1101.00MeV.root", "READ");
    TFile* f2 = new TFile("tmp_2143.00MeV.root", "READ");

    
    TGraphErrors* g1 = (TGraphErrors*)f1->Get("GE_Q2_Graph_1");
    
    TGraphErrors* g2 = (TGraphErrors*)f2->Get("GE_Q2_Graph_1");

    LoadTable();
    
    Double_t x1[100];
    Double_t y1[100];
    Double_t ex1[100] = {0};
    Double_t ey1[100];
    Double_t sy1[100];
    Int_t count1 = 0;
    
    for (int i=0; i<g1->GetN(); i++){

        double thisx, thisy;
        g1->GetPoint(i, thisx, thisy);
        //if (thisx < 1.4e-3 || thisx > 4.5e-3) continue;
        x1[count1] = thisx;
        y1[count1] = thisy/fitPara[2];// + (1.- fitPara[2]);
        ey1[count1] = g1->GetErrorY(i);
        
        count1++;
    }
    TGraphErrors* graph1 = new TGraphErrors(count1, x1, y1, ex1, ey1);
    graph1->SetMarkerColor(2);
    graph1->SetMarkerStyle(21);
    graph1->SetMarkerSize(1);
    
    Double_t x2[100];
    Double_t y2[100];
    Double_t ex2[100] = {0};
    Double_t ey2[100];
    Double_t sy2[100];
    Int_t count2 = 0;
    
    for (int i=0; i<g2->GetN(); i++){

        double thisx, thisy;
        g2->GetPoint(i, thisx, thisy);
        //if ( thisx > 0.017 || thisx < 2e-4) continue;
        x2[count2] = thisx;
        y2[count2] = thisy/fitPara[3];// + (1.-fitPara[3]);
        ey2[count2] = g2->GetErrorY(i);
        
        count2++;
    }
    TGraphErrors* graph2 = new TGraphErrors(count2, x2, y2, ex2, ey2);
    graph2->SetMarkerColor(2);
    graph2->SetMarkerStyle(20);
    graph2->SetMarkerSize(1);
    
    TGraphErrors* graphc = new TGraphErrors(carlCount, carlX, carlY, carlEX, carlEY);
    graphc->SetMarkerColor(1);
    graphc->SetMarkerStyle(26);
    graphc->SetMarkerSize(0.7);
    
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    /*c1->SetGridx();
    c1->SetGridy();
    c1->SetTickx();
    c1->SetTicky();
    */
    //c1->SetLogx();
    //graph2->GetXaxis()->SetMoreLogLabels();
    c1->SetMargin(0.12, 0.05, 0.12, 0.05);
    graph2->GetYaxis()->SetNdivisions(508);
    graph2->GetXaxis()->SetNdivisions(508);
    graph2->GetYaxis()->SetRangeUser(base-0.03, top);
    graph2->GetXaxis()->SetLimits(2e-4, 0.06);
    graph2->GetXaxis()->SetLabelSize(0.045);
    graph2->GetYaxis()->SetLabelSize(0.045);
    graph2->GetXaxis()->SetTitleSize(0.055);
    graph2->GetYaxis()->SetTitleSize(0.055);
    graph2->GetXaxis()->CenterTitle();
    graph2->GetYaxis()->CenterTitle();
    graph2->SetTitle("");
    graph2->GetXaxis()->SetTitle("Q^{2} (GeV^{2}/c^{2})");
    graph2->GetYaxis()->SetTitleOffset(1.);
    graph2->GetYaxis()->SetTitle("G_{E}^{p}");
    graph2->GetXaxis()->SetTitleFont(62);
    graph2->GetYaxis()->SetTitleFont(62);
    graph2->GetXaxis()->SetLabelFont(62);
    graph2->GetYaxis()->SetLabelFont(62);
    graph2->Draw("AP");
    //graph2->Fit("pol1");
    graphc->Draw("Psame");
    graph1->Draw("Psame");
    graph2->Draw("Psame");
    
    
    
    double SimonX[18] = {0.14, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.80, 0.85, 0.90, 1.00, 1.2, 1.4};
    double SimonEX[18] = {0};
    double SimonY[18] = {0.9858, 0.9767, 0.9722, 0.9619, 0.9612, 0.9511, 0.9463, 0.9428, 0.9353, 0.9320, 0.9246, 0.9165, 0.9064, 0.9043, 0.9023, 0.8839, 0.8637, 0.8466};
    double SimonEY[18] = {0.22, 0.25, 0.22, 0.43, 0.22, 0.50, 0.16, 0.31, 0.16, 0.30, 0.18, 0.35, 0.27, 0.32, 0.37, 0.80, 0.60, 0.40};
    int SimonCount = 18;
    
    for (int i=0; i<SimonCount; i++){
        SimonX[i] /= 25.68189504;
        SimonEY[i] /= 100;
        SimonEY[i] = SimonY[i]*SimonEY[i];
    }
    
    TGraphErrors* graph3 = new TGraphErrors(SimonCount, SimonX, SimonY, SimonEX, SimonEY);
    graph3->SetMarkerColor(4);
    graph3->SetMarkerStyle(23);
    graph3->SetMarkerSize(1);
    graph3->Draw("Psame");
    
    double MurphyX[11] = {0.150, 0.295, 0.297, 0.347, 0.390, 0.396, 0.440, 0.493, 0.530, 0.678, 0.794};
    double MurphyY[11] = {0.981, 0.969, 0.966, 0.961, 0.962, 0.955, 0.951, 0.942, 0.939, 0.922, 0.914};
    double MurphyEX[11] = {0};
    double MurphyEY[11] = {0.005, 0.005, 0.003, 0.003, 0.003, 0.003, 0.003, 0.004, 0.009, 0.004, 0.004};
    int MurphyCount = 11;
    
    for (int i=0; i<MurphyCount; i++){
        MurphyX[i] /= 25.68189504;
    }
    TGraphErrors* graph4 = new TGraphErrors(MurphyCount, MurphyX, MurphyY, MurphyEX, MurphyEY);
    graph4->SetMarkerColor(kGreen+2);
    graph4->SetMarkerStyle(22);
    graph4->SetMarkerSize(1);
    graph4->Draw("Psame");
     
     
    double HandX[15] = {0.28, 0.30, 0.30, 0.36, 0.49, 0.57, 0.60, 0.62, 0.79, 0.93, 1.00, 1.05, 1.30, 1.38, 1.60};
    double HandY[15] = {0.973, 0.959, 0.974, 0.967, 0.933, 0.915, 0.940, 0.922, 0.920, 0.848, 0.881, 0.884, 0.867, 0.873, 0.849};
    double HandEX[15] = {0};
    double HandEY[15] = {0.014, 0.010, 0.006, 0.040, 0.009, 0.037, 0.007, 0.010, 0.037, 0.034, 0.009, 0.009, 0.025, 0.036, 0.004};
    int HandCount = 15;
    for (int i=0; i<HandCount; i++){
        HandX[i] /= 25.68189504;
    }
    TGraphErrors* graph5 = new TGraphErrors(HandCount, HandX, HandY, HandEX, HandEY);
    graph5->SetMarkerColor(kViolet);
    graph5->SetMarkerStyle(25);
    graph5->SetMarkerSize(1);
    //graph5->Draw("Psame");
    
     
     
    TLegend* leg = new TLegend(0.15, 0.18, 0.6, 0.38);
    
    leg->AddEntry(graph1, "PRad 1.1 GeV data", "p");
    leg->AddEntry(graph2, "PRad 2.2 GeV data", "p");
    leg->AddEntry(graphc, "Mainz 2010", "p");
    leg->AddEntry(graph3, "Simon 1980", "p");
    leg->AddEntry(graph4, "Murphy 1974", "p");
    //leg->AddEntry(graph5, "Hand 1963", "p");
    leg->SetTextSize(0.04);
    leg->Draw("same");
    
}
