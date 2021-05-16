void Draw_bias_rms()
{
    gStyle->SetFrameLineWidth(2);
    gStyle->SetEndErrorSize(5);
    const int NPOINT = 10;
    //Double_t x[NPOINT] = {0.8751, 0.879, 0.8759, 0.84087, 0.84184, 0.8335, 0.877, 0.831, 0.879, 0.8414, 0.870, 0.833, 0.831, 0.831}; //CODATA, CODATA_scattering, CODATA H spectroscopic, mu-p 2013, mu-p 2010, H spectroscopic 2017, H 2018, PRad, Bernauer, CODATA-2018, ISR, York, PRad2-proj1, PRad2-proj2
    //Double_t ex[NPOINT] = {0.0061, 0.011, 0.0077, 0.00039, 0.00067, 0.0095, 0.013, 0.014, 0.0078, 0.0019, 0.028, 0.010, 0.0061, 0.0057};
    //Pohl-2010, mu-p 2013, *Beyer 2017, CODATA-2018, *Bezginov 2019, *PRad, PRad-II, Grinin, CODATA H spectroscopic, Fleurbaey, 
    Double_t x[NPOINT] = {0.84184, 0.84087, 0.8335, 0.8414, 0.833, 0.831, 0.831, 0.8482, 0.8759, 0.877};
    Double_t ex[NPOINT] = {0.00067, 0.00039, 0.0095, 0.0019, 0.010, 0.014, 0.0036, 0.0038, 0.0077, 0.013};
    //Double_t y[NPOINT] = {3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.0, 3.5};
    Double_t y[NPOINT] = {2.5, 2.0, 1.5, 2.0, 0.5, 1.0, 0.5, 0.0, 1.2, 0.7};
    Double_t ey[NPOINT] = {0};
    //electronic weighted average
    Double_t center = (x[2]/pow(ex[2],2)+x[4]/pow(ex[4],2)+x[5]/pow(ex[5],2))/(1/pow(ex[2],2)+1/pow(ex[4],2)+1/pow(ex[5],2));
    Double_t uncertainty = TMath::Sqrt(1/(pow(ex[2],-2)+pow(ex[4],-2)+pow(ex[5],-2))); 
    Double_t center_2 = (x[2]/pow(ex[2],2)+x[4]/pow(ex[4],2)+x[5]/pow(ex[5],2)+x[7]/pow(ex[7],2))/(1/pow(ex[2],2)+1/pow(ex[4],2)+1/pow(ex[5],2)+1/pow(ex[7],2));
    Double_t uncertainty_2 = TMath::Sqrt(1/(pow(ex[2],-2)+pow(ex[4],-2)+pow(ex[5],-2)+pow(ex[7],-2))); 
    cout<<center<<","<<center_2<<endl; 
    cout<<uncertainty<<","<<uncertainty_2<<endl; 
    Double_t ax = 0.;
    Double_t ay = 0;
    Double_t left = 0.78;
    Double_t right = 0.935;
    Double_t low = -0.5;
    Double_t up = 3.0;
    TGraphErrors* graph = new TGraphErrors(1, &ax, &ay, &ax, &ay);
    graph->SetTitle("");
    graph->GetXaxis()->SetTitle("Proton charge radius r_{p} [fm]");
    graph->GetXaxis()->SetTitleSize(0.06);
    graph->GetXaxis()->SetTitleOffset(1.1);
    graph->GetXaxis()->SetLabelSize(0.05);
    graph->GetXaxis()->SetLimits(left, right);
    graph->GetYaxis()->SetRangeUser(low, up);
    //graph->GetXaxis()->SetLimits(0.78, 0.88);
    //graph->GetYaxis()->SetRangeUser(-0.5, 4);
    graph->GetYaxis()->SetNdivisions(0);
    
    TCanvas* c1 = new TCanvas("c1", "c1", 1200, 600);
    //TCanvas* c1 = new TCanvas("c1", "c1", 1200, 800);
    c1->SetBottomMargin(0.15);
    c1->SetLeftMargin(0.05);
    c1->SetRightMargin(0.05);

    graph->GetXaxis()->SetTitleFont(62);
    graph->GetXaxis()->SetLabelFont(62);
    graph->GetXaxis()->CenterTitle();
    graph->Draw("AP");

    
    TBox* muBox = new TBox(x[1]-ex[1], low, x[1]+ex[1], up);
    muBox->SetFillStyle(3001);
    muBox->SetFillColorAlpha(807, 0.3);

    TBox* eBox = new TBox(center-uncertainty, low, center+uncertainty, up);
    eBox->SetFillStyle(3001);
    eBox->SetFillColorAlpha(kAzure-9, 0.3);

    TBox* eBox2 = new TBox(center_2-uncertainty_2, low, center_2+uncertainty_2, up);
    eBox2->SetFillStyle(3001);
    eBox2->SetFillColorAlpha(16, 0.3);
   
    TBox* HBox = new TBox(x[8]-ex[8], low, x[8]+ex[8], up);
    HBox->SetFillStyle(3001);
    HBox->SetFillColorAlpha(kAzure-9, 0.3);

    //eBox->Draw("same");
    //eBox2->Draw("same");
    muBox->Draw("same");
    HBox->Draw("same");

    TGraphErrors* g[NPOINT];
    for (int i=0; i<NPOINT; i++){
        g[i] = new TGraphErrors(1, &x[i], &y[i], &ex[i], &ey[i]);
        g[i]->SetMarkerSize(1.5);
        g[i]->SetLineWidth(2);
    }
    TLatex tBeam;
    tBeam.SetTextAlign(12);
    tBeam.SetTextSize(0.043);
    
    
    tBeam.SetTextSize(0.04);
    g[0]->SetMarkerColor(46);
    g[0]->SetLineColor(46);
    g[0]->SetMarkerStyle(20);
    g[0]->SetMarkerSize(0.8);
    g[0]->Draw("Psame");
    //tBeam.SetTextAlign(12);
    tBeam.DrawLatex(left+0.001, y[0], "#color[46]{Pohl 2010 (#muH spect.)}");
    
    g[1]->SetMarkerColor(807);
    g[1]->SetLineColor(807);
    g[1]->SetMarkerStyle(20);
    g[1]->SetMarkerSize(0.8);
    g[1]->Draw("Psame");
    tBeam.DrawLatex(left+0.001, y[1], "#color[807]{Antognini 2013 (#muH spect.)}");

    g[2]->SetMarkerColor(6);
    g[2]->SetLineColor(6);
    g[2]->SetMarkerSize(1.5);
    g[2]->SetMarkerStyle(20);
    g[2]->Draw("Psame");
    tBeam.DrawLatex(left+0.001, y[2], "#color[6]{Beyer 2017 (H spect.)}");


    g[3]->SetMarkerColor(9);
    g[3]->SetLineColor(9);
    g[3]->SetMarkerSize(1.5);
    g[3]->SetMarkerStyle(20);
    //g[3]->Draw("Psame");
    //tBeam.DrawLatex(left+0.001, y[3], "#color[9]{CODATA-2018}");
    
    g[4]->SetMarkerColor(403);
    g[4]->SetLineColor(403);
    g[4]->SetMarkerSize(1.5);
    g[4]->SetMarkerStyle(20);
    g[4]->Draw("Psame");
    tBeam.DrawLatex(left+0.001, y[4], "#color[403]{Bezginov 2019 (H spect.)}");
    
    g[5]->SetMarkerColor(1);
    g[5]->SetMarkerSize(2.0);
    g[5]->SetMarkerStyle(21);
    g[5]->SetLineColor(1);
    g[5]->SetLineWidth(3);
    g[5]->Draw("Psame");
    //tBeam.SetTextSize(0.03);
    tBeam.DrawLatex(left+0.001, y[5], "#color[1]{PRad exp. (ep scatt.)}");

    
    g[6]->SetMarkerColor(2);
    g[6]->SetMarkerSize(2.0);
    g[6]->SetMarkerStyle(21);
    g[6]->SetLineColor(2);
    g[6]->SetLineWidth(3);
    //g[6]->Draw("Psame");
    //tBeam.SetTextSize(0.03);
    //tBeam.DrawLatex(left+0.001, y[6], "#color[2]{PRad-II proj.}");
    
    g[7]->SetMarkerColor(36);
    g[7]->SetMarkerSize(1.5);
    g[7]->SetMarkerStyle(20);
    g[7]->SetLineColor(36);
    g[7]->SetLineWidth(3);
    g[7]->Draw("Psame");
    tBeam.DrawLatex(left+0.001, y[7], "#color[36]{Grinin 2020 (H Spect.)}");
    
    g[8]->SetMarkerColor(861);
    g[8]->SetLineColor(861);
    g[8]->SetMarkerStyle(20);
    g[8]->Draw("Psame");
    tBeam.SetTextAlign(32);
    tBeam.DrawLatex(right-0.001, y[8], "#color[861]{CODATA-2014 (H spect.)}");

    g[9]->SetMarkerColor(4);
    g[9]->SetLineColor(4);
    g[9]->SetMarkerSize(1.5);
    g[9]->SetMarkerStyle(20);
    g[9]->Draw("Psame");
    tBeam.DrawLatex(right-0.001, y[9], "#color[4]{Fleurbaey 2018 (H spect.)}");



   //TArrow *arrow = new TArrow(x[0],3.5,x[3],3.5,0.005,"<|>");
   //arrow->SetLineWidth(3);
   // arrow->SetLineColor(16);
   //arrow->Draw();
   
   
   //tBeam.SetTextSize(0.02);
   //tBeam.DrawLatex(0.855, 3.7, "#color[16]{5.6 #sigma}");
   
   //tBeam.SetTextSize(0.045);
   //tBeam.DrawLatex(0.852, 0.3, "#color[16]{2.7 #sigma}");

   TLine* line[4];
   line[0] = new TLine(center, low, center, up);
   line[1] = new TLine(x[1], low, x[1], up);
   line[2] = new TLine(center_2, low, center_2, up);
   line[3] = new TLine(x[8], low, x[8], up);
   
   line[0]->SetLineWidth(1);
   line[0]->SetLineColor(861);
   //line[0]->Draw("same");
   
   line[1]->SetLineWidth(1);
   line[1]->SetLineColor(807);
   line[1]->Draw("same");
   
   line[2]->SetLineWidth(1);
   line[2]->SetLineColor(14);
   //line[2]->Draw("same");


   line[3]->SetLineWidth(1);
   line[3]->SetLineColor(861);
   line[3]->Draw("same");


   //tBeam.SetTextSize(0.05);
   //tBeam.DrawLatex(0.808, -1.2, "#color[1]{PRad R_{p} (current) = 0.831 #pm 0.007 (stat.) #pm 0.012 (syst.) fm}");
   //tBeam.SetTextSize(0.05);
   //tBeam.DrawLatex(0.88, 3, "#font[62]{b}");
}
