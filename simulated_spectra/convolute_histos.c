// get energy resolution
TFile* file_res = TFile::Open("calibration_WT-20_ch000.dat.root_spectrum_calibrated_0-1124s.root_resolution_function.root");

TCanvas* c2 = (TCanvas*) file_res->Get("c2");
TGraphErrors* graph = (TGraphErrors*) c2->GetPrimitive("Graph");
TF1* fresolution = graph->GetFunction("fitFunction");

TString ffolder="original/";
TString ffolder_conv="convoluted/";

void convolute_histos() {
    ConvHist("simulated_bkg_muons.root");
    ConvHist("simulated_bkg_Rn.root");
    ConvHist("simulated_bkg_cryostat_Co56.root");
    ConvHist("simulated_bkg_cryostat_Co57.root");
    ConvHist("simulated_bkg_cryostat_Co58.root");
    ConvHist("simulated_bkg_cryostat_Co60.root");
    ConvHist("simulated_bkg_cryostat_K40.root");
    ConvHist("simulated_bkg_cryostat_Mn54.root");
    ConvHist("simulated_bkg_cryostat_Th232.root");
    ConvHist("simulated_bkg_cryostat_U238.root");
    ConvHist("simulated_bkg_CuShielding_combined_Co56.root");
    ConvHist("simulated_bkg_CuShielding_combined_Co57.root");
    ConvHist("simulated_bkg_CuShielding_combined_Co58.root");
    ConvHist("simulated_bkg_CuShielding_combined_Co60.root");
    ConvHist("simulated_bkg_CuShielding_combined_K40.root");
    ConvHist("simulated_bkg_CuShielding_combined_Mn54.root");
    ConvHist("simulated_bkg_CuShielding_combined_Th232.root");
    ConvHist("simulated_bkg_CuShielding_combined_U238.root");
    ConvHist("simulated_bkg_Ge_crystal_Ge68.root");
    ConvHist("simulated_bkg_Ge_crystal_Zn65.root");
    ConvHist("simulated_bkg_Ge_crystal_Co56.root");
    ConvHist("simulated_bkg_Ge_crystal_Co57.root");
    ConvHist("simulated_bkg_Ge_crystal_Co58.root");
    ConvHist("simulated_bkg_Ge_crystal_Co60.root");
    ConvHist("simulated_bkg_Ge_crystal_Mn54.root");
    ConvHist("simulated_bkg_innerPbShieldingCombined_Pb210.root");
    ConvHist("simulated_bkg_innerPbShieldingCombinedFitted_Pb210.root");
    
}

void ConvHist(TString fileName) {
    

    // open root file
    TFile* dataFile = TFile::Open(ffolder+fileName);
    
    if (!dataFile)
    {
        cout << "no root file found!" << endl;
        abort();
        
    }
    TH1D* hist = (TH1D*) dataFile->Get("hist");
    hist->SetDirectory(0);
    dataFile->Close();
    
    int nBins = hist->GetNbinsX();
    int Bin_min = hist->GetBinLowEdge(1);
    int Bin_max = hist->GetBinLowEdge(nBins+1);
    
    
    //// convolute with resolution
    // gaussian histogram
    TH1D* hgaus = new TH1D("hgaus","",nBins,Bin_min,Bin_max);
    
    // convoluted histogram
    TH1D* hist_conv = new TH1D("hist","",nBins,Bin_min,Bin_max);
    
    // convolute simulated histogram
    for (int j=1; j<nBins; ++j) {
        
        double energy = hist->GetBinCenter(j);
        double counts = hist->GetBinContent(j);
        double xBin,yBin;
        
        for (int i=1; i<=nBins; i++) {
            
            xBin = hgaus->GetBinCenter(i);
            yBin = TMath::Gaus(xBin, energy, fresolution->Eval(energy));
            
            hgaus->SetBinContent(i,yBin);
            
        }
        
        hgaus->Scale(counts/hgaus->Integral());
        
        hist_conv->Add(hgaus);
        
        
    }
    
    TFile* file = new TFile(ffolder_conv+fileName,"recreate");
    hist_conv->Write();
    file->Close();
    
}






