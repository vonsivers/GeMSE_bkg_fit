void variable_binning() {
    
    
    TString FileName = "/Users/sivers/Desktop/GeMSE/background_simulations/background_fit/simulated_spectra/convoluted/simulated_bkg_CuShielding_combined_Th232.root";
    
    // get original histo
    
    // check for root file
    TFile* File = TFile::Open(FileName);
    
    if (!File) {
        
        std::cout << "###### ERROR: could not open " << FileName << std::endl;
        return 0;
    }
    
    if (!File->GetListOfKeys()->Contains("hist")) {
        std::cout << "###### ERROR: no histogram in file " << FileName << std::endl;
        return 0;
    }
    
    TH1D* hist = (TH1D*) File->Get("hist");
    hist->SetDirectory(0);
    File->Close();
    
    // scale to counts/keV
    hist->Scale(1,"width");

    cout << "Counts in original histo: " << hist->Integral(hist->FindBin(100),hist->FindBin(3000),"width") << endl;
    
    
    // get histo with variable binning
    
    TH1D* hist_var = GetSimHist(FileName);
    
    // scale to counts/keV
    hist_var->Scale(1,"width");
    
    cout << "Counts in variable histo: " << hist_var->Integral(hist_var->FindBin(100),hist_var->FindBin(3000),"width") << endl;
    
    // Draw
    TCanvas* c = new TCanvas("c");
    hist->SetLineColor(1);
    hist_var->SetLineColor(2);
    hist->Draw();
    hist_var->Draw("same");
    
}

// ----------------------------------------------------
// get simulated histo
// ----------------------------------------------------


TH1D* GetSimHist(TString FileName) {
    
    
    double frange_low = 100.;
    double frange_high = 3000.;
    int fnBins = 300;
    double fBinSlope = 1.5;
    
    vector<double> fxBins;
    
    // lower edge of bins
    for (int i=0; i<=fnBins; ++i) {
        
        //fxBins[i] = frange_low+pow(i,fBinSlope)*(frange_high-frange_low)/pow(fnBins,fBinSlope);
        fxBins.push_back(frange_low+pow(i,fBinSlope)*(frange_high-frange_low)/pow(fnBins,fBinSlope));
        
    }

    
    // check for root file
    TFile* File = TFile::Open(FileName);
    
    if (!File) {
        
        std::cout << "###### ERROR: could not open " << FileName << std::endl;
        return 0;
    }
    
    if (!File->GetListOfKeys()->Contains("hist")) {
        std::cout << "###### ERROR: no histogram in file " << FileName << std::endl;
        return 0;
    }
    
    TH1D* hist = (TH1D*) File->Get("hist");
    hist->SetDirectory(0);
    File->Close();
    
    TH1D* hist_range = new TH1D("hist_var","",fnBins,&fxBins.at(0));
    //TH1D* hist_range = new TH1D("hist","",fnBins,frange_low,frange_high);
    
    for (int i=1; i<=fnBins; ++i) {
        
        // Low/High edge of new bin
        double LowEdge_new = hist_range->GetBinLowEdge(i);
        double HighEdge_new = hist_range->GetBinLowEdge(i+1);
        
        // Bin number range of old bins
        int BinL_old = hist->FindBin(LowEdge_new);
        int BinR_old = hist->FindBin(HighEdge_new);
        
        // content of new bin
        double content=0.;
        
        // fraction of old bin in new
        double fraction;
        
        // loop over old bins
        for (int j=BinL_old; j<=BinR_old; ++j) {
            
            // Low/High edge of old bin
            double LowEdge_old = hist->GetBinLowEdge(j);
            double HighEdge_old = hist->GetBinLowEdge(j+1);
            
            // new bin is contained in old bin
            if (LowEdge_old<=LowEdge_new && HighEdge_new<=HighEdge_old) {
                
                // calculate fraction of new bin in old bin
                fraction = hist_range->GetBinWidth(i)/hist->GetBinWidth(j);
                
            }
            else if (LowEdge_old<=LowEdge_new && HighEdge_new>HighEdge_old){ // new bin extends over several old bins
                
                // calculate fraction of old bin in new bin
                fraction = (HighEdge_old-LowEdge_new)/hist->GetBinWidth(j);
            }
            else if (LowEdge_old>LowEdge_new && HighEdge_new<=HighEdge_old){ // new bin extends over several old bins
                
                // calculate fraction of old bin in new bin
                fraction = (HighEdge_new-LowEdge_old)/hist->GetBinWidth(j);
            }
            else if (LowEdge_old>LowEdge_new && HighEdge_new>HighEdge_old){ // old bin contained in new bin
                
                // calculate fraction of old bin in new bin
                fraction = 1.;
            }
            
            // calculate content of new bin
            content += hist->GetBinContent(j)*fraction;
            
            
        }
        
        hist_range->SetBinContent(i,content);
    }
    
    return hist_range;
}




