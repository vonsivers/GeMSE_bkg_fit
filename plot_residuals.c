void plot_residuals() {
    
    
    TString FileName = "/Users/sivers/Desktop/GeMSE/background_simulations/background_fit/results/noRn_noCryostat/test_fit_All_stack.root";

        
        TFile* file_stack = TFile::Open(FileName);
        
        if (!file_stack) {
            
            std::cout << "###### ERROR: could not open " << FileName << std::endl;
            return 0;
        }
        
        if (!file_stack->GetListOfKeys()->Contains("c1")) {
            std::cout << "###### ERROR: no canvas named c1 in file " << FileName << std::endl;
            return 0;
        }
        
        TCanvas* c = (TCanvas*) file_stack->Get("c1");
        
        TH1D* hmeas = (TH1D*) c->FindObject("measurement");
        
        THStack* stack = (THStack*) c->GetPrimitive("");
        
        TH1D* hsim = (TH1D*) stack->GetStack()->Last();
    
    int nBins = hmeas->GetNbinsX();
    
    // get xbins
    double* xbins = hmeas->GetXaxis()->GetXbins()->fArray;
    
    // histo for residuals
    TH1D* hresidual = new TH1D("hres",";Energy (keV);(MC-Exp)/Exp",nBins,xbins);
    
    double y_exp, y_exp_err, y_sim;
    
    for (int i=1; i<=nBins; ++i) {
        
        y_exp = hmeas->GetBinContent(i);
        y_exp_err = hmeas->GetBinError(i);
        
        y_sim = hsim->GetBinContent(i);
        
        hresidual->SetBinContent(i,(y_exp-y_sim)/y_exp);
        hresidual->SetBinError(i,(y_exp_err/y_exp)*(y_exp-y_sim)/y_exp);
        
    }
    
    
    TLine* line = new TLine(100.,0.,3000.,0.);
    
        // Draw
        TCanvas* c2 = new TCanvas("c2");
        gStyle->SetOptStat(0);
    //graph_ratio->SetMarkerStyle(20);
    line->SetLineStyle(2);
    
    hresidual->Draw("");
    line->Draw("same");
    
    }



