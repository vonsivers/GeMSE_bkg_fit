
#include "GeMSE_bkg_fit.h"

#include <BAT/BCModelManager.h>
#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCModelOutput.h>
#include <BAT/BCSummaryTool.h>
#include <BAT/BCH1D.h>
#include <BAT/BCMTFChannel.h>
#include <BAT/BCParameter.h>

#include <TFile.h>
#include <TVectorD.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <sstream>



// ---------------------------------------------------------
GeMSE_bkg_fit::GeMSE_bkg_fit()
{}

// ---------------------------------------------------------
GeMSE_bkg_fit::~GeMSE_bkg_fit()
{}


// ----------------------------------------------------
// run main fit
// ----------------------------------------------------

int GeMSE_bkg_fit::RunFit() {
    
    if (!fParSet) {
        std::cout << "#### ERROR: parameters have not been set!" << std::endl;
        return 1;
    }
    
    
    // lower edge of bins
    for (int i=0; i<=fnBins; ++i) {
        
        //fxBins[i] = frange_low+pow(i,fBinSlope)*(frange_high-frange_low)/pow(fnBins,fBinSlope);
        fxBins.push_back(frange_low+pow(i,fBinSlope)*(frange_high-frange_low)/pow(fnBins,fBinSlope));
        
    }
    
    // measured spectrum
    if (!GetMeasuredHist(fMeas_name)) {
        return 0;
    }
    
    // add simulated spectra (name, filename, limit_low, limit_up, prior, prior mean, prior sigma)
    
    int j=0;

    if (!AddHistMuons("Muons", fSim_folder+fMuon_name, fparlimit_low[j], fparlimit_up[j],fmuon_flux_err)) {
        return 0;
    }
    j++;
    
    if (!AddHistPb("innerPbShielding_Pb210", fSim_folder+fPb210_name, fparlimit_low[j], fparlimit_up[j])) {
        return 0;
    }
    j++;
    
    if (!AddHistSim("SampleCavity_Rn222", fSim_folder+fRn_name, fparlimit_low[j], fparlimit_up[j])) {
        return 0;
    }
    j++;
    
    for (int i=0; i<fCuShielding_name.size(); ++i) {
        if (!AddHistSim("CuShielding_"+fCuShielding_isotope[i], fSim_folder+fCuShielding_name[i], fparlimit_low[j], fparlimit_up[j])) {
            return 0;
        }
        j++;
    }
    
    for (int i=0; i<fCryostat_name.size(); ++i) {
        if (!AddHistSim("Cryostat_"+fCryostat_isotope[i], fSim_folder+fCryostat_name[i], fparlimit_low[j], fparlimit_up[j])) {
            return 0;
        }
        j++;
    }
 
    for (int i=0; i<fGeCrystal_name.size(); ++i) {
        if (!AddHistSim("Ge_"+fGeCrystal_isotope[i], fSim_folder+fGeCrystal_name[i], fparlimit_low[j], fparlimit_up[j])) {
            return 0;
        }
        j++;
    }
    

    // ------------------------------
    // run fit with all spectra
    // ------------------------------


    // number of simulated spectra
    int Nprocesses = fname_process.size();
    
    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();
    
    // open log file
    BCLog::OpenLog(fresults_name+"_log.txt");
    BCLog::SetLogLevel(BCLog::warning);
    
    std::cout << "###### Running fit for all processes..." << std::endl;
    
	// create new fitter object
	BCMTF_HPGe* m = new BCMTF_HPGe();
    
    
    // set seed for reproducible results
    m->MCMCSetRandomSeed(21340);
    
    // set constant for integration
    m->SetConstant(fint_const);
    
    // set Metropolis as marginalization method
    m->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
    
    // set number of bins for marginalized distributions
    //m->SetNbins(500);
    
    // add channels (the order is important!)
    m->AddChannel("measured_spectrum");
    
    // set data
    m->SetData("measured_spectrum", fhist_meas);
    
    // loop over all simulated spectra
    for (int i=0; i<Nprocesses; ++i) {
        
        // add process
        m->AddProcess(fname_process[i], fparlimit_low[i], fparlimit_up[i]);
        
        // set template and histograms
        m->SetTemplate("measured_spectrum", fname_process[i], *(TH1D*)(fhist_process->At(i)));
        
        // set prior
        if (fprior_name[i]=="Const") {
            m->SetPriorConstant(fname_process[i]);
        }
        else if (fprior_name[i]=="Gauss") {
            m->SetPriorGauss(fname_process[i],fprior_mean[i],fprior_sigma[i]);
        }
        else {
            std::cout << "##### unknown prior " << fprior_name[i] << std::endl;
            return 0;
        }

    }
    
    
    // run marginalization to find right parameter limits
    std::cout << "#### pre-run fit to set parameter limits ..." << std::endl;
    
    m->MCMCSetPrecision(BCEngineMCMC::kLow);
    
    // do it 2 times
    for (int i=0; i<2; ++i) {
        
        m->MarginalizeAll();
        
        for (int i=0; i<Nprocesses; ++i) {
            SetNewLimits(m, fname_process[i]);
        }
    }

    
    // set precision
    if (fprecision=="low") {
        m->MCMCSetPrecision(BCEngineMCMC::kLow);
    }
    else if (fprecision=="medium") {
        m->MCMCSetPrecision(BCEngineMCMC::kMedium);
    }
    else if (fprecision=="high") {
        m->MCMCSetPrecision(BCEngineMCMC::kHigh);
    }
    else {
        std::cout << "##### unknown precision " << fprecision << std::endl;
        return 0;
    }
    

    
    TString results_name = fresults_name+"_All";
    
    // create new output object
    BCModelOutput* mout = new BCModelOutput(m, results_name+"_output.root");
    
    // create a new summary tool object
    BCSummaryTool * summary = new BCSummaryTool(m);
    
    std::cout << "#### running actual fit ..." << std::endl;
    
    // marginalize
    m->MarginalizeAll();
    
    // find global mode
    m->FindMode( m->GetBestFitParameters() );
    
    // calculate normalization
    m->Integrate();
    
    // get p-value
    m->CalculatePValue( m->GetBestFitParameters() );
    
    // print marginalized distributions
    m->PrintAllMarginalized(results_name+"_distributions.pdf");
    
    // print results file
    m->PrintResults(results_name+"_results.txt");
    
    // print summary results
    summary->PrintKnowledgeUpdatePlots(results_name+"_summary_update.pdf");
    
    // print templates and stacks
    BCMTFChannel * channel = m->GetChannel(0);
    channel->PrintTemplates(Form("%s_templates.pdf",results_name.Data()));
    m->PrintStack(0, m->GetBestFitParameters(), Form("%s_stack.pdf",results_name.Data()));
    m->PrintStack(0, m->GetBestFitParameters(), Form("%s_stack.root",results_name.Data()));

    // write marginalized distributions to output file
    mout->WriteMarginalizedDistributions();

    // calculate activities
    for (int i=0; i<Nprocesses; ++i) {
        
        // get activity
        double counts = m->GetMarginalized(fname_process[i])->GetMode();
        double counts_err_low = counts - m->GetMarginalized(fname_process[i])->GetQuantile(0.16);
        double counts_err_up = m->GetMarginalized(fname_process[i])->GetQuantile(0.84) - counts;
        
        // get upper limit
        double counts_limit = m->GetMarginalized(fname_process[i])->GetQuantile(fCL);

        
        if (fname_process[i]=="Muons") {
            factivity.push_back(counts/fcounts_muons*fmuon_flux);
            factivity_err_low.push_back(counts_err_low/fcounts_muons*fmuon_flux);
            factivity_err_up.push_back(counts_err_up/fcounts_muons*fmuon_flux);
            factivity_lim.push_back(counts_limit/fcounts_muons*fmuon_flux);
        }
        else {
            factivity.push_back(counts / fintegral[i] / ft_live);
            factivity_err_low.push_back(counts_err_low / fintegral[i] / ft_live);
            factivity_err_up.push_back(counts_err_up / fintegral[i] / ft_live);
            factivity_lim.push_back(counts_limit / fintegral[i] / ft_live);
        }
        
    }

    
    // make nice stack plot
    std::cout << "#### make nice stack plot ..." << std::endl;
    plot_stack();
    
    // close log file
    BCLog::CloseLog();
    
    
    // close output file
    mout->Close();
    
    
    // free memory
    //delete m;
    delete mout;
    delete summary;
    
    
    
    // ------------------------------
    // run bkg fit
    // ------------------------------
    
    // loop over all simulated spectra
    
    for (int j=0; j<Nprocesses; ++j) {
        TString results_name_bkg = fresults_name+"_"+fname_process[j];
        
        std::cout << "###### Running fit without " << fname_process[j] << "..." << std::endl;
        
        // create new fitter object
        BCMTF_HPGe* m_bkg = new BCMTF_HPGe();
        
        // set seed for reproducible results
        m_bkg->MCMCSetRandomSeed(21340);
        
        // set constant for integration
        m_bkg->SetConstant(fint_const);
        
        // set Metropolis as marginalization method
        m_bkg->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
        
        // set number of bins for marginalized distributions
        //m->SetNbins(500);
        
        // add channels (the order is important!)
        m_bkg->AddChannel("measured_spectrum");
        
        // set data
        m_bkg->SetData("measured_spectrum", fhist_meas);
        
        // loop over all simulated spectra
        for (int i=0; i<Nprocesses; ++i) {
            
            // leave out one process
            if (i!=j) {
                
                // add process
                m_bkg->AddProcess(fname_process[i], fparlimit_low[i], fparlimit_up[i]);
                
                // set template and histograms
                m_bkg->SetTemplate("measured_spectrum", fname_process[i], *(TH1D*)(fhist_process->At(i)));
                
                // set prior
                if (fprior_name[i]=="Const") {
                    m_bkg->SetPriorConstant(fname_process[i]);
                }
                else if (fprior_name[i]=="Gauss") {
                    m_bkg->SetPriorGauss(fname_process[i],fprior_mean[i],fprior_sigma[i]);
                }
                else {
                    std::cout << "##### unknown prior " << fprior_name[i] << std::endl;
                    return 0;
                }
            }
            
            
        }
        
        // run marginalization to find right parameter limits
        std::cout << "#### pre-run fit to set parameter limits ..." << std::endl;
        
        m_bkg->MCMCSetPrecision(BCEngineMCMC::kLow);
        
        // do it 2 times
        for (int i=0; i<2; ++i) {
            
            m_bkg->MarginalizeAll();
            
            for (int i=0; i<Nprocesses; ++i) {
                
                // leave out one process
                if (i!=j) {
                    SetNewLimits(m_bkg, fname_process[i]);
                }
            }
        }
        
        
        // set precision
        if (fprecision=="low") {
            m_bkg->MCMCSetPrecision(BCEngineMCMC::kLow);
        }
        else if (fprecision=="medium") {
            m_bkg->MCMCSetPrecision(BCEngineMCMC::kMedium);
        }
        else if (fprecision=="high") {
            m_bkg->MCMCSetPrecision(BCEngineMCMC::kHigh);
        }
        else {
            std::cout << "##### unknown precision " << fprecision << std::endl;
            return 0;
        }

        
        std::cout << "#### running actual fit ..." << std::endl;

        // marginalize
        m_bkg->MarginalizeAll();
        
        // find global mode
        m_bkg->FindMode( m_bkg->GetBestFitParameters() );
        
        // calculate normalization
        m_bkg->Integrate();
        
        // get p-value
        m_bkg->CalculatePValue( m_bkg->GetBestFitParameters() );
        
        // print marginalized distributions
        m_bkg->PrintAllMarginalized(results_name_bkg+"_distributions.pdf");
        
        // print results file
        m_bkg->PrintResults(results_name_bkg+"_results.txt");
        
        // set up model manager
        // create new BCTemplateFitterManager
        BCModelManager * smm = new BCModelManager();
        smm->AddModel(m_bkg);
        smm->AddModel(m);
        
        // compare models
        double BayesFactor = smm->BayesFactor(0,1);
        
        // check for positive signal and calculate activity
        if (BayesFactor<fBF_limit) {
            fActivityStr.push_back(TString::Format("%1.2e - %1.2e + %1.2e \t %1.2e",factivity[j],factivity_err_low[j],factivity_err_up[j],BayesFactor));
        }
        else {
            fActivityStr.push_back(TString::Format("< %1.2e \t %1.2e",factivity_lim[j],BayesFactor));
        }
        
        delete m_bkg;
        delete smm;
        
        
    }
    
    
    std::cout << "###### Writing summary file ..." << std::endl;
    
    //// write activities to file
    // open results file
    ofstream results_file;
    results_file.open(fresults_name+"_activities_summary.txt");
    results_file << "Process \t Activity/Flux (Bq/kg)/(Bq/m3)/(s^-1 cm^-2) \t Bayes Factor" << "\n";
    
    // loop over all processes
    for (int i=0; i<Nprocesses; ++i) {
        
        // print results to txt file
        results_file << fname_process[i] << "\t" << fActivityStr[i] << "\n";
        
    }
    results_file.close();
    
    
    
    return 1;
}


// ----------------------------------------------------
// add simulated histo Pb210
// ----------------------------------------------------


int GeMSE_bkg_fit::AddHistPb(TString name, TString file_name, double lim_low, double lim_up) {
    
    TH1D* hPb210 = GetSimHist(file_name,name);
    
    if (hPb210==0) {
        return 0;
    }
    
    double counts_Pb210 = factivity_Pb210*ft_live*fintegral.back();
    double counts_Pb210_err = factivity_Pb210_err*ft_live*fintegral.back();
    
    fname_process.push_back(name);
    fhist_process->Add(hPb210);
    fprior_name.push_back("Gauss");
    fprior_mean.push_back(counts_Pb210);
    fprior_sigma.push_back(counts_Pb210_err);
    
    return 1;
    
}

// ----------------------------------------------------
// add simulated muons histo
// ----------------------------------------------------


int GeMSE_bkg_fit::AddHistMuons(TString name, TString file_name, double lim_low, double lim_up, double prior_sigma) {
    
    TH1D* hmuons = GetSimHist(file_name,name);
    
    if (hmuons==0) {
        return 0;
    }
    
    fcounts_muons = ft_live*fintegral.back();
    
    fname_process.push_back(name);
    fhist_process->Add(hmuons);
    fprior_name.push_back("Gauss");
    fprior_mean.push_back(fcounts_muons);
    fprior_sigma.push_back(prior_sigma*fcounts_muons);
    
    return 1;
    
}


// ----------------------------------------------------
// add simulated histo
// ----------------------------------------------------


int GeMSE_bkg_fit::AddHistSim(TString name, TString file_name, double lim_low, double lim_up) {
    
    fname_process.push_back(name);
    
    TH1D* hist = GetSimHist(file_name,name);
    
    if (hist==0) {
        return 0;
    }
    
    fhist_process->Add(hist);
    fprior_name.push_back("Const");
    fprior_mean.push_back(0);
    fprior_sigma.push_back(0);
    
    return 1;
    
}

// ----------------------------------------------------
// get measured histo
// ----------------------------------------------------


int GeMSE_bkg_fit::GetMeasuredHist(TString FileName) {
    
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
    
    if (!File->GetListOfKeys()->Contains("t_live")) {
        std::cout << "###### ERROR: no live time in file " << FileName << std::endl;
        return 0;
    }
    
    TVectorD* v_live = (TVectorD*) File->Get("t_live");
    ft_live = ((*v_live)[0]);
    
    TH1D* hist = (TH1D*) File->Get("hist");
    hist->SetDirectory(0);
    File->Close();
    
    TH1D* hist_range = new TH1D("measurement","",fnBins,&fxBins.at(0));
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
    
    
    fhist_meas = (*hist_range);
    
    return 1;
    
}

// ----------------------------------------------------
// get simulated histo
// ----------------------------------------------------


TH1D* GeMSE_bkg_fit::GetSimHist(TString FileName, TString name) {
    
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
    
    TH1D* hist_range = new TH1D(name,"",fnBins,&fxBins.at(0));
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
    
    fintegral.push_back(hist_range->Integral());
    
    return hist_range;
}

// ----------------------------------------------------
// read global parameters from user specified file
// ----------------------------------------------------


int GeMSE_bkg_fit::ReadPar(TString FileName) {
    
    
    // clear all vectors
    fparlimit_low.clear();
    fparlimit_up.clear();
    fCuShielding_isotope.clear();
    fCryostat_isotope.clear();
    fGeCrystal_isotope.clear();
    fCuShielding_name.clear();
    fCryostat_name.clear();
    fGeCrystal_name.clear();
    
    ifstream File;
    File.open(FileName);
    
    if (!File.is_open()) {
        std::cout << "##### ERROR: could not open " << FileName << std::endl;
        return 0;
    }
    
    std::string headerline;
    TString spectrum_name;
    TString isotope_name;
    double parlimit_low, parlimit_up;
    
    getline(File, headerline);
    File >> fresults_name;
    getline(File, headerline);
    getline(File, headerline);
    File >> fprecision;
    getline(File, headerline);
    getline(File, headerline);
    File >> fBF_limit;
    getline(File, headerline);
    getline(File, headerline);
    File >> fCL;
    getline(File, headerline);
    getline(File, headerline);
    File >> fint_const;
    getline(File, headerline);
    getline(File, headerline);
    File >> fnBins;
    getline(File, headerline);
    getline(File, headerline);
    File >> fBinSlope;
    getline(File, headerline);
    getline(File, headerline);
    File >> frange_low >> frange_high;
    getline(File, headerline);
    getline(File, headerline);
    File >> fmuon_flux >> fmuon_flux_err;
    getline(File, headerline);
    getline(File, headerline);
    File >> factivity_Pb210 >> factivity_Pb210_err;
    getline(File, headerline);
    getline(File, headerline);
    File >> fMeas_name;
    getline(File, headerline);
    getline(File, headerline);
    File >> fSim_folder;
    getline(File, headerline);
    getline(File, headerline);
    File >> fMuon_name >> parlimit_low >> parlimit_up;
    fparlimit_low.push_back(parlimit_low);
    fparlimit_up.push_back(parlimit_up);
    getline(File, headerline);
    getline(File, headerline);
    File >> fPb210_name >> parlimit_low >> parlimit_up;
    fparlimit_low.push_back(parlimit_low);
    fparlimit_up.push_back(parlimit_up);
    getline(File, headerline);
    getline(File, headerline);
    File >> fRn_name >> parlimit_low >> parlimit_up;
    fparlimit_low.push_back(parlimit_low);
    fparlimit_up.push_back(parlimit_up);
    getline(File, headerline);
    getline(File, headerline);
    while (true)
    {
        getline(File, headerline);
        std::stringstream ss(headerline);
        ss >> isotope_name >> spectrum_name >> parlimit_low >> parlimit_up;
        
        if(isotope_name=="#") break;
        
        fCuShielding_isotope.push_back(isotope_name);
        fCuShielding_name.push_back(spectrum_name);
        fparlimit_low.push_back(parlimit_low);
        fparlimit_up.push_back(parlimit_up);
        
    }
    
    while (true)
    {
        getline(File, headerline);
        std::stringstream ss(headerline);
        ss >> isotope_name >> spectrum_name >> parlimit_low >> parlimit_up;
        
        if(isotope_name=="#") break;
        
        fCryostat_isotope.push_back(isotope_name);
        fCryostat_name.push_back(spectrum_name);
        fparlimit_low.push_back(parlimit_low);
        fparlimit_up.push_back(parlimit_up);
        
    }

    while (!File.eof())
    {
        getline(File, headerline);
        std::stringstream ss(headerline);
        ss >> isotope_name >> spectrum_name >> parlimit_low >> parlimit_up;
        
        fGeCrystal_isotope.push_back(isotope_name);
        fGeCrystal_name.push_back(spectrum_name);
        fparlimit_low.push_back(parlimit_low);
        fparlimit_up.push_back(parlimit_up);
        
    }
    
    File.close();
    
    // print information
    std::cout << "######################################" << std::endl;
    std::cout << "#### Reading "<< FileName << " ..." << std::endl;
    std::cout << "results name: " << fresults_name << std::endl;
    std::cout << "accuracy of MCMC: " << fprecision << std::endl;
    std::cout << "threshold for signal detection: " << fBF_limit << std::endl;
    std::cout << "CL for limit on activity: " << fCL << std::endl;
    std::cout << "integration constant: " << fint_const << std::endl;
    std::cout << "number of bins: " << fnBins << std::endl;
    std::cout << "bin slope: " << fBinSlope << std::endl;
    std::cout << "fit range (keV): " << frange_low << " - " << frange_high << std::endl;
    std::cout << "muon flux (cm^-2 s-^1): " << fmuon_flux << std::endl;
    std::cout << "muon flux uncertainty: " << fmuon_flux_err << std::endl;
    std::cout << "Pb210 activity (Bq/kg): " << factivity_Pb210 << " +- " << factivity_Pb210_err << std::endl;
    std::cout << "measured spectrum: " << fMeas_name << std::endl;
    std::cout << "folder simulated spectra: " << fSim_folder << std::endl;
    std::cout << "muon spectrum: " << fMuon_name << std::endl;
    std::cout << "Pb210 spectrum: " << fPb210_name << std::endl;
    std::cout << "Rn spectrum: " << fRn_name << std::endl;
    std::cout << "Cu shielding spectra:" << std::endl;
    for (int i=0; i<fCuShielding_name.size(); ++i)
    {
        std::cout << fCuShielding_name[i] << std::endl;
    }
    std::cout << "cryostat spectra:" << std::endl;
    for (int i=0; i<fCryostat_name.size(); ++i)
    {
        std::cout << fCryostat_name[i] << std::endl;
    }
    std::cout << "Ge crystal spectra:" << std::endl;
    for (int i=0; i<fGeCrystal_name.size(); ++i)
    {
        std::cout << fGeCrystal_name[i] << std::endl;
    }


    std::cout << "######################################" << std::endl;
    
    fParSet=true;
    
    return 1;
}


// ----------------------------------------------------
// make nice stack plot
// ----------------------------------------------------


int GeMSE_bkg_fit::plot_stack()
{
    
    TFile* file_stack = TFile::Open(fresults_name+"_All_stack.root");
    
    if (!file_stack) {
        
        std::cout << "###### ERROR: could not open " << fresults_name+"_All_stack.root" << std::endl;
        return 0;
    }
    
    if (!file_stack->GetListOfKeys()->Contains("c1")) {
        std::cout << "###### ERROR: no canvas named c1 in file " << fresults_name+"_All_stack.root" << std::endl;
        return 0;
    }
    
    TCanvas* c = (TCanvas*) file_stack->Get("c1");
    
    TH1D* hmeas = (TH1D*) c->FindObject("measurement");
    
    THStack* stack = (THStack*) c->GetPrimitive("");
    TList* list = stack->GetHists();
    
    TH1D* hMuons = (TH1D*) list->FindObject("Muons");
    TH1D* hRn = (TH1D*) list->FindObject("SampleCavity_Rn222");
    TH1D* hPb210 = (TH1D*) list->FindObject("innerPbShielding_Pb210");
    
    // combine CuShielding
    THStack* stack_CuShielding = new THStack("stack_CuShielding",";Energy (keV);Counts");
    
    for (int i=0; i<fCuShielding_name.size(); ++i) {

        TH1D* hCuShielding_isotope = (TH1D*) list->FindObject("CuShielding_"+fCuShielding_isotope[i]);
        stack_CuShielding->Add(hCuShielding_isotope);
        //delete hCuShielding_isotope;
    }
    TH1D* hCuShielding = (TH1D*) stack_CuShielding->GetStack()->Last();


    // combine cryostat
    THStack* stack_Cryostat = new THStack("stack_Cryostat",";Energy (keV);Counts");
    
    for (int i=0; i<fCryostat_name.size(); ++i) {
        
        TH1D* hCryostat_isotope = (TH1D*) list->FindObject("Cryostat_"+fCryostat_isotope[i]);
        stack_Cryostat->Add(hCryostat_isotope);
        //delete hCryostat_isotope;
    }
    TH1D* hCryostat = (TH1D*) stack_Cryostat->GetStack()->Last();

    // combine Ge crystal
    THStack* stack_Ge = new THStack("stack_Ge",";Energy (keV);Counts");
    
    for (int i=0; i<fGeCrystal_name.size(); ++i) {
        
        TH1D* hGeCrystal_isotope = (TH1D*) list->FindObject("Ge_"+fGeCrystal_isotope[i]);
        stack_Ge->Add(hGeCrystal_isotope);
        //delete hGeCrystal_isotope;
    }
    TH1D* hGe = (TH1D*) stack_Ge->GetStack()->Last();
    

    // stack all histos
    THStack* stack_all = new THStack("stack_all",";Energy (keV);Counts");
    stack_all->Add(hMuons);
    stack_all->Add(hRn);
    stack_all->Add(hPb210);
    stack_all->Add(hCuShielding);
    stack_all->Add(hCryostat);
    stack_all->Add(hGe);
    
    TH1D* hsum = (TH1D*) stack_all->GetStack()->Last();
    
    
    // scale all histos
    hmeas->Sumw2();
    hmeas->Scale(1./ft_live*86400.,"width");
    hsum->Scale(1./ft_live*86400.,"width");
    hMuons->Scale(1./ft_live*86400.,"width");
    hRn->Scale(1./ft_live*86400.,"width");
    hPb210->Scale(1./ft_live*86400.,"width");
    hCuShielding->Scale(1./ft_live*86400.,"width");
    hCryostat->Scale(1./ft_live*86400.,"width");
    hGe->Scale(1./ft_live*86400.,"width");
    
    //hmeas->GetXaxis()->SetRangeUser(100.,500.);
    
    // Draw
    TCanvas* c2 = new TCanvas("c2");
    gStyle->SetOptStat(0);
    c2->SetLogy();
    
    hmeas->SetMinimum(1.e-4);
    hmeas->GetXaxis()->SetTitle("Energy (keV)");
    hmeas->GetYaxis()->SetTitle("Counts (keV^{-1} day^{-1})");
    
    hmeas->SetLineColor(1);
    hmeas->SetMarkerColor(1);
    hsum->SetLineColor(2);
    hsum->SetFillStyle(0);
    hMuons->SetLineColor(3);
    hMuons->SetFillStyle(0);
    hRn->SetLineColor(4);
    hRn->SetFillStyle(0);
    hPb210->SetLineColor(5);
    hPb210->SetFillStyle(0);
    hCuShielding->SetLineColor(6);
    hCuShielding->SetFillStyle(0);
    hCryostat->SetLineColor(7);
    hCryostat->SetFillStyle(0);
    hGe->SetLineColor(8);
    hGe->SetFillStyle(0);
    
    hmeas->Draw("pe");
    stack_all->Draw("nostacksame");
    hsum->Draw("same");
    
    TLegend* leg1 = new TLegend(0.7,0.5,0.9,0.9);
    leg1->AddEntry(hmeas,"Data","p");
    leg1->AddEntry(hsum,"Simulation","l");
    leg1->AddEntry(hMuons,"Muons","l");
    leg1->AddEntry(hRn,"Rn222","l");
    leg1->AddEntry(hPb210,"Pb210","l");
    leg1->AddEntry(hCuShielding,"CuShielding","l");
    leg1->AddEntry(hCryostat,"Cryostat","l");
    leg1->AddEntry(hGe,"Ge","l");
    
    leg1->Draw();
    
    c2->SaveAs(fresults_name+"_All_stack_nice.pdf");
    c2->SaveAs(fresults_name+"_All_stack_nice.root");
    
    return 1;
    
}

// ----------------------------------------------------
// Set new parameter limits (eliminate empty bins)
// ----------------------------------------------------

void GeMSE_bkg_fit::SetNewLimits(BCMTF_HPGe* m, TString parameter) {
    
    TH1D* hist = m->GetMarginalized(parameter)->GetHistogram();
    
    int nBins = hist->GetNbinsX();
    
    // find first non-empty bin
    int firstBin = GetEmptyBin(hist, "first");
    
    // find last non-empty bin
    int lastBin = GetEmptyBin(hist, "last");
    
    // set new parameter limits
    double limit_low, limit_up;
    
    if (firstBin>1) {
        limit_low = hist->GetBinLowEdge(firstBin-1);
    }
    else {
        limit_low = hist->GetBinLowEdge(1);
    }
    if (lastBin<nBins-1) {
        limit_up = hist->GetBinLowEdge(lastBin+2);
    }
    else {
        limit_up = hist->GetBinLowEdge(lastBin)+hist->GetBinWidth(0);
    }
    
    m->GetParameter(parameter)->SetLimits(limit_low,limit_up);
    
    
}

// get first or last empty bin
int GeMSE_bkg_fit::GetEmptyBin(TH1D* hist, TString option) {
    
    int nBins = hist->GetNbinsX();
    int bin;
    double counts;
    
    if (hist->GetEntries()==0) {
        std::cout << "###### ERROR in GetEmptyBin(): histogram is empty" << std::endl;
        return 0;
    }
    
    if (option=="first") {
        bin = 1;
        for (int i=1; i<=nBins; ++i) {
            
            counts = hist->GetBinContent(i);
            
            if (counts>0) {
                bin = i;
                break;
            }
        }
    }
    
    else if (option=="last") {
        bin = nBins;
        for (int i=nBins; i>0; --i) {
            
            counts = hist->GetBinContent(i);
            
            if (counts>0) {
                bin = i;
                break;
            }
        }
    }
    else {
        std::cout << "###### ERROR in GetEmptyBin(): unknown option" << std::endl;
        return 0;
    }
    
    return bin;
    
}



