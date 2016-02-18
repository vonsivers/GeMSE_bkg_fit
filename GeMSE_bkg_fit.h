
#ifndef __GeMSE_bkg_fit__H
#define __GeMSE_bkg_fit__H

// ---------------------------------------------------------

#include <TH1D.h>
#include <TString.h>
#include <TObjArray.h>

#include <BAT/BCMTF.h>


// ---------------------------------------------------------
class GeMSE_bkg_fit
{
    
public:
    
    /**
     * The default constructor. */
    GeMSE_bkg_fit();
    
    /**
     * The default destructor. */
    ~GeMSE_bkg_fit();
    
    /**
     * Read parameters from file
     */
    int ReadPar(TString FileName);
    
    /**
     * Run the fit
     */
    int RunFit();


private:
    
    // add Pb210 histo
    int AddHistPb(TString name, TString file_name, double lim_low, double lim_up);
    
    // add muon histo
    int AddHistMuons(TString name, TString file_name, double lim_low, double lim_up, double prior_sigma);
    
    // add simulated histo
    int AddHistSim(TString name, TString file_name, double lim_low, double lim_up);
    
    // get measured histo
    int GetMeasuredHist(TString FileName);
    
    // get simulated histo
    TH1D* GetSimHist(TString FileName, TString name);
    
    // plot stack
    int plot_stack();

    // set new parameter limits (remove empty bins)
    void SetNewLimits(BCMTF* m, TString parameter);
    
    // get first/last empty bin
    int GetEmptyBin(TH1D* hist, TString option);
    
    // check if parameters have been set
    bool fParSet;
    
    // simulated spectra
    TObjArray *fhist_process = new TObjArray();
    
    // name of simulated spectra
    std::vector<TString> fname_process;
    
    // fit parameter limits of simulated spectra
    std::vector<double> fparlimit_low;
    std::vector<double> fparlimit_up;
    
    // priors of simulated spectra
    std::vector<TString> fprior_name;
    std::vector<double> fprior_mean;
    std::vector<double> fprior_sigma;
    
    // folder siwth simulated spectra
    TString fSim_folder;
    
    // name of results file
    TString fresults_name;
    
    // name of measurement file
    TString fMeas_name;
    
    // name of muon file
    TString fMuon_name;
    
    // name of Pb210 file
    TString fPb210_name;
    
    // name of Rn222 file
    TString fRn_name;
    
    // name of CuShielding files/isotopes
    std::vector<TString> fCuShielding_name;
    std::vector<TString> fCuShielding_isotope;
    
    // name of Cryostat files/isotopes
    std::vector<TString> fCryostat_name;
    std::vector<TString> fCryostat_isotope;

    // name of Ge files/isotopes
    std::vector<TString> fGeCrystal_name;
    std::vector<TString> fGeCrystal_isotope;

    // measured spectrum
    TH1D fhist_meas;
    
    // precision of MCMC
    TString fprecision;
    
    // CL of limit
    double fCL;
    
    // threshold of BF
    double fBF_limit;
    
    // measurement time
    double ft_live;
    
    // integral of simulated histo
    std::vector<double> fintegral;
    
    // number of bins
    int fnBins;
    
    // fit range (keV)
    double frange_low;
    double frange_high;
    
    // increase of bins
    double fBinSlope;
    
    // lower edge of bins
    //double *fxBins = new double;
    std::vector<double> fxBins;
    
    // simulated integral muon flux (cm-2 s-1)
    double fmuon_flux;
    
    // uncertainty of muon flux
    double fmuon_flux_err;
    
    // counts in simulated muon spectrum
    double fcounts_muons;
    
    // Pb210 acitivty (Bq/kg)
    double factivity_Pb210;
    double factivity_Pb210_err;
    
    // activity, limits, errors
    std::vector<double> factivity;
    std::vector<double> factivity_err_low;
    std::vector<double> factivity_err_up;
    std::vector<double> factivity_lim;
    
    // activity string for summary file
    std::vector<TString> fActivityStr;


};
// ---------------------------------------------------------

#endif


