/*
 *  W180.c
 *  
 *
 *  Created by sivers on 20.10.13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */
#include "GeMSE_bkg_fit.h"
#include <iostream>
using namespace std;

int main(int argc, char *argv[]) {
    
    // argc should be 2 for correct execution
    if ( argc != 2 ) {
        
        // We print argv[0] assuming it is the program name
        std::cout<<"usage: "<< argv[0] <<" <parameters.txt>" << std::endl;
        return 0;
    }

    TString FileName = argv[1];
    
    // create a new fitter
    GeMSE_bkg_fit* fit = new GeMSE_bkg_fit();
    
    // read parameters from file
    fit->ReadPar(FileName);
    
    // run the fit
    fit->RunFit();
    
    return 1;
}






