//c++ vissani.cpp -o vissani `root-config --glibs --cflags` 

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>
#include <ctime>
#include <time.h>
#include <vector>
#include <random>

#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TMarker.h"
#include "TEllipse.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TApplication.h"

using namespace std;

double f_rng (double min, double max) {

    double val = ( (double)rand() / RAND_MAX ) * ( max - min ) + min;
    return val;
}

double diff_spec (double *E_nu , double *par) {

    double N_targets = 4.55e32;
    double cross_sec = 7.5e-44 * E_nu[0];
    double sn_energy = 5e52;
    double sn_energy_mev = sn_energy * 6.24e5;
    double distance = 50;
    double distance_cm = distance * 3.086e21;
    double Temp = 4;

    double val = par[0]*N_targets * cross_sec * ( sn_energy_mev * pow(E_nu[0],2) * exp(-E_nu[0]/Temp)) / (4. * M_PI * pow(distance_cm,2) * 6 * pow(Temp,4) );
    //double val = par[0]*E_nu[0]*2;
    return val;
}

//repeated because i dont know how to take values from the function above -> it's with function->eval(XX)
double diff_spec_test (double E_nu, double norm, double Temp, double SN_E, double N_targ, double cs_0, double d) {

    double cross_sec = cs_0 * E_nu;

    double val = norm * N_targ * cross_sec * ( SN_E * pow(E_nu,2) * exp(-E_nu/Temp)) / (4. * M_PI * pow(d,2) * 6 * pow(Temp,4) );
    return val;
}

double chi_square_calc (vector<double> list_rgn, double Temp, double SN_E, double N_targ, double cs_0, double d) {

    double summation = 0;
    for ( int k = 0; k < list_rgn.size(); k++){

        summation +=  log ( N_targ * cs_0 * SN_E * pow(list_rgn[k],4) * exp ( - list_rgn[k] / Temp ) / (4*M_PI * pow(d,2) * 6 * pow(Temp,4) ) );     
    }

    double val = ( 2 * (N_targ * cs_0 * SN_E * Temp / M_PI / pow(d,2)) ) - 2 * summation;
    return val;
}




int main() {
    TApplication *myapp=new TApplication("myapp",0,0);

    clock_t c_start0 = clock();

    ///////////// Constants //////////////////////////
    double temperature = 4; //MeV
    double sn_energy_erg = 5e52; //erg
    double sn_energy_mev = sn_energy_erg * 6.24e5; // MeV
    double N_targets = 4.55e32;
    double cross_sec_0 = 7.5e-44; //cm^2/MeV^-2
    double distance = 50; // kpc
    double distance_cm = distance * 3.086e21; //cm



    ///////////////////// Function /////////////////////
    double xMin = 0.;
    double xMax = 65.;

    TF1 *spectrum = new TF1 ( "spectrum" , diff_spec, xMin , xMax, 1);
    spectrum->SetParameters(1,1);

    //cout << spectrum->Eval(10) << endl; // to obtain value from the function

    double maximum = spectrum->GetMaximum();
    cout << maximum << endl;
    //spectrum->SetParameters(1./maximum,1);

    ////////////////////// Random points ////////////////////

    TH1F * hRandDist = new TH1F("","",200,0,65);
    TGraphErrors *gRandom = new TGraphErrors();
    TGraphErrors *gExperimental = new TGraphErrors();

    vector<double> exp_data = {19, 22, 28, 36, 36, 37, 38, 39};
    vector<double> exp_dataerror = {5, 5, 6, 6, 9, 7, 7, 7};

    for ( int i = 0; i < exp_data.size(); i++){
        gExperimental->SetPoint(i,i+1,exp_data[i]);
        gExperimental->SetPointError(i,0,exp_dataerror[i]);
    } 

    TLine ** lPoints;
    lPoints = new TLine*[ 10 ];

    vector<double> rngs;
    double x;
    double y;

    int counter=0;
    int counter_N=0;
    while (rngs.size() < 10){ 
        counter_N ++;
        x = f_rng ( 20 , xMax);
        y = f_rng ( 0 , maximum);

        if ( y < diff_spec_test (x , 1., temperature, sn_energy_mev, N_targets, cross_sec_0, distance_cm) ){
            
            rngs.push_back( x );
            lPoints[ counter ] = new TLine( x , 0. , x , maximum);
            hRandDist->Fill(x);
            counter++;
        }
    }

    cout << "counter: " << counter_N << endl;

    for ( int k = 0; k < 1e6; k++){    
    
        x = f_rng ( 0 , xMax);
        y = f_rng ( 0 , maximum);

        if ( y < diff_spec_test (x , 1., temperature, sn_energy_mev, N_targets, cross_sec_0, distance_cm) ) hRandDist->Fill(x);
    }

    sort ( rngs.begin(), rngs.end());
    for (int  i=0; i<rngs.size();i++){
        gRandom -> SetPoint(i, i+1, rngs[i]);
        gRandom -> SetPointError(i, 0, 6);
        //gRandom -> SetPointError( i, 0, sqrt(counter_N));
        cout << i<<": "<< rngs[i] << endl;
    }


    /////////////////// Integration ////////////////////////////
    double N = 1e6;
    double middle_point;
    double delta_x = (xMax - xMin) / N;
    double threshold = 20.;
    
    double integral_full = 0;
    double integral_partial = 0;

    for (int step = 0; step < N ; step++){

        middle_point = delta_x*step + delta_x/2;

        integral_full += diff_spec_test ( middle_point, 1., temperature ,sn_energy_mev, N_targets, cross_sec_0, distance_cm) * delta_x;

        if ( middle_point > threshold){

            integral_partial += diff_spec_test ( middle_point, 1., temperature, sn_energy_mev, N_targets, cross_sec_0, distance_cm) * delta_x;
        }
    }
    cout << "Total integral = " << integral_full << endl;
    cout << "partial integral = " << integral_partial << endl;

    cout << "Percentage above 20 MeV: " << (double) integral_partial/integral_full * 100. << endl;


    ///////////////////// Chi-square analysis //////////

    double chi_square = chi_square_calc(rngs, temperature, sn_energy_mev, N_targets, cross_sec_0, distance_cm);
    cout << "chi_square: " << chi_square << endl;

    ///////// minimum /////

    double E_sum = 0;
    for ( int t = 0; t < rngs.size(); t++){

        E_sum += rngs[t];    
    }
    
    double temp_min = 1. / ( 5 * rngs.size() ) * E_sum;
    double sn_e_min = ( 5 * pow(rngs.size(),2) * M_PI * pow(distance_cm,2) ) / ( N_targets * cross_sec_0 * E_sum ); // in mev
    double sn_e_min_erg = sn_e_min * 1.6e-6; // in ergs

    cout << " T min: " << temp_min << endl;
    cout << " SN e min: " << sn_e_min_erg << endl;

    TMarker *minimum_point = new TMarker(temp_min,sn_e_min_erg / 1.e52,0);

    ///////// allowed regions ////

    double N_region = 1e3;
    double bin_t_min = 3;
    double bin_t_max = 8.5;
    double dt = (bin_t_max - bin_t_min) / N_region;
    double bin_energy_min = 0.8e57;
    double bin_energy_max = 10e57;
    double dE = (bin_energy_max - bin_energy_min) / N_region;

    double sweep_temp;
    double sweep_energy;

    double diff_chi_square;

    TGraph ** gAllowedRegions;
    gAllowedRegions = new TGraph*[3];

    vector<double> CL = { 0.95, 0.90, 0.683};
    double counter_region;
    for ( int j = 0; j < CL.size(); j++){

        counter_region = 0;
        gAllowedRegions[j] = new TGraph();

        for ( int i = 0; i < N_region; i++ ){

            sweep_temp = bin_t_min + i * dt;

            for ( int k = 0; k < N_region; k++ ){

                sweep_energy = bin_energy_min + k * dE;

                diff_chi_square= chi_square_calc ( rngs, sweep_temp, sweep_energy, N_targets,cross_sec_0,distance_cm ) 
                    - chi_square_calc ( rngs, temp_min, sn_e_min, N_targets,cross_sec_0,distance_cm );

                if( diff_chi_square < (-2 * log ( 1 - CL[j] ) ) ){

                    gAllowedRegions[j]->SetPoint(counter_region,sweep_temp,sweep_energy*1.6e-6 / 1.e52);
                    counter_region++;
                } 
            }
        }
    }
    


    //////////////////// Plots /////////////////////
    clock_t c_end0 = clock();
    cout << "\nTime needed: " << (double) (c_end0-c_start0) / CLOCKS_PER_SEC << " seconds." << endl;


    TCanvas * c1 = new TCanvas("c1","c1",1000,700);
	c1->cd();
    c1->SetGridx();
    c1->SetGridy();
    spectrum->SetLineColor(kAzure-5);
    spectrum->SetLineWidth(4);
    spectrum->SetTitle("");
    spectrum->GetYaxis()->SetTitle("dN/dE [MeV^{-1}]");
    spectrum->GetYaxis()->SetTitleOffset(1.0);
    spectrum->GetYaxis()->SetTitleSize(0.045);
    spectrum->GetXaxis()->SetTitleOffset(1.0);
    spectrum->GetXaxis()->SetTitleSize(0.045);
    spectrum->GetXaxis()->SetTitle("Energy [MeV]");
    spectrum->Draw();
    c1 -> SaveAs("Differential_spectrum.pdf"); 

    TCanvas * c2 = new TCanvas("c2","c2",1000,700);
	c2->cd();
    c2->SetGridx();
    c2->SetGridy();
    spectrum->SetTitle("Random points generated");
    spectrum->Draw();
    for ( int k = 0; k < 10; k++){
        lPoints[k]->SetLineColor(kGray+2);
        //lPoints[k]->SetLineStyle(2);
        lPoints[k]->SetLineWidth(2);
        lPoints[k]->Draw("same");
    }
    c2 -> SaveAs("distr_rgns.pdf"); 

    TCanvas * c3 = new TCanvas("c3","c3",1000,700);
	c3->cd();
    c3->SetGridx();
    c3->SetGridy();
    gRandom->SetMarkerColor(kAzure-5);
    gRandom->SetMarkerStyle(20);
    gRandom->SetLineColor(kAzure-5);
    gRandom->SetLineWidth(2);
    gRandom->SetTitle("");
    gRandom->GetYaxis()->SetTitle("Energy [MeV]");
    gRandom->GetYaxis()->SetTitleOffset(1.0);
    gRandom->GetYaxis()->SetTitleSize(0.045);
    gRandom->GetYaxis()->SetRangeUser(0.,60.);
    gRandom->GetXaxis()->SetTitleOffset(1.0);
    gRandom->GetXaxis()->SetTitleSize(0.045);
    gRandom->Draw("ap");
    gExperimental->SetMarkerColor(kGray+2);
    gExperimental->SetMarkerStyle(20);
    gExperimental->SetLineColor(kGray+2);
    gExperimental->SetLineWidth(2);
    gExperimental->Draw("same p");

    auto legend = new TLegend(0.6,0.14,0.88,0.4);
    legend->AddEntry(gRandom,"Monte Carlo","pl");
    legend->AddEntry(gExperimental,"Experimental data","pl");
    legend->Draw(); 
    c3->SaveAs("Comparison_of_data.pdf");

    TCanvas * c4 = new TCanvas("c4","c4",1000,700);
    c4->cd();
    hRandDist->SetLineColor(kAzure-5);
    hRandDist->SetTitle("Distribution of random points generated");
    hRandDist->SetLineWidth(2);
    hRandDist->SetFillStyle(3003);
    hRandDist->SetFillColor(kAzure+5);
    hRandDist->GetYaxis()->SetTitle("Counts");
    hRandDist->GetYaxis()->SetTitleOffset(1.0);
    hRandDist->GetYaxis()->SetTitleSize(0.045);
    hRandDist->GetXaxis()->SetTitleOffset(1.0);
    hRandDist->GetXaxis()->SetTitleSize(0.045);
    hRandDist->GetXaxis()->SetTitle("Energy [MeV]");
    hRandDist->Draw();  
    c4->SaveAs("Distribution_of_random_points_generated.pdf");

    TCanvas * c5 = new TCanvas("c5","c5",1000,700);
    c5->cd();
    gAllowedRegions[0]->SetTitle("");
    gAllowedRegions[0]->GetYaxis()->SetTitle("#varepsilon [10^{52} erg]");
    gAllowedRegions[0]->GetYaxis()->SetTitleOffset(1);
    gAllowedRegions[0]->GetYaxis()->SetRangeUser(0,2);
    gAllowedRegions[0]->GetYaxis()->SetTitleSize(0.045);
    gAllowedRegions[0]->GetXaxis()->SetTitleOffset(1.0);
    gAllowedRegions[0]->GetXaxis()->SetTitleSize(0.045);
    gAllowedRegions[0]->GetXaxis()->SetLimits(bin_t_min,bin_t_max);
    gAllowedRegions[0]->GetXaxis()->SetTitle("Temperature [MeV]");

    gAllowedRegions[0]->SetMarkerColor(kAzure-9);
    gAllowedRegions[0]->SetMarkerStyle(7);
    //gAllowedRegions[0]->SetMarkerSize(1.5);
    gAllowedRegions[0]->Draw("AP");  

    gAllowedRegions[1]->SetMarkerColor(kAzure-8);
    gAllowedRegions[1]->SetMarkerStyle(7);
    //gAllowedRegions[1]->SetMarkerSize(3);
    gAllowedRegions[1]->Draw("sameP");  

    gAllowedRegions[2]->SetMarkerColor(kAzure-7);
    gAllowedRegions[2]->SetMarkerStyle(7);
    //gAllowedRegions[2]->SetMarkerSize(1.5);
    gAllowedRegions[2]->Draw("sameP");  

    auto legend1 = new TLegend(0.65,0.65,0.88,0.88);
    gAllowedRegions[0]->SetFillColor(kAzure-9);
    gAllowedRegions[1]->SetFillColor(kAzure-8);
    gAllowedRegions[2]->SetFillColor(kAzure-7);
    legend1->AddEntry(gAllowedRegions[0],"3 #sigma","f");
    legend1->AddEntry(gAllowedRegions[1],"2 #sigma","f");
    legend1->AddEntry(gAllowedRegions[2],"1 #sigma","f");
    legend1->Draw("");

    minimum_point->SetMarkerColor(kWhite);
    minimum_point->SetMarkerStyle(20);
    minimum_point->SetMarkerSize(1.5);
    minimum_point->Draw("same");

    c5->SaveAs("Allowed_regions.pdf");

    // TF1 * fit = new TF1("fit", "[0]*exp(-[1]*x)", xMin,xMax);
    // fit->SetParName ( 0, "Normalization" );
    // fit->SetParName ( 1, "#mu" );
    // h_s_distr -> Fit ( fit );

    // fit->Draw("same");
    

    myapp->Run();
    return 0;
}

