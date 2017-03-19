//###################################################
//
// Create a latex table for glauber paper. 
// You will also need CentralityHelperFunctions.C 
//
// For any entry in the table, you need to enter in the: 
// - centrality range
// - sqrt(sNN)
// - inelastic xsec_NN
// - filename of glauber output
//
// Make sure you add them to the vectors in the right order :)
//
// run with: 
//
// ~$ root -l MakeLatexTables.C
//
// note: still need to add the systematic uncertainties... 
//
// -- J Kamin, Jan 2016
//
//###################################################


#include <vector>

static const int nCentralities = 100;
Double_t binedges[nCentralities+1] = {0};

std::vector<Int_t> Centrality_lo;
std::vector<Int_t> Centrality_hi;
std::vector<Double_t> Sqrts;
std::vector<Double_t> XSec;
std::vector<TString>  FileName;

std::vector<Double_t> NColl_lo;
std::vector<Double_t> NColl_hi;
std::vector<Double_t> NColl_mean;
std::vector<Double_t> NColl_rms;
std::vector<Double_t> NPart_lo;
std::vector<Double_t> NPart_hi;
std::vector<Double_t> NPart_mean;
std::vector<Double_t> NPart_rms;
std::vector<Double_t> b_lo;
std::vector<Double_t> b_hi;
std::vector<Double_t> b_mean;
std::vector<Double_t> b_rms;
std::vector<Double_t> TAA;
std::vector<Double_t> TAA_rms;


TH1D *h_b;
TH1D *h_ncoll;
TH1D *h_npart;
char saythis[500];
char saythis2[500];

#include "CentralityHelperFunctions.C"
int InitializeVectors();

void MakeLatexTables()
{

  TCanvas *c1 = new TCanvas("c1","c1");
  TFile *f1;
  TNtuple *myTuple;


  int nEntries = InitializeVectors();
  cout << "We have " << nEntries << " entries in our table." << endl;


  int plot_pbpb_or_ppb = 1;
  int xbins=200; float xlo=0; float xhi=30; int ybins=3000; float ylo=0.5; float yhi=3000.5; int rebinx=4; int rebiny=20;
  if(plot_pbpb_or_ppb==1){
    xbins=500; xlo=0; xhi=30; ybins=3000; ylo=0.5; yhi=3000.5; rebinx=1; rebiny=1;//rebinx=4; rebiny=20;
  }
  else if(plot_pbpb_or_ppb==2){
    xbins=1000; xlo=0; xhi=16; ybins=40; ylo=0.5; yhi=40.5; rebinx=1; rebiny=1;
  }
  int xbins_for_b = 1000*xbins;

  for(int i=0; i<nEntries; i++){

    cout << "Starting table entry " << i << endl;
    f1 = TFile::Open(FileName[i]);
    myTuple = (TNtuple*)f1->Get("nt_Pb_Pb");

    int FindNewCentralities = 0;
    if(i==0)
      FindNewCentralities = 1;
    else if(FileName[i]!=FileName[i-1])
      FindNewCentralities = 1;

    //################# Find the Centralities ######################
    c1->cd()->SetLogy();
    if(FindNewCentralities){
      sprintf(saythis,"B>>h_b(%d,%2.5f,%2.5f)",xbins_for_b,xlo,xhi);
      myTuple->Draw(saythis);
      sprintf(saythis,"h_b");
      h_b = (TH1D*)gPad->GetPrimitive(saythis);

      double HFefficiency = 1.00;
      double fullintegral_hb = h_b->Integral(-1,-1);
      double extra_chunk = fullintegral_hb - HFefficiency*fullintegral_hb;
      cout << "we're calculating with an HF eff = " << HFefficiency << "   extra_chunk = " << extra_chunk << endl;
      int lastfilledbin = h_b->GetNbinsX();
      for(int j=h_b->GetNbinsX(); j>0; j--){
        if(h_b->GetBinContent(j)!=0){
          lastfilledbin = j;
          break;
        }
      }
      if(extra_chunk>0)
        h_b->SetBinContent(lastfilledbin+1,extra_chunk);
      else{
        //h_b->SetBinContent(lastfilledbin,h_b->GetBinContent(lastfilledbin)+extra_chunk);
        //if(h_b->GetBinContent(lastfilledbin)<0){
        //  lastfilledbin--;
        //  h_b->SetBinContent(lastfilledbin,h_b->GetBinContent(lastfilledbin)+h_b->GetBinContent(lastfilledbin+1));
        //  h_b->GetBinContent(lastfilledbin+1,0);
        //}
        int killbin = h_b->GetNbinsX();
        double quickintegral = 0;
        for(int k=h_b->GetNbinsX(); k>0; k--){
          quickintegral+=h_b->GetBinContent(k);
          if(quickintegral > -1.0*extra_chunk){
            killbin = k;
            break;
          }
        }
        for(int k=killbin; k<h_b->GetNbinsX(); k++){
          h_b->SetBinContent(k,0);
          //h_b->SetBinError  (k,0);
        }
      }

      Double_t precision_quick = FindCentralities(nCentralities, binedges, h_b);
      cout << "precision:  " << precision_quick << endl;
    }
    //###############################################################


    Double_t b_lo_temp   = 0.0;
    Double_t b_hi_temp   = 100.0;
    Double_t b_mean_temp = -1.0;
    int binlo = int(float(nCentralities)/100.0)*Centrality_lo[i];
    int binhi = int(float(nCentralities)/100.0)*Centrality_hi[i];
    sprintf(saythis2,"B>%2.5f && B<%2.5f",binedges[binlo],binedges[binhi]);

    h_b->Reset();
    sprintf(saythis,"B>>h_b(%d,%2.5f,%2.5f)",xbins_for_b,xlo,xhi);
    //sprintf(saythis2,"B>%2.5f && B<%2.5f",binedges[Centrality_lo[i]],binedges[Centrality_hi[i]]);
    cout << saythis << endl;
    cout << saythis2 << endl;
    myTuple->Draw(saythis,saythis2);
    sprintf(saythis,"h_b");
    h_b = (TH1D*)gPad->GetPrimitive(saythis);
    b_mean.push_back(h_b->GetMean());
    b_rms.push_back(h_b->GetStdDev());
    b_lo.push_back(binedges[Centrality_lo[i]]);
    b_hi.push_back(binedges[Centrality_hi[i]]);


    sprintf(saythis,"Ncoll>>h_ncoll(%d,%2.5f,%2.5f)",ybins,ylo,yhi);
    //sprintf(saythis2,"B>%2.5f && B<%2.5f",binedges[Centrality_lo[i]],binedges[Centrality_hi[i]]);
    //cout << saythis << endl;
    //cout << saythis2 << endl;
    myTuple->Draw(saythis,saythis2);
    sprintf(saythis,"h_ncoll");
    h_ncoll = (TH1D*)gPad->GetPrimitive(saythis);
    NColl_mean.push_back(h_ncoll->GetMean());
    NColl_rms.push_back(h_ncoll->GetStdDev());


    sprintf(saythis,"Npart>>h_npart(%d,%2.5f,%2.5f)",ybins,ylo,yhi);
    //sprintf(saythis2,"B>%2.5f && B<%2.5f",binedges[Centrality_lo[i]],binedges[Centrality_hi[i]]);
    //cout << saythis << endl;
    //cout << saythis2 << endl;
    myTuple->Draw(saythis,saythis2);
    sprintf(saythis,"h_npart");
    h_npart = (TH1D*)gPad->GetPrimitive(saythis);
    NPart_mean.push_back(h_npart->GetMean());
    NPart_rms.push_back(h_npart->GetStdDev());

    TAA.push_back(NColl_mean[i]/XSec[i]);
    TAA_rms.push_back(NColl_rms[i]/XSec[i]);

    //################# Make the Histos #######################
    //sprintf(saythis,"Ncoll_bcent_Nom");
    //NColl_bcent_Nom = new TH1D(saythis,saythis,nCentralities,binedges);
    //sprintf(saythis,"Ncoll_cent_Nom");
    //NColl_cent_Nom  = new TH1D(saythis,saythis,nCentralities,0,100);
    //sprintf(saythis,"Ncoll"); sprintf(saythis2,"");
    //MakeCentHistos(myTuple,saythis,saythis2,nCentralities,binedges,NColl_bcent_Nom,NColl_cent_Nom);
    //###############################################################

  }

  cout << endl << endl;
  cout << "\\begin{table*}[t]" << endl;
  cout << "\\begin{tabular}{ccccccccc}\\hline" << endl;

  cout << "Centrality & $b_{min}$ (fm) & $b_{max}$ (fm) & $\\langle N_{Coll} \\rangle$ & $N_{Coll}$ RMS & $\\langle N_{Part} \\rangle$ & $N_{Part}$ RMS & $\\langle T_{AA} \\rangle$ (1/mb) & $T_{AA}$ RMS (1/mb) \\\\ \\hline" << endl;

  for(int i=0; i<Centrality_lo.size(); i++){
  
    cout << std::setprecision(4) << Centrality_lo[i] << "-" << Centrality_hi[i] << "\\\%  & "
      << b_lo[i]       << " & " << b_hi[i]      << " & " 
      << NColl_mean[i] << " & " << NColl_rms[i] << " & " 
      << NPart_mean[i] << " & " << NPart_rms[i] << " & " 
      << TAA[i]        << " & " << TAA_rms[i]   << "  \\\\" << endl;

  }

  cout << "\\end{tabular}"                   << endl;
  cout << "\\caption{\\label{tab:cent_pars}" << endl;
  cout << "  my table caption.}"             << endl;
  cout << "\\end{table*}"                    << endl << endl;

  cout << endl << endl;
  cout << "## NColl ## " << endl;
  for(int i=0; i<Centrality_lo.size(); i++)
    cout << std::setprecision(5) << NColl_mean[i] << ", ";
  cout << endl << "## NPart ## " << endl;
  for(int i=0; i<Centrality_lo.size(); i++)
    cout << std::setprecision(5) << NPart_mean[i] << ", ";
  cout << endl << "## TAA   ## " << endl;
  for(int i=0; i<Centrality_lo.size(); i++)
    cout << std::setprecision(5) << TAA[i]        << ", ";
  cout << endl;


  cout << endl << endl;
  cout << "## NColl ## " << endl;
  for(int i=0; i<Centrality_lo.size(); i++)
    cout << std::setprecision(5) << NColl_mean[i] << endl;
  cout << endl << "## NPart ## " << endl;
  for(int i=0; i<Centrality_lo.size(); i++)
    cout << std::setprecision(5) << NPart_mean[i] << endl;
  cout << endl << "## TAA   ## " << endl;
  for(int i=0; i<Centrality_lo.size(); i++)
    cout << std::setprecision(5) << TAA[i]        << endl;
  cout << endl;




}




int InitializeVectors()
{


  //Centrality_lo.push_back(0  ); Centrality_hi.push_back(5  ); Sqrts.push_back(2.76); XSec.push_back(64.0);
  //Centrality_lo.push_back(5  ); Centrality_hi.push_back(10 ); Sqrts.push_back(2.76); XSec.push_back(64.0);
  //Centrality_lo.push_back(10 ); Centrality_hi.push_back(20 ); Sqrts.push_back(2.76); XSec.push_back(64.0);
  //Centrality_lo.push_back(20 ); Centrality_hi.push_back(40 ); Sqrts.push_back(2.76); XSec.push_back(64.0);
  //Centrality_lo.push_back(40 ); Centrality_hi.push_back(60 ); Sqrts.push_back(2.76); XSec.push_back(64.0);
  //Centrality_lo.push_back(60 ); Centrality_hi.push_back(80 ); Sqrts.push_back(2.76); XSec.push_back(64.0);
  //Centrality_lo.push_back(80 ); Centrality_hi.push_back(100); Sqrts.push_back(2.76); XSec.push_back(64.0);

  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_64.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_64.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_64.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_64.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_64.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_64.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_64.0mb_0.4fm_1000000evt.root");



  //Centrality_lo.push_back(0  ); Centrality_hi.push_back(5  ); Sqrts.push_back(5.02); XSec.push_back(65.0);
  //Centrality_lo.push_back(5  ); Centrality_hi.push_back(10 ); Sqrts.push_back(5.02); XSec.push_back(65.0);
  //Centrality_lo.push_back(10 ); Centrality_hi.push_back(20 ); Sqrts.push_back(5.02); XSec.push_back(65.0);
  //Centrality_lo.push_back(20 ); Centrality_hi.push_back(40 ); Sqrts.push_back(5.02); XSec.push_back(65.0);
  //Centrality_lo.push_back(40 ); Centrality_hi.push_back(60 ); Sqrts.push_back(5.02); XSec.push_back(65.0);
  //Centrality_lo.push_back(60 ); Centrality_hi.push_back(80 ); Sqrts.push_back(5.02); XSec.push_back(65.0);
  //Centrality_lo.push_back(80 ); Centrality_hi.push_back(100); Sqrts.push_back(5.02); XSec.push_back(65.0);

  //Centrality_lo.push_back(0  ); Centrality_hi.push_back(5  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(5  ); Centrality_hi.push_back(10 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(10 ); Centrality_hi.push_back(15 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(15 ); Centrality_hi.push_back(20 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(20 ); Centrality_hi.push_back(25 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(25 ); Centrality_hi.push_back(30 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(30 ); Centrality_hi.push_back(35 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(35 ); Centrality_hi.push_back(40 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(40 ); Centrality_hi.push_back(45 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(45 ); Centrality_hi.push_back(50 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(50 ); Centrality_hi.push_back(55 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(55 ); Centrality_hi.push_back(60 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(60 ); Centrality_hi.push_back(65 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(65 ); Centrality_hi.push_back(70 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(70 ); Centrality_hi.push_back(75 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(75 ); Centrality_hi.push_back(80 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(80 ); Centrality_hi.push_back(85 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(85 ); Centrality_hi.push_back(90 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(90 ); Centrality_hi.push_back(95 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(95 ); Centrality_hi.push_back(100); Sqrts.push_back(5.02); XSec.push_back(70.0);

  //######################################################################################################
  //a whole bunch of new bins: March 14, 2017
  //Centrality_lo.push_back(0  ); Centrality_hi.push_back(5   ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(5  ); Centrality_hi.push_back(10  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(10 ); Centrality_hi.push_back(15  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(15 ); Centrality_hi.push_back(20  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(20 ); Centrality_hi.push_back(25  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(25 ); Centrality_hi.push_back(30  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(30 ); Centrality_hi.push_back(35  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(35 ); Centrality_hi.push_back(40  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(40 ); Centrality_hi.push_back(45  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(45 ); Centrality_hi.push_back(50  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(50 ); Centrality_hi.push_back(55  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(55 ); Centrality_hi.push_back(60  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(60 ); Centrality_hi.push_back(65  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(65 ); Centrality_hi.push_back(70  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(70 ); Centrality_hi.push_back(75  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(75 ); Centrality_hi.push_back(80  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(80 ); Centrality_hi.push_back(85  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(85 ); Centrality_hi.push_back(90  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(90 ); Centrality_hi.push_back(95  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(95 ); Centrality_hi.push_back(100 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(0  ); Centrality_hi.push_back(10  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(10 ); Centrality_hi.push_back(20  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(20 ); Centrality_hi.push_back(30  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(30 ); Centrality_hi.push_back(40  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(40 ); Centrality_hi.push_back(50  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(50 ); Centrality_hi.push_back(60  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(60 ); Centrality_hi.push_back(70  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(70 ); Centrality_hi.push_back(80  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(80 ); Centrality_hi.push_back(90  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(90 ); Centrality_hi.push_back(100 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(0  ); Centrality_hi.push_back(20  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(20 ); Centrality_hi.push_back(40  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(40 ); Centrality_hi.push_back(60  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(60 ); Centrality_hi.push_back(80  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(80 ); Centrality_hi.push_back(100 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(0  ); Centrality_hi.push_back(50  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(50 ); Centrality_hi.push_back(100 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(70 ); Centrality_hi.push_back(100 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(0  ); Centrality_hi.push_back(30  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  Centrality_lo.push_back(10 ); Centrality_hi.push_back(30  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(30 ); Centrality_hi.push_back(100 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(50 ); Centrality_hi.push_back(80  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(40 ); Centrality_hi.push_back(100 ); Sqrts.push_back(5.02); XSec.push_back(70.0);

  FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.8fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");

  //######################################################################################################

  //######################################################################################################
  //charged particle RAA centralites
  //  Centrality_lo.push_back(0  ); Centrality_hi.push_back(5  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //  Centrality_lo.push_back(5  ); Centrality_hi.push_back(10 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //  Centrality_lo.push_back(10 ); Centrality_hi.push_back(30 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //  Centrality_lo.push_back(30 ); Centrality_hi.push_back(50 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //  Centrality_lo.push_back(50 ); Centrality_hi.push_back(70 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //  Centrality_lo.push_back(70 ); Centrality_hi.push_back(90 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //  Centrality_lo.push_back(0  ); Centrality_hi.push_back(10 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //  Centrality_lo.push_back(0  ); Centrality_hi.push_back(100); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //
  //  FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //  FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //  FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //  FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //  FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //  FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //  FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //  FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //
  //######################################################################################################


  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  //FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");


  return Centrality_lo.size();
}
