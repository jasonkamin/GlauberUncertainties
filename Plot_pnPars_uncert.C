static const int nCentralities = 20;

Double_t sqrtsNN[8] = {2.76,  5.0,  5.5, 8.16,  8.8,   39,   63,   100};
Double_t sigmaNN[8] = {64.0, 65.0, 71.0, 73.0, 77.0, 95.0, 98.0, 100.0};
Double_t erSigNN[8] = { 2.0,  2.0,  2.0,  2.0,  2.0,  2.0,  3.0,   3.0};

char systems[4][2][5] = {{"Pb","Pb"}, {"Pbpn","Pbpn"}, {"Pb","p"}, {"Pbpn","p"}};

Double_t dmin[3] = {0.0,0.4,0.8};
int COLORS[8][3] = {
  {2, 1,2},
  {3, 1,3},
  {4, 1,4},
  {41,1,41},
  {6, 1,6},
  {7, 1,7},
  {8, 1,8},
  {9, 1,9}
};

TProfile *myProfNom;
TH1D *NColl_bNom;
TH2D *htempNom;
TH2D *htempNom2;
TProfile *myProf[100];
TH1D *NColl_b[100];
TH1D *NColl_br[100];
TH2D *htemp[100];
TH2D *htemp2[100];
TH2D *h_2Dvaris;

TH1D *h_b;

TH1D *NColl_bcent_Nom;
TH1D *NColl_cent_Nom;
TH1D *NColl_bcent[100];
TH1D *NColl_cent[100];
TH1D *NColl_bcentr[100];
TH1D *NColl_centr[100];

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


int Uniform_or_Custom_bins = 1;// Uniform:1  Custom:2 (note: for Custom, set nCentralities=100 !)
int plot_pbpb_or_ppb = 1;
int ncoll_or_npart = 1;// NColl:1  NPart:2

#include "CentralityHelperFunctions.C"
int InitializeVectors();

void Plot_pnPars_uncert()
{

  // first safety/sanity check. 
  if(Uniform_or_Custom_bins==2 and nCentralities!=100){
    cout << "##########  HEY !!  ##########" << endl;
    cout << " for custom bins ( fed via InitializeVectors() ), you must set nCentralities=100. " << endl;
    cout << "##########  BREAKING  ########" << endl << endl;
    return;
  }

  InitializeVectors();

  Double_t binedges[nCentralities+1] = {0};
  const int nCentVectors = Centrality_lo.size();

  Double_t binedge_lo[nCentVectors];
  Double_t binedge_hi[nCentVectors];
  char saythis[500];
  char saythis2[500];
  char pathtovariations[500];
  TString st_pathtovar;
  TNtuple *myTuple;

  TCanvas *c1 = new TCanvas("c1","c1");

  TFile *f0[100];
  TFile *f2;

  int xbins=200; float xlo=0; float xhi=30; int ybins=3000; float ylo=0.5; float yhi=3000.5; int rebinx=4; int rebiny=20;
  if(plot_pbpb_or_ppb==1){
    xbins=500; xlo=0; xhi=30; ybins=3000; ylo=0.5; yhi=3000.5; rebinx=1; rebiny=1;//rebinx=4; rebiny=20;
    f2 = TFile::Open("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
    //f2 = TFile::Open("/Users/jason/Dropbox/MyDocs/ALICE/Papers/glauber/glauber-2016/macros/mc/Pbpn_Rp6.68_Ap0.447_Rn6.70_An0.55_Pbpn_Rp6.68_Ap0.447_Rn6.70_An0.55_65.0mb_0.4fm_1000000evt.root");
    //sprintf(pathtovariations,"~/Desktop/rootfiles/Pb_Pb_nominalvariations/Nom_65mb/Pb_Pb_RAvariations_65.0mb_0.4fm_100100evt_p");
    sprintf(pathtovariations,"~/Desktop/rootfiles/Pb_Pb_nominalvariations/Nom_70mb/Pb_Pb_RAvariations_70.0mb_0.4fm_100100evt_p");
    //sprintf(pathtovariations,"~/Desktop/rootfiles/Pb_Pb_nominalvariations/Pbpn_Pbpn_noPerrors/Pbpn_Pbpn_noPerrors_RAvariations_65.0mb_0.4fm_100100evt_p");

    st_pathtovar = pathtovariations;
  }
  else if(plot_pbpb_or_ppb==2){
    xbins=1000; xlo=0; xhi=16; ybins=40; ylo=0.5; yhi=40.5; rebinx=4; rebiny=1;
    cout << "pPb currently non-functional..." << endl;
    cout << "BREAKING !" << endl;
    return;
    // must set path to pPb variations here ! 
    f2 = TFile::Open("/Users/jason/Dropbox/MyDocs/ALICE/Papers/glauber/glauber-2016/macros/mc/Pbpn_Rp6.68_Ap0.447_Rn6.70_An0.55_p_73.0mb_0.4fm_1000000evt.root");
  }
  int xbins_for_b = 100*xbins;


  h_2Dvaris = new TH2D("h_2Dvaris","h_2Dvaris",xbins/rebinx,xlo,xhi, ybins/rebiny,ylo,yhi);

  //##### do the nominal case #####
  cout << endl << "Starting the 'central values' case... " << endl;
  myTuple = (TNtuple*)f2->Get(f2->GetListOfKeys()->At(0)->GetName())->Clone("myTuple");
  if(ncoll_or_npart==1)      sprintf(saythis,"Ncoll:B>>htempNom(%d,%2.1f,%2.1f,%d,%2.1f,%2.1f)",xbins,xlo,xhi, ybins,ylo,yhi);
  else if(ncoll_or_npart==2) sprintf(saythis,"Npart:B>>htempNom(%d,%2.1f,%2.1f,%d,%2.1f,%2.1f)",xbins,xlo,xhi, ybins,ylo,yhi);
  myTuple->Draw(saythis,"","colz");
  htempNom = (TH2D*)gPad->GetPrimitive("htempNom")->Clone("htempNom");
  htempNom ->Draw("colz");
  htempNom2 = (TH2D*)htempNom->Rebin2D(rebinx,rebiny)->Clone("htempNom2");

  myProfNom = (TProfile*)(htempNom->ProfileX("_pfx",1,-1,"s")->Clone("myProfNom"));
  NColl_bNom = new TH1D("NColl_bNom","NColl_bNom",myProfNom->GetNbinsX(),myProfNom->GetBinLowEdge(1),myProfNom->GetBinLowEdge(myProfNom->GetNbinsX()+1));
  NColl_bCom = new TH1D("NColl_bCom","NColl_bCom",myProfNom->GetNbinsX(),myProfNom->GetBinLowEdge(1),myProfNom->GetBinLowEdge(myProfNom->GetNbinsX()+1));

  //################# Find the Centralities #####################ALLlxplus#
  c1->cd();
  sprintf(saythis,"B>>h_b(%d,%2.5f,%2.5f)",xbins_for_b,xlo,xhi);
  myTuple->Draw(saythis);
  sprintf(saythis,"h_b");
  h_b = (TH1D*)gPad->GetPrimitive(saythis);
  Double_t precision_quick = FindCentralities(nCentralities, binedges, h_b);
  cout << "precision:  " << precision_quick << endl;
  for(int i=0; i<nCentVectors; i++){
    binedge_lo[i] = binedges[Centrality_lo[i]];
    binedge_hi[i] = binedges[Centrality_hi[i]];
  }
  //###############################################################

  //################# Make the Histos #######################
  sprintf(saythis,"Ncoll_bcent_Nom");
  if(Uniform_or_Custom_bins==1)       NColl_bcent_Nom = new TH1D(saythis,saythis,nCentralities,binedges);
  else if(Uniform_or_Custom_bins==2)  NColl_bcent_Nom = new TH1D(saythis,saythis,nCentVectors,0,100);
  // not in use yet.  this is for when we will bin in HF energy. 
  //if(binedges[0]<binedges[nCentralities])
  //  NColl_bcent_Nom = new TH1D(saythis,saythis,nCentVectors,0,100);
  //  //NColl_bcent_Nom = new TH1D(saythis,saythis,nCentralities,binedges);
  //else{
  //  Double_t binedges2[nCentralities] = {0};
  //  int j = nCentralities;
  //  for(int i=0; i<nCentralities; i++){
  //    binedges2[i] = binedges[j--];
  //    cout << binedges2[i] << endl;
  //  }
  //  NColl_bcent_Nom = new TH1D(saythis,saythis,nCentralities,binedges2);
  //}
  sprintf(saythis,"Ncoll_cent_Nom");
  if(Uniform_or_Custom_bins==1)       NColl_cent_Nom  = new TH1D(saythis,saythis,nCentralities,0,100);
  else if(Uniform_or_Custom_bins==2)  NColl_cent_Nom  = new TH1D(saythis,saythis,nCentVectors,0,100);
  if(ncoll_or_npart==1)      {sprintf(saythis,"Ncoll"); sprintf(saythis2,"");}
  else if(ncoll_or_npart==2) {sprintf(saythis,"Npart"); sprintf(saythis2,"");}
  if(Uniform_or_Custom_bins==1)       MakeCentHistos(myTuple,saythis,saythis2,nCentralities,binedges,NColl_bcent_Nom,NColl_cent_Nom);
  else if(Uniform_or_Custom_bins==2)  MakeCentHistosCustomEdges(myTuple,saythis,saythis2,nCentralities,binedge_lo,binedge_hi,NColl_bcent_Nom,NColl_cent_Nom);
  //###############################################################

  for(int ibin=0; ibin<NColl_bNom->GetNbinsX()+1; ibin++){
    if(myProfNom ->GetBinContent(ibin)!=0){
      NColl_bNom ->SetBinContent(ibin,myProfNom->GetBinContent(ibin));
      NColl_bNom ->SetBinError  (ibin,myProfNom->GetBinError  (ibin));
    }
  }
  NColl_bNom->SetLineColor  (1);
  NColl_bNom->SetMarkerColor(1);
  NColl_bNom->SetMarkerStyle(20);
  NColl_bNom->SetMarkerSize (0.80);
  NColl_bNom->SetLineWidth  (3);

  NColl_bcent_Nom->SetLineColor  (1);
  NColl_bcent_Nom->SetMarkerColor(1);
  NColl_bcent_Nom->SetMarkerStyle(20);
  NColl_bcent_Nom->SetMarkerSize (0.80);
  NColl_bcent_Nom->SetLineWidth  (3);
  NColl_cent_Nom ->SetLineColor  (1);
  NColl_cent_Nom ->SetMarkerColor(1);
  NColl_cent_Nom ->SetMarkerStyle(20);
  NColl_cent_Nom ->SetMarkerSize (0.80);
  NColl_cent_Nom ->SetLineWidth  (3);

  cout << "... Finished the 'central values' case ! " << endl;

  //return;

  Double_t runAve[200] = {0};
  Double_t runStD[200] = {0};

  //##### do the variations ! #####
  cout << endl << "Starting the variations... " << endl;
  TList *mylist;
  int nKeys = 99;
  cout << "We found " << nKeys << " keys in file." << endl;
  for(int i=0; i<nKeys; i++){

    sprintf(saythis,"%s%d.root",st_pathtovar.Data(),i);
    cout << " try to open  " << saythis << endl;
    f0[i] = TFile::Open(saythis);
    mylist = f0[i]->GetListOfKeys();
    myTuple = (TNtuple*)f0[i]->Get(mylist->At(0)->GetName())->Clone("myTuple");
    cout << "KEY[" << i << "]   name: " << mylist->At(0)->GetName() << "    title: " << ((TNtuple*)f0[i]->Get(mylist->At(0)->GetName()))->GetTitle() << "   at: " << myTuple << endl;    

    if(ncoll_or_npart==1)      sprintf(saythis,"Ncoll:B>>htemp_%d(%d,%2.1f,%2.1f,%d,%2.1f,%2.1f)",i,xbins,xlo,xhi, ybins,ylo,yhi);
    else if(ncoll_or_npart==2) sprintf(saythis,"Npart:B>>htemp_%d(%d,%2.1f,%2.1f,%d,%2.1f,%2.1f)",i,xbins,xlo,xhi, ybins,ylo,yhi);
    cout << "Draw --> " << saythis << endl;
    myTuple->Draw(saythis,"","colz");
    sprintf(saythis,"htemp_%d",i);
    htemp[i] = (TH2D*)gPad->GetPrimitive(saythis)->Clone(saythis);
    htemp[i] ->Draw("colz");
    sprintf(saythis,"htemp2_%d",i);
    htemp2[i] = (TH2D*)htemp[i]->Rebin2D(rebinx,rebiny)->Clone(saythis);
    h_2Dvaris->Add(htemp2[i]);

    sprintf(saythis,"myProf_%d",i);
    myProf[i] = (TProfile*)(htemp[i]->ProfileX("_pfx",1,-1,"s")->Clone(saythis));



    //################# Find the Centralities #######################
    c1->cd();
    sprintf(saythis,"B>>h_b(%d,%2.3f,%2.3f)",xbins_for_b,xlo,xhi);
    myTuple->Draw(saythis);
    sprintf(saythis,"h_b");
    h_b = (TH1D*)gPad->GetPrimitive(saythis);
    precision_quick = FindCentralities(nCentralities, binedges, h_b);
    cout << "precision:  " << precision_quick << endl;
    for(int i=0; i<nCentVectors; i++){
      binedge_lo[i] = binedges[Centrality_lo[i]];
      binedge_hi[i] = binedges[Centrality_hi[i]];
    }
    //###############################################################

    //################# Make the Histos #######################
    sprintf(saythis,"Ncoll_bcent_%d",i);
    if(Uniform_or_Custom_bins==1)       NColl_bcent[i] = new TH1D(saythis,saythis,nCentralities,binedges);
    else if(Uniform_or_Custom_bins==2)  NColl_bcent[i] = new TH1D(saythis,saythis,nCentVectors,0,100);
    sprintf(saythis,"Ncoll_cent_%d",i);
    if(Uniform_or_Custom_bins==1)       NColl_cent[i]  = new TH1D(saythis,saythis,nCentralities,0,100);
    else if(Uniform_or_Custom_bins==2)  NColl_cent[i]  = new TH1D(saythis,saythis,nCentVectors,0,100);
    if(ncoll_or_npart==1)      {sprintf(saythis,"Ncoll"); sprintf(saythis2,"");}
    else if(ncoll_or_npart==2) {sprintf(saythis,"Npart"); sprintf(saythis2,"");}
    if(Uniform_or_Custom_bins==1)       MakeCentHistos(myTuple,saythis,saythis2,nCentralities,binedges,NColl_bcent[i],NColl_cent[i]);
    else if(Uniform_or_Custom_bins==2)  MakeCentHistosCustomEdges(myTuple,saythis,saythis2,nCentralities,binedge_lo,binedge_hi,NColl_bcent[i],NColl_cent[i]);
    //###############################################################


    sprintf(saythis,"NColl_b_%d",i);
    NColl_b[i] = new TH1D(saythis,saythis,myProf[i]->GetNbinsX(),myProf[i]->GetBinLowEdge(1),myProf[i]->GetBinLowEdge(myProf[i]->GetNbinsX()+1));
    for(int ibin=0; ibin<NColl_b[i]->GetNbinsX()+1; ibin++){
      if(myProf[i] ->GetBinContent(ibin)!=0){
         NColl_b[i] ->SetBinContent(ibin,myProf[i]->GetBinContent(ibin));
         NColl_b[i] ->SetBinError  (ibin,myProf[i]->GetBinError  (ibin));
         //h_2Dvaris->Fill(NColl_b[i]->GetBinCenter(ibin),NColl_b[i]->GetBinContent(ibin));
         runAve[ibin] += NColl_b[i] ->GetBinContent(ibin)/nKeys;
      }
    }
    NColl_b[i]->SetLineColor  (i);
    NColl_b[i]->SetMarkerColor(i);
    NColl_b[i]->SetMarkerStyle(24);
    NColl_b[i]->SetMarkerSize (0.40);

    NColl_bcent[i]->SetLineColor  (i);
    NColl_bcent[i]->SetMarkerColor(i);
    NColl_bcent[i]->SetMarkerStyle(20);
    NColl_bcent[i]->SetMarkerSize (0.80);

    NColl_cent[i]->SetLineColor  (i);
    NColl_cent[i]->SetMarkerColor(i);
    NColl_cent[i]->SetMarkerStyle(20);
    NColl_cent[i]->SetMarkerSize (0.80);


  }
  cout << "... Finished grabbing the parameter variations ! " << endl;


  for(int ibin=0; ibin<NColl_bNom->GetNbinsX()+1; ibin++){
    for(int i=0; i<nKeys; i++)
      runStD[ibin] += (NColl_b[i] ->GetBinContent(ibin)-runAve[ibin])*(NColl_b[i] ->GetBinContent(ibin)-runAve[ibin])/nKeys;
    runStD[ibin] = TMath::Sqrt(runStD[ibin]);
    NColl_bCom->SetBinContent(ibin,runAve[ibin]);
    NColl_bCom->SetBinError  (ibin,runStD[ibin]);
  }
  NColl_bCom->SetLineColor  (2);
  NColl_bCom->SetMarkerColor(2);
  NColl_bCom->SetMarkerStyle(24);
  NColl_bCom->SetMarkerSize (0.80);
  NColl_bCom->SetLineWidth  (2);

  TH1D *h_frame = new TH1D("h_frame","h_frame",10000,0,100);
  h_frame->GetXaxis()->SetRangeUser(0,xhi);
  h_frame->GetYaxis()->SetRangeUser(0,yhi);
  c1->cd();
  h_frame->DrawCopy();

  NColl_bNom->Draw("same");
  //NColl_bCom->Draw("same");
  for(int i=0; i<nKeys; i++){
    NColl_bcent[i]->Draw("same");
    //NColl_b[i]->Draw("same");
  }


  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);
  c2->cd(1)->SetLogz();
  h_2Dvaris->Draw("colz");//lego2z
  c2->cd(2)->SetLogz();
  htempNom2->Draw("colz");//lego2z


  TH2D *NColl_br_2D = new TH2D("NColl_br_2D","NColl_br_2D",NColl_bNom->GetNbinsX(),NColl_bNom->GetXaxis()->GetBinLowEdge(1),NColl_bNom->GetXaxis()->GetBinLowEdge(NColl_bNom->GetNbinsX()+1), 300,0.5,1.5);

  for(int i=0; i<nKeys; i++){
    sprintf(saythis,"NColl_br_%d",i);
    NColl_br[i] = (TH1D*)NColl_b[i]->Clone(saythis);
    NColl_br[i] ->Divide(NColl_bNom);
    NColl_br[i]->Draw("same");

    for(int j=0; j<NColl_bNom->GetNbinsX()+1; j++){
      NColl_br[i] ->SetBinError  (j,0.0);
      NColl_br_2D->Fill(NColl_br[i] ->GetBinCenter(j),NColl_br[i] ->GetBinContent(j));
    }
  }


  TH1D *proj_for_fit[500];
  TH1D *envelope = new TH1D("envelope","systematic uncertainty",NColl_bNom->GetNbinsX(),NColl_bNom->GetXaxis()->GetBinLowEdge(1),NColl_bNom->GetXaxis()->GetBinLowEdge(NColl_bNom->GetNbinsX()+1));
  for(int j=0; j<NColl_bNom->GetNbinsX(); j++){
    sprintf(saythis,"proj_for_fit_%d",j);
    proj_for_fit[j] = (TH1D*)NColl_br_2D->ProjectionY(saythis,j+1,j+1);
    envelope->SetBinContent(j+1,1.0);
    envelope->SetBinError  (j+1,proj_for_fit[j]->GetStdDev());
  }

  TCanvas *c4 = new TCanvas("c4","c4",1200,800);
  c4->Divide(5,5);
  //c4->cd(1);
  //NColl_br_2D->Draw("colz");
  int maxloopforExamples = 25;
  if(nCentralities<25 && Uniform_or_Custom_bins==1)      maxloopforExamples = nCentralities;
  else if(nCentVectors<25  && Uniform_or_Custom_bins==2) maxloopforExamples = nCentVectors;
  for(int j=0; j<maxloopforExamples; j++){
    c4->cd(j+1);
    proj_for_fit[j]->GetXaxis()->SetRangeUser(0.8,1.2);
    proj_for_fit[j]->DrawCopy();
  }


  //#############################################################
  TLegend *leg3 = new TLegend(0.6,0.6,0.89,0.89);
  TCanvas *c3   = new TCanvas("c3","c3",1050,750);

  c3->Clear();
  c3->Divide(2,1);
  c3->cd(1)->SetPad(0.01,0.4,0.99,0.99);
  c3->cd(1)->SetBottomMargin(0.0);
  if(ncoll_or_npart==1)      h_frame->GetYaxis()->SetTitle("N_{Coll}");
  else if(ncoll_or_npart==2) h_frame->GetYaxis()->SetTitle("N_{Part}");
  h_frame->GetXaxis()->SetTitle("b [fm]");
  h_frame->GetXaxis()->SetRangeUser(0,xhi);
  h_frame->GetYaxis()->SetRangeUser(0.9,yhi);

  NColl_bNom->Draw("same");
  for(int i=0; i<nKeys; i++){
    NColl_b[i]->Draw("same");
  }

  leg3->Draw();

  c3->cd(2)->SetPad(0.01,0.01,0.99,0.4);
  c3->cd(2)->SetTopMargin   (0.0);
  c3->cd(2)->SetBottomMargin(0.1);
  TLine *line1 = new TLine(0,1,80,1);
  line1->SetLineStyle(2);
  line1->SetLineColor(16);
  //h_frame->GetYaxis()->SetRangeUser(0.80,1.20);
  h_frame->GetYaxis()->SetRangeUser(0.90,1.10);
  h_frame->GetYaxis()->SetTitle("ratio");
  h_frame->SetTitle("ratio");
  h_frame->DrawCopy();

  envelope->SetFillColorAlpha(2,0.6);
  envelope->Draw("E3,same");

  for(int i=0; i<nKeys; i++)
    NColl_br[i]->Draw("pe,same");
  line1->Draw("same");

  //
  //#############################################################

  cTemp->cd();

  TH2D *NColl_bcentr_2D = new TH2D("NColl_bcentr_2D","NColl_bcentr_2D",NColl_bcent_Nom->GetNbinsX(),NColl_bcent_Nom->GetXaxis()->GetBinLowEdge(1),NColl_bcent_Nom->GetXaxis()->GetBinLowEdge(NColl_bcent_Nom->GetNbinsX()+1), 1000,0.0,2.0);
  TH2D *NColl_centr_2D = new TH2D("NColl_centr_2D","NColl_centr_2D",NColl_cent_Nom->GetNbinsX(),NColl_cent_Nom->GetXaxis()->GetBinLowEdge(1),NColl_cent_Nom->GetXaxis()->GetBinLowEdge(NColl_cent_Nom->GetNbinsX()+1), 1000,0.0,2.0);

  for(int i=0; i<nKeys; i++){
    sprintf(saythis,"NColl_bcentr_%d",i);
    NColl_bcentr[i] = (TH1D*)NColl_bcent[i]->Clone(saythis);
    //NColl_bcentr[i] ->Divide(NColl_bcent_Nom);
    //NColl_bcentr[i]->Draw("same");

    sprintf(saythis,"NColl_centr_%d",i);
    NColl_centr[i] = (TH1D*)NColl_cent[i]->Clone(saythis);
    //NColl_centr[i] ->Divide(NColl_cent_Nom);
    //NColl_centr[i]->Draw("same");

    for(int j=0; j<NColl_bcent_Nom->GetNbinsX()+1; j++){
      NColl_bcentr[i] ->SetBinContent(j,NColl_bcentr[i]->GetBinContent(j)/NColl_bcent_Nom->GetBinContent(j));
      NColl_bcentr[i] ->SetBinError  (j,0.0);
      NColl_bcentr_2D->Fill(NColl_bcentr[i] ->GetBinCenter(j),NColl_bcentr[i] ->GetBinContent(j));

      NColl_centr[i] ->SetBinContent(j,NColl_centr[i]->GetBinContent(j)/NColl_cent_Nom->GetBinContent(j));
      NColl_centr[i] ->SetBinError  (j,0.0);
      NColl_centr_2D->Fill(NColl_centr[i] ->GetBinCenter(j),NColl_centr[i] ->GetBinContent(j));
    }
  }

  TH1D *proj_for_fita[500];
  TH1D *envl_NColl_bcent = new TH1D("envl_NColl_bcent","systematic uncertainty",NColl_bcent_Nom->GetNbinsX(),NColl_bcent_Nom->GetXaxis()->GetBinLowEdge(1),NColl_bcent_Nom->GetXaxis()->GetBinLowEdge(NColl_bcent_Nom->GetNbinsX()+1));
  for(int j=0; j<NColl_bcent_Nom->GetNbinsX(); j++){
    sprintf(saythis,"proj_for_fit_%d",j);
    proj_for_fita[j] = (TH1D*)NColl_bcentr_2D->ProjectionY(saythis,j+1,j+1);
    envl_NColl_bcent->SetBinContent(j+1,1.0);
    envl_NColl_bcent->SetBinError  (j+1,proj_for_fita[j]->GetStdDev());
  }
  TCanvas *c4a = new TCanvas("c4a","c4a",1200,800);
  c4a->Divide(5,5);
  for(int j=0; j<maxloopforExamples; j++){
    c4a->cd(j+1);
    proj_for_fita[j]->GetXaxis()->SetRangeUser(0.8,1.2);
    proj_for_fita[j]->DrawCopy();
    proj_for_fita[j]->Reset();
  }

  if(ncoll_or_npart==1)      cout << endl << "   Starting the envelope calculation for NColl ..." << endl;
  else if(ncoll_or_npart==2) cout << endl << "   Starting the envelope calculation for NPart ..." << endl;
  TH1D *proj_for_fitb[500];
  TH1D *envl_NColl_cent = new TH1D("envl_NColl_cent","systematic uncertainty",NColl_cent_Nom->GetNbinsX(),NColl_cent_Nom->GetXaxis()->GetBinLowEdge(1),NColl_cent_Nom->GetXaxis()->GetBinLowEdge(NColl_cent_Nom->GetNbinsX()+1));
  for(int j=0; j<NColl_cent_Nom->GetNbinsX(); j++){
    sprintf(saythis,"proj_for_fit_%d",j);
    proj_for_fitb[j] = (TH1D*)NColl_centr_2D->ProjectionY(saythis,j+1,j+1);
    envl_NColl_cent->SetBinContent(j+1,1.0);
    envl_NColl_cent->SetBinError  (j+1,proj_for_fitb[j]->GetStdDev());
    cout << "  centrality " << j+1 << "  StdDev: " << proj_for_fitb[j]->GetStdDev() << endl;
  }
  TCanvas *c4b = new TCanvas("c4b","c4b",1200,800);
  c4b->Divide(5,5);
  for(int j=0; j<maxloopforExamples; j++){
    c4b->cd(j+1);
    proj_for_fitb[j]->GetXaxis()->SetRangeUser(0.8,1.2);
    proj_for_fitb[j]->DrawCopy();
    proj_for_fitb[j]->Reset();
  }



  //#############################################################
  TLegend *leg5 = new TLegend(0.6,0.6,0.89,0.89);
  TCanvas *c5   = new TCanvas("c5","c5",1050,750);

  c5->Clear();
  c5->Divide(2,1);
  c5->cd(1)->SetPad(0.01,0.4,0.99,0.99);
  c5->cd(1)->SetBottomMargin(0.0);
  h_frame->GetYaxis()->SetTitle("N_{Coll}");
  h_frame->GetXaxis()->SetTitle("b [fm]");
  h_frame->GetXaxis()->SetRangeUser(0,xhi);
  h_frame->GetYaxis()->SetRangeUser(0.9,yhi);

  NColl_bcent_Nom->Draw("same");
  for(int i=0; i<nKeys; i++){
   NColl_bcent[i]->Draw("same");
  }

  leg5->Draw();

  c5->cd(2)->SetPad(0.01,0.01,0.99,0.4);
  c5->cd(2)->SetTopMargin   (0.0);
  c5->cd(2)->SetBottomMargin(0.1);
  //h_frame->GetYaxis()->SetRangeUser(0.80,1.20);
  h_frame->GetYaxis()->SetRangeUser(0.90,1.10);
  h_frame->GetYaxis()->SetTitle("ratio");
  h_frame->SetTitle("ratio");
  h_frame->DrawCopy();

  envl_NColl_bcent->SetFillColorAlpha(2,0.6);
  envl_NColl_bcent->Draw("E3,same");

  for(int i=0; i<nKeys; i++)
    NColl_bcentr[i]->Draw("pe,same");
  line1->Draw("same");

  //
  //#############################################################



  //#############################################################
  TLegend *leg6 = new TLegend(0.6,0.6,0.89,0.89);
  TCanvas *c6   = new TCanvas("c6","c6",1050,750);

  c6->Clear();
  c6->Divide(2,1);
  c6->cd(1)->SetPad(0.01,0.4,0.99,0.99);
  c6->cd(1)->SetBottomMargin(0.0);
  h_frame->GetYaxis()->SetTitle("N_{Coll}");
  h_frame->GetXaxis()->SetTitle("b [fm]");
  h_frame->GetXaxis()->SetRangeUser(0,100);
  h_frame->GetYaxis()->SetRangeUser(0.9,yhi);

  NColl_cent_Nom->Draw("same");
  for(int i=0; i<nKeys; i++){
   NColl_cent[i]->Draw("same");
  }

  leg6->Draw();

  c6->cd(2)->SetPad(0.01,0.01,0.99,0.4);
  c6->cd(2)->SetTopMargin   (0.0);
  c6->cd(2)->SetBottomMargin(0.1);
  //h_frame->GetYaxis()->SetRangeUser(0.80,1.20);
  h_frame->GetYaxis()->SetRangeUser(0.90,1.10);
  h_frame->GetYaxis()->SetTitle("ratio");
  h_frame->GetXaxis()->SetTitle("Centrality [%]");
  h_frame->SetTitle("ratio");
  h_frame->DrawCopy();

  envl_NColl_cent->SetFillColorAlpha(2,0.6);
  envl_NColl_cent->Draw("E3,same");

  for(int i=0; i<nKeys; i++)
    NColl_centr[i]->Draw("pe,same");
  line1->Draw("same");

  //
  //#############################################################




  //#############################################################
  TLegend *leg7 = new TLegend(0.6,0.6,0.89,0.89);
  TCanvas *c7   = new TCanvas("c7","c7",1050,750);

  c7->Clear();
  c7->Divide(2,1);
  c7->cd(1)->SetPad(0.01,0.4,0.99,0.99);
  c7->cd(1)->SetBottomMargin(0.0);
  h_frame->GetYaxis()->SetTitle("N_{Coll}");
  h_frame->GetXaxis()->SetTitle("b [fm]");
  h_frame->GetXaxis()->SetRangeUser(0,30);
  h_frame->GetYaxis()->SetRangeUser(0.9,yhi);

  NColl_bcent_Nom->SetMarkerColor(2);
  NColl_bNom     ->SetMarkerColor(4);
  NColl_bcent_Nom->Draw("same");
  NColl_bNom     ->Draw("same");
  //for(int i=0; i<nKeys; i++){
  // NColl_cent[i]->Draw("same");
  //}

  leg7->Draw();

  c7->cd(2)->SetPad(0.01,0.01,0.99,0.4);
  c7->cd(2)->SetTopMargin   (0.0);
  c7->cd(2)->SetBottomMargin(0.1);
  //h_frame->GetYaxis()->SetRangeUser(0.80,1.20);
  h_frame->GetYaxis()->SetRangeUser(0.90,1.10);
  h_frame->GetYaxis()->SetTitle("ratio");
  h_frame->GetXaxis()->SetTitle("Centrality [%]");
  h_frame->SetTitle("ratio");
  h_frame->DrawCopy();

  envl_NColl_bcent->SetFillColorAlpha(2,0.4);
  envl_NColl_bcent->Draw("E3,same");

  envelope->SetFillColorAlpha(4,0.4);
  envelope->Draw("E3,same");

  line1->Draw("same");

  //
  //#############################################################



}





int InitializeVectors()
{


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
  //Centrality_lo.push_back(30 ); Centrality_hi.push_back(100 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(50 ); Centrality_hi.push_back(80  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  //Centrality_lo.push_back(40 ); Centrality_hi.push_back(100 ); Sqrts.push_back(5.02); XSec.push_back(70.0);

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
  Centrality_lo.push_back(70 ); Centrality_hi.push_back(100 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  Centrality_lo.push_back(0  ); Centrality_hi.push_back(30  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  Centrality_lo.push_back(10 ); Centrality_hi.push_back(30  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  Centrality_lo.push_back(30 ); Centrality_hi.push_back(100 ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  Centrality_lo.push_back(50 ); Centrality_hi.push_back(80  ); Sqrts.push_back(5.02); XSec.push_back(70.0);
  Centrality_lo.push_back(40 ); Centrality_hi.push_back(100 ); Sqrts.push_back(5.02); XSec.push_back(70.0);

  FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");
  FileName.push_back("/Users/jason/Desktop/rootfiles/Pb_Pb_70.0mb_0.4fm_1000000evt.root");

  //######################################################################################################

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
