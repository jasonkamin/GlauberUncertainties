Double_t FindCentralities(int nCents, Double_t *bedges, TH1 *h_bdist);
Double_t FindCentralitiesFromHydjet(int nCents, Double_t *bedges, TH2 *h2_hf_npart);
void     MakeCentHistos(TNtuple *mytup, char *variable, char *cond, int nCent, Double_t *bedges, TH1D *h_vsbcent, TH1D *h_vscent);
void     MakeCentHistosCustomEdges(TNtuple *mytup, char *variable, char *cond, int nCent, Double_t *bedges_lo, Double_t *bedges_hi, TH1D *h_vsbcent, TH1D *h_vscent);

TCanvas *cTemp = new TCanvas("cTemp","cTemp");

Double_t FindCentralities(int nCents, Double_t *bedges, TH1 *h_bdist)
{

  cout << "______ Starting CentralityHelperFunctions::FindCentralities() _____" << endl;
  Double_t binprecision = 0.0;
  binprecision = double(nCents)/double(h_bdist->GetNbinsX());

  //##########################################################################
  // Get the edges for the centralities ! 
  Int_t b_lo = 1;
  Int_t b_hi = 1;

  Double_t binboundaryN[nCentralities+1] = {0};
  binboundaryN[0] = 1;
  bedges[0] = h_bdist->GetBinLowEdge(binboundaryN[0]);
  binboundaryN[nCentralities] = h_bdist->GetNbinsX();
  bedges[nCentralities] = h_bdist->GetBinLowEdge(binboundaryN[nCentralities]);
  Double_t fullintegral = h_bdist->Integral(-1,-1);
  for(int ibin=1; ibin<nCents+1; ibin++){

    b_lo = b_hi;
    for(Int_t iquick=b_lo; iquick<h_bdist->GetNbinsX()+1; iquick++){
      Double_t partintegral = h_bdist->Integral(b_lo,iquick);
      if(partintegral/fullintegral>1.0/nCentralities){
        b_hi = iquick;
        binboundaryN[ibin] = b_hi;
        bedges[ibin] = h_bdist->GetBinLowEdge(binboundaryN[ibin]);
        break;
      }
    }
    cout << ibin << " - b_lo,b_hi : " << binboundaryN[ibin-1] << "," << binboundaryN[ibin] << "   " 
      << bedges[ibin-1] << "," << bedges[ibin] << endl;

  }
  int maxbbin = -1;
  for(int ibin=h_bdist->GetNbinsX(); ibin>0; ibin--){
    if(h_bdist->GetBinContent(ibin)>0){
      maxbbin = ibin;
      break;
    }
  }

  binboundaryN[nCentralities] = maxbbin+1;
  //bedges[nCentralities] = h_bdist->GetBinLowEdge(binboundaryN[nCentralities]);
  //if(bedges[nCentralities] != h_bdist->GetBinLowEdge(binboundaryN[nCentralities]))
  //cout << "  !!!!!  CAREFUL -- you've messed with the HF efficiency !!!!! " << endl;
  cout << nCentralities+1 << " - b_lo,b_hi : " << binboundaryN[nCentralities-1] << "," << binboundaryN[nCentralities] << "   " << bedges[nCentralities-1] << "," << bedges[nCentralities] << endl;
  //##########################################################################

  //TCanvas *cc = new TCanvas("c_impact_dist","c_impact_dist");

  //h_bdist->DrawCopy();
  //TLine *lineDiv[nCents+1];
  //for(int i=0; i<nCents+1; i++){
  //  lineDiv[i] = new TLine(bedges[i],0,bedges[i],3500);
  //  lineDiv[i]->SetLineColor(2);
  //  lineDiv[i]->SetLineStyle(2);
  //  lineDiv[i]->DrawClone("same");
  //}

  if(binprecision>0.025)
    cout << "Warning... bin precision of centrality bins is greater than 2.5\% ! " << endl;

  cout << "______ Finishing CentralityHelperFunctions::FindCentralities() _____" << endl;
  return binprecision;
}




Double_t FindCentralitiesFromHydjet(int nCents, Double_t *bedges, TH1 *h_hfdist)
{
  // Currently a work in progress....  1 Mar 2017

  cout << "______ Starting CentralityHelperFunctions::FindCentralities() _____" << endl;
  Double_t binprecision = 0.0;
  binprecision = double(nCents)/double(h_hfdist->GetNbinsX());

  //##########################################################################
  // Get the edges for the centralities ! 
  Int_t b_lo = 1;
  Int_t b_hi = 1;

  Double_t binboundaryN[nCentralities+1] = {0};
  binboundaryN[0] = 1;
  bedges[0] = h_hfdist->GetBinLowEdge(binboundaryN[0]);
  binboundaryN[nCentralities] = h_hfdist->GetNbinsX();
  bedges[nCentralities] = h_hfdist->GetBinLowEdge(binboundaryN[nCentralities]);
  Double_t fullintegral = h_hfdist->Integral(-1,-1);
  for(int ibin=1; ibin<nCents+1; ibin++){

    b_lo = b_hi;
    for(Int_t iquick=b_lo; iquick<h_hfdist->GetNbinsX()+1; iquick++){
      Double_t partintegral = h_hfdist->Integral(b_lo,iquick);
      if(partintegral/fullintegral>1.0/nCentralities){
        b_hi = iquick;
        binboundaryN[ibin] = b_hi;
        bedges[ibin] = h_hfdist->GetBinLowEdge(binboundaryN[ibin]);
        break;
      }
    }
    cout << ibin << " - b_lo,b_hi : " << binboundaryN[ibin-1] << "," << binboundaryN[ibin] << "   " 
      << bedges[ibin-1] << "," << bedges[ibin] << endl;

  }
  int maxbbin = -1;
  for(int ibin=h_hfdist->GetNbinsX(); ibin>0; ibin--){
    if(h_hfdist->GetBinContent(ibin)>0){
      maxbbin = ibin;
      break;
    }
  }

  binboundaryN[nCentralities] = maxbbin+1;
  //bedges[nCentralities] = h_hfdist->GetBinLowEdge(binboundaryN[nCentralities]);
  //if(bedges[nCentralities] != h_hfdist->GetBinLowEdge(binboundaryN[nCentralities]))
  //cout << "  !!!!!  CAREFUL -- you've messed with the HF efficiency !!!!! " << endl;
  cout << nCentralities+1 << " - b_lo,b_hi : " << binboundaryN[nCentralities-1] << "," << binboundaryN[nCentralities] << "   " << bedges[nCentralities-1] << "," << bedges[nCentralities] << endl;
  //##########################################################################

  //TCanvas *cc = new TCanvas("c_impact_dist","c_impact_dist");

  //h_hfdist->DrawCopy();
  //TLine *lineDiv[nCents+1];
  //for(int i=0; i<nCents+1; i++){
  //  lineDiv[i] = new TLine(bedges[i],0,bedges[i],3500);
  //  lineDiv[i]->SetLineColor(2);
  //  lineDiv[i]->SetLineStyle(2);
  //  lineDiv[i]->DrawClone("same");
  //}

  if(binprecision>0.025)
    cout << "Warning... bin precision of centrality bins is greater than 2.5\% ! " << endl;

  cout << "______ Finishing CentralityHelperFunctions::FindCentralities() _____" << endl;
  return binprecision;
}




void MakeCentHistos(TNtuple *mytup, char *variable, char *cond, int nCent, Double_t *edges, TH1D *h_vsbcent, TH1D *h_vscent)
{

  cout << "______ Starting CentralityHelperFunctions::MakeCentHistos() _____" << endl;

  char sayme[100];
  char sayme2[100];
  TH1D *htemp;
  //##########################################################################
  // make the centrality histograms...

  sprintf(sayme,"%s_bcent",variable);
  //h_vbcent = new TH1D(sayme,sayme,nCentralities,binedges);

  cout << h_vscent->GetName() << endl;

  for(int ibin=1; ibin<h_vscent->GetNbinsX()+1; ibin++){
    sprintf(sayme, "%s>>htemp(10000,0,8000)",variable);
    sprintf(sayme2,"B>%2.7f && B<%2.7f",edges[ibin-1],edges[ibin]);
    cout << sayme  << endl;
    cout << sayme2 << endl;
    //((TCanvas*)gROOT->Get("cTemp"))->cd();
    cTemp->cd();
    mytup->Draw(sayme,sayme2);
    htemp = (TH1D*)gPad->GetPrimitive("htemp");
    h_vscent ->SetBinContent(ibin,htemp->GetMean());
    h_vscent ->SetBinError  (ibin,htemp->GetStdDev());
    h_vsbcent->SetBinContent(ibin,htemp->GetMean());
    h_vsbcent->SetBinError  (ibin,htemp->GetStdDev());
    //NColl_cent[iFile][jVar]->SetBinContent(ibin,htempNColl->GetStdDev()*TMath::Sqrt(htempNColl->GetEntries()));
    cout << float(ibin)/nCentralities << "  " << htemp->GetMean() << " +/- " << htemp->GetStdDev() << endl;
    htemp->Reset();
  }

  //##########################################################################

  cout << "______ Finishing CentralityHelperFunctions::MakeCentHistos() _____" << endl;
  return;
}





void MakeCentHistosCustomEdges(TNtuple *mytup, char *variable, char *cond, int nCent, Double_t *edges_lo, Double_t *edges_hi, TH1D *h_vsbcent, TH1D *h_vscent)
{

  cout << "______ Starting CentralityHelperFunctions::MakeCentHistosCustomEdges() _____" << endl;
  cout << "                 " << h_vscent->GetNbinsX() << " bins." << endl;

  char sayme[500];
  char sayme2[500];
  TH1D *htemp;
  //##########################################################################
  // make the centrality histograms...

  sprintf(sayme,"%s_bcent",variable);
  //h_vbcent = new TH1D(sayme,sayme,nCentralities,binedges);

  cout << h_vscent->GetName() << endl;

  for(int ibin=0; ibin<h_vscent->GetNbinsX(); ibin++){
    sprintf(sayme, "%s>>htemp(10000,0,8000)",variable);
    sprintf(sayme2,"B>%2.7f && B<%2.7f",edges_lo[ibin],edges_hi[ibin]);
    cout << sayme  << endl;
    cout << sayme2 << endl;
    //((TCanvas*)gROOT->Get("cTemp"))->cd();
    cTemp->cd();
    mytup->Draw(sayme,sayme2);
    htemp = (TH1D*)gPad->GetPrimitive("htemp");
    h_vscent ->SetBinContent(ibin+1,htemp->GetMean());
    h_vscent ->SetBinError  (ibin+1,htemp->GetStdDev());
    h_vsbcent->SetBinContent(ibin+1,htemp->GetMean());
    h_vsbcent->SetBinError  (ibin+1,htemp->GetStdDev());
    //NColl_cent[iFile][jVar]->SetBinContent(ibin+1,htempNColl->GetStdDev()*TMath::Sqrt(htempNColl->GetEntries()));
    //cout << float(ibin+1)/nCentralities << "  " << htemp->GetMean() << " +/- " << htemp->GetStdDev() << endl;
    cout << float(ibin+1) << "  " << htemp->GetMean() << " +/- " << htemp->GetStdDev() << endl;
    htemp->Reset();
  }

  //##########################################################################

  cout << "______ Finishing CentralityHelperFunctions::MakeCentHistosCustomEdges() _____" << endl;
  return;
}

