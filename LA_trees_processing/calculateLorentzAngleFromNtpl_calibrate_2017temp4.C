{
//  gROOT->ProcessLine(".L ./style-CMSTDR.C");
//  gROOT->ProcessLine("setTDRStyle()");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0000);
  gStyle->SetPalette(1,0);


  TFile *f = new TFile("make_ntuple_output/pixelav_events_n28721.root");


  f->cd();
	
  double minfit = 5.;
  double maxfit = 280.;
  double minpq = minfit+10.;
  double maxpq = maxfit-10.;
  TF1 *f1 = new TF1("f1","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x + [7]*x*x*x*x*x*x*x",minfit, maxfit);
//  TF1 *f1 = new TF1("f1","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x",minfit, maxfit);  
  f1->SetParName(0,"offset");
  f1->SetParName(1,"tan#theta_{LA}");
  f1->SetParName(2,"quad term");
  f1->SetParName(3,"cubic term");
  f1->SetParName(4,"quartic term");
  f1->SetParName(5,"quintic term");
  f1->SetParName(6,"sextic term");
  f1->SetParName(7,"septic term");
  f1->SetParameter(0,0);
  f1->SetParameter(1,0.4);
  f1->SetParameter(2,0.0);
  f1->SetParameter(3,0.0);
  f1->SetParameter(4,0.0);
  f1->SetParameter(5,0.0);
  f1->SetParameter(6,0.0);
  f1->SetParameter(7,0.0);
    f1->SetLineColor(2);
  f1->SetLineWidth(3);
  
// Add Linear fit for comparison with previous results
  	
   TF1 *f0 = new TF1("f0","[0] + [1]*x",50., 235.);
   f0->SetParName(0,"offset");
   f0->SetParName(1,"tan#theta_{LA}");
   f0->SetParameter(0,0);
   f0->SetParameter(1,0.4);
   f0->SetLineColor(2);
   f0->SetLineWidth(3);
	
  int hist_drift_ = 285;
  int hist_depth_ = 50;
  double min_drift_ = -1000;
  double max_drift_ = 1000;
  double min_depth_ = 000;
  double max_depth_ = 285;
  double thick_ = 0.0285;
  double ypitch_ = 0.0150;
  double hypitch_ = ypitch_/2.;
  double cotbeta_min = 4.*ypitch_/thick_;

//  ofstream fAngles("Run2012A_modules1to4_layer3.txt", ios::trunc); 

  // 	fLorentzFit << "module" << "\t" << "layer" << "\t" << "offset" << "\t" << "error" << "\t" << "slope" << "\t" << "error" << "\t" "rel.err" << "\t" "pull" << "\t" << "chi2" << "\t" << "prob" << endl;
  TH2F * h_drift_depth_adc = new TH2F("h_drift_depth_adc", "h_drift_depth_adc",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
  TH2F * h_drift_depth_adc2 = new TH2F("h_drift_depth_adc2","h_drift_depth_adc2",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
  TH2F * h_drift_depth_noadc = new TH2F("h_drift_depth_noadc",";drift [#mum];depth [#mum]",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
  TProfile * pq_vs_depth = new TProfile("pq_vs_depth","Avg Pixel Charge;depth [#mum]",hist_depth_, minpq, maxpq,"");
  pq_vs_depth->SetLineColor(2);
  pq_vs_depth->SetMinimum(0.);

  TH1F *h_pt = new TH1F("h_pt","h_pt",100,0,5);
  TH1F *h_ndof = new TH1F("h_ndof","h_ndof",101,0,100);
  TH1I *h_trackQuality = new TH1I("h_trackQuality","h_trackQuality",6,0,6);
  TH1I *h_nHitsPerTrack = new TH1I("h_nHitsPerTrack","h_nHitsPerTrack",100,0,100);
  TH1F *h_qclus = new TH1F("h_qclus","h_qclus",200,0.,50.);
	
   Int_t run_;
   ULong64_t event_;
   Int_t module_;
   Int_t ladder_;
   Int_t layer_;
   Int_t isflipped_;
   Float_t pt_;
   Float_t p_;
   Float_t eta_;
   Float_t phi_;
   Double_t chi2_;
   Double_t ndof_;
   Int_t trackQuality_;
   Int_t nHitsPerTrack_;
   Int_t isHighPurity_;
   const Int_t maxpix = 100;
   
   struct {
      Int_t npix;
      Float_t row[maxpix];
      Float_t col[maxpix];
      Float_t adc[maxpix];
      Float_t x[maxpix];
      Float_t y[maxpix];
   } pixinfo_;
   
   struct {
      Int_t ncol;
      Int_t dcol[maxpix];
      Float_t adc[maxpix];
      Float_t depth[maxpix];
   } colinfo_;
   
   struct {
      Float_t x;
      Float_t y;
      Double_t alpha;
      Double_t beta;
      Double_t gamma;
   } simhit_, trackhit_;
   
   struct {
      Float_t x;
      Float_t y;
      Float_t charge;
      Int_t size_x;
      Int_t size_y;
      Int_t maxPixelCol;
      Int_t maxPixelRow;
      Int_t minPixelCol;
      Int_t minPixelRow;
   } clust_;
   
   struct {
      Float_t x;
      Float_t y;
   } rechit_;
 
	
  // fill the histrograms with the ntpl
  TTree * LATree = (TTree*)f->Get("SiPixelLorentzAngleTree_");
  int nentries = LATree->GetEntries();
  LATree->SetBranchAddress("run", &run_);
  LATree->SetBranchAddress("event", &event_);
  LATree->SetBranchAddress("module", &module_);
  LATree->SetBranchAddress("ladder", &ladder_);
  LATree->SetBranchAddress("layer", &layer_);
  LATree->SetBranchAddress("isflipped", &isflipped_);
  LATree->SetBranchAddress("pt", &pt_);
  LATree->SetBranchAddress("p", &p_);//M
  LATree->SetBranchAddress("eta", &eta_);
  LATree->SetBranchAddress("phi", &phi_);
  LATree->SetBranchAddress("chi2", &chi2_);
  LATree->SetBranchAddress("ndof", &ndof_);
  LATree->SetBranchAddress("trackhit", &trackhit_);
  LATree->SetBranchAddress("simhit", &simhit_);
  LATree->SetBranchAddress("npix", &pixinfo_.npix);
  LATree->SetBranchAddress("rowpix", pixinfo_.row);
  LATree->SetBranchAddress("colpix", pixinfo_.col);
  LATree->SetBranchAddress("adc", pixinfo_.adc);
  LATree->SetBranchAddress("xpix", pixinfo_.x);
  LATree->SetBranchAddress("ypix", pixinfo_.y);
  LATree->SetBranchAddress("clust", &clust_); // charge is given in 1000 e
  LATree->SetBranchAddress("rechit", &rechit_);
//  LATree->SetBranchAddress("trackQuality", &trackQuality_);
//  LATree->SetBranchAddress("isHighPurity", &isHighPurity_);
//  LATree->SetBranchAddress("nHitsPerTrack", &nHitsPerTrack_);
	
  cout << "Running over " << nentries << " hits" << endl;

  //cuts
  float pt_cut = 0.6;
  float clusterSizeY_cut = 2.9999;
  float residual_cut = 0.005;
  float normChi2_cut = 4.0;
  float clusterCharge_cut = 30.; //charge in ke per unit thickness
  float trackQuality_cut = 2.0;
  int highPurity_cut = 1;
	
  for(int ientrie = 0 ; ientrie < nentries; ientrie++){
  //  for(int ientrie = 0 ; ientrie < 20000000; ientrie++){
     LATree->GetEntry(ientrie);


//    if(module_ <=4) continue;
//    if(layer_ !=3) continue;
    //if(run_ != 175990)continue;
    //if(!(module_==1 && layer_==1)) continue;
    //cout<<"module = "<<module_<<endl;
    //cout<<"layer = "<<layer_<<endl;


  
     h_trackQuality->Fill(trackQuality_);
     bool large_pix = false;
     for (int j = 0; j <  pixinfo_.npix; j++){
        int colpos = static_cast<int>(pixinfo_.col[j]+0.01);
        if (pixinfo_.row[j] == 0 || pixinfo_.row[j] == 79 || pixinfo_.row[j] == 80 || pixinfo_.row[j] == 159 || colpos % 52 == 0 || colpos % 52 == 51 ){
	      large_pix = true;	}
     }

    //if(ndof_<10)continue;
    //if(ndof_==0)continue;//because of some crashes
    //if(trackQuality_<trackQuality_cut)continue;
//    if(isflipped_ != 0) continue;
//     if(pt_< pt_cut)continue;
     double residual = TMath::Sqrt( (trackhit_.x - rechit_.x) * (trackhit_.x - rechit_.x) + (trackhit_.y - rechit_.y) * (trackhit_.y - rechit_.y) );

//    if( (clust_.size_y >= clusterSizeY_cut) && (chi2_/ndof_) < normChi2_cut && !large_pix && residual < residual_cut && clust_.charge < clusterCharge_cut){
        
      if(clust_.size_y < clusterSizeY_cut) continue;

      //fAnles->open();
      //write in file (cotan(alpha) cotan(beta) p)
      //cout<<trackhit_.x<<"\t"<<rechit_.x<<endl;
      //cout<<"pt: "<<pt_<<endl;
      //cout<<"p: "<<p_<<endl;
//      fAngles << TMath::Tan(TMath::Pi()/2. - trackhit_.alpha) << "\t" << TMath::Tan(TMath::Pi()/2. - trackhit_.beta) << "\t"<<p_ <<"\t"<<isflipped_<<endl;
        h_pt->Fill(pt_);
        h_ndof->Fill(ndof_);
        h_nHitsPerTrack->Fill(nHitsPerTrack_);
      //h_trackQuality->Fill(trackQuality_);
      //fAngles.close();
       double cotbeta = 1./TMath::Tan(trackhit_.beta);
//       if(fabs(cotbeta) <= cotbeta_min) continue;
       double cotalpha = 1./TMath::Tan(trackhit_.alpha);
       double drdz = sqrt(1.+cotalpha*cotalpha+cotbeta*cotbeta);
       float ccc = clusterCharge_cut*drdz;
       double qclus = clust_.charge/drdz;
       h_qclus->Fill(qclus);
       if(clust_.charge > ccc) continue;
       double drcor = drdz/fabs(cotbeta);
       float ylim1 = trackhit_.y - thick_*cotbeta/2.;
       float ylim2 = trackhit_.y + thick_*cotbeta/2.;
       float xlim1 = trackhit_.x - thick_*cotalpha/2.;  
       colinfo_.ncol = 0;
       int k = 0;
       bool increment_depth = true;
       for (int j = 0; j <  pixinfo_.npix; j++){
           float ypixlow = pixinfo_.y[j]-hypitch_;
           float ypixhigh = pixinfo_.y[j]+hypitch_;
           if(cotbeta > 0.) {
              if(ylim1 > ypixlow) ypixlow = ylim1;
              if(ylim2 < ypixhigh) ypixhigh = ylim2;
           } else {
              if(ylim2 > ypixlow) ypixlow = ylim2;
              if(ylim1 < ypixhigh) ypixhigh = ylim1;
           }           
           float ypixavg = 0.5*(ypixlow+ypixhigh);
           float dypix = fabs(ypixhigh-ypixlow);
	       float dx = (pixinfo_.x[j] - xlim1) * 10000.;
	       float dy = (ypixavg - ylim1) * 10000.;
	       float depth = dy * tan(trackhit_.beta);
	       float drift = dx - dy * tan(trackhit_.gamma);
	       double drpix = 10000. * dypix*drcor;
	       if(drpix < 30.) increment_depth = false;
//	       double qpix = pixinfo_.adc[j]/drpix;
	       double qpix = pixinfo_.adc[j];
	       h_drift_depth_adc->Fill(drift, depth, qpix);
	       h_drift_depth_adc2->Fill(drift, depth, qpix*qpix);
	       h_drift_depth_noadc->Fill(drift, depth);
          int dcol = static_cast<int>(pixinfo_.col[j]-clust_.minPixelCol+0.01);
          for(int i = 0; i < colinfo_.ncol; ++i) {
             if(dcol == colinfo_.dcol[i]) {
                colinfo_.adc[i] += pixinfo_.adc[j]/drpix;
                goto skip;
             }
          }
          if(colinfo_.ncol < maxpix) {
             colinfo_.dcol[colinfo_.ncol] = dcol;
             colinfo_.depth[colinfo_.ncol] = depth;
             colinfo_.adc[colinfo_.ncol] = pixinfo_.adc[j]/drpix;
             ++colinfo_.ncol;
          }
    skip: ++k;
       }
       if(increment_depth) {
          for(int i = 0; i < colinfo_.ncol; ++i) {
             if(colinfo_.depth[i] > minpq && colinfo_.depth[i] < maxpq) pq_vs_depth->Fill(colinfo_.depth[i], colinfo_.adc[i]);
          }
        }
      
//     }
  }
	
  TH1F * h_mean = new TH1F("h_mean",";depth [#mum];drift [#mum]", hist_depth_, min_depth_, max_depth_);
  TH1F * h_drift_depth_adc_slice_ = new TH1F("h_slice","h_slice", hist_drift_, min_drift_, max_drift_);
  //loop over bins in depth (z-local-coordinate) (in order to fit slices)
  for( int i = 1; i <= hist_depth_; i++){
    // 				findMean(i, (i_module + (i_layer - 1) * 8));
    double npix = 0;

    h_drift_depth_adc_slice_->Reset("ICE");
		
    // determine sigma and sigma^2 of the adc counts and average adc counts
    //loop over bins in drift width
    for( int j = 1; j<= hist_drift_; j++){
      if(h_drift_depth_noadc->GetBinContent(j, i) >= 1){
	double adc_error2 = (h_drift_depth_adc2->GetBinContent(j,i) - h_drift_depth_adc->GetBinContent(j,i)*h_drift_depth_adc->GetBinContent(j, i) / h_drift_depth_noadc->GetBinContent(j, i)) /  h_drift_depth_noadc->GetBinContent(j, i);
	h_drift_depth_adc_slice_->SetBinContent(j, h_drift_depth_adc->GetBinContent(j,i));
	h_drift_depth_adc_slice_->SetBinError(j, sqrt(adc_error2));
	npix += h_drift_depth_noadc->GetBinContent(j,i);	
      }else{
	h_drift_depth_adc_slice_->SetBinContent(j, h_drift_depth_adc->GetBinContent(j,i));
	h_drift_depth_adc_slice_->SetBinError(j, 0);
      }
    } // end loop over bins in drift width
			
    double mean = h_drift_depth_adc_slice_->GetMean(1); 
    double error = 0;
    if(npix != 0){
      error = h_drift_depth_adc_slice_->GetRMS(1) / sqrt(npix);
    }
			
    h_mean->SetBinContent(i, mean);
    h_mean->SetBinError(i, error);	
  }// end loop over bins in depth 
	
  TCanvas * c1 = new TCanvas("c1", "c1", 800, 600);
  //c1->Divide(2,1);
  //c1->cd(1);
  // 	h_drift_depth_noadc->
  h_drift_depth_noadc->Draw("colz");
  // 	c1->cd(2);
  // 	h_drift_depth_adc->Draw();
  //c1->cd(2);
  c1->SaveAs("LA_from_Ntpl_output/c1.pdf");

  TCanvas * c11 = new TCanvas("c11", "c11", 800, 600);

  h_mean->Draw();
  // 	c1->cd(4);
	
  h_mean->Fit(f1,"ERQ");
  double p0 = f1->GetParameter(0);
  double e0 = f1->GetParError(0);
  double p1 = f1->GetParameter(1);
  double e1 = f1->GetParError(1);
  double p2 = f1->GetParameter(2);
  double e2 = f1->GetParError(2);
  double p3 = f1->GetParameter(3);
  double e3 = f1->GetParError(3);
  double p4 = f1->GetParameter(4);
  double e4 = f1->GetParError(4);
  double p5 = f1->GetParameter(5);
  double e5 = f1->GetParError(5);
  double p6 = f1->GetParameter(6);
  double e6 = f1->GetParError(6);
  double p7 = f1->GetParameter(7);
  double e7 = f1->GetParError(7);
  double chi2 = f1->GetChisquare();
  double prob = f1->GetProb();	

  c11->SaveAs("LA_from_Ntpl_output/c11.pdf");
  c11->SaveAs("LA_from_Ntpl_output/c11.C");

  TCanvas * c2 = new TCanvas("c2", "c2", 800, 600);
  h_pt->Draw();
  c2->SaveAs("LA_from_Ntpl_output/c2.pdf");

  TCanvas * c3 = new TCanvas("c3", "3", 800, 600);
  h_ndof->Draw();
  c3->SaveAs("LA_from_Ntpl_output/c3.pdf");

  TCanvas * c4 = new TCanvas("c4", "c4", 800, 600);
  h_nHitsPerTrack->Draw();
  c4->SaveAs("LA_from_Ntpl_output/c4.pdf");

  TCanvas * c5 = new TCanvas("c5", "c5", 800, 600);
  h_trackQuality->Draw();
  c5->SaveAs("LA_from_Ntpl_output/c5.pdf");
   
  TCanvas * c6 = new TCanvas("c6", "c6", 800, 600);
   pq_vs_depth->Draw("HIST");
   c6->SaveAs("LA_from_Ntpl_output/c6.C");
   c6->SaveAs("LA_from_Ntpl_output/c6.pdf");
   
  TCanvas * c7 = new TCanvas("c7", "c7", 800, 600);
   h_qclus->Draw("HIST");
   c7->SaveAs("LA_from_Ntpl_output/c7.C");
   c7->SaveAs("LA_from_Ntpl_output/c7.pdf");
  
	
  // 	delete h_mean;
  // 	delete h_drift_depth_adc_slice_;
  cout << "offset" << "\t" << "error" << "\t" << "slope" << "\t" << "error" << "\t" "rel.err" << "\t" << "chi2" << "\t" << "prob" << endl;
  cout  << p0 << "\t" << e0 << "\t" << p1 << "\t" << e1 << "\t" << e1 / p1 *100. << "\t" << chi2 << "\t" << prob << endl;
  cout  << p2 << "\t" << e2 << "\t" << p3 << "\t" << e3 << "\t" << p4 << "\t" << e4  << "\t" << p5 << "\t" << e5 << endl;
  cout  << p6 << "\t" << e6  << "\t" << p7 << "\t" << e7 << endl;

  double cwidth = p1*max_depth_+p2*max_depth_*max_depth_+p3*max_depth_*max_depth_*max_depth_
            +p4*max_depth_*max_depth_*max_depth_*max_depth_
            +p5*max_depth_*max_depth_*max_depth_*max_depth_*max_depth_;
  cout  << "charge width = " << cwidth << endl;
       
  double temp = 263.;
  double vs = 1.53e9*pow(temp, -0.87);
  double ec = 1.01*pow(temp, 1.55);
  double beta = 2.57e-2*pow(temp, 0.66);
  double ibeta = 1./beta;
  double rH = 1.02;
  double arg0 = vs*rH*3.8*1.e-4/ec;
  
   	ofstream fLorentzFit( "LA_from_Ntpl_output/c_lorentzFit.txt", ios::trunc );
  
// make electric field profile from the slope

  double minplot = minfit+5.;
  double maxplot = maxfit - 4.95;

     for(double dep = minplot; dep < maxplot ; dep += 5.) {
         double slope = p1+2.*p2*dep+3.*p3*dep*dep+4.*p4*dep*dep*dep+5.*p5*dep*dep*dep*dep
             +6.*p6*dep*dep*dep*dep*dep+7.*p7*dep*dep*dep*dep*dep*dep;
         double arg1 = pow(arg0/slope, beta);
         double efield = ec*pow(arg1-1., ibeta);
         fLorentzFit << dep << "\t" << efield << endl;  
     }
     
     fLorentzFit.close();
     

  TCanvas * c12 = new TCanvas("c12", "c12", 800, 600);

  h_mean->Draw();
	
  h_mean->Fit(f0,"ERQ");
  double fp0 = f0->GetParameter(0);
  double fe0 = f0->GetParError(0);
  double fp1 = f0->GetParameter(1);
  double fe1 = f0->GetParError(1);
  double fchi2 = f0->GetChisquare();
  double fprob = f0->GetProb();	


  c12->SaveAs("LA_from_Ntpl_output/c12.pdf");
  c12->SaveAs("LA_from_Ntpl_output/c12.C");
  
    cout << "offset" << "\t" << "error" << "\t" << "slope" << "\t" << "error" << "\t" "rel.err" << "\t" << "chi2" << "\t" << "prob" << endl;
  cout  << fp0 << "\t" << fe0 << "\t" << fp1 << "\t" << fe1 << "\t" << fe1 / fp1 *100. << "\t" << fchi2 << "\t" << fprob << endl;


}
