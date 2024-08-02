//template_utils.h

        auto vav_pars1 = get_vavilov_pars(hp);


//vavilov distribution (to fit to)
Double_t vavilov(Double_t *v, Double_t *par)
{
    Float_t arg = 0.;
    if (par[2] != 0) arg = (float)(v[0] - par[1])/fabs(par[2]);

    float beta2 = 1.f;   
    float kappa = (float)par[3];
    VVIObjF vvidist(kappa, beta2, 0);
    float xl, xu;
    vvidist.limits(xl,xu);
    if(arg < xl) arg = xl;
    if(arg > xu) arg = xu;
    Double_t fitval = par[0]*(double)vvidist.fcn(arg);
    return fitval;
}



std::vector<float> get_vavilov_pars(TH1F *h){

    TF1 *vfunc = new TF1("vavilov",vavilov,-10.,50.,4);
    vfunc->SetParNames("norm","mean","sigma","kappa");

    float par_guess [4];
    float mean = h->GetMean();
    par_guess[0] = h->GetEntries()*20000./mean; 
    par_guess[1] = mean;
    par_guess[2] = 0.1*mean;
    par_guess[3] = 0.02*mean/20000.;
    vfunc->SetParameters(par_guess[0], par_guess[1], par_guess[2], par_guess[3]);
    vfunc->SetParLimits(1, 0., 1e9);
    vfunc->SetParLimits(2, 0., 1e9);
    vfunc->SetParLimits(3, 0.01, 10.);
    h->Fit("vavilov");
    std::vector<float> pars;
    //we don't care about norm
    for(int i=1; i<4; i++){
        float par = vfunc->GetParameter(i);
        if((i==1 || i == 2) && par < 1e-4) par = par_guess[i];
        pars.push_back(par);
    }
    return pars;

}

