#include <TGraph2D.h>
#include <vector>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <iostream>
#include "RooRealVar.h"
#include "conditions.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "RooExponential.h"
#include "RooNLLVar.h"
#include "RooMinimizer.h"
#include "RooLandau.h"
//#include <RooCrystalBall.h>
#include "RooPolynomial.h"
#include <RooAbsCollection.h>
#include "RooGenericPdf.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include <TMath.h>


using namespace RooFit;
using namespace std;

	char* Round(float f, int d){
		char buf[16];
	      	sprintf(buf, "%.*g", d, f);
		return strdup(buf);
	}

vector<double_t>* fit(TTree* tree[4][6], RooDataSet* data,RooRealVar x, vector<float_t>* conditions, vector<float_t> efficiency){
	
	
	//takes data, plots fit and delivers fitted parameters
	//conditions and efficiency are for plot title 
	//
	//Fitted parameters are in vector<float_t> = [nsig, mean,sigma, nbkg, a0]
		
	RooRealVar mean("mean","mean",5.1,5.4);
	RooRealVar sigma("sigma","sigma",0.005,0.1);	
	RooRealVar width("width","width",0,2);
	RooVoigtian sig("signal","signal",x,mean,width,sigma);
	//RooGaussian sig("signal","signal",x,mean,sigma);
	//RooBreitWigner sig("bw","bw",x,mean,width);

	RooRealVar mean1("mean1","mean1",5.0,5.1);
	RooRealVar sigma1("sigma1","sigma1",0.005,0.05);	
	RooRealVar b0("b0","b0",-10,10);
	RooRealVar b1("b1","b1",-10,10);
	RooRealVar b2("b2","b2",-10,10);
	//RooGaussian bkg_left("bkg_gauss","bkg_gauss",x,mean1,sigma1);
	

	RooRealVar mean2("mean2","mean2",4.9,5.2);
	RooRealVar sigmaL("sigmaL","sigmaL",0.005,0.5);
	RooRealVar sigmaR("sigmaR","sigmaR",0.005,0.1);
	RooRealVar alphaL("alphaL","alphaL",0,5);
	RooRealVar nL("nL","nL",0,10);
	RooRealVar alphaR("alphaR","alphaR",0,5);
	RooRealVar nR("nR","nR",0,10);
	RooCBShape bkg_left("cb_bkg","cb_bkg",x,mean2,sigmaL,alphaL,nL);

	RooRealVar erf_xshift("erf_xshift","erf_xshift",-6,6);
	RooRealVar erf_yshift("erf_yshift","erf_yshift",-6,6);
	RooRealVar erf_xfactor("erf_xfactor","erf_xfactor",-10,10);
	//RooGenericPdf bkg_left("erf", "erf", "TMath::Erf(erf_xfactor*x-erf_xshift)+erf_yshift",RooArgList(x,erf_xshift,erf_yshift,erf_xfactor));

	RooRealVar a0("a0", "a0", -1,1);
	RooRealVar a1("a1", "a1",-10, 10 );
	RooRealVar a2("a2", "a2", -5, 5);
	RooRealVar a3("a3", "a3", 0.01, 2.0);
	RooRealVar a4("a4", "a4", -1, 1.0);
	RooGenericPdf bkg("bkg","Background","exp(a0*x)",RooArgSet(x,a0)) ;

	RooRealVar nsig("N_{sig}","number of signal events", 10000,300000);
	RooRealVar nbkg("N_{bkg}","number of background events",0,1600000) ;
	RooRealVar nbkgleft("N_{bkgleft}","number of gauss background events",100,1600000) ;

	RooAddPdf*  model = new RooAddPdf("model","model",RooArgList(sig,bkg, bkg_left), RooArgList(nsig,nbkg, nbkgleft)) ;

	model->fitTo(*data);

/*
	string title = "cos2d>" + to_string(conditions->at(0)) + " , prob>" + to_string(conditions->at(1)) +
	" ,|dimu_mass-literature|<" + to_string(conditions->at(2)) + " eff. " + to_string(efficiency[1]/efficiency[0]);
	char* p;
	p = &title[0];
*/
	 
	// plot fit and pull
	TCanvas *can = new TCanvas("canv","canv",800,700);

	RooPlot* frame5 = x.frame(Title("Fitted Data")) ;
	RooPlot* Data = data->plotOn(frame5, XErrorSize(1), Binning(80), Name("Data"), MarkerSize(0.6), DrawOption("PZ["));
	RooPlot* modelplot = model->plotOn(frame5, Name("model"),Components("model"));
	RooPlot* sigplot = model->plotOn(frame5, Name("sig"),Components(sig), DrawOption("F"),FillColor(3), FillStyle(3001), LineColor(0));
	RooPlot* bkgplot = model->plotOn(frame5, Name("bkg"),Components(bkg), LineStyle(kDashed), LineColor(6));
	RooPlot* bkg_leftplot = model->plotOn(frame5, Name("bkg_left"),Components(bkg_left), LineStyle(7), LineColor(7));
	data->plotOn(frame5, XErrorSize(1), Binning(80), Name("Data"), MarkerSize(0.6), DrawOption("PZ"));
	

	//sigplot->getAttFill()->SetFillColorAlpha(2, 0.5);
	//sigplot->Draw("SAMEF");
		
	//TH1F* raw = new TH1F("raw","raw",100, 4.7, 6 );
	//tree[0][0]->Draw("b_mass>>raw");
	//raw->SetOption("C"); //in case need raw data on it
	//raw->SetLineStyle(1);
	//raw->SetLineColor(3);
	//raw->Draw("SAMEC");


	// compute the chisquare
	
	RooAbsCollection *flparams = model->getParameters(data)->selectByAttrib("Constant",kFALSE);
	Int_t nflparams = flparams->getSize();
	cout<<nflparams<<"\n";
	Double_t chisquare = -1;
	chisquare = frame5->chiSquare("model", "Data",nflparams);
	cout<<chisquare<<"\n";
        can->Divide(1,2);

        // plot of the curve and the fit
	can->cd(1);
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.03); 
	gPad->SetPad(0.01,0.2,0.99,0.99);
	frame5->GetYaxis()->SetTitleSize(0.04);
	frame5->GetYaxis()->SetTitleOffset(1.2);
	frame5->GetXaxis()->SetLabelSize(0.);
	frame5->GetXaxis()->SetTitleSize(0.);
	frame5->Draw();
	
	
	TLegend* leg = new TLegend(0.5, 0.6, 0.9,0.85);
	leg->AddEntry("Data", "Data points", "p");
	leg->AddEntry("model", "Fit to data", "l");
	leg->AddEntry("sig","Signal","F");
	leg->AddEntry("bkg", "Combinatorial background", "l");
	leg->AddEntry("bkg_left", "CB background", "l");
	leg->SetTextSize(0.04);
	leg->SetTextFont(42);
	leg->Draw();
	
	TPaveText* label_2 = new TPaveText(0.55,0.23,0.8,0.50, "NDC");
	label_2->SetBorderSize(0);
	label_2->SetFillColor(0);
	label_2->SetTextSize(0.041);
	label_2->SetTextFont(42);
	gStyle->SetStripDecimals(kTRUE);
	label_2->SetTextAlign(11);
	Int_t total_pre_mc = efficiency[2];
	Int_t total_post_mc = efficiency[3];
	TString Nsig = to_string(int(round(nsig.getValV())));
	label_2->AddText("N_{sig}^{sel} = "+Nsig);
	Int_t Neff = int(round(nsig.getValV()*13787971./total_post_mc));
	//TString Neffstr = to_string(Neff); 
	TString gammabkg = string(Round(1 - efficiency[5]/efficiency[4],3));
	//label_2->AddText("#gamma_{bkg} = "+gammabkg); 
	TString epsilonMC = string(Round(total_post_mc/13787971,5));
	label_2->AddText("#varepsilon_{MC} = 9.36 x 10^{-4}");
	TString chisquarestr = string(Round(chisquare,3));
	Float_t lum = 0.774; //updated with trigger efficiency
	Float_t BRB = 0.001020;
	Float_t BRjpsi = 0.05961;
	Float_t sigmaB = Neff/BRB/BRjpsi/lum/1000000000;
	TString sigmaBstr = string(Round(sigmaB,5));
	label_2->AddText("#sigma(B^{#pm}) = (518 #pm 12) #mub");
	label_2->AddText("#chi^{2} = "+chisquarestr);
	label_2->AddText("#zeta = 0.386");
	label_2->Draw();

	cout<<"\n Total post mc: " << total_post_mc << "\n";		
	
	TPaveText* cms = new TPaveText(0.14,0.922,0.3,0.93, "NDC");
	cms->AddText("CMS private work");
	cms->SetBorderSize(0);
	cms->SetFillColor(0);
	cms->SetTextSize(0.04);
	cms->SetTextFont(42);
	cms->SetTextAlign(11);
	//cms->SetTextStyle(2);
	cms->Draw();

	TPaveText* lum_g = new TPaveText(0.7,0.92,1,0.93, "NDC");
	lum_g->SetBorderSize(0);
	lum_g->SetFillColor(0);
	lum_g->SetTextSize(0.04);
	lum_g->SetTextFont(42);
	lum_g->SetTextAlign(11);
	char* lumchar = Round(lum,3);
	TString lumstr = string(lumchar);
	lum_g->AddText(lumstr+" fb^{-1}      13 TeV");
	lum_g->Draw();

	
	// plot of the residuals
	
	can->cd(2);
	
	RooHist* hpull = frame5->pullHist("Data", "model");
	RooPlot* frame2 = x.frame(Title("Pull"));
	frame2->addPlotable(hpull,"P");//,"E3");
	frame2->SetMarkerStyle(2);
	frame2->SetMarkerSize(0.01);

	gPad->SetLeftMargin(0.15); 
	gPad->SetPad(0.01,0.01,0.99,0.2);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.5);
	frame2->GetYaxis()->SetNdivisions(202);
	frame2->GetYaxis()->SetRangeUser(-4,4);
	frame2->GetXaxis()->SetTitle("m_{B^{#pm}} [GeV]");
	frame2->GetYaxis()->SetTitle("Pulls");	
	frame2->GetXaxis()->SetTitleSize(0.2);
	frame2->GetYaxis()->SetTitleSize(0.2);
	frame2->GetXaxis()->SetLabelSize(0.15);
	frame2->GetYaxis()->SetLabelSize(0.15);
	frame2->GetXaxis()->SetLabelOffset(0.01);
	frame2->GetYaxis()->SetLabelOffset(0.01);
	frame2->GetYaxis()->SetTitleOffset(0.2);
	frame2->GetXaxis()->SetTickLength(0.1);
	gPad->SetFrameFillColor(0);
	gPad->SetFrameBorderMode(0);
	gPad->SetFrameFillColor(0);
	gPad->SetFrameBorderMode(0);	
	frame2->Draw();
	
	//TLine *line = new TLine();
	//line->DrawLine(rangeMin,0,rangeMax,0);
	//line->SetLineColor(2);
	//line->DrawLine(rangeMin,-3,rangeMax,-3);
	//line->DrawLine(rangeMin,3,rangeMax,3);
				  
	 //save output
	
	string titlestr = "fitvoig.png";
	char* title = &titlestr[0];
	can->SaveAs(title);


	//pull distribution	
	//
	//
	TCanvas* canv_pull = new TCanvas("canv_pull", "canv_pull", 800,700);
	TH1D* hist_pull = new TH1D("hist_pull","hist_pull",130,-5,5);
	
	for(Int_t i = 0; i<hpull->GetN();i++){
		Double_t x, point;
		hpull->GetPoint(i,x,point);
		//if(x<5. || x>6. continue);//boh	
		hist_pull->Fill(point);
	}
	hist_pull->Rebin(4);
	//hist_pull->SetMarkerStyle(4);
	//hist_pull->SetMarkerSize(2);
	TAxis* Xaxis = hist_pull->GetXaxis();
	TAxis* Yaxis = hist_pull->GetYaxis();
	Xaxis->SetTitle("Residual value");
	Xaxis->SetTitleSize(0.045);
	Xaxis->SetLabelSize(0.045);
	Xaxis->SetTitleOffset(1.1);
	Yaxis->SetTitle("Frequency");
	Yaxis->SetTitleSize(0.045);
	Yaxis->SetLabelSize(0.045);
	Yaxis->SetTitleOffset(0.8);
	//gStyle->SetOptStat(0);
	hist_pull->Draw();

	//fit it
	
	TF1* fgauss = new TF1("fgauss","gaus",-5,5);
	fgauss->SetLineColor(4);
	hist_pull->Fit("fgauss");
	fgauss->Draw("same");
	gStyle->SetOptFit(1);
	cout<< fgauss->GetChisquare() / fgauss->GetNDF() ;

	
	TPaveText* cms1 = new TPaveText(0.14,0.922,0.3,0.93, "NDC");
	cms1->AddText("CMS private work");
	cms1->SetBorderSize(0);
	cms1->SetFillColor(0);
	cms1->SetTextSize(0.04);
	cms1->SetTextFont(42);
	cms1->SetTextAlign(11);
	//cms->SetTextStyle(2);
	cms->Draw();
	
	TLegend* leg3 = new TLegend(0.65, 0.7, 0.9,0.85);
	leg3->AddEntry("hist_pull", "Entries", "l");
	leg3->AddEntry("fgauss", "Gaussian fit", "l");
	leg3->SetTextSize(0.04);
	leg3->SetTextFont(42);
//	leg3->Draw();
	

	canv_pull->SaveAs("pull_distr_reb4.png");


	vector<double_t>* params = new vector<double_t> {nsig.getValV(),mean.getValV(),sigma.getValV(),nbkg.getValV(),a0.getValV()};
	return params;
}
	

void plot(){

	Float_t Nvoigcb[1]= {22802};
	Float_t errvoigcb[1] = {158};
	Float_t ndegvoigcb[1] = {11};
	Float_t chivoigcb[1] = {1.04};

	Float_t Nvoiggaus[1] = {20759};
	Float_t errvoiggaus[1] = {387};
	Float_t ndegvoiggaus[1] = {9};
	Float_t chivoiggaus[1] = {4.14};

	Float_t Ngauscb[1] = {18489};
	Float_t errgauscb[1] = {146};
	Float_t ndeggauscb[1] = {10};
	Float_t chigauscb[1] = {2.05};
	
	Float_t Ngausgaus[1] = {19126};
	Float_t errgausgaus[1] = {172};
	Float_t ndeggausgaus[1] = {8};
	Float_t chigausgaus[1] = {4.48};

	Float_t Nbwcb[1] = {25121};
	Float_t errbwcb[1] =  {164}; 
	Float_t ndegbwcb[1] = {10};
	Float_t chibwcb[1] =  {6.93};

	Float_t Nbwgaus[1] = {26211};
	Float_t errbwgaus[1] = {185};
	Float_t ndegbwgaus[1] = {8};
	Float_t chibwgaus[1] = {14.83};
	

	TCanvas* c1 = new TCanvas("c1","c1",700,500);
	
	// create the multigraph
	TGraph *mg1 = new TGraphErrors(1,Nvoigcb,ndegvoigcb,errvoigcb,0);
	mg1->SetMarkerColor(1);
	mg1->GetXaxis()->SetTitle("N_{sig}");
	TGraph *mg2 = new TGraphErrors(1,Nvoiggaus,ndegvoiggaus,errvoiggaus,0);
	mg2->SetMarkerColor(2);
	
	TGraph *mg3 = new TGraphErrors(1,Ngauscb,ndeggauscb,errgauscb,0);
	mg3->SetMarkerColor(3);

	TGraph *mg4 = new TGraphErrors(1,Ngausgaus,ndeggausgaus,errgausgaus,0);
	mg4->SetMarkerColor(4);

	TGraph *mg5 = new TGraphErrors(1,Nbwcb,ndegbwcb,errbwcb,0);
	mg5->SetMarkerColor(5);
	
	TGraph *mg6 = new TGraphErrors(1,Nbwgaus,ndegbwgaus,errbwgaus,0);
	mg6->SetMarkerColor(6);

	mg1->Draw("CP");
	mg2->Draw("SAMEP");
	mg3->Draw("SAMEP");
	mg4->Draw("SAMEP");
	mg5->Draw("SAMEP");
	mg6->Draw("SAMEP");
	
	c1->Update();
	c1->SaveAs("uncertainties.png");

}



RooDataSet* select(TTree* tree[4][6], TTree* treemc, RooRealVar x, vector<float_t>* conditions, vector<float_t>& efficiency){

	// inputs: 
	// -  tree to select from
	// - mc tree for pure signal
	// - mass variable x
	// - vector containing selection conditions
	// - reference to vector in main to set efficiency parameters [total input, selected]
	
	RooDataSet* mass =  new RooDataSet("mass","mass", RooArgSet(x));
	
        float_t b_mass = 0;
	float_t b_cos2d = 0;
	float_t sv_prob = 0;
	float_t dimu_mass = 0;
	float_t k_pt = 0;
        float_t l1_pt = 0;
	float_t l2_pt = 0;
	float_t k_eta = 0;
        float_t l1_eta = 0;
        float_t l2_eta = 0;
	float_t sv_lxy = 0;
	UChar_t hlt_mu9_ip6;

	Int_t total_pre = 0; // initialize parameters for efficiency evaluation
	Int_t total_post = 0;
	Int_t total_pre_mc = 0;
	Int_t total_post_mc = 0;
	Int_t bkg_pre = 0;
	Int_t bkg_post = 0;

	// efficiency: compare amount of events through selection  and bkg rejection = 1 - bkg_post/bkg_pre

	//get conditions
	float_t cond_cos2d = conditions->at(0);
	float_t cond_prob = conditions->at(1);
	float_t cond_dimu_mass = conditions->at(2);	
	float_t cond_k_pt =  conditions->at(3);
	float_t cond_l1_pt = conditions->at(4);
	float_t cond_l2_pt = conditions->at(5);
	float_t cond_sv_lxy = conditions->at(6);
	float_t cond_k_eta = conditions->at(7); // for now use same condition on strongly correlated etas
	float_t cond_l1_eta = conditions->at(8);
	float_t cond_l2_eta = conditions->at(9);
	

	for(int k=0;k<1;k++){
		for(int l=0;l<1;l++){//loop over all trees
			int nEntries = tree[k][l]->GetEntries();
			total_pre +=tree[k][l]->GetEntries("b_mass>5.0 && hlt_mu9_ip6 == 1");//add up total entires in tree [k,l]
			bkg_pre += tree[k][l]->GetEntries("b_mass>5.0 && hlt_mu9_ip6 && (b_mass<5.179 || b_mass> 5.379)");
			
			tree[k][l]->SetBranchAddress("b_mass",&b_mass);
			tree[k][l]->SetBranchAddress("b_cos2d", &b_cos2d);
			tree[k][l]->SetBranchAddress("sv_prob", &sv_prob);
			tree[k][l]->SetBranchAddress("dimu_mass",&dimu_mass);
			tree[k][l]->SetBranchAddress("k_pt",&k_pt);
			tree[k][l]->SetBranchAddress("l1_pt",&l1_pt);
			tree[k][l]->SetBranchAddress("l2_pt",&l2_pt);
			tree[k][l]->SetBranchAddress("sv_lxy",&sv_lxy);
			tree[k][l]->SetBranchAddress("k_eta",&k_eta);
			tree[k][l]->SetBranchAddress("l1_eta",&l1_eta);
			tree[k][l]->SetBranchAddress("l2_eta",&l2_eta);
			tree[k][l]->SetBranchAddress("hlt_mu9_ip6",&hlt_mu9_ip6);
	
			for(int i=0; i<nEntries; i++){
				tree[k][l]->GetEntry(i);//set entries in corresponding variable
				if (    b_cos2d>cond_cos2d && 
					hlt_mu9_ip6 == 1 &&
					sv_prob>cond_prob &&
				 	abs(3.097- dimu_mass)<cond_dimu_mass && 
					b_mass>5.0 &&
					k_pt > cond_k_pt &&
					//l1_pt < cond_l1_pt &&
					l2_pt > cond_l2_pt &&
					sv_lxy > cond_sv_lxy &&//needs to be > apparently
					abs(k_eta) < cond_k_eta &&
					abs(l1_eta) < cond_l1_eta &&
					abs(l2_eta) < cond_l2_eta ){

					//conditions to add entry to mass set
					x.setVal(b_mass);
					mass->add(RooArgSet(x));
					total_post++;
					if(b_mass<5.179 || b_mass>5.379) {
						bkg_post++;
					}	
				}
			}
			if (k>1 && l==4){break;}
		}
	}


	// measure efficiency on mc
	
	Int_t ismatched = 0;
	treemc->SetBranchAddress("ismatched",&ismatched);
	
	treemc->SetBranchAddress("b_mass",&b_mass);
	treemc->SetBranchAddress("b_cos2d", &b_cos2d);
	treemc->SetBranchAddress("sv_prob", &sv_prob);
	treemc->SetBranchAddress("dimu_mass",&dimu_mass);
	treemc->SetBranchAddress("k_pt",&k_pt);
	treemc->SetBranchAddress("l1_pt",&l1_pt);
	treemc->SetBranchAddress("l2_pt",&l2_pt);
	treemc->SetBranchAddress("sv_lxy",&sv_lxy);
	treemc->SetBranchAddress("k_eta",&k_eta);
	treemc->SetBranchAddress("l1_eta",&l1_eta);
	treemc->SetBranchAddress("l2_eta",&l2_eta);
	treemc->SetBranchAddress("hlt_mu9_ip6",&hlt_mu9_ip6);

	Int_t mc_entries = treemc->GetEntries();
	total_pre_mc = treemc->GetEntries("b_mass>5.0 && hlt_mu9_ip6 && ismatched==1");
        
	for(int i=0; i<mc_entries; i++){
        	treemc->GetEntry(i);//set entries in corresponding variable
		if (ismatched &&
		    hlt_mu9_ip6 == 1 &&
		    b_mass>5.0 &&
		    b_cos2d>cond_cos2d &&
		    sv_prob>cond_prob &&
		    abs(3.097- dimu_mass)<cond_dimu_mass &&
		    k_pt > cond_k_pt &&
		   // l1_pt < cond_l1_pt &&
		    l2_pt > cond_l2_pt &&
		    sv_lxy > cond_sv_lxy &&//needs to be > apparently
		    abs(k_eta) < cond_k_eta &&
		    abs(l1_eta) < cond_l1_eta &&
		    abs(l2_eta) < cond_l2_eta ){
	                    //conditions to mc signal passing the cut
		            total_post_mc++;
		}
	}
	
	gROOT->LoadMacro("/work/mratti/CMS_style/tdrstyle.C");
	gROOT->ProcessLine("setTDRStyle()");
	

	bool sig_bkg = 0;

	if(sig_bkg){

	//Plot sig_bkg		
	TCanvas *can = new TCanvas("canv","canv",700,500);

	TH1F* bkg = new TH1F("bkg","bkg",100,0,0.7);
	tree[0][0]->Draw("sv_lxy>>bkg", "(b_mass<5.179||b_mass>5.379) &&  hlt_mu9_ip6 == 1");
	bkg->SetLineColor(1);
	bkg->SetLineStyle(0);
	bkg->SetLineWidth(2);

	TH1F* bkg_filled = new TH1F(*bkg);
	bkg_filled->SetAxisRange(0,0.035);//update with range to be excluded
	bkg_filled->SetLineWidth(0);
	bkg_filled->SetFillColorAlpha(3,0.5);
/*
	TH1F* bkg_filled2 = new TH1F(*bkg);
	bkg_filled2->SetAxisRange(3.247,5);//update with range to be excluded
	bkg_filled2->SetLineWidth(0);
	bkg_filled2->SetFillColorAlpha(3,0.5);
*/
	TH1F* signal = new TH1F("signal","signal",100, 0,0.7);
	treemc->Draw("sv_lxy>>signal", "hlt_mu9_ip6 == 1");
	signal->SetLineColor(4);
	signal->SetLineStyle(0);
	signal->SetLineWidth(2);
	signal->GetXaxis()->SetLabelSize(0.05);
	signal->GetYaxis()->SetLabelSize(0.05);
	signal->GetXaxis()->SetTitleSize(0.05);
	signal->GetYaxis()->SetTitleSize(0.05);
	signal->GetYaxis()->SetTitle("Entries - Normalized scale");	
	signal->GetXaxis()->SetTitle("l_{xy} [cm]");	
	signal->GetXaxis()->SetLabelOffset(0.015);	
	signal->GetXaxis()->SetTitleOffset(1.1);
	signal->GetYaxis()->SetLabelOffset(0.01);
	signal->GetYaxis()->SetTitleOffset(1.3);
	//signal->GetXaxis()->SetNdivisions(505);
	//signal->GetXaxis()->SetRangeUser(0.999,1);
	signal->GetYaxis()->SetRangeUser(0,1700);
	//gPad->SetTicky(0);

	TH1F* sig_filled = new TH1F(*signal);
	sig_filled->SetAxisRange(0.,0.035);//update with range to be excluded
	sig_filled->SetFillColorAlpha(1,0.3);
	sig_filled->SetLineWidth(0);
/*
	TH1F* sig_filled2 = new TH1F(*signal);
	sig_filled2->SetAxisRange(3.247,5);//update with range to be excluded
	sig_filled2->SetFillColorAlpha(1,0.3);
	sig_filled2->SetLineWidth(0);
*/
	signal->DrawNormalized("",1);
	bkg->DrawNormalized("SAME",1);
	bkg_filled->DrawNormalized("SAME",1);
	sig_filled->DrawNormalized("SAME",1);
	//bkg_filled2->DrawNormalized("SAME",1);
	//sig_filled2->DrawNormalized("SAME",1);
	
	TLegend* leg = new TLegend(0.5, 0.7, 0.76,0.9);
	leg->AddEntry(bkg,"Background","l");
	leg->AddEntry(signal,"MC signal","l");
	leg->AddEntry(bkg_filled, "Cut background", "f");
	leg->AddEntry(sig_filled, "Cut signal", "f");
	leg->SetTextFont(42);
	leg->SetTextSize(0.04);
	leg->Draw();   

	TPaveText* label_2 = new TPaveText(0.5,0.5,0.65,0.65, "NDC");
	label_2->SetBorderSize(0);
	label_2->SetFillColor(0);
	label_2->SetTextSize(0.041);
	label_2->SetTextFont(42);
	label_2->SetTextAlign(11);
	label_2->AddText("Cut: l_{xy} > 0.035 GeV");
	char* se = Round((float)total_post_mc/total_pre_mc, 3);
	TString sig_eff = TString(se);
	label_2->AddText("Signal efficiency = "+sig_eff);
	char* br = Round(1.0 - (float)bkg_post/bkg_pre,3);
	TString bkg_rej = TString(br);
	label_2->AddText("Background rejection = "+bkg_rej);
	label_2->Draw();
			
	
	// AGGIORNA CONDITION SOPRA; 
	// nome file!!!
	// cambia label axis;
	// label cut
	// range di hist e di fill
	//  location label/legend; binning?; variabile da importare in hists
	
	TPaveText* cms = new TPaveText(0.14,0.92,0.3,0.93, "NDC");
	cms->AddText("CMS");
	cms->SetBorderSize(0);
	cms->SetFillColor(0);
	cms->SetTextSize(0.04);
	cms->SetTextFont(42);
	cms->SetTextAlign(11);
	//cms->SetTextStyle(2);
	cms->Draw();

	TPaveText* lum_g = new TPaveText(0.7,0.92,1,0.93, "NDC");
	lum_g->SetBorderSize(0);
	lum_g->SetFillColor(0);
	lum_g->SetTextSize(0.04);
	lum_g->SetTextFont(42);
	lum_g->SetTextAlign(11);
	lum_g->AddText("0.774 fb^{-1}     13 TeV");
	lum_g->Draw();

	
	string sig_bkg1 = "new_sig_bkg_lxy.png";
	char* sb1 = &sig_bkg1[0];
	can->SaveAs(sb1);
	}
	
	efficiency[0] = total_pre;// alternative output in main
	efficiency[1] = total_post;
	efficiency[2] = total_pre_mc;
	efficiency[3] = total_post_mc;
	efficiency[4] = bkg_pre;
	efficiency[5] = bkg_post;
	

	cout<< "Total amount of entries: " << total_pre<<"\n";
	cout << "Selected amount of entries: " << total_post <<"\n";
	cout << "bkg_pre: "<< bkg_pre<<"\n";
	cout << "bkg_post: " << bkg_post << "\n";
	cout << "Efficiency: " << 100.*total_post_mc/total_pre_mc <<"% \n \n \n ";


	return mass;
}


void plot_efficiency(TTree* t[4][6], TTree* tmc, RooRealVar x, vector<float_t>& efficiency){
	// set standard cuts. then vary each cut separately, and extract two essential values: Nsig/Nbkg as in the fit, and bkg efficiency
	// plot all in a single plot. 	
	// input: set of cuts
	// output: generate plot
	
	int n_cond = 10; //amount of different values in each condition, e.g. we set cos2d = [0.90,0.95,0.99,0.995,0.999]

	float_t signal_efficiencies[10][n_cond];
	float_t bkg_rejection[10][n_cond];

	for(int i=0;i<10;i++){
		if (i==4) continue;//skip l1
		for(int j=0;j<n_cond;j++){
			vector<float_t>* cond = cond_func();
			cond->at(i) = conditions[i][j]; // change current condition to evaluate
			RooDataSet* mass = select(t,tmc, x,cond,efficiency);//performs selection and overwrites efficiency pars
			//vector<double_t>* params = fit(mass,x,cond,efficiency,title);
			// note that bkg efficiency does not depend on the fit, but signal eff does (not anymore)
			
			signal_efficiencies[i][j] = efficiency[3] /efficiency[2];//Nsig/Nbkg
			bkg_rejection[i][j] = 1-efficiency[5]/efficiency[4];//1 - bkg_post/bkg_pre
		}
	}

	TCanvas* eff = new TCanvas("eff","eff",800,600);

	TGraph* cos2d = new TGraph(n_cond,signal_efficiencies[0],bkg_rejection[0]);
	cos2d->SetMarkerStyle(22);
//	cos2d->GetXaxis()->SetRangeUser(0.9,1);
	cos2d->SetMarkerColor(8);
	cos2d->SetLineColor(8);
	cos2d->Draw("ACP");

	TGraph* prob = new TGraph(n_cond,signal_efficiencies[1],bkg_rejection[1] );
	prob->SetMarkerStyle(22);
	prob->SetLineColor(1);
	prob->Draw("SAMECP");
	
	TGraph* dimu_mass = new TGraph(n_cond,signal_efficiencies[2],bkg_rejection[2]);
	dimu_mass->SetMarkerStyle(22);
	dimu_mass->SetMarkerColor(3);
	dimu_mass->SetLineColor(3);
	dimu_mass->Draw("SAMECP");
	
	TGraph* k_pt = new TGraph(n_cond,signal_efficiencies[3],bkg_rejection[3]);
	k_pt->SetMarkerStyle(22);
	k_pt->SetMarkerColor(2);
	k_pt->SetLineColor(2);
	k_pt->Draw("SAMECP");

	TGraph* l1_pt = new TGraph(n_cond,signal_efficiencies[4],bkg_rejection[4] );
	l1_pt->SetMarkerStyle(22);
	l1_pt->SetMarkerColor(4);
	l1_pt->SetLineColor(4);
	//l1_pt->Draw("SAMECP");

	TGraph* l2_pt = new TGraph(n_cond,signal_efficiencies[5], bkg_rejection[5]);
	l2_pt->SetMarkerStyle(22);
	l2_pt->SetMarkerColor(5);
	l2_pt->SetLineColor(5);
	l2_pt->Draw("SAMECP");
	
	TGraph* sv_lxy  = new TGraph(n_cond,signal_efficiencies[6],bkg_rejection[6]);
	sv_lxy->SetMarkerStyle(22);
	sv_lxy->SetMarkerColor(6);	
	sv_lxy->SetLineColor(6);	
	sv_lxy->Draw("SAMECP");
	
	TGraph* k_eta = new TGraph(n_cond,signal_efficiencies[7],bkg_rejection[7]);
	k_eta->SetMarkerStyle(22);
	k_eta->SetMarkerColor(7);
	k_eta->SetLineColor(7);
	k_eta->Draw("SAMECP");

	TGraph* l1_eta = new TGraph(n_cond,signal_efficiencies[8],bkg_rejection[8]);
	l1_eta->SetMarkerStyle(22);
	l1_eta->SetMarkerColor(9);
	l1_eta->SetLineColor(9);
	l1_eta->Draw("SAMECP");
	
	TGraph* l2_eta = new TGraph(n_cond,signal_efficiencies[9],bkg_rejection[9]);
	l2_eta->SetMarkerStyle(22);
	l2_eta->SetMarkerColor(12);
	l2_eta->SetLineColor(12);
	l2_eta->Draw("SAMECP");

	//gPad->SetLogx();
	//gPad->SetLogy();	
	
	cos2d->GetXaxis()->SetTitle("Signal efficiency (N_{post}^{MC}/ N_{pre}^{MC})");
	cos2d->GetYaxis()->SetTitle("Background rejection (1 - N_{post}^{bkg}/ N_{pre}^{bkg})");
	cos2d->GetXaxis()->SetTitleOffset(1.2);
	cos2d->GetXaxis()->SetLabelOffset(0.01);
	cos2d->GetYaxis()->SetTitleOffset(1.2);
	cos2d->GetXaxis()->SetTitleSize(0.04);
	cos2d->GetYaxis()->SetTitleSize(0.04);
	cos2d->GetXaxis()->SetLabelSize(0.04);
	cos2d->GetYaxis()->SetLabelSize(0.04);
	
	//gStyle->SetPadLeftMargin(0.25);

	TLegend* leg = new TLegend(0.2,0.17,0.32,0.67);
	leg->AddEntry(cos2d,"cos(#theta)","lp");
	leg->AddEntry(prob,"#rho","lp");
	leg->AddEntry(dimu_mass,"m_{#mu#mu}","lp");
	leg->AddEntry(k_pt,"p_{T,K}","lp");
	//leg->AddEntry(l1_pt,"l1_pt","lp");
	leg->AddEntry(l2_pt,"p_{T,2}","lp");
	leg->AddEntry(sv_lxy,"l_{xy}","lp");
	leg->AddEntry(k_eta,"#eta_{K}","lp");
	leg->AddEntry(l1_eta,"#eta_{1}","lp");
	leg->AddEntry(l2_eta,"#eta_{2}","lp");
        leg->SetTextFont(42);
	leg->SetTextSize(0.04);
	leg->Draw();
	
	TPaveText* cms = new TPaveText(0.14,0.922,0.3,0.93, "NDC");
	cms->AddText("CMS private work");
	cms->SetBorderSize(0);
	cms->SetFillColor(0);
	cms->SetTextSize(0.04);
	cms->SetTextFont(42);
	cms->SetTextAlign(11);
	cms->Draw();
	
	TPaveText* lum_g = new TPaveText(0.7,0.92,1,0.93, "NDC");
	lum_g->SetBorderSize(0);
	lum_g->SetFillColor(0);
	lum_g->SetTextSize(0.04);
	lum_g->SetTextFont(42);
	lum_g->SetTextAlign(11);
	lum_g->AddText("0.774 fb^{-1}     13 TeV");
	lum_g->Draw();
								
		
	eff->SaveAs("efficiency_zoom_nopt1_final.png");


}		


void experiment(){
	


	TFile* f[4][6];
	TTree* t[4][6];
	
	f[0][0] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_A1.root");
	f[0][1] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_A2.root");
	f[0][2] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_A3.root");
	f[0][3] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_A4.root");
	f[0][4] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_A5.root");
	f[0][5] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_A6.root");
	
	
	f[1][0] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_B1.root");
	f[1][1] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_B2.root");
	f[1][2] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_B3.root");
	f[1][3] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_B4.root");
	f[1][4] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_B5.root");
	f[1][5] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_B6.root");

	f[2][0] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_C1.root");
	f[2][1] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_C2.root");
	f[2][2] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_C3.root");
	f[2][3] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_C4.root");
	f[2][4] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_C5.root");
	f[2][5] = nullptr;
	
	f[3][0] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_D1.root");
	f[3][1] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_D2.root");
	f[3][2] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_D3.root");
	f[3][3] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_D4.root");
	f[3][4] = TFile::Open("/scratch/anlyon/BHNL/samples/data/flat_bparknano_D5.root");
	f[3][5] = nullptr;
	
	for(int k = 0; k<4 ; k++){
  		for(int l=0; l<6; l++){
    			t[k][l] = (TTree*)f[k][l]->Get("control_tree");
			if (k>1 && l==4){break;}
  		}
 	}

	TFile* fmc;
	TTree* tmc;
	
	fmc = TFile::Open("/scratch/anlyon/BHNL/samples/mc/flat_bparknano.root");
	tmc = (TTree*)fmc->Get("control_tree");	

	vector<float_t>* conditions = cond_func_fine();
	RooRealVar x("x","x",5.0,6.0);
	
	vector<float_t> efficiency(6);//alternative output. (Ntotal before cuts, N after cuts,bkg before cuts, bkg after cuts)
	
	
	//to generate efficiency plot
	//plot_efficiency(t,tmc, x,efficiency);
	
	//to generate standard fit
	
	string sig_bkg = "sig_bkg_dimumass.png";
	char* sb = &sig_bkg[0];
	
	RooDataSet* mass = select(t,tmc,x,conditions, efficiency);
	vector<double_t>* params = fit(t,mass,x,conditions,efficiency);
	//plot();
}
