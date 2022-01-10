#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"


// For cuts
enum { kAll = 1, kET, kSiBunch, kTrigger, kOppTrack, kTheta, kx0y0, kFiducial, kMax };

// For triggers
enum { t_CPT2 = 1, t_CPT2noBBCL, t_SDT, t_ET, t_CPTnoBBCL, t_Max };

// Fiducial cuts                                                                                                                                                                     
const Int_t NFB = 50 ; // no. of fiducial bins                                                                                                                                       
Double_t fidcut[8][NFB] ;
// scanned from the other direction ( vertical or horizontal )                                    
Double_t fidneg[8][NFB] ;
Double_t fidpos[8][NFB] ;


const Int_t nsigma = 3 ;
//const Int_t nsigma = 2 ;

const Int_t ERRCODE = -999 ;

const Int_t MaxBunch = 120 ;

short BlueSpin[MaxBunch], YellowSpin[MaxBunch] ;

Double_t beam_momentum =  254.865736997 ;

const short orientations[32] = {1,1,1,1,  -1,-1,-1,-1,  1,1,1,1, -1,-1,-1,-1,  1,-1,1,-1,  -1,1,-1,1,  1,-1,1,-1, -1,1,-1,1 };


// 2017 : added by K. Yip (Mar. 5, 2018)
const double RP_POS_Z[ 32 ] = { 
  -15.76648, -15.77318, -15.77988, -15.78658,   
  -15.76819, -15.77489, -15.78159, -15.78829,  
  -17.56672, -17.57342, -17.58012, -17.58682,  
  -17.56635, -17.57305, -17.57975, -17.58645,   
  15.76805,  15.77475,  15.78145,   15.78815,  
  15.76899,  15.77569,  15.78239,   15.78909,
  17.56906,  17.57576,  17.58246,   17.58916,
  17.56841,  17.57511,  17.58181,   17.58851
};


const double Emin[8][4] = {  {20, 28, 50, 66}, {20, 28, 50, 60}, {20, 28, 50, 60}, {20, 28, 45, 60}, 
			     {20, 28, 50, 60}, {20, 26, 50, 60}, {20, 28, 50, 60}, {20, 28, 50, 60} } ;

const int MaxClSize = 5 ;

const double CutDiff2 = 0.0015 ;
//const double CutDiff = 0.0005 ;
const double CutDiff = 0.0003 ;
//  const double CutDiff = 0.00025 ;
// const double CutDiff = 0.0002 ;

/*
class myTrackPt {

  int rpid ; // 0 ... 7 ( E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D
  int planeid ; // 0 ... 3 ( A, B, C, D )
  

} ;
*/


TH1D *hx0, *hy0, *hz0 ;
TH1D *hx02, *hy02, *hz02 ;
TH1D *diff_thetax, *diff_thetay ;
TH1D *hthetay1, *hthetay2 ;

TH1D *hy5, *hy7, *hdiff_y5_y7 ;

TH1D *hm_t, *hphi, *hm_t_EUWD, *hm_t_EDWU ;
TH1D *hm_tv, *hm_tv_EUWD, *hm_tv_EDWU ; // variable bins

TH2D *ht_vs_phi ;

TH1D *hdiff_from_fit_thetax,  *hdiff_from_fit_thetay ;

TH2D *E1U_W1D_xy, *E2U_W2D_xy, *E1D_W1U_xy, *E2D_W2U_xy ;
TH2D *E1U_E1D_xy, *E2U_E2D_xy, *W1U_W1D_xy, *W2U_W2D_xy ;
TH2D *E1U_E1D, *E2U_E2D, *W1U_W1D, *W2U_W2D ;

TH1D *HowManyPlanes ;

enum { intv = 6, deg_delta = 10, dim_phi = 90/deg_delta } ; // 10 degrees

TH1D *hUU[intv], *hDD[intv] ; /// for counting Up-Up and Down-Down entries
TH1D *hUD[intv], *hDU[intv] ; /// for counting Up-Down and Down-Up entries

TH1D *raw_asym[intv] ; /// raw_asymmetries 
TH1D *false_asym[intv] ; /// false raw_asymmetries 

TH1I *hAnalysisFlow ; // Cuts Flow

TH1I *hTriggers ; // Triggers

TH1D *htrig_ADC, *htrig_TAC ; // ADC/TAC in RP trigger PMT

TH1D *hTOFHits_C ; // No. of TOF hits for CPTnoBBCL
TH1D *hTOFHits_O ; // No. of TOF hits for non-CPTnoBBCL events

TH1D *hTOFMult_C ; // No. of TOF Mult (trigger) for CPTnoBBCL
TH1D *hTOFMult_O ; // No. of TOF Mult (trigger) for non-CPTnoBBCL events


TH2D *hTOFHits_vs_Mult ; // TOF Hits vs TOF Trigger Multiplicity


const double div1 = 0.2, div2 = 0.35, div3 = 1.0, div4 = 10. ;

const Double_t LowerBoundY = -0.077 ;
const Double_t UpperBoundY =  0.079 ;
const Double_t FidScaleY = 50/0.156 ;  

const Double_t LowerBoundX = -0.041 ;
const Double_t UpperBoundX =  0.065 ;
const Double_t FidScaleX = 50/0.106 ;  

// To be used for E1D(1), E2D(3), W1D(5), W2D(7)
Bool_t FiducialNeg(Int_t RPId, Double_t findbin, Double_t findbound) {

  //  return kTRUE ;
  //
  if ( findbin < LowerBoundX || findbin > UpperBoundX ) 
    return kFALSE ;
  else if ( findbound < LowerBoundY || findbound > UpperBoundY ) 
    return kFALSE ;
  else { 
    return ( ( findbound >= fidcut[RPId][ int((findbin-LowerBoundX)*FidScaleX) ] ) &&
	     ( findbin >= fidneg[RPId/2][ int((findbound-LowerBoundY)*FidScaleY) ] ) &&  // for the other direction, x & y (ie. findbin 
	     ( findbin <= fidpos[RPId/2][ int((findbound-LowerBoundY)*FidScaleY) ] ) ) ; // & findbound ) relationship interchanged
  }
  //

}

// To be used for E1U(0), E2U(2), W1U(4), W2U(6)
Bool_t FiducialPos(Int_t RPId, Double_t findbin, Double_t findbound) {

  //  return kTRUE ;
  //
  if ( findbin < LowerBoundX || findbin > UpperBoundX ) 
    return kFALSE ;
  else if ( findbound < LowerBoundY || findbound > UpperBoundY ) 
    return kFALSE ;
  else {
    return ( ( findbound <= fidcut[RPId][ int((findbin-LowerBoundX)*FidScaleX) ] ) &&
	     ( findbin >= fidneg[RPId/2][ int((findbound-LowerBoundY)*FidScaleY) ] ) &&  // for the other direction, x & y (ie. findbin 
	     ( findbin <= fidpos[RPId/2][ int((findbound-LowerBoundY)*FidScaleY) ] ) ) ; // & findbound ) relationship interchanged
  }
  //

}


Int_t read_fiducial() {

  ifstream icuts ;
  Int_t isec, lines ;
  string linestring ;

  // fiducial cuts

  const Int_t NPairs = 4 ;
  //  const string filenames[NPairs] = { ".edge.EHI_EHO", ".edge.EVD_EVU", ".edge.WHI_WHO" , ".edge.WVD_WVU" } ; 
  //  const Int_t NegSide[NPairs] = { 0, 3, 4, 6 } ;
  //  const Int_t PosSide[NPairs] = { 1, 2, 5, 7 } ;
  const string filenames[NPairs] = { ".edge.E1U_E1D", ".edge.E2U_E2D", ".edge.W1U_W1D" , ".edge.W2U_W2D" } ; 
  const Int_t NegSide[NPairs] = { 1, 3, 5, 7 } ;
  const Int_t PosSide[NPairs] = { 0, 2, 4, 6 } ;

  // Just recycling the variable "isec"
  for ( isec = 0 ; isec<NPairs ; isec++ ) {

    icuts.open( (filenames[isec]).c_str(), ifstream::in );

    if ( ! icuts.is_open() ) {
      cerr << "Error opening file " << filenames[isec] << " ! " << endl ;
      return -1 ;
    }

    lines = 0 ;
    while ( icuts.good() && !icuts.eof() && icuts.peek() != EOF ) {
      getline(icuts, linestring,'\n');
      istringstream linestore(linestring);

      linestore >> fidcut[NegSide[isec]][lines] >> fidcut[PosSide[isec]][lines] ;

      cout << filenames[isec] << " : " << lines << " : " << fidcut[NegSide[isec]][lines] << " , " << fidcut[PosSide[isec]][lines] << endl ;
    
      lines++ ;

      if ( lines > NFB ) { 
        cerr << filenames[isec] << " has more than " << NFB << " lines ! " << endl ;
        return -1 ;
      }
    }

    icuts.close();

  }

  // Instead of scanning over y, we scan over x here.
  const string filenames2[NPairs] = { ".edge.E1U_E1D.projx", ".edge.E2U_E2D.projx", ".edge.W1U_W1D.projx" , ".edge.W2U_W2D.projx" } ; 

  // Just recycling the variable "isec"
  for ( isec = 0 ; isec<NPairs ; isec++ ) {

    icuts.open( (filenames2[isec]).c_str(), ifstream::in );

    if ( ! icuts.is_open() ) {
      cerr << "Error opening file " << filenames2[isec] << " ! " << endl ;
      return -1 ;
    }

    lines = 0 ;
    while ( icuts.good() && !icuts.eof() && icuts.peek() != EOF ) {
      getline(icuts, linestring,'\n');
      istringstream linestore(linestring);

      linestore >> fidneg[isec][lines] >> fidpos[isec][lines] ;

      cout << filenames2[isec] << " : " << lines << " : " << fidneg[isec][lines] << " , " << fidpos[isec][lines] << endl ;
    
      lines++ ;

      if ( lines > NFB ) { 
	cerr << filenames2[isec] << " has more than " << NFB << " lines ! " << endl ;
	return -1 ;
      }
    }

    icuts.close();

  }

  return 0 ;

}


void make_histograms() {

  //  char RPPlane[4][2] = { "A", "B", "C", "D" } ;

  hx0 = new TH1D("hx0", "x_{0}(W) - x_{0}(E)", 100, -0.02, 0.02);
  hy0 = new TH1D("hy0", "y_{0}(W) - y_{0}(E)", 100, -0.02, 0.02);
  hz0 = new TH1D("hz0", "Vertices of z_{IP}", 100, -10.0, 10.0);

  hx02 = new TH1D("hx02", "Vertices of x_{IP}", 100, -0.01, 0.01);
  hy02 = new TH1D("hy02", "Vertices of y_{IP}", 100, -0.01, 0.01);
  hz02 = new TH1D("hz02", "Vertices of z_{IP}", 100, -10.0, 10.0);

  diff_thetax = new TH1D("diff_thetax", "Difference of theta_x", 100, -0.002, 0.002);
  diff_thetay = new TH1D("diff_thetay", "Difference of theta_y", 100, -0.002, 0.002);

  hthetay1 = new TH1D("hthetay1", "theta_y1", 100, -0.01, 0.01);
  hthetay2 = new TH1D("hthetay2", "theta_y2", 100, -0.01, 0.01);

  hy5 = new TH1D("hy5", "y coordinates of W1D", 60, -0.08, -0.02); 
  hy7 = new TH1D("hy7", "y coordinates of W2D", 60, -0.075, -0.015); 
  hdiff_y5_y7 = new TH1D("hdiff_y5_y7", "diff. of y coordinates of W1D & W2D", 100, -0.01, 0.01); 

  hm_t = new TH1D("hm_t","-t", 50, 0., 1.) ;
  hm_t_EUWD = new TH1D("hm_t_EUWD","-t for EU-WD", 50, 0., 1.) ;
  hm_t_EDWU = new TH1D("hm_t_EDWU","-t for ED-WU", 50, 0., 1.) ;

  // Using UA4 variable -t bins
  const Int_t NBINS = 36;
  Double_t UA4edges[NBINS + 1] = {
    0.19,0.21,0.23,0.25,0.27,0.29,0.31,0.33,0.35,0.37,0.39,0.41,0.43,
    0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.65,0.67,0.69,
    0.71,0.75,0.79,0.83,0.87,0.91,0.95,0.99,1.03,1.07,1.11
  };
  // Bin 1 corresponds to range [0.19, 0.21]
  // Bin 2 corresponds to range [0.21, 0.23] etc.
  // The last ten bins are size 0.04

  hm_tv = new TH1D("hm_tv","-t with UA4 binning", NBINS, UA4edges) ;
  hm_tv_EUWD = new TH1D("hm_tv_EUWD","-t with UA4 binning for EU-WD", NBINS, UA4edges) ;
  hm_tv_EDWU = new TH1D("hm_tv_EDWU","-t with UA4 binning for ED-WU", NBINS, UA4edges) ;

  hphi = new TH1D("hphi","#phi", 100., -180., 180.) ;
  ht_vs_phi = new TH2D("ht_vs_phi","-t vs #phi", 100, -180., 180., 100, 0., 1.) ;

  hdiff_from_fit_thetax = new TH1D("hdiff_from_fit_thetax", "theta_x(fit) - (thetax1+thetax2)/2", 100, -0.001, 0.001);
  hdiff_from_fit_thetay = new TH1D("hdiff_from_fit_thetay", "theta_y(fit) - (thetay1+thetay2)/2", 100, -0.001, 0.001);

  E1U_W1D_xy = new TH2D("E1U_W1D_xy","xy positions of E1U and W1D", 100, -0.041, 0.065, 100, -0.077, 0.079) ;
  E2U_W2D_xy = new TH2D("E2U_W2D_xy","xy positions of E2U and W2D", 100, -0.041, 0.065, 100, -0.077, 0.079) ;
  E1D_W1U_xy = new TH2D("E1D_W1U_xy","xy positions of E1D and W1U", 100, -0.041, 0.065, 100, -0.077, 0.079) ;
  E2D_W2U_xy = new TH2D("E2D_W2U_xy","xy positions of E2D and W2U", 100, -0.041, 0.065, 100, -0.077, 0.079) ;

  E1U_E1D_xy = new TH2D("E1U_E1D_xy","xy positions of E1U and E1D", 100, -0.041, 0.065, 100, -0.077, 0.079) ;
  E2U_E2D_xy = new TH2D("E2U_E2D_xy","xy positions of E2U and E2D", 100, -0.041, 0.065, 100, -0.077, 0.079) ;
  W1U_W1D_xy = new TH2D("W1U_W1D_xy","xy positions of W1U and W1D", 100, -0.041, 0.065, 100, -0.077, 0.079) ;
  W2U_W2D_xy = new TH2D("W2U_W2D_xy","xy positions of W2U and W2D", 100, -0.041, 0.065, 100, -0.077, 0.079) ;

  // Before the Fiducial cuts
  E1U_E1D = new TH2D("E1U_E1D","xy positions of E1U and E1D", 100, -0.041, 0.065, 100, -0.077, 0.079) ;
  E2U_E2D = new TH2D("E2U_E2D","xy positions of E2U and E2D", 100, -0.041, 0.065, 100, -0.077, 0.079) ;
  W1U_W1D = new TH2D("W1U_W1D","xy positions of W1U and W1D", 100, -0.041, 0.065, 100, -0.077, 0.079) ;
  W2U_W2D = new TH2D("W2U_W2D","xy positions of W2U and W2D", 100, -0.041, 0.065, 100, -0.077, 0.079) ;

  HowManyPlanes = new TH1D("HowManyPlanes","No. of planes in a RP point", 7, 0., 7.) ;

  char strs[19], strl[100];

  // All -t bins
  hUU[0] = new TH1D("hUU", "Spins Up   Up  ", 360/deg_delta, -180., 180.);
  hDD[0] = new TH1D("hDD", "Spins Down Down", 360/deg_delta, -180., 180.);
  hUD[0] = new TH1D("hUD", "Spins Up   Down", 360/deg_delta, -180., 180.);
  hDU[0] = new TH1D("hDU", "Spins Down Up", 360/deg_delta, -180., 180.);
  raw_asym[0] = new TH1D("raw_asym", "Raw Asymmetry", 180/deg_delta, -90., 90.);
  false_asym[0] = new TH1D("false_asym", "Raw Asymmetry", 180/deg_delta, -90., 90.);


  // 5-bins (for historical reason: the smallest is at the last)
  //
  //  char CRanges[intv-1][25] = {  "-t < 0.1", "0.1 #leq -t < 0.3", "0.3 #leq -t < 4.", 
  //				"4. #leq -t < 40000.", "40000. #leq -t"  } ;
  char CRanges[intv-1][25] = {  "-t < 0.2", "0.2 #leq -t < 0.35", "0.35 #leq -t < 1", 
  				"1. #leq -t", "Nothing"  } ;

  for ( Int_t i = 1 ; i<intv; i++ ) {

    sprintf(strs, "hUU%d", i);
    sprintf(strl, "Spins Up Up (%s)", CRanges[i-1]);
    hUU[i] = new TH1D(strs, strl, 360/deg_delta, -180., 180.);     

    sprintf(strs, "hDD%d", i);
    sprintf(strl, "Spins Down Down (%s)", CRanges[i-1]);
    hDD[i] = new TH1D(strs, strl, 360/deg_delta, -180., 180.);     

    sprintf(strs, "hUD%d", i);
    sprintf(strl, "Spins Up Down (%s)", CRanges[i-1]);
    hUD[i] = new TH1D(strs, strl, 360/deg_delta, -180., 180.);     

    sprintf(strs, "hDU%d", i);
    sprintf(strl, "Spins Down Up (%s)", CRanges[i-1]);
    hDU[i] = new TH1D(strs, strl, 360/deg_delta, -180., 180.);     

    sprintf(strs, "raw_asym%d", i);
    sprintf(strl, "Raw Asymmetry (%s)", CRanges[i-1]);
    raw_asym[i] = new TH1D(strs, strl, 180/deg_delta, -90., 90.); 

    sprintf(strs, "false_asym%d", i);
    sprintf(strl, "False Asymmetry (%s)", CRanges[i-1]);
    false_asym[i] = new TH1D(strs, strl, 180/deg_delta, -90., 90.); 

  }


  hAnalysisFlow = new TH1I("AnalysisFlow", "CutsFlow", kMax-1, 1, kMax);
  TString summaryLabels[] = { TString("All"), TString("ET"), TString("Si Bunch"), TString("Trig & RP"), 
			      TString("Opp Track"), TString("3#sigma(#Delta#theta)"), TString("x0/y0"), TString("Fiducial") };
  for(int tb=1; tb<kMax; ++tb)
    hAnalysisFlow->GetXaxis()->SetBinLabel(tb, summaryLabels[tb-1]);

  hTriggers = new TH1I("Triggers", "Which Triggers", t_Max-1, 1, t_Max);
  TString triggerLabels[] = { TString("CPT2"), TString("CPT2noBBCL"), TString("SDT"), TString("ET"), TString("CPTnoBBCL") };
  for(int tb=1; tb<t_Max; ++tb)
    hTriggers->GetXaxis()->SetBinLabel(tb, triggerLabels[tb-1]);

  htrig_ADC = new TH1D("trig_ADC", "Trigger PMT ADC", 100, 0., 2000.);
  htrig_TAC = new TH1D("trig_TAC", "Trigger PMT TAC", 100, 0., 2000.);

  hTOFHits_C = new TH1D("TOFHits_C", "No. of TOF Hits for CPTnoBBCL triggers", 50, 0., 50);
  hTOFHits_O = new TH1D("TOFHits_O", "No. of TOF Hits for non-CPTnoBBCL triggers", 50, 0., 50);

  hTOFMult_C = new TH1D("TOFMult_C", "No. of TOF Multiplicity for CPTnoBBCL triggers", 50, 0., 50);
  hTOFMult_O = new TH1D("TOFMult_O", "No. of TOF Multiplicity for non-CPTnoBBCL triggers", 50, 0., 50);

  hTOFHits_vs_Mult = new TH2D("TOFHits_vs_Mult","TOF Hits vs Trigger Multiplicity", 50, 0., 50., 50, 0., 50.) ;

  // Remember to put in ->Write() in finish_histogram()

  return ;

}

/*
void finish_histograms() {

  hx0->Write();
  hy0->Write();
  hz0->Write();
  hx02->Write();
  hy02->Write();
  hz02->Write();
  diff_thetax->Write() ;
  diff_thetay->Write() ;

  hthetay1->Write();
  hthetay2->Write();

  hy5->Write();
  hy7->Write();

  hdiff_y5_y7->Write();

  hm_t->Write();
  hphi->Write();
  ht_vs_phi->Write();

  return ;

}
*/

std::pair<int,int> search_matched(vector<double> xy1, vector<double> energy1, vector<Short_t> length1, 
				  vector<double> xy2, vector<double> energy2, vector<Short_t> length2, int RPID) {

  double Diff, min_Diff = ERRCODE  ;

  int  size1 = xy1.size() ;
  int  size2 = xy2.size() ;

  int  index1 = ERRCODE ;
  int  index2 = ERRCODE ;

  for ( int i = 0 ; i < size1 ; ++i ) {

    if ( ( length1[i] <= MaxClSize ) && (energy1[i] >= Emin[RPID][ length1[i]<4 ? length1[i]-1 : 3 ]) ) {

      for ( int j = 0 ; j < size2 ; ++j ) {

	if ( ( length2[j] <= MaxClSize ) && (energy2[j] >= Emin[RPID][ length2[j]<4 ? length2[j]-1 : 3 ]) ) {

	  Diff = xy1[i] - xy2[j] ;

	  if ( TMath::Abs(Diff) < TMath::Abs(min_Diff) ) {
	    min_Diff = Diff ;
	    index1 = i ;
	    index2 = j ;
	  }

	} // 	if ( length2[i] <= MaxClSize && (energy2[i] >= Emin[RPID][ length2[i]<4 ? length2[i]-1 : 3 ]) ) {

      } //       for ( int j = 0 ; j < size2 ; ++j ) {

    } //     if ( length1[i] <= MaxClSize && (energy1[i] >= Emin[RPID][ length1[i]<4 ? length1[i]-1 : 3 ]) ) {

  } //   for ( int i = 0 ; i < size1 ; ++i ) {

  //
  if ( ( size1 == 0 ) && ( size2 == 1 ) ) {
    min_Diff = 0 ;
    index1 = -1 ;
    index2 =  0 ;   
  }
  else if ( ( size1 == 1 ) && ( size2 == 0 ) ) {
    min_Diff = 0 ;
    index1 = 0 ;
    index2 = -1 ;
  }
  //

  if ( TMath::Abs(min_Diff) <= CutDiff )
    return std::make_pair( index1, index2 ) ;
  else
    return std::make_pair( ERRCODE, ERRCODE ) ;


}

void fill_phi(double phi_deg, double m_t, int BunchNumber) {

	if (  ( BlueSpin[BunchNumber] > 0 ) &&  ( YellowSpin[BunchNumber] > 0 ) ) {

	  hUU[0]->Fill( phi_deg ) ;

	  if ( m_t<div1 ) 
	    hUU[1]->Fill(phi_deg) ;
	  else if ( m_t>=div1 && m_t<div2 ) 
	    hUU[2]->Fill(phi_deg) ;
	  else if ( m_t>=div2 && m_t<div3 ) 
	    hUU[3]->Fill(phi_deg) ;
	  else if ( m_t>=div3 && m_t<div4 ) 
	    hUU[4]->Fill(phi_deg) ;
	  else if ( m_t>=div4 ) 
	    hUU[5]->Fill(phi_deg) ;

	}
	else if (  ( BlueSpin[BunchNumber] < 0 ) &&  ( YellowSpin[BunchNumber] < 0 ) ) {

	  hDD[0]->Fill( phi_deg ) ;

	  if ( m_t<div1 ) 
	    hDD[1]->Fill(phi_deg) ;
	  else if ( m_t>=div1 && m_t<div2 ) 
	    hDD[2]->Fill(phi_deg) ;
	  else if ( m_t>=div2 && m_t<div3 ) 
	    hDD[3]->Fill(phi_deg) ;
	  else if ( m_t>=div3 && m_t<div4 ) 
	    hDD[4]->Fill(phi_deg) ;
	  else if ( m_t>=div4 ) 
	    hDD[5]->Fill(phi_deg) ;

	}
	else if (  ( BlueSpin[BunchNumber] > 0 ) &&  ( YellowSpin[BunchNumber] < 0 ) ) {

	  hUD[0]->Fill( phi_deg ) ;

	  if ( m_t<div1 ) 
	    hUD[1]->Fill(phi_deg) ;
	  else if ( m_t>=div1 && m_t<div2 ) 
	    hUD[2]->Fill(phi_deg) ;
	  else if ( m_t>=div2 && m_t<div3 ) 
	    hUD[3]->Fill(phi_deg) ;
	  else if ( m_t>=div3 && m_t<div4 ) 
	    hUD[4]->Fill(phi_deg) ;
	  else if ( m_t>=div4 ) 
	    hUD[5]->Fill(phi_deg) ;

	}
	else if (  ( BlueSpin[BunchNumber] < 0 ) &&  ( YellowSpin[BunchNumber] > 0 ) ) {

	  hDU[0]->Fill( phi_deg ) ;

	  if ( m_t<div1 ) 
	    hDU[1]->Fill(phi_deg) ;
	  else if ( m_t>=div1 && m_t<div2 ) 
	    hDU[2]->Fill(phi_deg) ;
	  else if ( m_t>=div2 && m_t<div3 ) 
	    hDU[3]->Fill(phi_deg) ;
	  else if ( m_t>=div3 && m_t<div4 ) 
	    hDU[4]->Fill(phi_deg) ;
	  else if ( m_t>=div4 ) 
	    hDU[5]->Fill(phi_deg) ;

	}
	else
	  cout << "Blue spin : " << (int) BlueSpin[BunchNumber] << " and Yellow spin : " << (int) YellowSpin[BunchNumber] << endl ;

  return ;

}


void Calculate_Asymmetry(TH1D *h1, TH1D *h2, TH1D *hasymmetry) {

  Double_t rup, rdn, lup, ldn, denominator;

  for ( Int_t i = 0; i<dim_phi ; i++ ) {

    // -ve phi

    // raw asymmetry
    rup = h1->GetBinContent(i+1) ;
    lup = h1->GetBinContent(180/deg_delta-i) ;
    rdn = h2->GetBinContent(i+1) ;
    ldn = h2->GetBinContent(180/deg_delta-i) ;

    denominator = TMath::Sqrt(rup * ldn) + TMath::Sqrt(rdn * lup) ;

    if ( denominator > 0. ) {
      hasymmetry->SetBinContent(dim_phi-i,(TMath::Sqrt(rup * ldn) - TMath::Sqrt(rdn * lup)) / denominator ) ;
      hasymmetry->SetBinError(dim_phi-i, TMath::Sqrt(rup*ldn*(rdn + lup) + rdn*lup*(rup + ldn)) / denominator / denominator ) ;
    }
    else {
      hasymmetry->SetBinContent(dim_phi-i, -10.) ;
      hasymmetry->SetBinError(dim_phi-i, 100000.) ;
    }

    // +ve phi

    // raw asymmetry
    rup = h1->GetBinContent(360/deg_delta-i) ;
    lup = h1->GetBinContent(180/deg_delta+i+1) ;
    rdn = h2->GetBinContent(360/deg_delta-i) ;
    ldn = h2->GetBinContent(180/deg_delta+i+1) ;
    
    denominator = TMath::Sqrt(rup * ldn) + TMath::Sqrt(rdn * lup) ;

    if ( denominator > 0. ) {
      hasymmetry->SetBinContent(i+dim_phi+1,(TMath::Sqrt(rup * ldn) - TMath::Sqrt(rdn * lup)) / denominator ) ;
      hasymmetry->SetBinError(i+dim_phi+1, TMath::Sqrt(rup*ldn*(rdn + lup) + rdn*lup*(rup + ldn)) / denominator / denominator ) ;
    }
    else {
      hasymmetry->SetBinContent(i+dim_phi+1, -10.) ;
      hasymmetry->SetBinError(i+dim_phi+1, 100000.) ;
    }

  }

}



UInt_t countSetBits (UInt_t n) {

  UInt_t count = 0;

  while (n) {
    count += n & 1;
    n >>= 1;
  }

  return count;

}
