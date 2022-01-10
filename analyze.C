/*

Aug. 2, 2021: today, I've implemented cluster-size cut and cluster-energy (Emin) cuts.  This has made the histograms
              ClsEnergy and ClsLength not very meaningful.

	      In order to see all energy and all lengths, one has to take away all those energy/cluster-size cuts.

 */

// c++ headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <utility>

// ROOT headers
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraph.h"

// picoDst headers
#include "StUPCEvent.h"
#include "StRPEvent.h"
#include "StUPCRpsCluster.h"

#include "analyze.h"

using namespace std;

// 2017
const Int_t trigger_ADC_thres = 30 ; 
const Int_t trigger_TAC_thres = 200 ;
const Int_t trigger_TAC_thres_up = 1750 ;


const double Angle_fit_agree = 0.000125 ;




Double_t fitf(Double_t *x, Double_t *par) {
  Double_t fitval = par[0] + par[1]*(x[0]-par[2]);
  return fitval;
}

double RP_xy[64] ; 

Int_t Read_Spin(TString filename) {

  Int_t i, lines=0 ;

  //ifstream ispin(".spin", ifstream::in);
  ifstream ispin(filename, ifstream::in);
  if ( ! ispin.is_open() ) {
    cerr << "The file " << filename << " does NOT exist. " << endl ;
    return -1 ;
  }

  string linestring ;

  while ( ispin.good() && !ispin.eof() && ispin.peek() != EOF ) {
    getline(ispin, linestring,'\n');
    istringstream linestore(linestring);
    
    //    linestore >> istore >> i;
    linestore >> i; // No store
    linestore >> BlueSpin[i] >> YellowSpin[i] ;

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Akio said: "CDEV report pattern they are putting out from source. You need to take into account all spin flips from linac to booster to ags to RHIC to IP12 to IP6. 
    // I think IP12 is seeing even# of flips from source, and IP6 odd# for years." => +ve (CDEV) => -ve (IP6)
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    BlueSpin[i] = -1 * BlueSpin[i] ;
    YellowSpin[i] = -1 * YellowSpin[i] ;
    cout << i << " " << BlueSpin[i] << " " <<  YellowSpin[i] << endl ;
    lines++ ;

  } ;

  cout << lines << " lines (of spin patterns) read. " << endl ;

  ispin.close();

  //
  return 0 ;

}



void read_RP_xy(string searchpattern) {

  string line ;

  ifstream infile("pp2ppOffsets_v2017.0.2.txt") ;

  size_t findrun ;

  if ( infile.is_open() ) {
    while ( infile.good() ) {

      getline(infile, line) ;
      findrun = line.find(searchpattern) ;
      if ( findrun != string::npos ) {
	cout << " Found Bogdan's alignment ...! " << endl ;
	getline(infile,line) ;
	stringstream sline(line);
	for ( int i =0; i<16 ; i++) {
	  sline >> RP_xy[i] ;
	  cout << RP_xy[i] << " " ;
	}
	cout << endl ;

	sline.str(string());
	getline(infile,line) ;
	sline.clear();
	sline.str(line) ;

	for ( int i =16; i<32 ; i++) {
	  sline >> RP_xy[i] ;
	  cout << RP_xy[i] << " " ;
	}
	cout << endl ;

	sline.str(string());
	getline(infile,line) ;
	sline.clear();
	sline.str(line) ;

	for ( int i =32; i<48 ; i++) {
	  sline >> RP_xy[i] ;
	  cout << RP_xy[i] << " " ;
	}
	cout << endl ;

	sline.str(string());
	getline(infile,line) ;
	sline.clear();
	sline.str(line) ;

	for ( int i =48; i<64 ; i++) {
	  sline >> RP_xy[i] ;
	  cout << RP_xy[i] << " " ;
	}
	cout << endl ;

	break ;
      }
    }

  }
  else {
    cerr << "Bogdan's alignment does NOT exist !" << endl ;
    exit(EXIT_FAILURE);
  }

  infile.close() ;

  return ;

}

//_____________________________________________________________________________
void analyze(UInt_t RunNumber=18176018) {

  //open input file
  //  TFile *infile = TFile::Open("/gpfs01/star/pwg/jaroslav/star-upcDst-data/test_productions/mc/StUPC_slight14e1x1.root", "read");
  //  TFile *infile = TFile::Open("/star/data01/pwg_tasks/upc02/PartRHICf/18176017/18176017.root", "read");
  //  TFile *infile = TFile::Open("/star/data01/pwg_tasks/upc02/PartRHICf/18174055/18174055.root", "read");
  //  TFile *infile = TFile::Open("/star/data01/pwg_tasks/upc02/18174055/18174055/18174055.root","read");

  ostringstream stm ;
  //  ofstream outdata("event.list") ; // writing out events

  // Can't use string ; otherwise it'd crash
  TString filename ;

  stm << RunNumber ;

  filename = ".spin." + stm.str() ;

  if ( Read_Spin(filename) < 0 ) {
    cerr << "Read_Spin failed ?!" << endl ;
    exit(EXIT_FAILURE);
  }

  read_RP_xy(stm.str()) ;

  if ( read_fiducial() != 0) {
    cerr << "Something is wrong while reading fiducial cuts ..." << endl ;
    return ;
  }

  filename = "/star/data01/pwg_tasks/upc02/PartRHICf/" + stm.str() + "/" + stm.str() + ".root" ;

  cout << filename << endl ;

  TFile *infile = TFile::Open(filename, "read");
  //

  //get picoDst tree in file
  TTree *upcTree = dynamic_cast<TTree*>( infile->Get("mUPCTree") );

  filename = "t" + stm.str() + ".root" ;
  //open output file
  //  TFile *outfile = TFile::Open("output.root", "recreate");
  TFile *outfile = TFile::Open(filename, "recreate");

  //pT histogram
  //  TH1D *hPt = new TH1D("hPt", "hPt", 100, 0.0, 0.01);

  //connect upc event to the tree
  static StUPCEvent *upcEvt = 0x0;
  upcTree->SetBranchAddress("mUPCEvent", &upcEvt);

  static StRPEvent *rpEvt = 0x0;
  upcTree->SetBranchAddress("mRPEvent", &rpEvt);
  
  UInt_t NRP = rpEvt->mNumberOfRomanPots ;

  //  assert( NRP == 8 ) ;

  Int_t index, this_index ;
  UInt_t i, j ;
  //  Double_t x, y ;

  double mslope, thetax1, thetay1, thetax2, thetay2 ;
  double thetax, thetay, m_t, phi ;

  char hname[10], htitle[100];

  char RPNAME[8][5] = {"E1U", "E1D", "E2U", "E2D", "W1U", "W1D", "W2U", "W2D" } ;
  char RPPlane[4][2] = { "A", "B", "C", "D" } ;

  //  sprintf(txt,"Run: #%d with Silicon Clusters",RunNumber);  
  // Create the 8*2 histograms for finding 2D efficiencies 
  make_histograms() ;

  Int_t triggerID, BunchNumber ;

  UInt_t  Silicon_Bunch, NCls, NTrgCnt ;

  //ask for number of events
  Long64_t nev = upcTree->GetEntries();
  cout << "Number of events: " << nev << endl;

  //  nev= 200 ;

  UShort_t ADC[NRP][2], TAC[NRP][2] ;

  UInt_t LastPlane, LastRP ;

  pair<int,int> matched_indices ;

  int ix1, ix2, iy1, iy2 ;

  int ncount = 0 ;

  Bool_t BCounted = kFALSE, BSelected ;
  Bool_t E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, EU, ED, WU, WD, EA, EB ;
  int iboth = 0 ;

  Int_t nafter_theta3sigma = 0, ngoodfit = 0, nfiducial = 0, nV0EW = 0, nET = 0 ;

  int NCls_RP[NRP] [4];

  int planebits[NRP] ;

  //  double Diff, min_Diff  ;
  double sigma_thetax, sigma_thetay, mean_thetax, mean_thetay ;
  double sigma_x0, sigma_y0, mean_x0, mean_y0 ;

  const int NFIT = 4 ;
  double xfit[NFIT], yfit[NFIT] ;

  double a1, b1, c1, a2, b2, c2, x[NRP], y[NRP], z[NRP], x1, x2, y1, y2, z1, z2 ;
  int Good[NRP] ;

  double x0E, y0E, x0W, y0W ;


  // Reading in the mean and sigma of the angular comparison ===========================================
  filename = "." + stm.str() + ".mean_sigma" ;

  cout << "Mean and Sigma file is : " << filename << endl ;

  ifstream msfile(filename) ;

  if ( msfile.fail() ) {
    cerr << "Failure in reading : " << filename << endl ;
    exit(EXIT_FAILURE);
  }

  double temp ; // to read the uncertainties of mean during the fit --- not used so far
  msfile >> mean_thetax >> temp >> sigma_thetax >> temp ;
  msfile >> mean_thetay >> temp >> sigma_thetay >> temp ;

  /* 18176034 (2nd peak)
  mean_thetay = -5.34409E-04 ;
  sigma_thetay = 5.45644E-5 ;
  */

  cout << " mean_thetax : " << mean_thetax  << ";  sigma_thetax : " << sigma_thetax  << endl ;
  cout << " mean_thetay : " << mean_thetay  << ";  sigma_thetay : " << sigma_thetay  << endl ;

  msfile.close();
  // ==================================================================================


  // Reading in the mean and sigma of the angular comparison ===========================================
  filename = "." + stm.str() + ".0_mean_sigma" ;

  cout << "hx0/hy0 : Mean and Sigma file is : " << filename << endl ;

  msfile.open(filename) ;

  if ( msfile.fail() ) {
    cerr << "Failure in reading : " << filename << endl ;
    exit(EXIT_FAILURE);
  }

  msfile >> mean_x0 >> temp >> sigma_x0 >> temp ;
  msfile >> mean_y0 >> temp >> sigma_y0 >> temp ;

  cout << " mean_x0 : " << mean_x0  << ";  sigma_x0 : " << sigma_x0  << endl ;
  cout << " mean_y0 : " << mean_y0  << ";  sigma_y0 : " << sigma_y0  << endl ;

  msfile.close();
  // ==================================================================================

  //event loop
  for(Long64_t iev=0; iev<nev; iev++) {

    //    if ( ( iev % 10000) == 0 ) cout << "No. of events comes to : " << iev << " : =================== " << endl ;

    //get the event
    upcTree->GetEntry(iev);

    hAnalysisFlow->Fill(kAll);

    BunchNumber = upcEvt->getBunchCrossId7bit() ;
    //    cout << "Bunch number: " << upcEvt->getBunchCrossId7bit() << " ; "  << upcEvt->getBunchCrossId() << endl ;
 
    // cout << "Run Number " << upcEvt->getRunNumber() << endl ;

    //    if ( upcEvt->getRunNumber() > 18174055 ) {
    if ( RunNumber > 18174055 ) {
      triggerID = 590709 ;
    }
    else {  // for 18174055 
      triggerID = 4 ;
    }

    //   cout << "triggerID : " << triggerID << endl ;
    
    // RP_ET cut

    // RP_ET cut
    // Due to prescale, ET trigger ID may not be present even when it should
    //	      if ( ! upcEvt->isTrigger(triggerID) )
    //		continue ;
    // 3rd bit in TCU output, the top table at the last page of https://www.star.bnl.gov/public/trg/TSL/Software/TOF_MTD_PP2PP.pdf
    if ( ( upcEvt->getLastDSM0() & 4 ) != 4 )
      continue ;

    nET++ ;

    hAnalysisFlow->Fill(kET);


    vector<double> planeA[ NRP ], energyA[ NRP ] ;
    vector<double> planeB[ NRP ], energyB[ NRP ] ;
    vector<double> planeC[ NRP ], energyC[ NRP ] ;
    vector<double> planeD[ NRP ], energyD[ NRP ] ;

    vector<Short_t> lengthA[ NRP ],  lengthB[ NRP ],  lengthC[ NRP ],  lengthD[ NRP ] ;

    //    vector<double>::iterator itl, its ;
    //    vector<Short_t>::iterator sitl ;

    BCounted = kFALSE ;
    LastPlane = TMath::Abs(ERRCODE) ;
    LastRP = TMath::Abs(ERRCODE) ;

    Silicon_Bunch = rpEvt->siliconBunch() ;
    NCls = rpEvt->getNumberOfClusters() ;

    //    NCls_vs_siliconBunch->Fill( double(Silicon_Bunch), double(NCls) ) ;

    // First cut on Silicon-bunch
    if ( Silicon_Bunch > 0 && Silicon_Bunch <=8 )
      continue ;

    hAnalysisFlow->Fill(kSiBunch);

    for ( i =0; i<NRP; i++ ) {

      planebits[i] = 0 ;

      for ( j=0; j<4; j++ ) {
	NCls_RP[i][j] = 0 ;
      }

      //      assert( planeA[i].size() == 0 ) ;
      //      assert( lengthD[i].size() == 0 ) ;

      ADC[i][0] = rpEvt->adc(i, 0); 
      ADC[i][1] = rpEvt->adc(i, 1); 
      TAC[i][0] = rpEvt->tac(i, 0); 
      TAC[i][1] = rpEvt->tac(i, 1); 

      if ( ! upcEvt->isTrigger(triggerID) ) {
	htrig_ADC->Fill( ADC[i][0] ) ;
	htrig_ADC->Fill( ADC[i][1] ) ;
	htrig_TAC->Fill( TAC[i][0] ) ;
	htrig_TAC->Fill( TAC[i][1] ) ;
      }

    }

    BSelected = kFALSE ;

    /* Trigger verification:
    E1U = ( (TAC[0][0] > trigger_TAC_thres) && (TAC[0][0] < trigger_TAC_thres_up) && (ADC[0][0] > trigger_ADC_thres) ) ||
      ( (TAC[0][1] > trigger_TAC_thres) && (TAC[0][1] < trigger_TAC_thres_up) && (ADC[0][1] > trigger_ADC_thres) ) ;

    E1D = ( (TAC[1][0] > trigger_TAC_thres) && (TAC[1][0] < trigger_TAC_thres_up) && (ADC[1][0] > trigger_ADC_thres) ) ||
          ( (TAC[1][1] > trigger_TAC_thres) && (TAC[1][1] < trigger_TAC_thres_up) && (ADC[1][1] > trigger_ADC_thres) ) ;

    E2U = ( (TAC[2][0] > trigger_TAC_thres) && (TAC[2][0] < trigger_TAC_thres_up) && (ADC[2][0] > trigger_ADC_thres) ) ||
          ( (TAC[2][1] > trigger_TAC_thres) && (TAC[2][1] < trigger_TAC_thres_up) && (ADC[2][1] > trigger_ADC_thres) ) ;

    E2D = ( (TAC[3][0] > trigger_TAC_thres) && (TAC[3][0] < trigger_TAC_thres_up) && (ADC[3][0] > trigger_ADC_thres) ) ||
          ( (TAC[3][1] > trigger_TAC_thres) && (TAC[3][1] < trigger_TAC_thres_up) && (ADC[3][1] > trigger_ADC_thres) ) ;

    W1U = ( (TAC[4][0] > trigger_TAC_thres) && (TAC[4][0] < trigger_TAC_thres_up) && (ADC[4][0] > trigger_ADC_thres) ) ||
          ( (TAC[4][1] > trigger_TAC_thres) && (TAC[4][1] < trigger_TAC_thres_up) && (ADC[4][1] > trigger_ADC_thres) ) ;

    W1D = ( (TAC[5][0] > trigger_TAC_thres) && (TAC[5][0] < trigger_TAC_thres_up) && (ADC[5][0] > trigger_ADC_thres) ) ||
          ( (TAC[5][1] > trigger_TAC_thres) && (TAC[5][1] < trigger_TAC_thres_up) && (ADC[5][1] > trigger_ADC_thres) ) ;

    W2U = ( (TAC[6][0] > trigger_TAC_thres) && (TAC[6][0] < trigger_TAC_thres_up) && (ADC[6][0] > trigger_ADC_thres) ) ||
          ( (TAC[6][1] > trigger_TAC_thres) && (TAC[6][1] < trigger_TAC_thres_up) && (ADC[6][1] > trigger_ADC_thres) ) ;

    W2D = ( (TAC[7][0] > trigger_TAC_thres) && (TAC[7][0] < trigger_TAC_thres_up) && (ADC[7][0] > trigger_ADC_thres) ) ||
          ( (TAC[7][1] > trigger_TAC_thres) && (TAC[7][1] < trigger_TAC_thres_up) && (ADC[7][1] > trigger_ADC_thres) )  ;

    //  Logical deductions :
    EU = E1U || E2U ;
    ED = E1D || E2D ;
    WU = W1U || W2U ;
    WD = W1D || W2D ;

    EA = EU && WD ;
    EB = ED && WU ;

    if ( EA || EB ) {// elastic trigger

      if ( RunNumber > 18174055 ) {

	if ( upcEvt->isTrigger(590701) ) { // CPT2 
	  hTrigORNot->Fill(1) ;
	}
	else {
	  hTrigORNot->Fill(0) ;
	}

	if ( upcEvt->isTrigger(590705) ) { // CPT2noBBCL
	  hTrigORNot->Fill(3) ;
	}
	else {
	  hTrigORNot->Fill(2) ;
	}

	if ( upcEvt->isTrigger(590703) ) { // SDT
	  hTrigORNot->Fill(5) ;
	}
	else {
	  hTrigORNot->Fill(4) ;
	}

	if ( upcEvt->isTrigger(590709) ) { // ET
	  hTrigORNot->Fill(7) ;
	}
	else {
	  hTrigORNot->Fill(6) ;
	}

	if ( upcEvt->isTrigger(590708) ) { // CPTnoBBCL
	  hTrigORNot->Fill(9) ;
	}
	else {
	  hTrigORNot->Fill(8) ;
	}


      }
      else {  // for 18174055 

	if ( upcEvt->isTrigger(1) ) {	// CPT2
	  hTrigORNot->Fill(1) ;
	}
	else {
	  hTrigORNot->Fill(0) ;
	}

	if ( upcEvt->isTrigger(2) ) {	// CPT2noBBCL
	  hTrigORNot->Fill(3) ;
	}
	else {
	  hTrigORNot->Fill(2) ;
	}

	if ( upcEvt->isTrigger(3) ) {	// SDT
	  hTrigORNot->Fill(5) ;
	}
	else {
	  hTrigORNot->Fill(4) ;
	}

	if ( upcEvt->isTrigger(4) ) {	// ET
	  hTrigORNot->Fill(7) ;
	}
	else {
	  hTrigORNot->Fill(6) ;
	}

	if ( upcEvt->isTrigger(5) ) {   // CPTnoBBCL
	  hTrigORNot->Fill(9) ;
	}
	else {
	  hTrigORNot->Fill(8) ;
	}

      }

    } //     if ( EA || EB ) {// elastic trigger

    */

    //Cluster loop
    for( i=0; i<NCls; i++) {

      StUPCRpsCluster *cls = rpEvt->getCluster(i);

      /* I've checked that the romanPotId and planeId is non-decreasing
      if ( (cls->romanPotId() < LastRP) && (LastRP != TMath::Abs(ERRCODE)) ) {
	cout << "Shout : LastRP = " << LastRP << " and this RP is = " << cls->romanPotId() << endl ;
      }


      if ( cls->romanPotId() == LastRP )
	if ( cls->planeId() < LastPlane ) 
	  cout << "Shout : for RP : " << LastRP << " and LastPlane = " << LastPlane << " and this plane is = " << cls->planeId() << endl ;

      */

      LastPlane = cls->planeId() ;
      LastRP = cls->romanPotId() ;

      /*
      if ( ! ( (TAC[0] > trigger_TAC_thres) && (TAC[0] < trigger_TAC_thres_up) && (ADC[0] > trigger_ADC_thres) ) &&
	   ! ( (TAC[1] > trigger_TAC_thres) && (TAC[1] < trigger_TAC_thres_up) && (ADC[1] > trigger_ADC_thres) ) ) {
      //      if ( ! ( (TAC[0] > trigger_TAC_thres) && (TAC[0] < trigger_TAC_thres_up) && (ADC[0] > trigger_ADC_thres) ) ||
      //	   ! ( (TAC[1] > trigger_TAC_thres) && (TAC[1] < trigger_TAC_thres_up) && (ADC[1] > trigger_ADC_thres) ) )
      	continue ;
      }
      */

      //      if ( ( (TAC[0] > trigger_TAC_thres) && (TAC[0] < trigger_TAC_thres_up) && (ADC[0] > trigger_ADC_thres) ) ||
      //	   ( (TAC[1] > trigger_TAC_thres) && (TAC[1] < trigger_TAC_thres_up) && (ADC[1] > trigger_ADC_thres) ) ) {

      if ( ( (TAC[LastRP][0] > trigger_TAC_thres) && (TAC[LastRP][0] < trigger_TAC_thres_up) && (ADC[LastRP][0] > trigger_ADC_thres) ) ||
	   ( (TAC[LastRP][1] > trigger_TAC_thres) && (TAC[LastRP][1] < trigger_TAC_thres_up) && (ADC[LastRP][1] > trigger_ADC_thres) ) ) {

	(NCls_RP[LastRP][LastPlane])++ ;

	if ( LastPlane == 0 ) {

	  //	cout << LastPlane << " : " << cls->xy() -  << " " << RP_xy[8*LastRP+LastPlane*2]/1000. << " " 
	  //	     << (RP_xy[8*LastRP+LastPlane*2] + RP_xy[8*LastRP+LastPlane*2 + 1])/1000. << endl ;

	  //	planeA[LastRP].push_back( cls->xy() ) ;
	  planeA[LastRP].push_back( (RP_xy[8*LastRP+LastPlane*2] + RP_xy[8*LastRP+LastPlane*2 + 1])/1000. + orientations[4*LastRP+LastPlane]*cls->position() ) ;
	  energyA[LastRP].push_back( cls->energy() ) ;
	  lengthA[LastRP].push_back( cls->length() ) ;
	}
	else if ( LastPlane == 1 ) {
	  //	cout << LastPlane << " : " << cls->xy() - orientations[4*LastRP+LastPlane]*cls->position() << " " << RP_xy[8*LastRP+LastPlane*2]/1000. << " " 
	  //	     << (RP_xy[8*LastRP+LastPlane*2] + RP_xy[8*LastRP+LastPlane*2 + 1])/1000. << endl ;

	  //	planeB[LastRP].push_back( cls->xy() ) ;
	  planeB[LastRP].push_back( (RP_xy[8*LastRP+LastPlane*2] + RP_xy[8*LastRP+LastPlane*2 + 1])/1000. + orientations[4*LastRP+LastPlane]*cls->position() ) ;
	  energyB[LastRP].push_back( cls->energy() ) ;
	  lengthB[LastRP].push_back( cls->length() ) ;
	}
	else if ( LastPlane == 2 ) {
	  //	cout << LastPlane << " : " << cls->xy() - orientations[4*LastRP+LastPlane]*cls->position() << " " << RP_xy[8*LastRP+LastPlane*2]/1000. << " " 
	  //	     << (RP_xy[8*LastRP+LastPlane*2] + RP_xy[8*LastRP+LastPlane*2 + 1])/1000. << endl ;

	  //	planeC[LastRP].push_back( cls->xy() ) ;
	  planeC[LastRP].push_back( (RP_xy[8*LastRP+LastPlane*2] + RP_xy[8*LastRP+LastPlane*2 + 1])/1000. + orientations[4*LastRP+LastPlane]*cls->position() ) ;
	  energyC[LastRP].push_back( cls->energy() ) ;
	  lengthC[LastRP].push_back( cls->length() ) ;
	}
	else if ( LastPlane == 3 ) {
	  //	cout << LastPlane << " : " << cls->xy() - orientations[4*LastRP+LastPlane]*cls->position() << " " << RP_xy[8*LastRP+LastPlane*2]/1000. << " " 
	  //	     << (RP_xy[8*LastRP+LastPlane*2] + RP_xy[8*LastRP+LastPlane*2 + 1])/1000. << endl ;

	  //	planeD[LastRP].push_back( cls->xy() ) ;
	  planeD[LastRP].push_back( (RP_xy[8*LastRP+LastPlane*2] + RP_xy[8*LastRP+LastPlane*2 + 1])/1000. + orientations[4*LastRP+LastPlane]*cls->position() ) ;
	  energyD[LastRP].push_back( cls->energy() ) ;
	  lengthD[LastRP].push_back( cls->length() ) ;
	}
	else
	  cerr << "planeId = " << LastPlane << " for romanPotId = " << LastRP << endl ;

	//      cout << cls->romanPotId() << " , " << cls->planeId() << " : " << cls->position() << " ; " << cls->xy() << endl;
	//      cout << trk->getPt() << endl;
	//      hPt->Fill(trk->getPt());

      }

    } //     for( i=0; i<NCls; i++) {

    NTrgCnt = 0 ;
    for ( i =0; i<NRP; i++ ) {

      if ( ( NCls_RP[i][0] + NCls_RP[i][1] + NCls_RP[i][2] + NCls_RP[i][3] ) > 0 )
	NTrgCnt++ ;
 
    }
    
    if ( NTrgCnt >= 2 ) hAnalysisFlow->Fill(kTrigger);


    /*
    mean_thetax = 6.56538E-5 ;
    sigma_thetax = 6.86404E-5 ;

    mean_thetay = -5.15183E-5 ;
    sigma_thetay = 7.23879E-5 ;

    if (  ( (NCls_RP[7][0] >1) || (NCls_RP[7][1] >1) || (NCls_RP[7][2] >1) || (NCls_RP[7][3] >1) )  ||  // avoiding too many clusters
	  ( (NCls_RP[5][0] >1) || (NCls_RP[5][1] >1) || (NCls_RP[5][2] >1) || (NCls_RP[5][3] >1) )  ||
	  ( (NCls_RP[0][0] >1) || (NCls_RP[0][1] >1) || (NCls_RP[0][2] >1) || (NCls_RP[0][3] >1) )  ||
	  ( (NCls_RP[2][0] >1) || (NCls_RP[2][1] >1) || (NCls_RP[2][2] >1) || (NCls_RP[2][3] >1) )  )
      continue ;
    */

    for ( i = 0 ; i<NRP ; i++ ) {

      Good[i] = -1 ;

      matched_indices = search_matched(planeA[i], energyA[i], lengthA[i],
				       planeC[i], energyC[i], lengthC[i], i);
      iy1 = matched_indices.first ;
      iy2 = matched_indices.second ;

      matched_indices = search_matched(planeB[i], energyB[i], lengthB[i],
				       planeD[i], energyD[i], lengthD[i], i);
      ix1 = matched_indices.first ;
      ix2 = matched_indices.second ;

      //      cout << "indices : " << ix1 << " " << ix2 << " ; " << iy1 << " " << iy2 << " " << endl ;
      //      if (  (ix1 >= 0) && (ix2 >= 0) && (iy1 >= 0) && (iy2 >= 0) ) {
      if (  ( (ix1 >= 0) || (ix2 >= 0) ) &&  ( (iy1 >= 0) || (iy2 >= 0) ) ) {

	//	cout << i << " : " << rpEvt->z(i, 0) << " " << rpEvt->z(i, 1) << " " << rpEvt->z(i, 2) << " " << rpEvt->z(i, 3) << endl ;
	//	thetay = ( planeC[i][iy2] - planeA[i][iy1] ) / ( rpEvt->z(i,2) - rpEvt->z(i,0) ) ;

	Good[i] = 1 ;

	if ( (ix1 >= 0) && (ix2 >= 0) ) {
	  x[i] = (planeB[i][ix1] + planeD[i][ix2])/2. ;
	  planebits[i] = 3 ;
	}
	else if ( ix1 >= 0 ) {
	  x[i] = planeB[i][ix1] ;
	  planebits[i] = 1 ;
	}
	else if ( ix2 >= 0 ) {
	  x[i] = planeD[i][ix2] ;
	  planebits[i] = 2 ;
	}

	if ( (iy1 >= 0) && (iy2 >= 0) ) {
	  y[i] = (planeA[i][iy1] + planeC[i][iy2])/2. ;
	  planebits[i] += 12 ;
	}
	else if ( iy1 >= 0 ) {
	  y[i] = planeA[i][iy1] ;
	  planebits[i] += 4 ;
	}
	else if ( iy2 >= 0 ) {
	  y[i] = planeC[i][iy2] ;
	  planebits[i] += 8 ;
	}

	// 2021-8-27: Somehow, the z positions of 2017 in DST (from MuDst) were using the 2015 values (a bit wrong !).
	//	z[i] = (rpEvt->z(i,1) + rpEvt->z(i,2))/2. ;
	z[i] = (RP_POS_Z[i*4 + 1] + RP_POS_Z[i*4 + 2])/2. ;

	//	if ( i==5  || i==7 ) {
	  //	  cout << i << " : " << x[i] << " " << y[i] << " " << z[i] << " --- " << rpEvt->z(i,1) << " " << rpEvt->z(i,2) << endl ;
	//	  cout << i << " : " << x[i] << " " << y[i] << " " << z[i] << " --- " << RP_POS_Z[i*4 + 1] << " " << RP_POS_Z[i*4 + 2] << endl ;
	//	}

      }

    }

    if ( (Good[0] > 0)  && (Good[2] > 0) && (Good[5] > 0) && (Good[7] > 0) ) {       

      ncount++ ;
      BCounted = kTRUE ;

      hAnalysisFlow->Fill(kOppTrack);

      /*  This is for the intersection of East-Track and West-Track

      // (y - y1)/(x-x1) = (y2-y1)/(x2-x1) = m => mx - y + (y1- mx1) = 0
      // z is "y" in the above equation, x/y are "x" in the above equation.


      //      For x-z ( z is "y" here )
      // E1U-E2U (x):
      mslope = (z[2] - z[0])/(x[2] - x[0]) ;
      a1 = mslope ;
      b1 = -1. ;
      c1 = z[0] - mslope*x[0] ;

      // W1D-W2D (x):
      mslope = (z[7] - z[5])/(x[7] - x[5]) ;
      a2 = mslope ;
      b2 = -1. ;
      c2 = z[5] - mslope*x[5] ;

      if ( (a1*b2-a2*b1) != 0 ) {
	hx0->Fill( (b1*c2-b2*c1)/(a1*b2-a2*b1) ) ;
	hz0->Fill( (c1*a2-c2*a1)/(a1*b2-a2*b1) ) ;
      }


      //      For y-z ( z is "y" here )
      // E1U-E2U (y):
      mslope = (z[2] - z[0])/(y[2] - y[0]) ;
      a1 = mslope ;
      b1 = -1. ;
      c1 = z[0] - mslope*y[0] ;

      // W1D-W2D (y):
      mslope = (z[7] - z[5])/(y[7] - y[5]) ;
      a2 = mslope ;
      b2 = -1. ;
      c2 = z[5] - mslope*y[5] ;

      if ( (a1*b2-a2*b1) != 0 ) {
	hy0->Fill( (b1*c2-b2*c1)/(a1*b2-a2*b1) ) ;
	hz0->Fill( (c1*a2-c2*a1)/(a1*b2-a2*b1) ) ;
      }

      */

      // (y-y1)/(z-z1) = (y2-y1)/(z2-z1) => when z=0, y = y1 - z1*(y2-y1)/(z2-z1) 

      x0E = x[0] - z[0]*(x[2]-x[0])/(z[2]-z[0]) ;
      y0E = y[0] - z[0]*(y[2]-y[0])/(z[2]-z[0]) ;

      x0W = x[5] - z[5]*(x[7]-x[5])/(z[7]-z[5]) ;
      y0W = y[5] - z[5]*(y[7]-y[5])/(z[7]-z[5]) ;

      hx0->Fill( x0W - x0E );
      hy0->Fill( y0W - y0E );

      // Anoter method ===> turns out to be the same as above --- intersection of two lines EU - WD 
      // XRP1 = XIP + theta1*(ZRP1 - ZIP)
      // XRP2 = XIP + theta2*(ZRP2 - ZIP)

      // => ZIP = [ XRP1 - XRP2  - theta1*ZRP1 + theta2*ZRP2 ] / (theta2 - theta1)
      // => XIP = [ theta2*XRP1 - theta1*XRP2 - theta1*theta2*(ZRP1-ZRP2) ] / (theta2 - theta1)

      // E1U-E2U : W1D-W2D
      thetax1 = ( x[2] - x[0] ) / ( z[2] - z[0] ) ;
      thetax2 = ( x[7] - x[5] ) / ( z[7] - z[5] ) ;

      thetay1 = ( y[2] - y[0] ) / ( z[2] - z[0] ) ;
      thetay2 = ( y[7] - y[5] ) / ( z[7] - z[5] ) ;

      hy5->Fill( y[5] ) ;
      hy7->Fill( y[7] ) ;

      hdiff_y5_y7->Fill( y[7] - y[5] ) ;

      /*

      hx02->Fill( ( thetax2*x[2] - thetax1*x[7] - thetax1*thetax2*(z[2]-z[7]) ) / (thetax2- thetax1) ) ;
      hz02->Fill( ( x[2] - x[7] - thetax1*z[2] + thetax2*z[7] ) / (thetax2- thetax1) ) ;

      hy02->Fill( ( thetay2*y[2] - thetay1*y[7] - thetay1*thetay2*(z[2]-z[7]) ) / (thetay2- thetay1) ) ;
      hz02->Fill( ( y[2] - y[7] - thetay1*z[2] + thetay2*z[7] ) / (thetay2- thetay1) ) ;

      */

      hthetay1->Fill(thetay1);
      hthetay2->Fill(thetay2);
      diff_thetax->Fill(thetax1 - thetax2);
      diff_thetay->Fill(thetay1 - thetay2);


      thetax = ERRCODE ;
      if ( TMath::Abs( thetax1-thetax2 - mean_thetax ) <= nsigma*sigma_thetax ) {

	xfit[0] = z[0] ;
	yfit[0] = x[0] ;
	xfit[1] = z[2] ;
	yfit[1] = x[2] ;
	xfit[2] = z[5] ;
	yfit[2] = x[5] ;
	xfit[3] = z[7] ;
	yfit[3] = x[7] ;

	TF1 *func = new TF1("fitf",fitf,-17.6,17.6,3);
	TGraph* gr = new TGraph(NFIT,xfit,yfit) ;
	func->SetParameter(1, (thetax1+thetax2)/2.) ;
	gr->Fit("fitf","QRM");
	//	cout << func->GetParameter(0) << " " << func->GetParameter(1) << " " << func->GetParameter(2) << endl ;

	hx02->Fill( func->GetParameter(0) ) ;
	thetax = func->GetParameter(1) ;
	hz02->Fill( func->GetParameter(2) ) ;

	/*
	if ( TMath::Abs(thetax - (thetax1+thetax2)/2.) > 0.001 ) {
	  cout << "thetax - (thetax1+thetax2)/2.  too big => thetax = " << thetax << " and " << " (thetax1+thetax2)/2. = "
	       << (thetax1+thetax2)/2. << " and thetax1 = " << thetax1 << ", thetax2 = " << thetax2 
	       << " ChiSquare/NDF = " << func->GetChisquare()/func->GetNDF() << endl ;
	  for ( Int_t itest = 0 ; itest<4 ; itest++) {
	    cout << "xfit[" << itest << "] = " << xfit[itest] << " ;" << endl ;
	    cout << "yfit[" << itest << "] = " << yfit[itest] << " ;" << endl ;
	  }
	}
	*/
	hdiff_from_fit_thetax->Fill(thetax - (thetax1+thetax2)/2.) ;

	// The fit sometimes give very small theta
	//	thetax = (thetax1+thetax2)/2. ;

	delete gr ;
	delete func ;

      }

      thetay = ERRCODE ;
      if ( ( TMath::Abs( thetay1-thetay2 - mean_thetay ) <= nsigma*sigma_thetay ) ||
	   ( ( RunNumber==18176034 ) && ( TMath::Abs( thetay1-thetay2 + 5.34409e-04 ) <= nsigma*5.45644e-05 ) ) ) {

	xfit[0] = z[0] ;
	yfit[0] = y[0] ;
	xfit[1] = z[2] ;
	yfit[1] = y[2] ;
	xfit[2] = z[5] ;
	yfit[2] = y[5] ;
	xfit[3] = z[7] ;
	yfit[3] = y[7] ;

	/* hack
	cout << "EUWD: " << endl ;
	for ( Int_t itest = 0 ; itest<4 ; itest++) {
	  cout << "xfit[" << itest << "] = " << xfit[itest] << " ;" << endl ;
	  cout << "yfit[" << itest << "] = " << yfit[itest] << " ;" << endl ;
	}
	cout << " ================== " << endl ;
	*/

	TF1 *func = new TF1("fitf",fitf,-17.6,17.6,3);
	TGraph* gr = new TGraph(NFIT,xfit,yfit) ;
	func->SetParameter(1, (thetay1+thetay2)/2.) ;
	gr->Fit("fitf","QRM");
	//	cout << "theta input : " << (thetay1+thetay2)/2. << endl ;

	//	cout << func->GetParameter(0) << " " << func->GetParameter(1) << " " << func->GetParameter(2) << endl ;

	hy02->Fill( func->GetParameter(0) ) ;
	thetay = func->GetParameter(1) ;
	hz02->Fill( func->GetParameter(2) ) ;

	/*
	if ( TMath::Abs(thetay - (thetay1+thetay2)/2.) > 0.001 )
	  cout << "thetay - (thetay1+thetay2)/2.  too big => thetay = " << thetay << " and " << " (thetay1+thetay2)/2. = "
	       << (thetay1+thetay2)/2. << " and thetay1 = " << thetay1 << ", thetay2 = " << thetay2 
	       << " ChiSquare/NDF = " << func->GetChisquare()/func->GetNDF() << endl ;
	*/

	hdiff_from_fit_thetay->Fill(thetay - (thetay1+thetay2)/2.) ;

	// The fit sometimes give very small theta
	//	thetay = (thetay1+thetay2)/2. ;

	delete gr ;
	delete func ;

      }

      if ( ( thetax > ERRCODE ) && (thetay > ERRCODE ) ) {

	nafter_theta3sigma++ ;
	hAnalysisFlow->Fill(kTheta);

	if ( ( TMath::Abs(x0W - x0E - mean_x0) <  nsigma*sigma_x0 ) && 
	     ( ( TMath::Abs(y0W - y0E - mean_y0) < nsigma*sigma_y0 ) || 
	       ( ( RunNumber==18176034 ) && ( TMath::Abs(y0W - y0E + 8.47275E-03 ) < nsigma*4.88228E-4 ) ) ) ) {

	  nV0EW++ ;
	  hAnalysisFlow->Fill(kx0y0);

	  //	  if ( ( TMath::Abs(thetax - (thetax1+thetax2)/2.) < Angle_fit_agree ) && 
	  //	       ( ( TMath::Abs(thetay - (thetay1+thetay2)/2.) < Angle_fit_agree ) || 
	  //		 ( ( RunNumber==18176034 ) && ( TMath::Abs(thetay - (thetay1+thetay2)/2. + 2.51424E-4 ) < Angle_fit_agree ) ) ) ) {

	  //	    ngoodfit++ ;

	    thetax = (thetax1+thetax2)/2. ;
	    thetay = (thetay1+thetay2)/2. ;

	    E1U_E1D->Fill(x[0], y[0]) ;
	    W1U_W1D->Fill(x[5], y[5]) ;

	    E2U_E2D->Fill(x[2], y[2]) ;
	    W2U_W2D->Fill(x[7], y[7]) ;

	    if ( FiducialPos(0, x[0], y[0]) &&  FiducialNeg(5, x[5], y[5]) && 
		 FiducialPos(2, x[2], y[2]) &&  FiducialNeg(7, x[7], y[7]) ) {

	      nfiducial++ ;
	      hAnalysisFlow->Fill(kFiducial);

	      BSelected = kTRUE ;

	      m_t = beam_momentum*beam_momentum*(thetax*thetax + thetay*thetay) ;
	      if ( m_t > 0.64 ) cout << "Big: -t = " << m_t << endl ;
	      phi = TMath::ATan2(TMath::Tan(thetay), TMath::Tan(thetax))*180./TMath::Pi();

	      hphi->Fill(phi) ;
	      ht_vs_phi->Fill(phi, m_t) ;

	      // For radial spin
	      //	      phi = phi - 90. ; 
	      //	      if ( phi < -180. ) phi = phi + 360. ;
	      phi = phi + 90. ; 
	      if ( phi > 180. ) phi = phi - 360. ;

	      //	      assert( (phi <= 180.) && (phi>=-180.) ) ;


	      hm_t_EUWD->Fill(m_t) ;
	      hm_t->Fill(m_t) ;

	      hm_tv_EUWD->Fill(m_t) ;
	      hm_tv->Fill(m_t) ;

	      fill_phi(phi, m_t, BunchNumber) ; 

	      E1U_W1D_xy->Fill(x[0], y[0]) ;
	      E1U_W1D_xy->Fill(x[5], y[5]) ;

	      E2U_W2D_xy->Fill(x[2], y[2]) ;
	      E2U_W2D_xy->Fill(x[7], y[7]) ;

	      E1U_E1D_xy->Fill(x[0], y[0]) ;
	      W1U_W1D_xy->Fill(x[5], y[5]) ;

	      E2U_E2D_xy->Fill(x[2], y[2]) ;
	      W2U_W2D_xy->Fill(x[7], y[7]) ;

	      // Count how many planes
	      HowManyPlanes->Fill( countSetBits(planebits[0]) + 0.1 ) ;
	      HowManyPlanes->Fill( countSetBits(planebits[2]) + 0.1 ) ;
	      HowManyPlanes->Fill( countSetBits(planebits[5]) + 0.1 ) ;
	      HowManyPlanes->Fill( countSetBits(planebits[7]) + 0.1 ) ;

	      //	      outdata << RunNumber << " " << upcEvt->getEventNumber() << endl ;

	    } //	  if ( FiducialPos(0, x[0], y[0]) &&  FiducialNeg(5, x[5], y[5]) && FiducialPos(2, x[2], y[2]) &&  FiducialNeg(7, x[7], y[7]) ) {

	    //	  } // 	  if ( ( TMath::Abs(thetax - (thetax1+thetax2)/2.) < Angle_fit_agree ) && 

	}

      }

    }



    if ( (Good[1] > 0)  && (Good[3] > 0) && (Good[4] > 0) && (Good[6] > 0) ) {

      ncount++ ;
      if ( BCounted ) {
	iboth++ ;
	//      	cout << "Both EU-WD  and ED-WU in the same event ! "<< iboth << endl ;
      }

      hAnalysisFlow->Fill(kOppTrack);

      x0E = x[1] - z[1]*(x[3]-x[1])/(z[3]-z[1]) ;
      y0E = y[1] - z[1]*(y[3]-y[1])/(z[3]-z[1]) ;

      x0W = x[4] - z[4]*(x[6]-x[4])/(z[6]-z[4]) ;
      y0W = y[4] - z[4]*(y[6]-y[4])/(z[6]-z[4]) ;

      hx0->Fill( x0W - x0E );
      hy0->Fill( y0W - y0E );

      //      cout << ncount << " : E1D E2D -- W1U W2U ! " << endl ;

      // E1D-E2D : W1U-W2U
      thetax1 = ( x[3] - x[1] ) / ( z[3] - z[1] ) ;
      thetax2 = ( x[6] - x[4] ) / ( z[6] - z[4] ) ;

      thetay1 = ( y[3] - y[1] ) / ( z[3] - z[1] ) ;
      thetay2 = ( y[6] - y[4] ) / ( z[6] - z[4] ) ;


      diff_thetax->Fill(thetax1 - thetax2);
      diff_thetay->Fill(thetay1 - thetay2);

      thetax = ERRCODE ;
      if ( TMath::Abs( thetax1-thetax2 - mean_thetax ) <= nsigma*sigma_thetax ) {

	xfit[0] = z[1] ;
	yfit[0] = x[1] ;
	xfit[1] = z[3] ;
	yfit[1] = x[3] ;
	xfit[2] = z[4] ;
	yfit[2] = x[4] ;
	xfit[3] = z[6] ;
	yfit[3] = x[6] ;

	TF1 *func = new TF1("fitf",fitf,-17.6,17.6,3);
	TGraph* gr = new TGraph(NFIT,xfit,yfit) ;
	func->SetParameter(1, (thetax1+thetax2)/2.) ;
	gr->Fit("fitf","QRM");
	//	cout << func->GetParameter(0) << " " << func->GetParameter(1) << " " << func->GetParameter(2) << endl ;

	hx02->Fill( func->GetParameter(0) ) ;
	thetax = func->GetParameter(1) ;
	hz02->Fill( func->GetParameter(2) ) ;

	/*
	if ( TMath::Abs(thetax - (thetax1+thetax2)/2.) > 0.001 )
	  cout << "thetax - (thetax1+thetax2)/2.  too big => thetax = " << thetax << " and " << " (thetax1+thetax2)/2. = "
	       << (thetax1+thetax2)/2. << " and thetax1 = " << thetax1 << ", thetax2 = " << thetax2 
	       << " ChiSquare/NDF = " << func->GetChisquare()/func->GetNDF() << endl ;
	*/

	hdiff_from_fit_thetax->Fill(thetax - (thetax1+thetax2)/2.) ;

	// The fit sometimes give very small theta
	//	thetax = (thetax1+thetax2)/2. ;

	delete gr ;
	delete func ;

      }

      thetay = ERRCODE ;
      if ( TMath::Abs( thetay1-thetay2 - mean_thetay ) <= nsigma*sigma_thetay ) {

	xfit[0] = z[1] ;
	yfit[0] = y[1] ;
	xfit[1] = z[3] ;
	yfit[1] = y[3] ;
	xfit[2] = z[4] ;
	yfit[2] = y[4] ;
	xfit[3] = z[6] ;
	yfit[3] = y[6] ;

	/* hack
	cout << "EDWU: " << endl ;
	for ( Int_t itest = 0 ; itest<4 ; itest++) {
	  cout << "xfit[" << itest << "] = " << xfit[itest] << " ;" << endl ;
	  cout << "yfit[" << itest << "] = " << yfit[itest] << " ;" << endl ;
	}
	cout << " ================== " << endl ;
	*/

	TF1 *func = new TF1("fitf",fitf,-17.6,17.6,3);
	TGraph* gr = new TGraph(NFIT,xfit,yfit) ;
	func->SetParameter(1, (thetay1+thetay2)/2.) ;
	gr->Fit("fitf","QRM");
	//	cout << "theta input : " << (thetay1+thetay2)/2. << endl ;

	//	cout << func->GetParameter(0) << " " << func->GetParameter(1) << " " << func->GetParameter(2) << endl ;

	hy02->Fill( func->GetParameter(0) ) ;
	thetay = func->GetParameter(1) ;
	hz02->Fill( func->GetParameter(2) ) ;

	/*
	if ( TMath::Abs(thetay - (thetay1+thetay2)/2.) > 0.001 )
	  cout << "thetay - (thetay1+thetay2)/2.  too big => thetay = " << thetay << " and " << " (thetay1+thetay2)/2. = "
	       << (thetay1+thetay2)/2. << " and thetay1 = " << thetay1 << ", thetay2 = " << thetay2 
	       << " ChiSquare/NDF = " << func->GetChisquare()/func->GetNDF() << endl ;
	*/

	hdiff_from_fit_thetay->Fill(thetay - (thetay1+thetay2)/2.) ;

	// The fit sometimes give very small theta
	//	thetay = (thetay1+thetay2)/2. ;

	delete gr ;
	delete func ;

      }

      if ( ( thetax > ERRCODE ) && (thetay > ERRCODE ) ) {

	nafter_theta3sigma++ ;
	hAnalysisFlow->Fill(kTheta);

	if ( ( TMath::Abs(x0W - x0E - mean_x0) < nsigma*sigma_x0 ) && 
	     ( TMath::Abs(y0W - y0E - mean_y0) < nsigma*sigma_y0 ) ) {

	  nV0EW++ ;
	  hAnalysisFlow->Fill(kx0y0);


	  //	  if ( ( TMath::Abs(thetax - (thetax1+thetax2)/2.) < Angle_fit_agree ) && 
	  //	       ( TMath::Abs(thetay - (thetay1+thetay2)/2.) < Angle_fit_agree )  ) {

	  //	    ngoodfit++ ;

	    thetax = (thetax1+thetax2)/2. ;
	    thetay = (thetay1+thetay2)/2. ;

	    E1U_E1D->Fill(x[1], y[1]) ;
	    W1U_W1D->Fill(x[4], y[4]) ;

	    E2U_E2D->Fill(x[3], y[3]) ;
	    W2U_W2D->Fill(x[6], y[6]) ;

	    if ( FiducialPos(4, x[4], y[4]) &&  FiducialNeg(1, x[1], y[1]) && 
		 FiducialPos(6, x[6], y[6]) &&  FiducialNeg(3, x[3], y[3]) ) {

	      nfiducial++ ;
	      hAnalysisFlow->Fill(kFiducial);

	      BSelected = kTRUE ;

	      m_t = beam_momentum*beam_momentum*(thetax*thetax + thetay*thetay) ;
	      if ( m_t > 0.64 ) cout << "Big: -t = " << m_t << endl ;
	      phi = TMath::ATan2(TMath::Tan(thetay), TMath::Tan(thetax))*180./TMath::Pi() ;

	      hphi->Fill(phi) ;
	      ht_vs_phi->Fill(phi, m_t) ;

	      // For radial spin
	      //	      phi = phi - 90. ; 
	      //	      if ( phi < -180. ) phi = phi + 360. ;

	      phi = phi + 90. ; 
	      if ( phi > 180. ) phi = phi - 360. ;

	      //	      assert( (phi <= 180.) && (phi>=-180.) ) ;


	      hm_t_EDWU->Fill(m_t) ;
	      hm_t->Fill(m_t) ;

	      hm_tv_EDWU->Fill(m_t) ;
	      hm_tv->Fill(m_t) ;

	      fill_phi(phi, m_t, BunchNumber) ;

	      E1D_W1U_xy->Fill(x[1], y[1]) ;
	      E1D_W1U_xy->Fill(x[4], y[4]) ;

	      E2D_W2U_xy->Fill(x[3], y[3]) ;
	      E2D_W2U_xy->Fill(x[6], y[6]) ;

	      E1U_E1D_xy->Fill(x[1], y[1]) ;
	      W1U_W1D_xy->Fill(x[4], y[4]) ;

	      E2U_E2D_xy->Fill(x[3], y[3]) ;
	      W2U_W2D_xy->Fill(x[6], y[6]) ;

	      // Count how many planes
	      HowManyPlanes->Fill( countSetBits(planebits[1]) + 0.1 ) ;
	      HowManyPlanes->Fill( countSetBits(planebits[3]) + 0.1 ) ;
	      HowManyPlanes->Fill( countSetBits(planebits[4]) + 0.1 ) ;
	      HowManyPlanes->Fill( countSetBits(planebits[6]) + 0.1 ) ;

	      //	      outdata << RunNumber << " " << upcEvt->getEventNumber() << endl ;

	    } // 	  if ( FiducialPos(4, x[4], y[4]) &&  FiducialNeg(1, x[1], y[1]) && FiducialPos(6, x[6], y[6]) &&  FiducialNeg(3, x[3], y[3]) ) {

	    //	  } // 	  if ( ( TMath::Abs(thetax - (thetax1+thetax2)/2.) < Angle_fit_agree ) && 

	}

      }

    }

    if ( BSelected ) {

      // TOF Hits vs TOF Trigger Multiplicity
      hTOFHits_vs_Mult->Fill( upcEvt->getTOFMultiplicity(), upcEvt->getNumberOfHits() ) ; 

      //      if ( upcEvt->getRunNumber() > 18174055 ) {
      if ( RunNumber > 18174055 ) {

	if ( upcEvt->isTrigger(590701) )	hTriggers->Fill(t_CPT2);
	if ( upcEvt->isTrigger(590703) )	hTriggers->Fill(t_SDT);
	if ( upcEvt->isTrigger(590705) )	hTriggers->Fill(t_CPT2noBBCL);
	if ( upcEvt->isTrigger(590708) )	{
	  hTriggers->Fill(t_CPTnoBBCL);
	  hTOFHits_C ->Fill( upcEvt->getNumberOfHits() )  ;
	  hTOFMult_C ->Fill( upcEvt->getTOFMultiplicity() )  ;
	}
	else {
	  hTOFHits_O ->Fill( upcEvt->getNumberOfHits() )  ;
	  hTOFMult_O ->Fill( upcEvt->getTOFMultiplicity() )  ;
	}	  

	if ( upcEvt->isTrigger(590709) )	hTriggers->Fill(t_ET);      

      }
      else {  // for 18174055 

	if ( upcEvt->isTrigger(1) )	hTriggers->Fill(t_CPT2);
	if ( upcEvt->isTrigger(2) )	hTriggers->Fill(t_CPT2noBBCL);
	if ( upcEvt->isTrigger(3) )	hTriggers->Fill(t_SDT);
	if ( upcEvt->isTrigger(4) )	hTriggers->Fill(t_ET);      
	if ( upcEvt->isTrigger(5) ) {
	  hTriggers->Fill(t_CPTnoBBCL);
	  hTOFHits_C ->Fill( upcEvt->getNumberOfHits() )  ;
	  hTOFMult_C ->Fill( upcEvt->getTOFMultiplicity() )  ;
	}
	else {
	  hTOFHits_O ->Fill( upcEvt->getNumberOfHits() )  ;
	  hTOFMult_O ->Fill( upcEvt->getTOFMultiplicity() )  ;
	}	  

      }

    } // if ( BSelected ) {


  } // event loop

  infile->Close();
  //  outdata.close();

  cout << "No. of tracks after ET Cuts : " << nET << endl ;
  cout << "No. of opposite tracks in the East and West : " << ncount << endl ;
  cout << "No. of events with both EU-WD and ED-WU track pairs : "<< iboth << endl ;
  cout << "No. of opposite tracks passing correlation " << nsigma << "sigma cut : " << nafter_theta3sigma << endl ;

  cout << "No. of opposite tracks with very close x0/y0's extrapolated from E & W at z=0 : " <<  nV0EW << endl ;

  // Commented out as this cut doesn't cut any more events after the above cut
  //  cout << "No. of opposite tracks with a good fit from East to West : " <<  ngoodfit << endl ;

  cout << "No. of tracks after Fiducial Cuts : " << nfiducial << endl ;


  for ( j = 0; j<intv ; j++ ) {
    Calculate_Asymmetry(hUU[j], hDD[j], raw_asym[j]) ;
    Calculate_Asymmetry(hUD[j], hDU[j], false_asym[j]) ;
  }


  outfile->Write(); // this will write all the histograms instead of writing each histogram one by one
  outfile->Close();

  //  cout <<

  return ;
}


