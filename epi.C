
// Daniel Pitzl, DESY, Aug 2019

// root -l epi.C

{
  // set styles:

  gStyle->SetTextFont(62);//62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 2.0, "y" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.99 ); // global title

  gStyle->SetLineWidth(1);// frames
  gStyle->SetHistLineColor(4); // 4=blau
  gStyle->SetHistLineWidth(3);

  gStyle->SetHistFillStyle(0); // no fill

  gStyle->SetFrameLineWidth(2);

  gStyle->SetHistMinimumZero(); // no zero suppression

  //gStyle->SetOptDate();
 
  gROOT->ForceStyle();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cout << "try to open eesig-e5000.dat";

  ifstream sFile( "eesig-e5000.dat" );

  if( !sFile ) {
    cout << " : failed " << endl;
    return;
  }

  cout << " : success " << endl;

  string START {" START"};
  string hd;

  while( hd != START ) {
    getline( sFile, hd ); // read one line into string
    cout << "  " << hd << endl;
  }

  vector <double> vde;
  vde.reserve(1250);

  vector <double> veesig;
  veesig.reserve(1250);

  while( ! sFile.eof() ) {

    int j;
    double de;
    double eesig;
    sFile >> j;
    sFile >> de;
    sFile >> eesig; // this is dE^2 dsig/ddE
    vde.push_back(de);
    veesig.push_back( eesig );

  } // while

  int n = vde.size();
  cout << "  read " << n << " values" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cout << "try to open eesig-pi45000.dat";

  ifstream pFile( "eesig-pi45000.dat" );

  if( !pFile ) {
    cout << " : failed " << endl;
    return;
  }

  cout << " : success " << endl;

  hd.clear();

  while( hd != START ) {
    getline( pFile, hd ); // read one line into string
    cout << "  " << hd << endl;
  }

  vector <double> wde;
  wde.reserve(1250);

  vector <double> weesig;
  weesig.reserve(1250);

  while( ! pFile.eof() ) {

    int j;
    double de;
    double eesig;
    pFile >> j;
    pFile >> de;
    pFile >> eesig; // this is dE^2 dsig/ddE
    wde.push_back(de);
    weesig.push_back( eesig );

  } // while

  int m = wde.size();
  cout << "  read " << m << " values" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // square canvas:
  //                              topleft x, y, width x, y
  TCanvas * c1 = new TCanvas( "c1", "c1", 0, 0, 813, 837 );
  //                                 to get fCw 800 fCh 800

  c1->Print( "epi.pdf[", "pdf" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.05);

  c1->SetTitle("dE2 sigma" );
  gPad->SetLogx(1);
  gPad->SetLogy(1);

  TGraph * g1 = new TGraph( n, &vde[0], &veesig[0] );
  g1->SetTitle( "#DeltaE^{2} #sigma(#DeltaE)" );
  g1->SetLineWidth(2);
  g1->SetLineColor(2);
  g1->Draw("ac");

  TGraph * g2 = new TGraph( m, &wde[0], &weesig[0] );
  g2->SetLineWidth(2);
  g2->SetLineColor(1);
  g2->Draw("c");

  TLegend * lgnd = new TLegend( 0.60, 0.4, 0.94, 0.5 );
  lgnd->AddEntry( g1, "  5 GeV e", "l" );
  lgnd->AddEntry( g2, " 45 GeV pi", "l" );
  lgnd->Draw( "same" );

  c1->Print( "epi.pdf" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1->Print( "epi.pdf]" ); // ] closes file
  cout << "evince epi.pdf" << endl;

}
