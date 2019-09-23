
// Daniel Pitzl, Sep 2019
// simulate multiple electron inelastic scattering (ionization energy loss)
// small compared to elastic nuclear Coulomb scattering
// rlq simscat.C+

#include <iostream> // cout
#include <fstream> // files
#include <time.h> // clock_gettime
#include <stack>
#include <map>

#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TRandom1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

//------------------------------------------------------------------------------
void simscat()
{
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
  vde.reserve(1699);
  vector <double> vlogde;
  vlogde.reserve(1699);

  vector <double> vesig;
  vesig.reserve(1699);

  while( ! sFile.eof() ) {

    int j;
    double de;
    double eesig;
    sFile >> j;
    sFile >> de;
    sFile >> eesig; // this is dE^2 dsig/ddE
    vde.push_back(de);
    vlogde.push_back(log(de));
    // dsig/dlndE = dE dsig/ddE
    vesig.push_back( eesig/de );

  } // while

  int n = vde.size();
  cout << "  read " << n << " values" << endl;

  gPad->SetLogx(0);
  gPad->SetLogy(0);
  TGraph * g4 = new TGraph( n, &vlogde[0], &vesig[0] );
  g4->SetTitle( "d#sigma/dlog(#DeltaE)" );
  g4->Draw("ac"); // M-shell dominates

  // inversion method for dsigma/dlog(dE) = dE dsigma/ddE
  // equidstant bins in log(dE)

  double x0 = vlogde[0]; // first bin is half-bin
  double sum = 0;
  vector <double> vsum(n);
  vector <double> vx(n);
  for( int i = 0; i < n-1; ++i ) {
    double x = vlogde[i];
    double xl = 0.5 * ( x0 +x );
    double xr = 0.5 * ( x + vlogde[i+1] );
    vx[i] = x;
    sum += vesig[i]*(xr-xl);
    vsum[i] = sum;
    x0 = x;
  }
  int l = n-1; // last half-bin
  double dx = 0.5*(vlogde[l-1]+vlogde[l]);
  sum  += vesig[l]*dx;
  vsum[l] = sum;
  double path = 1E4/sum; // [um]
  cout << "sum " << sum << ", path = " << path << " um" << endl; // 0.26 um

  // normalize integral, put into look-up map:

  static map < double, double > xmap;

  for( int i = 0; i < n; ++i )
    xmap[vsum[i]/sum] = vlogde[i]; // norm, insert in map

  double r0 = vsum[0]/sum; // first bin

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TFile * hFile = new TFile( "simscat.root", "RECREATE" );

  double thick = 900; // [um]

  TH1I honede( "onede", "single E loss;log_{10}(single energy loss [eV]);scatters",
	       120, 0, 6 );

  TH1I hsumde( "sumde", Form( "E loss in %d um;energy loss [keV] in %d um Si;tracks",
			      int(thick+0.1), int(thick+0.1) ),
	       500, 0, 3*thick );

  TH1I hangmax( "angmax", "max delta ray angle;max delta ray angle [deg]);tracks",
		90, 0, 90 );

  // Highland: 13.6 mrad/GeV*sqrt(d/X0)

  double sx0 = sqrt( thick / 9.36e4 );

  TH1I hsumpx( "sumpx", Form( "px in %d um;px [eV] in %d um Si;tracks",
			      int(thick+0.1), int(thick+0.1) ),
	       200, -2000*sx0, 2000*sx0 );
  TH1I hsumpy( "sumpy", Form( "py in %d um;py [eV] in %d um Si;tracks",
			      int(thick+0.1), int(thick+0.1) ),
	       200, -2000*sx0, 2000*sx0 );

  double log10 = log(10);
  double twome = 1.022E6; // [eV]
  double pi = 4*atan(1);
  double wt = 180/pi;

  TRandom1 * myRandom = new TRandom1(); // 1 = Ranlux

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  timespec ts;
  clock_gettime( CLOCK_REALTIME, &ts );
  long sc0 = ts.tv_sec;
  long ns0 = ts.tv_nsec;

  for( int i = 0; i < 10*1000; ++i ) { // tracks

    if( i%10000 == 0 ) cout << " " << i << flush;

    double z = 0;
    double sumde = 0;
    double maxsin = 0;
    double sumpx = 0;
    double sumpy = 0;

    while( z < thick ) {

      double r0 = myRandom->Rndm();
      if( r0 < 1e-25 ) // avoid zero
	r0 = 1e-25;

      z += -path * log(r0); // exponential path distribution

      double r1 = myRandom->Rndm(); // [0,1]
      auto it = xmap.upper_bound(r1);
      double de = vde[0];
      if( it != xmap.begin() ) {
	--it;
	de = exp( it->second );
      }

      if( i < 10*1000 )
	honede.Fill( log(de)/log10 );

      sumde += de; // [eV]

      double cost = sqrt( de / (twome + de) ); // elastic kinematics
      double sint = sqrt( 1 - cost*cost ); // mostly 90 deg
      double phi = 2*pi*myRandom->Rndm();

      sumpx += sint * cos(phi);
      sumpy += sint * sin(phi);

    } // while thick

    hsumde.Fill( sumde*1E-3 ); // [keV]
    hangmax.Fill( asin(maxsin)*wt );
    hsumpx.Fill( sumpx ); // 400 rad*eV * sqrt(t/X0) = 0.4 mrad*MeV * sqrt(t/X0)
    hsumpy.Fill( sumpy ); // elastic Coulomb scattering dominates

  } // tracks i

  cout << endl;

  clock_gettime( CLOCK_REALTIME, &ts );
  cout << "duration " << ts.tv_sec - sc0 + ( ts.tv_nsec - ns0 )*1E-9 << endl; // 12.8692, 13.3

  hFile->Write();
  hFile->ls();
  cout << "size " << hFile->GetSize() << " B"
       << " = " << hFile->GetSize()/1024 << " kB"
       << " = " << hFile->GetSize()/1024/1024 << " MB"
       << endl;
  cout << "rlq " << hFile->GetName() << endl;

}
