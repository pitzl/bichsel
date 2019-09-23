
// Daniel Pitzl, Aug 2019
// from pixelav.c
// simulate delta rays

// simdelta -n 100100 -t 150

#include <cstdlib> // atoi
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

using namespace std;

struct track {
  double x;
  double y;
  double z;
  double tet;
  double phi;
  double E;
  double rng;
};

//------------------------------------------------------------------------------
double range( double E ) // [um]([keV])
{
  // "practical" range of e in Si, with multiple scattering along the path
  // W. Blum, L. Rolandi: Particle Detection with Drift Chambers, Springer 1994 p 7
  // M. Brigida et al, NIMA533 (2004) 322-43

  double A = 1.9841; // um/keV
  double B = 0.9957;
  double C = 4.603E-3; // 1/keV
  return A*E * ( 1 - B / ( 1 + C*E ) ); // um
}

//------------------------------------------------------------------------------
//void simdelta()
int main( int argc, char* argv[] )
{
  int ntr = 100*1000;

  double thck = 150; // [um]

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-n" ) )
      ntr = atoi( argv[++i] ); // last event

    if( !strcmp( argv[i], "-t" ) )
      thck = atof( argv[++i] ); // [um]

  } // argc

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cout << "try to open eesig-e5000.dat";
  ifstream sFile( "eesig-e5000.dat" );

  if( !sFile ) {
    cout << " : failed " << endl;
    return 1;
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

  vector <double> veesig;
  veesig.reserve(1699);

  vector <double> velsig;
  velsig.reserve(1699);

  vector <double> vesig;
  vesig.reserve(1699);

  vector <double> vsig;
  vsig.reserve(1699);

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
    veesig.push_back( eesig );
    vesig.push_back( eesig/de );
    vsig.push_back( eesig/de/de );
    velsig.push_back( log(de)*eesig/de );

  } // while
  int n = vde.size();
  cout << "  read " << n << " values" << endl;

  TCanvas * c0 = new TCanvas;
  c0->SetTitle("dE2 sigma" );
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  TGraph * g0 = new TGraph( n, &vde[0], &veesig[0] );
  g0->SetTitle( "#DeltaE^{2} #sigma(#DeltaE)" );
  g0->Draw("ac");

  TCanvas * c1 = new TCanvas;
  c1->SetTitle("dE sigma" );
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  TGraph * g1 = new TGraph( n, &vde[0], &vesig[0] );
  g1->SetTitle( "#DeltaE #sigma(#DeltaE)" );
  g1->Draw("ac");

  TCanvas * c2 = new TCanvas;
  c2->SetTitle("sigma" );
  c2->SetName("sigma" );
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  TGraph * g2 = new TGraph( n, &vde[0], &vsig[0] );
  g2->SetTitle( "#sigma(#DeltaE)" );
  g2->Draw("ac");

  TCanvas * c3 = new TCanvas;
  c3->SetTitle("ln(dE) dE sigma" );
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  TGraph * g3 = new TGraph( n, &vde[0], &velsig[0] );
  g3->SetTitle( "d#sigma/dlog(ln(#DeltaE)) = ln(#DeltaE) #DeltaE d#sigma/d#DeltaE" );
  g3->Draw("ac");

  TCanvas * c4 = new TCanvas;
  c4->SetTitle("dsigma/dlog(dE)" );
  gPad->SetLogx(0);
  gPad->SetLogy(0);
  TGraph * g4 = new TGraph( n, &vlogde[0], &vesig[0] );
  g4->SetTitle( "d#sigma/dlog(#DeltaE)" );
  g4->Draw("ac"); // M-shell dominates

  // inversion method for dsigma/dlog(dE) = dE dsigma/ddE
  // non-equidstant bins in dE and log(dE)

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

  TFile * hFile = new TFile( "simdelta.root", "RECREATE" );

  TH1I honede( "onede", "single E loss;log_{10}(single energy loss [eV]);scatters",
	       120, 0, 6 );

  TH1I hrng( "range", "delta ray range;log_{10}(delta ray range [um]);rays", 80, 0, 4 );

  TH1I hsumde( "sumde", Form( "E loss in %d um;energy loss [keV] in %d um Si;tracks",
			      int(thck+0.1), int(thck+0.1) ),
	       500, 0, 3*thck );

  TH1I hdemax( "demax", "max E loss;log_{10}(max energy loss [eV]);tracks", 120, 2, 8 );
  TH1I hrngmax( "rngmax", "max delta ray range;log_{10}(max delta ray range [um]);tracks",
		80, 0, 4 );
  TH1I hdxymax( "dxymax", "max delta ray radius;log_{10}(max delta ray radius [um]);tracks",
		80, 0, 4 );
  TH1I hangmax( "angmax", "max delta ray angle;max delta ray angle [deg]);tracks",
		90, 0, 90 );
  TH1I hndel( "ndel", "delta rays per track;delta rays;tracks", 20, -0.5, 19.5 );
  TH2I * hdxy = new
    TH2I( "dxy", "delta ray hits;x [um];y [um];delta ray hits",
	  400, -200, 200, 400, -200, 200 );

  TH1I hnx( "nx", Form( "ionizations in %d um;ionizations / %d um;tracks",
			int(thck+0.1), int(thck+0.1) ),
	    thck/path, 0, 5*int(thck/path) );
  TH1I hneh( "neh", "e-h pairs per track;e-h pairs;tracks", 500, 0, 500*thck );
  TH1I hmeh( "meh", "drifted e-h pairs per track;drifted e-h pairs;tracks",
	     500, 0, 50*thck );

  double log10 = log(10);
  double twome = 1.022E6; // [eV]
  double pi = 4*atan(1);
  double wt = 180/pi;
  double Weh = 3.645; // [ev] per eh pair

  TRandom1 * myRandom = new TRandom1(); // 1 = Ranlux

  bool Fill = 0; // 37% faster

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  timespec ts;
  clock_gettime( CLOCK_REALTIME, &ts );
  long sc0 = ts.tv_sec;
  long ns0 = ts.tv_nsec;

  for( int i = 0; i < ntr; ++i ) { // tracks

    if( i%10000 == 0 ) cout << " " << i << flush;

    double z = 0;
    int nx = 0;
    double sumde = 0;
    int neh = 0;
    int meh = 0;
    double maxde = 0;
    double maxrng = 0;
    double maxdxy = 0;
    double maxsin = 0;
    stack <track> deltas;

    while( z < thck ) {

      ++nx;

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

      if( Fill )
	honede.Fill( log(de)/log10 );

      sumde += de;

      double rng = range( de*1E-3 ); // [um]
      if( Fill )
	hrng.Fill( log(rng)/log10 );

      if( de > maxde ) {
	maxde = de;
	maxrng = rng;
      }

      double cost = sqrt( de / (twome + de) ); // elastic kinematics
      double sint = sqrt( 1 - cost*cost ); // mostly 90 deg
      double dxy = rng*sint;
      if( dxy > maxdxy ) {
	maxdxy = dxy;
	maxsin = sint;
      }

      if( rng < 1 ) { // [um] cutoff
	neh += de/Weh; // local ionization
      }
      else {
	// put delta ray on stack
	track t;
	t.x = 0;
	t.y = 0;
	t.z = z;
	t.tet = acos(cost); // 0..pi
	t.phi = 2*pi*myRandom->Rndm();
	t.E = de; // [eV]
	t.rng = rng; // [um]
	deltas.push(t);
      }

    } // while thck

    meh = neh/10; // sparse sampling for drift

    hsumde.Fill( sumde*1E-3 ); // [ke]
    hdemax.Fill( log(maxde)/log10 );
    hrngmax.Fill( log(maxrng)/log10 );
    hdxymax.Fill( log(maxdxy)/log10 );
    hangmax.Fill( asin(maxsin)*wt );
    hndel.Fill( deltas.size() ); // mean 0.3, max 4

    // deal with delta-e on the stack:
    // we need sigma(dE) for slow e
    // until then: equidistant e-h pairs along range

    while( ! deltas.empty() ) {

      track t = deltas.top();
      deltas.pop();
      int zeh = t.E / Weh;
      neh += zeh;
      double dl = t.rng/zeh;
      double dz = dl*cos(t.tet);
      double dr = dl*sin(t.tet);
      double dx = cos(t.phi)*dr;
      double dy = sin(t.phi)*dr;

      for( int k = 0; k < zeh; ++k ) { // propagate delta
	++nx;
	if( nx%10 == 1 ) ++meh;
	if( Fill )
	  hdxy->Fill( t.x, t.y );
	t.x += dx;
	t.y += dy;
	t.z += dz;
	if( t.z > thck ) break;
      }

    } // deltas

    hnx.Fill( nx );
    hneh.Fill( neh ); // peak 10.8 ke, FWHM 44%
    hmeh.Fill( meh );

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
