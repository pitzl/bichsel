
// extend sigmade table
{
  double e1249 = 1331586.625;
  double e1250 = 1346086.625; // [eV]
  double f = e1250/e1249; // step facor
  double sig = 178346.578; // asymptotic constant

  for( int j = 1; j < 450; ++j )
    cout << " " << 1250+j
	 << "  " << e1250 * pow(f,j)
	 << "  " << sig
	 << endl;
}
