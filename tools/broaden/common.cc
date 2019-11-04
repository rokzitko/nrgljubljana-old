double **buffers; // binary data buffers
int *sizes; // sizes of buffers

typedef map<double, double> mapdd;
typedef vector<double> vec;

mapdd spec; // Spectrum
unsigned int nr_spec; // Number of raw spectrum points
vec vfreq, vspec; // Same info as spectrum, but in vector<double> form

string tostring(int i)
{
  ostringstream S;
  S << i;
  return S.str();
}

// Load a file containing binary representation of raw spectral density.
// The grid is not assumed to be uniform.
void load(int i)
{
  // i-th z-value defines the name of the directory where the results of
  // the NRG calculation are contained.
  string filename;
  if (one && Nz == 1) {
    filename = name;
  } else {
    filename = tostring(i) + "/" + name;
  }
  ifstream f(filename.c_str(), ios::in | ios::binary);
  if (!f.good() || f.eof() || !f.is_open()) {
    cerr << "Error opening file " << filename << endl;
    exit(1);
  }
  if (verbose) {
    cout << "Reading " << filename << endl;
  }

  const int rows = 1+nrcol; // number of elements in a line

  // Determine the number of records
  f.seekg(0, ios::beg);
  const ios::pos_type begin_pos = f.tellg();
  f.seekg(0, ios::end);
  const ios::pos_type end_pos = f.tellg();
  const long len = end_pos-begin_pos;
  assert(len % (rows*sizeof(double)) == 0);
  const int nr = len / (rows*sizeof(double)); // number of lines
  if (verbose) {
    cout << "len=" << len
	 << " nr=" << nr << " data points" << endl;
  }

  // Allocate the read buffer. The data will be kept in memory for the
  // duration of the calculation!
  double *buffer = new double[rows*nr];
  f.seekg(0, ios::beg); // Return to the beginning of the file.
  f.read((char*)buffer, len);
  if (f.fail()) {
    cerr << "Error reading " << filename << endl;
    exit(1);
  }
  f.close();

  // Keep record of the the buffer and its size.
  buffers[i] = buffer;
  sizes[i] = nr;

  if (verbose) {
    // Check normalization.
    double sum = 0.0;
    for (int j = 0; j < nr; j++)
      sum += buffer[rows*j+col];
    cout << "Weight=" << sum << endl;
  }
}

// Load all the input data.
void read_files()
{
  buffers = new double*[Nz+1];
  sizes = new int[Nz+1];

  for (int i = 1; i <= Nz; i++) {
    load(i);
  }
}

// Combine data from all NRG runs (z-averaging).
void merge()
{
  const int rows = 1+nrcol; // number of elements in a line

  // Sum weight corresponding to equal frequencies.  Map of
  // (frequency,weight) pairs is used for this purpose.
  for (int i = 1; i <= Nz; i++) {
    for (int l = 0; l < sizes[i]; l++) {
      double & freq = buffers[i][rows*l];
      double & value = buffers[i][rows*l+col];
      mapdd::iterator I = spec.find(freq);
      if (I == spec.end()) {
	spec[freq] = value;
      } else {
	I->second += value;
      }
    }
  }

  nr_spec = spec.size();
  if (verbose) {
    cout << nr_spec << " unique frequencies." << endl;
  }

  // Normalize weight by 1/Nz, determine total weight, and store the
  // (frequency,weight) data in the form of linear vectors for faster
  // access in the ensuing calculations.
  double sum = 0.0;
  for (mapdd::iterator I = spec.begin(); I != spec.end(); I++) {
    const double weight = (I->second /= Nz); // Normalize weight on the fly
    const double freq = I->first;
    vfreq.push_back(freq);
    vspec.push_back(weight);
    sum += weight;
  }
  if (verbose) {
    cout << "Total weight=" << sum << endl;
  }
  assert(vfreq.size() == nr_spec && vspec.size() == nr_spec);
}

