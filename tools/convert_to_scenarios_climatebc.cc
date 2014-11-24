#include "common.h"
#include <map>

#define TIMESOFYEAR 17

/* Purpose of this program:
   To map a set of files onto one NetCDF file in the format desired
   What it will do:
   - Rename variables
   - Convert variables from float to double
   - Scale and multiply the input
   - Convert orography to sea-land mask file
   - Convert 1st file's lat/lon to latitude and longitude fields
   - Merge it all up into one big data file
   What it won't do:
   - Extract fields from 4D files
   - Combine fields to create a single field
*/



  /* Need:
     - Diurnal temperature range
     - Wind speed
     - Relative humidity 
     - Vapour pressure
     - Sea ice (need to download sit not sic)
      - Downloading
     - Evap (not possible without surface details)

     For DITR:
      - Calculate Tmax - Tmin
    
     For wind speed:
      - Combine u and v wind components using pythagorean theorem

     For RH:
      - Nice research, but a better idea might be to just retrieve the raw data...

      - es0 = 6.1078hPa (reference saturation vapour pressure at 0C)
        - As per Smithsonian Meteorological Tables, Sixth Revised Edition (pp. 351)
      - T0 = 273.15 (reference temperature, in K)
      - Td = dewpoint (in K)
      - T = temperature (in K)
      - lv = latent heat of vapourization (in j/kg)
        - 2.5008 E6 - 2.3E3 * (T - 273.15) (T in K)
	- Atmosphere-Ocean Dynamics, Adrian E. Gill (pp. 607)
      - Rv = gas constant for water vapour (j / (K * kg)
        - 461.5 (checked from Smithsonian Meteorological Tables. Sixth Revised Edition (pp. 289)
      Clausius-Clapyron equation (rearranged)
      - e = es0 * exp( lv/Rv * (1/T0 - 1/T))
      Arden L. Buck (1981): http://ams.allenpress.com/perlserv/?request=get-pdf&doi=10.1175%2F1520-0450%281981%29020%3C1527%3ANEFCVP%3E2.0.CO%3B2
      - New equation for saturation vapour pressure
      - p is pressure in hPa (millibars)
      - T is in degrees C
      - es = (.0007 + (3.46E-6 * p)) * 6.1121 * exp(17.502T / 240.97 + T)
      Relative Humidity:
      - RH = 100% * e(T)/es(T)
      Vapour pressure
      - Vapour pressure is e(T)

     as per http://www.faqs.org/faqs/meteorology/temp-dewpoint/

  */

enum TRANS_TYPE{AMOUNT,AVERAGE};

#define LAND_LEVEL 0.01

class VarTrans {
public:
  VarTrans() { add_factor = 0; multiplier = 0; name = "ERROR"; }
  VarTrans(string name, double add_factor, double multiplier): name(name), add_factor(add_factor), multiplier(multiplier) { }
  string name;
  double add_factor;
  double multiplier;
};

void add_var_trans_entries(map<string, VarTrans>& var_trans, map<string, VarTrans>& var_trans_future) {
  // Unit conversions
  // Total cloud cover fraction
  var_trans["clt"] = VarTrans("tcld", 0, 1);

  // Total cloud cover fraction (percent of baseline)
  var_trans_future["clt"] = VarTrans("tcld", 0, 1);

  // Diurnal temperature range (derived; tasmax - tasmin)
  var_trans["ditr"] = VarTrans("ditr", 0, 1);
  var_trans_future["ditr"] = VarTrans("ditr", 0, 1);

  // 500hPa geopotential height (extracted)
  var_trans["zg500"] = VarTrans("h500", 0, 1);
  var_trans_future["zg500"] = VarTrans("h500", 0, 1);

  // Surface relative humidity (extracted)
  var_trans["hur1000"] = VarTrans("rhum", 0, 1);
  var_trans_future["hur1000"] = VarTrans("rhum", 0, 1);

  // Surface specific humidity
  var_trans["huss"] = VarTrans("shum", 0, 1);
  var_trans_future["huss"] = VarTrans("shum", 0, 1);

  // Soil moisture content
  var_trans["mrso"] = VarTrans("somm", 0, 1); 
  var_trans_future["mrso"] = VarTrans("somm", 0, 1); 

  // Precipitation; convert from mm/month to mm/day by dividing by 30 (approximate)
  var_trans["pr"] = VarTrans("prec", 0, /*1.0/30.0*/ 86400);

  // Precipitation (percent of baseline)
  var_trans_future["pr"] = VarTrans("prec", 0, /*1*/ 86400);

  // Mean sea level pressure (convert to hPa by dividing by 100)
  var_trans["psl"] = VarTrans("mslp", 0, 0.01);
  var_trans_future["psl"] = VarTrans("mslp", 0, 0.01);

  // Surface shortware radiation
  var_trans["rsds"] = VarTrans("irad", 0, 1);
  var_trans_future["rsds"] = VarTrans("irad", 0, 1);

  // Sea ice thickness (convert from m to kg/m^2 by multiplying by mean density of sea ice (0.91MT/m^3 == 910kg/m^3)
  var_trans["sit"] = VarTrans("sice", 0, 910);
  var_trans_future["sit"] = VarTrans("sice", 0, 910);

  // Snow depth
  var_trans["snd"] = VarTrans("snod", 0, 1);
  var_trans_future["snd"] = VarTrans("snod", 0, 1);

  // Snow melt (convert from mm/s to mm/day by multiplying by 86400)
  var_trans["snm"] = VarTrans("melt", 0, 86400);
  var_trans_future["snm"] = VarTrans("melt", 0, 86400);

  // Snow amount
  var_trans["snw"] = VarTrans("snow", 0, 1);
  var_trans_future["snw"] = VarTrans("snow", 0, 1);

  // This needs to be calculated
  var_trans["soil"] = VarTrans("soil", 0, 1);
  var_trans_future["soil"] = VarTrans("soil", 0, 1);

  // Mean Temperature (convert from C*10 to C by dividing by 10)
  var_trans["tas"] = VarTrans("temp", 0, /*0.1*/ 1);
  var_trans_future["tas"] = VarTrans("temp", 0, /*0.1*/ 1);

  // Mean Daily Maximum Temperature (convert from C*10 to C by dividing by 10)
  var_trans["tasmax"] = VarTrans("tmax", 0, /*0.1*/ 1);
  var_trans_future["tasmax"] = VarTrans("tmax", 0, /*0.1*/ 1);

  // Mean Daily Minimum Temperature (convert from C*10 to C by dividing by 10)
  var_trans["tasmin"] = VarTrans("tmin", 0, /*0.1*/ 1);
  var_trans_future["tasmin"] = VarTrans("tmin", 0, /*0.1*/ 1);

  // Mean Surface Temperature (convert from K to C by subtracting 273.15)
  var_trans["ts"] = VarTrans("surt", -273.15, 1);
  var_trans_future["ts"] = VarTrans("surt", 0, 1);

  // Vapor Pressure (of water; compute using Buck eqn referenced below)
  var_trans["vapp"] = VarTrans("vapp", 0, 1);
  var_trans_future["vapp"] = VarTrans("vapp", 0, 1);

  // Wind speed (needs to be calculated)
  var_trans["wind"] = VarTrans("wind", 0, 1);

  // Wind speed (needs to be calculated) (percent of baseline)
  var_trans_future["wind"] = VarTrans("wind", 0, 1);

  // Can't have evap (no data)
}

template<typename T>
void transmogrify_ll(T* dat, const int num_lat, const int num_long) {
  const int max_lat = num_lat - 1;
  const int mid_lat = num_lat / 2;
  const int mid_long = num_long / 2;
  for(int top = 0, bottom = max_lat; top < mid_lat; top++, bottom--) {
    const int top_off = top * num_long;
    const int bottom_off = bottom * num_long;
    for(int left = 0, right = mid_long; left < mid_long; left++, right++) {
      swap(dat[top_off + left], dat[bottom_off + right]);
      swap(dat[bottom_off + left], dat[top_off + right]);
    }
    if(num_long % 2 == 1) {
      // Odd # of longs
      swap(dat[top_off + mid_long], dat[bottom_off + mid_long]);
    }
  }
  if(num_lat % 2 == 1) {
    // Odd # of lats
    const int mid_off = mid_lat * num_long;
    for(int left = 0, right = mid_long; left < mid_long; left++, right++) {
      swap(dat[mid_off + left], dat[mid_off + right]);
    }
  }
}

double u_lon_to_s(double lon) {
  return fmod((lon + 180.0f), 360.0f) - 180.0f;
}

string create_daterange_string(FileRecord& fr) {
  string filesub = fr.filename.substr(0, fr.filename.find_last_of('.'));
  boost::char_separator<char> sep("-");
  tokenizer<char_separator<char> > tok(filesub, sep);
  vector<string> tokens(tok.begin(), tok.end());
  assert(tokens.size() >= 6);
  return tokens[YEARSTART] + "_" + tokens[YEAREND];
}

void copy_lats_longs(NcFile& in, NcFile& out) {
  NcVar* lats = copy_var(in.get_var("lat"), out);
  NcVar* longs = copy_var(in.get_var("lon"), out);
  double* latdat = new double[lats->num_vals()];
  double* longdat = new double[longs->num_vals()];
  assert(lats->get(latdat, lats->num_vals()));
  assert(longs->get(longdat, longs->num_vals()));
  
  // Make sure that the lats are actually reversed,
  // that there are an even # of lats and longs,
  // and that the latitude origin is at 0 (greenwich)
  assert(latdat[0] < latdat[lats->num_vals() - 1]);
  
  // Reverse lats so north is at the top
  const int mid_lat = lats->num_vals() / 2;
  const int max_lat = lats->num_vals() - 1;
  for(int i = 0; i < mid_lat; i++)
    swap(latdat[i], latdat[max_lat - i]);
  
  // Move longs around so 180W is on the left edge
  const int mid_long = longs->num_vals() / 2;
  for(int i = 0; i < mid_long; i++)
    swap(longdat[i], longdat[i + mid_long]);

  // Convert unsigned (0 to 360) to signed (-180 to 180)
  for(int i = 0; i < longs->num_vals(); i++)
    longdat[i] = u_lon_to_s(longdat[i]);

  assert(lats->put(latdat, lats->num_vals()));
  assert(longs->put(longdat, longs->num_vals()));
  
  lats->rename("lats");
  longs->rename("longs");
  
  delete[] latdat;
  delete[] longdat;
}

void copy_lats_longs_bnds(NcFile& in, NcFile& out) {
  NcVar* lats = copy_var(in.get_var("lat_bnds"), out);
  NcVar* longs = copy_var(in.get_var("lon_bnds"), out);
  assert(lats && longs);

  long* lat_edges = lats->edges();
  int lat_recsize = get_recsize_and_edges(lats, lat_edges);
  long* long_edges = longs->edges();
  int long_recsize = get_recsize_and_edges(longs, long_edges);

  double* latdat = new double[lat_recsize];
  double* longdat = new double[long_recsize];
  assert(lats->get(latdat, lat_edges));
  assert(longs->get(longdat, long_edges));

  
  // Reverse lats so north is at the top
  const int mid_lat = lat_recsize / 2;
  const int max_lat = lat_recsize - 2;
  for(int i = 0; i < mid_lat; i += 2) {
    printf("%f, %f\n", latdat[i], latdat[max_lat - i]);
    swap(latdat[i], latdat[max_lat - i]);
    printf("%f, %f\n", latdat[i], latdat[max_lat - i]);
    swap(latdat[i + 1], latdat[max_lat - i + 1]);
  }
  
  // Move longs around so 180W is on the left edge
  const int mid_long = long_recsize / 2;
  for(int i = 0; i < mid_long; i += 2) {
    printf("%f, %f\n", longdat[i], longdat[i + mid_long]);
    swap(longdat[i], longdat[i + mid_long]);
    printf("%f, %f\n", longdat[i], longdat[i + mid_long]);
    swap(longdat[i + 1],longdat[i + mid_long + 1]);
  }

  // Convert unsigned (0 to 360) to signed (-180 to 180)
  for(int i = 0; i < long_recsize; i++)
    longdat[i] = u_lon_to_s(longdat[i]);

  // Dirty
  if(longdat[long_recsize - 1] == -180 && longdat[long_recsize - 2] > 0) {
    longdat[long_recsize - 1] = 180;
  }
  
  assert(lats->put(latdat, lat_edges));
  assert(longs->put(longdat, long_edges));
  
  delete[] latdat;
  delete[] longdat;
  delete[] lat_edges;
  delete[] long_edges;
}

void add_slmask(NcVar* invar, NcFile& out) {
  long* edges = invar->edges();
  int recsize = get_recsize_and_edges(invar, edges);
  int framesize = recsize / TIMESOFYEAR;
  float* data = new float[recsize];
  int* slmask = new int[framesize];

  // Create slmask variable
  NcDim* lat = out.get_dim("yc");
  NcDim* lon = out.get_dim("xc");
  assert(lat && lon);
  NcVar* slmaskvar = out.add_var("slmask", ncInt, lat, lon);
  assert(slmaskvar);
  float min = 1E20;
  float max = -1;

  // Compute sea-land mask
  assert(invar->get(data, edges));
  for(int i = 0; i < framesize; i++) {
    min = (min > data[i]) ? data[i] : min;
    max = (max < data[i]) ? data[i] : max;
  }
  float threshold = (min + max) / 2;
  for(int i = 0; i < framesize; i++) {
    slmask[i] = (data[i] >= threshold);
  }

  /*
  const int num_lat = lat->size();
  const int num_long = lon->size();
  transmogrify_ll(slmask, num_lat, num_long);
  */

  // Store sea-land mask
  long* slmask_edges = slmaskvar->edges();
  assert(slmaskvar->put(slmask, slmask_edges));

  delete[] edges;
  delete[] slmask_edges;
  delete[] data;
  delete[] slmask;
}

class VarFile {
public:
  VarFile(string file, NcVar* var, VarTrans& vt): vt(vt) { this->file = file; this->var = var; }
  string file;
  NcVar* var;
  VarTrans& vt;
};

int main(int argc, char** argv) {
  char buf[1024];
  list<FileRecord> l;
  bool first = true;
  bool seen_sftlf = false;
  map<string, VarTrans> var_trans;
  map<string, VarTrans> var_trans_future;
  map<string, string> date_trans;
  map<string, string> expt_trans;
  list<VarFile> varfile;
  list<VarFile>::const_iterator vfit;

  NcDim* rows = 0;
  NcDim* columns = 0;
  NcDim* timesofyear = 0;

  expt_trans["sresa1b"] = "A1B";
  expt_trans["sresa2"] = "A2";
  expt_trans["sresb1"] = "B1";
  expt_trans["hist"] = "HIST";

  // Translate real dates into climatology dates
  date_trans["1961_1990"] = "1961_1990";
  date_trans["1971_2000"] = "1971_2000";
  date_trans["1981_2000"] = "1981_2000";
  date_trans["2010_2039"] = "2020";
  date_trans["2040_2069"] = "2050";
  date_trans["2070_2099"] = "2080";

  add_var_trans_entries(var_trans, var_trans_future);

  ncopts = NC_VERBOSE;

  if(argc < 2) {
    printf("Usage: convert_to_scenarios <output_file>\n");
  }

  // Open new NetCDF file for writing
  NcFile out(argv[1], NcFile::Replace, NULL, 0, NcFile::Offset64Bits);
  assert(out.is_valid());
  out.set_fill(NcFile::NoFill);

  while(fgets(buf, 1024, stdin)) {
    // Chomp a la perl
    *(strchr(buf, '\n')) = '\0';

    string filename = buf;
    FileRecord fr(buf, false);
    if(!fr.is_ok) {
      printf("Failed to open file %s\n", buf);
      continue;
    }

    NcFile& in = *(fr.f);
    string daterange = create_daterange_string(fr);
    
    // From the first file...
    if(first) {
      first = false;

      // Copy all global attributes
      copy_atts(&in, &out);

      // Copy dimensions
      copy_dims(in, out);

      // Get dimensions
      timesofyear = out.get_dim("timeofyear");
      rows = out.get_dim("yc");
      columns = out.get_dim("xc");

      assert(timesofyear && rows && columns);

      NcVar* lats = copy_var(in.get_var("xc"), out);
      NcVar* longs = copy_var(in.get_var("yc"), out);
      //copy_lats_longs(in, out);
      //copy_lats_longs_bnds(in, out);

      NcVar* albers = in.get_var("albers_conical_equal_area");
      assert(albers);
      copy_var(albers, out);
    }

    // Parse out filename
    boost::char_separator<char> sep("-");
    string f = fr.file.substr(0, fr.file.rfind("."));
    tokenizer<char_separator<char> > tok(f, sep);
    vector<string> tokens(tok.begin(), tok.end());

    // Compare variable name...
    NcVar* invar = in.get_var(fr.var.c_str());
    assert(invar);
    NcDim* inlats = in.get_dim("yc");
    NcDim* inlongs = in.get_dim("xc");

    if(!inlats || !inlongs) {
      continue;
    }

    if(inlats->size() != rows->size() || inlongs->size() != columns->size()) {
      printf("Dimensions of input and output data do not match; file: %s\n", buf);
      continue;
    }
    
    if(fr.var == "sftlf") {
    } else {
      if(!seen_sftlf) {
	add_slmask(invar, out);
	seen_sftlf = true;
      } 

      // Create new var
      list<string> expts;
      list<string>::const_iterator expt;
      if(var_trans.find(fr.var) == var_trans.end()) {
	continue;
      }

      VarTrans& var = (fr.expt == "20c3m" || fr.expt == "hist") ? var_trans[fr.var] : var_trans_future[fr.var];
      if(fr.expt == "20c3m") {
	expts.push_back(expt_trans["sresa1b"]);
	expts.push_back(expt_trans["sresa2"]);
	expts.push_back(expt_trans["sresb1"]);
      } else {
        assert(expt_trans.find(fr.expt) != expt_trans.end());
	expts.push_back(expt_trans[fr.expt]);
      }

      for(expt = expts.begin(); expt != expts.end(); ++expt) {
	string date_translated;
	if(date_trans.find(daterange) == date_trans.end()) {
	  date_translated = daterange;
	} else {
	  date_translated = date_trans[daterange];
	}

	// Create variable name
	string varname = *expt + "-" + fr.run + "_" + date_translated + "_" + var.name;

	printf("Old var: %s, Var: %s\n", fr.var.c_str(), varname.c_str());

	// Get dimensions to use
	vector<NcDim*> dims;
	populate_dimvec(invar, out, dims);
      
	// Create output variable
	NcVar* outvar = out.add_var(varname.c_str(), ncFloat, (int)dims.size(), (const NcDim**)&dims[0]);
	copy_atts(invar, outvar);

	varfile.push_back(VarFile(filename, outvar, var));
      }
    }
  }

  if(!first) {
    rows->rename("rows");
    columns->rename("columns");
    timesofyear->rename("timesofyear");
  }

  printf("Step 1 complete\n");

  for(vfit = varfile.begin(); vfit != varfile.end(); vfit++) {
    NcVar* outvar = (*vfit).var;
    VarTrans& var = (*vfit).vt;
    FileRecord fr((*vfit).file, false);
    if(!fr.is_ok) {
      printf("Failed to open file %s\n", buf);
      continue;
    }

    printf("foo\n");

    NcFile& in = *(fr.f);
    NcVar* invar = in.get_var(fr.var.c_str());
    long* edges = invar->edges();
    int recsize = get_recsize_and_edges(invar, edges);

    // Allocate space for data
    float* indat = new float[recsize];
    assert(invar->get(indat, edges));
    
    const double m = var.multiplier;
    const double a = var.add_factor;
    for(int i = 0; i < recsize; i++) {
      indat[i] = indat[i] * m + a;
    }

    // Output old data
    assert(outvar->put(indat, edges));
	
    delete[] indat;
    delete[] edges;
  }

  out.close();

  printf("Step 2 complete\n");
}
