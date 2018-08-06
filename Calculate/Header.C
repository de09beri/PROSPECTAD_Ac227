//Macro containing variables used in RnPoVsCell.C and RnPoVsTime.C


using namespace std;
const int NUMCELLS = 154;
const double POLIFETIME = 2.57;   //[ms]
const double TIMEWINDOW = 5.0*POLIFETIME, TIMEOFFSET = 10.0*POLIFETIME;    //[ms]

const int NUMEXCLUDECELLS = 32;
int ExcludeCellArr[NUMEXCLUDECELLS] = {0, 1, 2, 3, 5, 6, 9, 10, 11, 12, 13, 18, 21, 23, 27, 31, 32, 34, 41, 44, 48, 52, 56, 63, 69, 79, 86, 87, 115, 122, 127, 139};


//Set up bins and ranges for histograms
double dtMin = 0.0, dtMax = TIMEWINDOW;	//[ms]
int numDtBins = (dtMax - dtMin)/0.1;
	
double PSDMin = 0.15, PSDMax = 0.37;
int numPSDBins = 150; 

double EnMin = 0.49, EnMax = 1.16;		//[MeVee]
int numEnBins = (EnMax - EnMin)/0.005;	// 5 keV bins
	
double dzMin = -250, dzMax = 250;	//[mm]
int numDzBins = (dzMax - dzMin)/2.5;	//0.25 cm bins

double posMin = -800, posMax = 800;	//[mm]
int numPosBins = (posMax - posMin)/10;	//1 cm bins 






