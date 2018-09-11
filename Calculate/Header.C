//Macro containing variables used in RnPoVsCell.C and RnPoVsTime.C


using namespace std;
const int NUMCELLS = 154;
const double POLIFETIME = 2.57;   //[ms]
const double TIMEWINDOW = 5.0*POLIFETIME, TIMEOFFSET = 10.0*POLIFETIME;    //[ms]

//const int NUMEXCLUDECELLS = 32;
//int ExcludeCellArr[NUMEXCLUDECELLS] = {0, 1, 2, 3, 5, 6, 9, 10, 11, 12, 13, 18, 21, 23, 27, 31, 32, 34, 41, 44, 48, 52, 56, 63, 69, 79, 86, 87, 115, 122, 127, 139};

// Phys_Neutrino_v2
//const int NUMEXCLUDECELLS = 31;
//int ExcludeCellArr[NUMEXCLUDECELLS] = {0,1,2,3,4,5,6,9,10,11,12,13,18,21,23,24,27,32,34,40,44,52,68,79,86,102,115,122,127,130,139};

// Phys_NuFact all runs
const int NUMEXCLUDECELLS = 40;
int ExcludeCellArr[NUMEXCLUDECELLS] = {0,1,2,3,4,5,6,9,10,11,12,13,18,21,23,24,27,31,32,34,40,41,44,48,52,56,63,68,69,79,86,87,102,111,115,122,127,130,136,139};

// Phys_NuFact all runs + ET's
//const int NUMEXCLUDECELLS = 68;
//int ExcludeCellArr[NUMEXCLUDECELLS] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,18,21,23,24,27,28,31,32,34,40,41,42,44,48,52,55,56,63,68,69,70,79,83,84,86,87,97,98,102,111,112,115,122,125,126,127,130,136,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153};


//const int NUMEXCLUDECELLS = 153;
//int ExcludeCellArr[NUMEXCLUDECELLS] = {0,1,2,3,4,5,6,9,10,11,12,13,18,21,23,24,27,31,32,34,40,41,44,48,52,56,63,68,69,79,86,87,102,111,115,122,127,130,136,139,7,   8,   14,  15,  16,  17,  19,  20,  22,  25,  26,  28,  29,  30,  33,  35,  36,  37,  38,  39,  42,  43,  45,  46,  47,  49,  50,  51,  53,  54,  55,  57,  58,59,  60,  61,  62,  64,  65,  66,  67,  70,  71,  72,  73,  74,  75,  76,  77,  78,  80,  81,  82,  83,  84,  85,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,98,  99,  100, 101, 103,105, 106, 107, 108, 109, 110, 112, 113, 114, 116, 117, 118, 119, 120, 121, 123, 124, 125, 126, 128, 129, 131, 132, 133, 134, 135, 137,138, 140, 141, 142, 143, 145, 146, 147, 148, 149, 150, 151, 152, 153};


// 7,   8,   14,  15,  16,  17,  19,  20,  22,  25,  26,  28,  29,  30,  33,  35,  36,  37,  38,  39,  42,  43,  45,  46,  47,  49,  50,  51,  53,  54,  55,  57,  58,  
// 59,  60,  61,  62,  64,  65,  66,  67,  70,  71,  72,  73,  74,  75,  76,  77,  78,  80,  81,  82,  83,  84,  85,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97, 
// 98,  99,  100, 101, 103, 104, 105, 106, 107, 108, 109, 110, 112, 113, 114, 116, 117, 118, 119, 120, 121, 123, 124, 125, 126, 128, 129, 131, 132, 133, 134, 135, 137, 
// 138, 140, 141, 142, 143, 145, 146, 147, 148, 149, 150, 151, 152, 153  


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

double pileupVetoT = 800*(1e-6);	//[ms]




