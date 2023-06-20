#define PI 3.14159265e0
#define Cvac 2.99792458e10   /* cm/s */ 
#define H 6.626196e-27    /* erg*sec */
#define k_B 1.380622e-16  /* erg/K */
#define KMAX  1200
#define TENLOG 2.30258509299405e0

/*             ODF-Kram           */
// Number of wavelength Intervals in the ODF
#define Nlam 328
// ODF sub-bins per wavelength interval,
#define Nbin 12
//  of frequency intervals used in the binning
#define Nnu 328
// Tau5000 bin, for a defined reference scale (tau_5000 continuum)
#define BIN5000 1
// FULL ODF, if this, skip sorting and just write out a file:
// Nnu - Nbin - NT - Np to be read in.
#define FULLODF 0
// Rosseland mean only
#define ROSS 1
// Planck mean only
#define PLANCK 0
// Grey bins (will be Ross only)
#define GREY 0
// Number of bands, set to 1 for grey, set to NBands = 328*12 = 3936
// for full ODF
#define Nbands 5
double Level_tau_bot[Nbands] = { -0.25, -1.5,-3.0,-99.0,-99};
double Level_tau_top[Nbands] = {99.0e0, -0.25,-1.5,-3.0,-3.0};
double Level_nu_bot[Nbands] = {  0, 0,   0,  0,175};
double Level_nu_top[Nbands] = {328, 328,328,175,328};
double TAU1_2[Nbands]       = {0.35, 0.35, 0.35, 0.35,0.35};
//
/* 12 bins, as in stagger
double Level_tau_bot[Nbands] = {- 1.46,- 3.81,-99.00,- 0.62,- 0.62,- 1.50,- 2.28,-99.00,- 0.62,- 1.50,- 2.28,-99.00};
double Level_tau_top[Nbands] = { 99.00,- 1.46,- 3.81, 99.00, 99.00,- 0.62,- 1.50,- 2.28, 99.00,- 0.62,- 1.50,- 2.28};
double Level_nu_bot[Nbands]  = {     0,     0,     0,   159,   176,   159,   159,   159,   256,   183,   191,   244};
double Level_nu_top[Nbands]  = {   159,   159,   159,   176,   256,   183,   191,   244,   328,   328,   328,   328};
double TAU1_2[Nbands]        = {  0.35, 0.35, 0.35,  0.35,  0.35, 0.35, 0.35, 0.35,  0.35, 0.35, 0.35, 0.35};
*/

// If full ODF define NBands = 328*12 = 3936
/* Some exampel bins, 2,4, and grey.
   Ross sets Rosseland mean only.
   Grey flag will output only one bin, using Rosseland mean. Grey 
   Over-rides both Nbands and Level.
   Nbands = levels output.
   Level = an array of Nbands-1 optical depths (relative to tau5000=1), which
   the bins are seperated */
/*
#define FULLODF 1
#define ROSS 0
#define PLANCK 0
#define GREY 0
#define Nbands 3936
*/
// 2 levels RT
/*
#define ROSS 0 
#define GREY 0
#define Nbands 2
double Level_tau_bot[Nbands] = {- 0.5,- 99.0};
double Level_tau_top[Nbands] = { 99.00,- 1.46};
double Level_nu_bot[Nbands]  = {     0,     0};
double Level_nu_top[Nbands]  = {   328, 328};
double TAU1_2[Nbands]        = {  0.35, 0.35};
*/

// 4 levels RT
/*
#define ROSS 0 
#define GREY 0
#define Nbands 4
double Level_tau_bot[Nbands] = {-0.5e0,-1.5e0,-3.0e0,-99.0};
double Level_tau_top[Nbands] = {99.0e0,-0.5e0,-1.5e0,-3.0e0};
double Level_nu_bot[Nbands] = {0,0,0,0};
double Level_nu_top[Nbands] = {328,328,328,328};
double TAU1_2[Nbands]       = {0.35, 0.35, 0.35,  0.35};
*/
// Splitting the last bin does the TR slightly better when scattering included
/*#define Nbands 5
double Level_tau_bot[Nbands] = { -0.25, -1.5,-3.0,-99.0,-99};
double Level_tau_top[Nbands] = {99.0e0, -0.25,-1.5,-3.0,-3.0};
double Level_nu_bot[Nbands] = {  0, 0,   0,  0,175};
double Level_nu_top[Nbands] = {328, 328,328,175,328};
double TAU1_2[Nbands]       = {0.35, 0.35, 0.35, 0.35,0.35};
*/

// Grey RT
/*
#define ROSS 1 
#define GREY 1
#define Nbands 4
double Level[Nbands-1] = {-0.1};
*/

// Temperature Axis - New ODF, remove [3.30] and 1 from NT for Old ODF
#define NT 57

double tab_T[NT]={3.30, 3.32, 3.34, 3.36, 3.38, 3.40, 3.42, 
                  3.44, 3.46, 3.48, 3.50, 3.52, 3.54, 3.56, 
                  3.58, 3.60, 3.62, 3.64, 3.66, 3.68, 3.70,
                  3.73, 3.76, 3.79, 3.82, 3.85, 3.88, 3.91, 
                  3.94, 3.97, 4.00, 4.05, 4.10, 4.15, 4.20, 
                  4.25, 4.30, 4.35, 4.40, 4.45, 4.50, 4.55, 
                  4.60, 4.65, 4.70, 4.75, 4.80, 4.85, 4.90, 
                  4.95, 5.00, 5.05, 5.10, 5.15, 5.20, 5.25, 
                  5.30};

// Pressure Axis, remove  [-4, -3.5, -3.0,-2.5] and 4 from Np for Old ODF
#define Np 25
double tab_p[Np]={-4., -3.5, -3., -2.5, -2., -1.5, -1., -.5, 0., .5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6. ,6.5, 7., 7.5, 8. };

// Weights for each ODF sub-bin
double wt[Nbin]={6./60.,6./60.,6./60.,6./60.,6./60.,6./60.,6./60.,
                  6./60.,6./60.,3./60.,2./60.,1./60.}; 

// ODF indice for the tau_5000 bin, arrays for kappa and B continuum.
int i_5000;

// Tau_5000 of the background atmosphere
double tau5000[KMAX];

// The minimum depth in the atmosphere (highest logtau5000)

int Kind;

/* Properties of the background atmosphere:
z - depth, dzi - z-spacing above [1] and below [0], kap0 - Reference opacity
p - pressure, tau0 - , rho - density, T - temperature */

double zax[KMAX],dzi[KMAX][2], kap0[KMAX], p[KMAX], tau0[KMAX], rho[KMAX], T[KMAX];
// B - blackbody emission, dBdT - derivative of blackbody with temperature 
double  B[Nnu][NT], dBdT[Nnu][NT];

// ODF, stored as log(kap)
short ODF[Nlam][NT][Np][Nbin] ;
/* acont_pT (wavelength,temperature,pressure) - Tabulated continuum opacity
   acont (height, frequency) - tabulated continuum opacity for the background
   atmosphere, scont = scattering, kcont = absorption */

double acont[KMAX][Nnu],scont[KMAX][Nnu],kcont[KMAX][Nnu];
float acont_pT[Nlam][NT][Np],scont_pT[Nlam][NT][Np],kcont_pT[Nlam][NT][Np];

// Begining and end wavelengths of each ODF sub-bin

double wbeg[Nlam] , wend[Nlam];

// middle frequency and frequency spacing of each sub-bin

double dnu[Nnu], numid[Nnu];

// ODF absorption Opacity for the background atmosphere, with frequency, sub-bin and depth.
double kap[Nnu][Nbin][KMAX],line_opacity[Nnu][Nbin][KMAX];
double ref_opacity[Nbands][KMAX];

// Tau grid for background atmosphere at each wavelength point and sub bin.
double tau_nu[Nnu][Nbin][KMAX];

/* Record what the z-height, z-index and what tau5000 is at the
point where the ODF gives tau = 1 */

int Index[Nnu][Nbin];

/* Record what the z-height, z-index and what tau5000 is at the
point where the ODF gives tau = 1 */
double z_at_tau_eq_1[Nnu][Nbin];
int    i_at_tau_eq_1[Nnu][Nbin];
double tref_at_tau_eq_1[Nnu][Nbin];

/*
Mean Binned results, tabulated with temperature and pressure:
kap_ro - Rosseland Opacity
kap_pl - Planck Opacity
kap_mean - Mean opacity of the Rosseland and Planck (weighted with optical depth)
B_band - Planck Function
dBdT - Derivative of Planck Function with temperature
*/

double kap_ro[Nbands][NT][Np];
double kap_pl[Nbands][NT][Np];
double kap_pp[Nbands][KMAX];

double abn_ro[Nbands][NT][Np];
double abn_pl[Nbands][NT][Np];
double abn_pp[Nbands][KMAX];

double sig_ro[Nbands][NT][Np];
double sig_pl[Nbands][NT][Np];
double sig_pp[Nbands][KMAX];


// outputs
float kap_ref_pp[Nbands][NT][Np];
float kap_mean[Nbands][NT][Np];
float abn_mean[Nbands][NT][Np];
float sig_mean[Nbands][NT][Np];

float * tau_pp_out;
float * sig_pp_out;
float * abn_pp_out;
float * kap_pp_out;
double * tau5000_pp;
float kap_full_ros[NT][Np];

float B_band[Nbands][NT];
float dBdT_band[Nbands][NT];
float kap_5000[NT][Np];
float B_5000[NT];

float nuout[Nnu];

// Create a structure to keep track of the sorted
// wavelengths, consisting of:
// nuix - the nu index
// binix - the bin index
// and *pNext - a pointer to the next wavelength in this band

struct member 
{
  int nuix;
  int binix;
  struct member *pNext;
};

// last wavelength point of each bin
struct member *pLast[Nbands];
// first wavelength point of each bin
struct member *pFirst[Nbands];

void input(void);
void initialize(void);
void sort(void);
void meanop(void);
void output(void);
void sort_full(void);
/* Kurucz Bins in Angstrom
  0  90.883297   1  93.499998   2  96.086501   3  97.666001   4  99.624547   5  102.04505   6  103.82255   7  105.57255   8  107.67725   9  110.42725
 10  114.00000  11  117.77810  12  121.27810  13  124.82825  14  127.06205  15  128.40900  16  130.47685  17  132.38780  18  133.92115  19  136.64765  
 20  139.83575  21  143.26220  22  147.23910  23  151.00000  24  155.24370  25  158.76100  26  162.01730  27  165.99999  28  170.25439  29  173.41730
 30  176.84590  31  180.18869  32  181.72760  33  186.11171  34  191.01090  35  193.87600  36  198.35424  37  201.77894  38  205.17959  39  210.49999
 40  216.23915  41  219.83830  42  223.00265  43  226.82235  44  229.91885  45  234.26105  46  239.76105  47  246.53720  48  252.36790  49  256.79265
 50  260.17774  51  263.71580  52  268.50000  53  273.50000  54  278.46315  55  283.75670  56  289.79355  57  296.29245  58  301.10545  59  307.81299
 60  317.99999  61  327.99999  62  337.99999  63  347.89653  64  357.10153  65  365.70499  66  375.07999  67  384.91209  68  394.83209  69  405.00000
 70  414.39960  71  421.89960  72  430.00000  73  441.16300  74  451.16300  75  460.00000  76  470.00000  77  480.00890  78  490.00890  79  499.62949
 80  506.51100  81  514.38150  82  530.00000  83  550.00000  84  567.03069  85  584.53069  86  605.00000  87  624.99999  88  645.18928  89  662.68928
 90  679.99999  91  699.99999  92  716.20056  93  731.20056  94  749.99999  95  769.99999  96  789.99999  97  809.99999  98  829.99999  99  849.99999 
100  869.99999 101  889.99999 102  905.87672 103  925.87672 104  959.99999 105  999.99999 106  1040.0000 107  1080.1528 108  1115.1528 109  1145.0000
110  1180.0000 111  1219.6464 112  1259.6464 113  1300.0000 114  1340.0000 115  1380.0000 116  1422.1696 117  1462.1696 118  1497.1742 119  1542.1742 
120  1595.7536 121  1647.7677 122  1692.0141 123  1730.0000 124  1775.0000 125  1820.0000 126  1860.0000 127  1905.0000 128  1952.3495 129  1997.3495
130  2045.6605 131  2085.6605 132  2125.0000 133  2175.0000 134  2225.0000 135  2275.0000 136  2325.0000 137  2375.0000 138  2425.0000 139  2480.6310 
140  2530.6310 141  2575.0000 142  2625.0000 143  2675.0000 144  2725.0000 145  2775.0000 146  2825.0000 147  2875.0000 148  2950.0000 149  3050.0000 
150  3150.0000 151  3250.0000 152  3350.0000 153  3450.0000 154  3550.0000 155  3623.5092 156  3673.5092 157  3750.0000 158  3850.0000 159  3950.0000 
160  4050.0000 161  4150.0000 162  4250.0000 163  4400.0000 164  4550.0000 165  4650.0000 166  4750.0000 167  4850.0000 168  4950.0000 169  5050.0000 
170  5150.0000 171  5250.0000 172  5350.0000 173  5450.0000 174  5550.0000 175  5650.0000 176  5750.0000 177  5850.0000 178  5950.0000 179  6050.0000 
180  6150.0000 181  6249.9999 182  6349.9999 183  6449.9999 184  6549.9999 185  6649.9999 186  6749.9999 187  6849.9999 188  6949.9999 189  7049.9999 
190  7149.9999 191  7249.9999 192  7349.9999 193  7449.9999 194  7549.9999 195  7649.9999 196  7749.9999 197  7849.9999 198  7949.9999 199  8049.9999 
200  8152.9134 201  8252.9134 202  8349.9999 203  8449.9999 204  8549.9999 205  8649.9999 206  8749.9999 207  8849.9999 208  8949.9999 209  9049.9999 
210  9149.9999 211  9249.9999 212  9349.9999 213  9449.9999 214  9549.9999 215  9649.9999 216  9749.9999 217  9849.9999 218  9949.9999 219  10125.000
220  10375.000 221  10625.000 222  10875.000 223  11125.000 224  11375.000 225  11625.000 226  11875.000 227  12125.000 228  12375.000 229  12625.000 
230  12875.000 231  13125.000 232  13375.000 233  13625.000 234  13875.000 235  14125.000 236  14419.083 237  14669.083 238  14875.000 239  15125.000 
240  15375.000 241  15625.000 242  15875.000 243  16200.000 244  16600.000 245  17000.000 246  17400.000 247  17800.000 248  18200.000 249  18600.000 
250  19000.000 251  19400.000 252  19800.000 253  20250.000 254  20750.000 255  21250.000 256  21750.000 257  22250.000 258  22647.016 259  22897.016 
260  23250.000 261  23750.000 262  24250.000 263  24750.000 264  25250.000 265  25750.000 266  26250.000 267  26750.000 268  27250.000 269  27750.000 
270  28250.000 271  28750.000 272  29250.000 273  29750.000 274  30250.000 275  30750.000 276  31250.000 277  31750.000 278  32411.716 279  33411.716 
280  34500.000 281  35500.000 282  36500.000 283  37500.000 284  38500.000 285  39500.000 286  40500.000 287  41500.000 288  42500.000 289  43500.000 
290  44500.000 291  45500.000 292  46500.000 293  47500.000 294  48500.000 295  49500.000 296  50500.000 297  51500.000 298  52500.000 299  53500.000
300  54500.000 301  55500.000 302  56500.000 303  57500.000 304  58500.000 305  59500.000 306  60500.000 307  61500.000 308  62499.999 309  63499.999 
310  64999.999 311  66999.999 312  68999.999 313  70999.999 314  72999.999 315  74999.999 316  76999.999 317  78999.999 318  80999.999 319  82999.999
320  84999.999 321  86999.999 322  88999.999 323  90999.999 324  92999.999 325  94999.999 326  96999.999 327  50049000. */
