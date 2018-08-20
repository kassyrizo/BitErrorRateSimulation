#ifndef COMMS_H
#define COMMS_H

using namespace arma;
using namespace std;


long genRandomSeed();

double q(double arg);

void readParityFile(umat &H, umat &G, int &n, int&k, string &fileName);

void getUserDecodingScheme(bool &SPA, bool &GAA);

bool getUserPuncturing(umat H, int n, int k, vector<int> & puncturedBits, vector<double> & ratesToPuncture, vector<int> & numBitsToPuncture, bool& randomPuncturing);

const char* createOutputFileName(int n, int k, bool puncture);

void printCheckDegreeDistribution(umat H);

void getSNRValues(vector<double> ratesToPuncture, vector< vector<double> > & dBSNR);

void checkForValidFilename();

bool checkForYesOrNo();

void readPunctureFile(vector<int> & puncturedBits);

void savePuncturingPattern(vector<int> puncturedBits, int n, int k);

void readCodeRates(vector<double> & ratesToPuncture, double codeRate);

void readSNRValues(vector<double> & dBSNR);

#endif
