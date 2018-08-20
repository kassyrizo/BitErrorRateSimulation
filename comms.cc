#include <armadillo>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <string>
#include <iostream>

#include "comms.h"
#include "Puncture.h"

long genRandomSeed()
{
  long randnum = 0;
  int fd = open ("/dev/urandom", O_RDONLY);
  if (fd != -1) {
    read( fd, &randnum, 4 );
    close(fd);
  }
  else{
  	printf("ERROR: Can't read /dev/urandom.\n");
  	exit(1);
  }

  if(randnum>0) randnum*=-1;

  if(!randnum){
  	printf("ERROR: Zero seed.\n");
  	exit(1);
  }

  return randnum;
}

double q(double arg)
{
	return 0.5 - 0.5*erf(arg/sqrt(2.0));
}

void readParityFile(umat &H, umat &G, int &n, int &k, string &fileName)
{
    // Fetch the name of the file that has the LDPC parameters
    cout << "Please enter the parity check matrix binary filename >> ";

    // check for valid filename
    ifstream inputFile;
    while(1)
    {
        cin >> fileName;
        inputFile.open(fileName.c_str(), ios::binary | ios::in);
        if(inputFile)
            break;
        cout << "You have entered an invalid filename. Please try again >> " << flush;
    }

    // Read the binary file data
    inputFile.read((char*)&n,sizeof(int));
    inputFile.read((char*)&k,sizeof(int));
    H.zeros(k,n);

    for(int i = 0; i < k; i++)
    {
        for(int j = 0; j < n; j++)
            inputFile.read((char*)&H(i,j),sizeof(char));
    }

    inputFile.read((char*)&k,sizeof(int));
    G.zeros(k,n);

    for(int i = 0; i < k; i++)
    {
        for(int j = 0; j < n; j++)
            inputFile.read((char*)&G(i,j),sizeof(char));
    }

    inputFile.read((char*)&k,sizeof(int));
    inputFile.close();

    return;
}

void getUserDecodingScheme(bool &SPA, bool &GAA)
{
    cout << "\nWhich decoding scheme would you like to use; SPA, or GAA? >> ";

    string decodingScheme;

    SPA = false;
    GAA = false;
    /*WBF = false;
    BF  = false;
    GDBF = false;*/

    while(1)
    {
        cin >> decodingScheme;

        if( decodingScheme.compare("SPA") == 0 )
            SPA = true;
        else if( decodingScheme.compare("GAA") == 0 )
        	GAA = true;
        /*else if( decodingScheme.compare("WBF") == 0 )
            WBF = true;
        else if( decodingScheme.compare("GDBF") == 0 )
        	GDBF = true;
        else if( decodingScheme.compare("BF") == 0 )
            BF = true;*/
        else
        {
            cout << "You have entered an invalid decoding scheme, please try again >> " << flush;
            continue;
        }

        return;
    }
}

bool getUserPuncturing(umat H, int n, int k, vector<int> &puncturedBits, vector<double> &ratesToPuncture, vector<int> &numBitsToPuncture, bool& randomPuncturing)
{
    cout << "\nWould you like to puncture this code? Y/y or N/n >> ";
    bool  puncture = checkForYesOrNo();

    if( !puncture )
    {
        numBitsToPuncture.push_back(0);
        ratesToPuncture.push_back( (double)k/(double)n );
        return false;
    }

    cout << "Would you like to use random puncturing? Y/y or N/n >> ";
    if(checkForYesOrNo())
    	randomPuncturing = true;
    else
    {
		// Ask the user for a file/option to save if there isn't one
		cout << "Would you like to load a puncturing file? Y/y or N/n >> ";
		bool loadFile = checkForYesOrNo();

		if( loadFile )
			readPunctureFile(puncturedBits);
		else
		{
			puncturedBits = PunctureCode(H);
			cout << "Would you like to save this puncturing pattern? Y/y or N/n >> ";

			if( checkForYesOrNo() )
				savePuncturingPattern(puncturedBits, n, k);
		}
    }

    // Ask the user for the code rates to puncture
    double coderate = (double)k/(double)n;

    cout << "\nWould you like to simulate the mothercode? Y/y or N/n >> ";

    if(checkForYesOrNo())
        ratesToPuncture.push_back(coderate);

    cout << "\nWhat rates would you like to simulate? Type D/d when finished >> ";
    readCodeRates(ratesToPuncture, coderate);

    cout << "\nThe mother code is rate " << coderate << " and you'd like to simulate rates ";
    for(vector<double>::iterator it = ratesToPuncture.begin(); it != ratesToPuncture.end(); ++it)
        cout << *it << " ";
    cout << endl;

    // Find the actual puncturing rates we can handle, and remove any duplicates
    for( int i = 0; i < (int) ratesToPuncture.size(); i++ )
    {
        int bitsToPuncture = n - floor( (double)k / ratesToPuncture.at(i) );

        if(bitsToPuncture > (int) puncturedBits.size() && !randomPuncturing)
            bitsToPuncture = puncturedBits.size();

        double codeRate = (double)k/(double)(n - bitsToPuncture);

        ratesToPuncture.at(i) = codeRate;
        numBitsToPuncture.push_back(bitsToPuncture);
    }

    // Sort and remove any duplicate rates
    sort( ratesToPuncture.begin(), ratesToPuncture.end() );
    ratesToPuncture.erase( unique( ratesToPuncture.begin(), ratesToPuncture.end() ), ratesToPuncture.end() );

    sort( numBitsToPuncture.begin(), numBitsToPuncture.end() );
    numBitsToPuncture.erase( unique( numBitsToPuncture.begin(), numBitsToPuncture.end() ), numBitsToPuncture.end() );

    cout << "The rates that this puncturing scheme can accommodate are ";
    for(vector<double>::iterator it = ratesToPuncture.begin(); it != ratesToPuncture.end(); ++it)
        cout << *it << " ";
    cout << endl;

    cout << "This corresponds to the following number of punctured bits: ";
    for(vector<int>::iterator it = numBitsToPuncture.begin(); it != numBitsToPuncture.end(); ++it)
    	cout << *it << " ";
    cout << endl;


    return puncture;
}

const char* createOutputFileName(int n, int k, bool puncture)
{
    string file;

    cout << "\nPlease enter the name of the file that you wish to save the output to >> ";

    cin >> file;
    const char* fileName;// = stream.str().c_str();

    fileName = file.c_str();

    return fileName;
}

void printCheckDegreeDistribution(umat H)
{
    umat checkDegreeDistribution = zeros<umat>(1,H.n_rows);

    for(int i = 0; i < (int)H.n_rows; i++)
    {
        uvec edges = find(H.row(i));
        checkDegreeDistribution.at(i) = edges.size();
    }
    vector<int> checkDegree = conv_to< vector<int> >::from(checkDegreeDistribution);

    sort(checkDegree.begin(), checkDegree.end());
    checkDegree.erase( unique( checkDegree.begin(), checkDegree.end() ), checkDegree.end() );

    cout << endl;
    for(int i = 0; i < (int)checkDegree.size(); i++)
    {
        uvec nodes = find(checkDegreeDistribution == checkDegree[i]);
        int numNodes = nodes.size();
        cout << "There are " << numNodes << " check nodes with " << checkDegree[i] << " edges." << endl;
    }
}

void getSNRValues(vector<double> ratesToPuncture, vector< vector<double> > & dBSNR)
{
    for(int i = 0; i < (int)ratesToPuncture.size(); i++)
    {
        cout << "\nPlease enter the SNRs you wish to simulate for rate " << ratesToPuncture.at(i) << ". When you are done, enter D/d. >> ";
        readSNRValues(dBSNR.at(i));
    }

    // sort the SNR's for each rate, remove any duplicates and echo back to user
    for(int i = 0; i < (int)dBSNR.size(); i++)
    {
        sort( dBSNR.at(i).begin(), dBSNR.at(i).end() );
        dBSNR.at(i).erase( unique( dBSNR.at(i).begin(), dBSNR.at(i).end() ), dBSNR.at(i).end() );
        cout << "\nSNR's to simulate ";
        for(vector<double>::iterator it = dBSNR[i].begin(); it != dBSNR[i].end(); ++it)
            cout << *it << " ";
        cout << endl;
    }
}

void checkForValidFilename(ifstream & inputFile)
{
    string fileName;

    while(1)
    {
        cin >> fileName;
        inputFile.open( fileName.c_str(), ios::binary | ios::in );
        if( inputFile )
            break;
        cout << "You have entered an invalid filename. Please try again >> " << flush;
    }

    return;
}

bool checkForYesOrNo()
{
    while(1)
    {
        string yesOrNo;
        cin >> yesOrNo;

        if( (yesOrNo.compare("Y") == 0) || (yesOrNo.compare("y") == 0) )
            return true;
        else if( (yesOrNo.compare("N") == 0) || (yesOrNo.compare("n") == 0) )
            return false;

        cout << "You have entered an invalid option. Please try again >> " << flush;
    }

    return false;
}

void readPunctureFile(vector<int> & puncturedBits)
{
    cout << "Please enter the puncturing pattern filename >> ";
    ifstream punctureFile;
    checkForValidFilename(punctureFile);

    // Read in the total number of punctured bits
    int numberPuncturedBits;
    punctureFile.read( (char*)&numberPuncturedBits, sizeof(int) );
    puncturedBits.resize(numberPuncturedBits);

    // Read in the punctured bits
    for( int i = 0; i < numberPuncturedBits; i++ )
        punctureFile.read( (char*)&puncturedBits[i], sizeof(int) );

    // Close the puncturing file
    punctureFile.close();
}

void savePuncturingPattern(vector<int> puncturedBits, int n, int k)
{
    stringstream outputFileName;
    outputFileName << k << "x" << n << "PuncturedBits.bin";

    // Open a binary file for output
    ofstream outputFile;
    outputFile.open( outputFileName.str().c_str(), ios::out | ios::binary );

    if( outputFile )
    {
        int size = puncturedBits.size();
        outputFile.write( (char*)&size, sizeof(int) );

        for( int i = 0; i < (int) puncturedBits.size(); i++ )
            outputFile.write( (char*)&puncturedBits[i], sizeof(int) );

        cout << "\nPuncturing pattern saved to file " << outputFileName.str().c_str() << "!";
    }
    else
        cout << "\nPuncturing pattern failed to write to file! ";

    return;
}

void readCodeRates(vector<double> & ratesToPuncture, double codeRate)
{
    while(1)
    {
        string input;
        cin >> input;

        if( (input.compare("D") == 0) || (input.compare("d") == 0) )
            break;

        double rate = atof(input.c_str());

        if(rate < codeRate || rate > 1.0)
        {
            cout << "\nYou have entered an invalid code rate. Please try again >> ";
            continue;
        }

        ratesToPuncture.push_back(rate);
    }

    return;
}

void readSNRValues(vector<double> & dBSNR)
{
    while(1)
    {
        string input;
        cin >> input;

        if( (input.compare("D") == 0) || (input.compare("d") == 0 ) )
            break;

        double snr = atof(input.c_str());

        // if atof fails, (i.e. no floating point number, it will return 0.0
        // check if the input was not 0.0 before putting out an error
        if( snr == 0.0 && !( (input.compare("0.0") == 0) || (input.compare("0") == 0) ) )
        {
            cout << "You have entered an invalid SNR value. Please try again >> ";
            continue;
        }

        dBSNR.push_back(snr);
    }
}

