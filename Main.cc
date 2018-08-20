////////////////////////////////////////////////////
////////////////////////////////////////////////////
// Sum Product Decoder with Puncturing            //
////////////////////////////////////////////////////
////////////////////////////////////////////////////

#include <armadillo>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <time.h>

#include "comms.h"
#include "decoders.h"

using namespace std;
using namespace arma;

const long maxWordsPerSNR = 100000000;     // The maximum number of words to run each simulation for
const long minErrsPerSNR = 100;          // The minimum number of errors to count for each simulation
const int  maxNumIterations = 100;    // Number of iterations for the Sum Product Algorithm

int main(int argc, char** argv)
{
    // Generate a random seed to ensure
    // all random numbers differ
    srand(genRandomSeed());

    int  n, k;
    umat H, G;

    // Get the user input binary file, and read H,G,n, and k from it
    string ldpcFileName;
    readParityFile(H, G, n, k, ldpcFileName);
    double codeRate = (double)k/(double)n;

    // Display the properties of the parity matrix
    cout << "\nThe code has the following properties >> " << endl;
    cout << "\nn = " << n << " k = " << k <<". The number of check nodes is " << H.n_rows << "." << endl;
    printCheckDegreeDistribution(H);

    // Get Decoding Algorithm
    bool SPA, GAA;
    getUserDecodingScheme(SPA, GAA);

    // Allow the user to select puncturing or no puncturing
    vector<int> puncturedBits;
    vector<double> ratesToPuncture;
    vector<int> numBitsToPuncture;
    bool randomPuncturing = false;
    bool puncture = getUserPuncturing(H,n,k, puncturedBits, ratesToPuncture, numBitsToPuncture, randomPuncturing);

    // Fetch the SNR values from the user for each rate that was indicated earlier
    vector< vector <double> > dBSNR(ratesToPuncture.size());
    getSNRValues(ratesToPuncture, dBSNR);

    //Create CSV file for output
    ofstream dataFile;
    const char* fileName = createOutputFileName(n,k,puncture);
    dataFile.open(fileName);

    // Display the start of the simulation
    time_t ctt = time(0);
    cout << "\n-- Bit Error Rate Simulation (" << n << ", " << k << ") LDPC Code --" << endl;
    cout << "Simulation started on " << asctime(localtime(&ctt)) << endl;

    /////////////////////////
	// Simulation Variables //
	//////////////////////////
	unsigned long int numWords;      // counter, running total of words per iteration (i.e. per SNR value)
	unsigned long int numBits;       // counter, running total of bits per iteration (numWords * k)
	unsigned long int numWordErrs;   // counter, running total of word errors per iteration
	unsigned long int numBitErrs;    // counter, running total of bit errors per iteration
	bool wordError;
	
	mat 	informationBits;       	// The information bits to be encoded (1 x k bit matrix)
	int*    codeWord = new int[n](); 	// The encoded bits (1 x n bit matrix)
	int*    modulatedCodeWord = new int[n](); // modulated code word (1 x n bit matrix)
	mat 	noiseVector; 			// The error vector (i.e. noise) (1 x n bit matrix)
	double* softReceivedVector = new double[n](); // The soft decision received (1 x n bit matrix)
	int*    estimatedCodeWord  = new int[n]();    		// The estimated codeword at the output of the decoder (1 x n bit matrix)
	int* 	estimatedInformationBits = new int[k]();	// The estimated information bits recovered from the estimated codeword (1 x k bit matrix)
	
	// Encoder/Decoder Objects
	BlockEncoder 		blockEncoder   		(&G);
	CodePuncturer 		codePuncturer 		(puncturedBits, n);
	SumProductDecoder 	sumProductDecoder 	(&H, maxNumIterations);
	GallagerADecoder 	gallagerADecoder 	(&H, maxNumIterations);

	/*WeightedBitFlipDecoder wbfDecoder   (&H, maxNumIterations);
	GDBFDecoder gdbfDecoder             (&H, maxNumIterations);*/
	BitFlipDecoder bitFlipDecoder       (&H, maxNumIterations);

	/////////////////////////////////////////////////////////////////////////////////////////
	// Begin the simulation, repeat for each desired code rate, and each desired SNR value //
	/////////////////////////////////////////////////////////////////////////////////////////
	for(int rate = 0; rate < (int)ratesToPuncture.size(); rate++)
	{
	    // Figure out the total number of SNR's to run for this instance of the simulation
	    int numSNR = dBSNR.at(rate).size();

	    // Find the code rate for the noise calculation
        codeRate = (double)k/(double)(n - numBitsToPuncture[rate]);
        cout << "\n-- Simulation for rate " << codeRate << " --" << endl;

        // Calculate the simulation noise power
        // SNR = Signal Power/Noise Power
        // It can be shown that (for white noise) SNR simplifies to Eb/N0
        // For coded modulation Eb = Es/R where R is the code rate of the code
        // Thus SNR simplifies to Eb/(RN0)
        // The power in white noise is also the variance of the noise, that is var = N0/2
        // We can assume that the energy per bit has been normalized, which results in
        // SNR = 1/(coderate * noisevariance), thus we can find the noise variance as follows
        mat linSNR = exp10(conv_to<mat>::from( dBSNR.at(rate) ) /10 );    // SNR values in linear units
        mat noiseVar = 1/(linSNR*codeRate*2);     // No/2 = sigma^2
        mat noiseStd = sqrt(noiseVar);

        mat wordErrRate = zeros(1, numSNR); // Store the word error rate for each iteration of the simulation
        mat bitErrRate = zeros(1, numSNR);  // Store the bit error rate for each iteration of the simulation


        for(int snr = 0; snr < numSNR; snr++)
        {
            printf("\n--SNR: %g dB -- \n", dBSNR.at(rate).at(snr));
            printf( " Number of word errors:        ");

            // Initialize all simulation variables, and set them all to zero
            numWords = 0;
            numBits = 0;
            numWordErrs = 0;
            numBitErrs = 0;

            // Run the simulation until the maximum number of words is reached or until the minimum number of errors is reached
            while(numWords < maxWordsPerSNR && numWordErrs < minErrsPerSNR)
            {
                // Initialize matrices for simulation
                informationBits.zeros(1, k); // information to be coded
                noiseVector.zeros(1, n); // error vector (from noise in channel)

                // Data Source
                informationBits.randu();
                informationBits = round(informationBits);

                /*cout << "\ninfobits: ";
                for(int i = 0; i < k; i++)
                	cout << informationBits[i] << " ";
                cout << endl;*/

                // Encoder (Encoding from the Generator Matrix)
                if(puncture)
                {
                	// need to resize, since puncturing removed some bits
                	delete [] codeWord;
                	codeWord = new int[n]();
                }
                blockEncoder.Encode( conv_to<umat>::from(informationBits), codeWord );

                /*cout << "\n" << std::setw(15) << "CodeWord: ";
				for(int i = 0; i < n; i++)
					cout << std::setw(3) << codeWord[i] * -2 + 1 << " ";
				cout << endl;*/


                // Puncture
                if( puncture )
                {
                	// Resize the vectors so that they are the size of the punctured code word
                	// This requires deleting their memory and reallocating.
                	delete [] softReceivedVector;
                	delete [] modulatedCodeWord;

                	softReceivedVector = new double[n-numBitsToPuncture[rate]];
                	modulatedCodeWord = new int[n-numBitsToPuncture[rate]];


                	if(randomPuncturing)
                		codePuncturer.ChangePuncturedBits(numBitsToPuncture[rate], n, k);
                	codePuncturer.Puncture(codeWord, numBitsToPuncture[rate]);
                }

                // BPSK Modulator, 0 -> 1, 1 -> -1
                for(int i = 0; i < n - numBitsToPuncture[rate]; i++)
                	modulatedCodeWord[i] = -2*codeWord[i] + 1;

                /*cout << "\nCodeWord:    ";
                for(int i = 0; i < n; i++)
                	cout << modulatedCodeWord[i] << " ";
                cout << endl;*/

                // Channel, AWGN
                noiseVector.randn(1, n-numBitsToPuncture[rate]);
                noiseVector = noiseVector*noiseStd(snr);

                for(int i = 0; i < n-numBitsToPuncture[rate]; i++)
                	softReceivedVector[i] = (double)modulatedCodeWord[i] + noiseVector[i];

                // Un-puncture
                if(puncture)
                	codePuncturer.UnPuncture(softReceivedVector, numBitsToPuncture[rate]);

                /*cout << "\nReceivedVector: ";
				for(int i = 0; i < n; i++)
					cout << softReceivedVector[i] << " ";
				cout << endl;*/

                // Decoder
                // Only use decoder if rate is not 1, since then we should have BPSK
                /*if(codeRate >= 1.0)
                {
            		for(int i = 0; i < n; i++)
            			estimatedCodeWord[i] = softReceivedVector[i] > 0 ? 1 : -1;
                }
                else*/
                {
                	if(SPA)
                		sumProductDecoder.Decode(softReceivedVector, noiseVar(snr), estimatedCodeWord);
                	else if(GAA && puncture)
                		gallagerADecoder.DecodeWithErasures(softReceivedVector, estimatedCodeWord, k);
                	else if(GAA)
                		gallagerADecoder.Decode(softReceivedVector, estimatedCodeWord);
                	/*else if(WBF)
                		wbfDecoder.Decode(softReceivedVector, estimatedCodeWord);
                	else if(GDBF)
                		gdbfDecoder.Decode(softReceivedVector, estimatedCodeWord);
                	else if(BF)
                		bitFlipDecoder.Decode(softReceivedVector, estimatedCodeWord);*/
                }
                // Demodulate
                for(int i = n-k; i < n; i++ )
                	estimatedInformationBits[i-n+k] = -1*estimatedCodeWord[i] > 0 ; // since the code is systematic

                wordError = false;
                for(int i = 0; i < k; i++)
                {
                	if(informationBits[i] != estimatedInformationBits[i])
                	{
                		numBitErrs++;
                		wordError = true;
                	}
                }

                if(wordError)
                	numWordErrs++;

                numWords++;
                numBits  += k;

                if(!(numBits % 5))
                {
                    printf("\b\b\b\b\b\b\b%7ld", numWordErrs);
                    fflush(0);
                }

            } // while(numWords < maxWordsPerSNR && numWordErrs < minErrsPerSNR), end of simulation for one SNR value

            wordErrRate(snr) = (double)numWordErrs/(numWords);
            bitErrRate(snr)  = (double)numBitErrs/(numBits);

            printf("\nWER: %g, BER: %g\n", wordErrRate(snr), bitErrRate(snr));
            dataFile << dBSNR.at(rate).at(snr) << "," << wordErrRate(snr) << "," << bitErrRate(snr) << "," << q(sqrt(2*pow(10, conv_to<mat>::from(dBSNR.at(rate))(snr)/10))) << "\n";

        } // for(snr < numSNR), end of whole simulation

        cout << "\n-- Results --" << endl;

        for(int snr = 0; snr < numSNR; snr++)
        {
            printf("SNR: %g dB, WER: %g, BER: %g\n", dBSNR.at(rate).at(snr), wordErrRate(snr), bitErrRate(snr));
        }
	}

	dataFile.close();

	ctt = time(0);

	cout << "Simulation completed on " << asctime(localtime(&ctt)) << endl;

	delete [] codeWord;
	delete [] modulatedCodeWord;
	delete [] softReceivedVector;
	delete [] estimatedCodeWord;
	delete [] estimatedInformationBits;

    return 0;
}
