#ifndef DECODERS_H
#define DECODERS_H

using namespace arma;
using namespace std;


class CodePuncturer
{
	public:
		CodePuncturer(vector<int> bits, int codeWordSize);
		~CodePuncturer();
		void Puncture(int* &codeWord, int numberOfBitsToPuncture);
		void UnPuncture(double* &puncturedword, int numOfBitsPunctured);
		void ChangePuncturedBits(int numBitsToPuncture, int n, int k);

	private:
		int sizeOfCodeWord;
		int totalNumberOfPuncturedBits;

		int * bitsThatCanBePunctured;

		void sortBits(int* bitsToSort, int numBitsToSort);
};

class BlockEncoder
{
	public:
		BlockEncoder(umat* Gptr); // constructor
		~BlockEncoder(); // destructor

		void Encode(umat u, int* &codeWord);

	private:
		int rows;
		int cols;

		int ** colSupport; // support of the column of the G matrix
		int *  colOnes;    // number of ones in each column of G matrix
};

class SumProductDecoder
{
	public:
		SumProductDecoder(umat* H, int maxNumIterations); // constructor
		~SumProductDecoder(); // destructor

		void Decode(double* receivedVectorPtr, double noiseVariance, int* &estimatedCodeWord);

	private:
		int totalIterations; // Total number of iterations for the SPA
		int numIterations;   // holds current number of iterations

		int rows; 			// Number of rows in the parity check matrix
		int cols; 			// Number of columns in the parity check matrix

		int** rowSupport;	// matrix, holds the indices of the support of each row of the parity check matrix, by row
		int** colSupport; 	// matrix, holds the indices of the support of each column of the parity check matrix, by column
		int*  rowOnes;		// array, holds the number of ones in each row of the parity check matrix
		int*  colOnes;		// array, holds the number of ones in each column of the parity check matrix

		double*  r;			// array, modified received vector for SPA
		double*  eps_tot;	// array, total extrinsic information provided to code bit by the rows with support for that codebit
		double** Y;			// matrix, channel information matrix
		double** Z;			// matrix, channel value matrix
		double** E;			// matrix, extrinsic matrix
		double** epsilon;	// matrix, extrinsic information provided to transmitted codebit by other code bits checked by h

		double* pSoftReceivedVector; // pointer, points to the soft decision received vector

		int*    hardReceivedVector; // array, hard decision vector at each iteration
		int*    syndrome;			// array, syndrome (vprime*Trans(H))

		bool calculateSyndrome();
};

class GallagerADecoder
{
	public:
		GallagerADecoder(umat* H, int maxNumIterations); // constructor
		~GallagerADecoder(); // destructor

		void Decode(double* receivedVectorPtr, int* &estimatedCodeWord);
		void DecodeWithErasures(double* receivedVectorPtr, int* &estimatedCodeWord, int k);

	private:
		int totalIterations;	// Total number of iteraitons for Gallager A algorithm
		int numIterations;	 	// holds the current number of iterations

		int rows; 				// Number of rows in the parity check matrix
		int cols;				// Number of columns in the parity check matrix

		int** rowSupport;		// matrix, holds the indices of the support of each row of the parity check matrix, by row
		int** colSupport; 		// matrix, holds the indices of the support of each column of the parity check amtrix, by column
		int*  rowOnes;			// array, holds the number of ones in each row of the parity check matrix
		int*  colOnes;			// array, holds the number of ones in each column of the parity check matrix

		int* initHardDecision;	// array, initial hard decision vector
		int* hardDecisionVector;	// array, hard decision vector at each iteration
		int* erasedBits;
		int* syndrome;			// array, syndrome (vprime*trans(H))

		int** varNodeVals;		// Messages sent from variable nodes to check nodes
		int** parityChecks;		// Messages sent from check nodes to variable nodes, parity checks involving all OTHER variable nodes in check equation

		int** peelingGuesses; 	// Matrix that holds the guesses from the peeling algorithm for each erased bit
		double*  peelingVector;

		int numErasedBits;
		int bitsPutBack;
		int positives;
		int negatives;
		int erasures;

		enum PeelingGuess
		{
			POSITIVE = 0,
			NEGATIVE = 1,
			ERASURE  = 2
		};

		bool calculateSyndrome();
};

class WeightedBitFlipDecoder
{
	public:
		WeightedBitFlipDecoder(umat* Hptr, int maxNumIterations); // constructor
		~WeightedBitFlipDecoder(); // destructor

		void Decode(double* receivedVectorPtr, int* &estimatedCodeWord);

	private:
		int totalIterations;

		int rows;
		int cols;

		int** rowSupport;
		int** colSupport;
		int* rowOnes;
		int* colOnes;

		double* syndromeReliability;
		double* bitReliability;

		double* pSoftReceivedVector;

		int* 	hardReceivedVector;
		int* 	syndrome;

		int index;

		bool calculateSyndrome();
};

class GDBFDecoder
{
	public:
		GDBFDecoder(umat* Hptr, int maxNumIterations); // constructor
		~GDBFDecoder(); // destructor

		void Decode(double* receivedVectorPtr, int* &estimatedCodeWord);

	private:
		int totalIterations;

		int index;

		int rows;
		int cols;

		int** rowSupport;
		int** colSupport;
		int* rowOnes;
		int* colOnes;

		int* decisionVector;
		int* syndromeComponents;
		double* inverseFunctions;

		double* pSoftReceivedVector;

		int* syndrome;
		int* hardDecisionVector;

		bool calculateSyndrome();
};

class BitFlipDecoder
{
	public:
		BitFlipDecoder(umat* Hptr, int maxNumIterations); // constructor
		~BitFlipDecoder(); // destructor

		void Decode(double* receivedVectorPtr, int* &estimatedCodeWord);

	private:
		int totalIterations;

		int rows;
		int cols;

		int** rowSupport;
		int** colSupport;
		int* rowOnes;
		int* colOnes;

		int* failedChecks;
		int  parityCheck;
		int  maxFailedChecks;

		double* pSoftReceivedVector;

		int* 	hardReceivedVector;
		int*    syndrome;

		int index;

		bool syndromeCheck;
};

double get_wall_time();

double get_cpu_time();

#endif
