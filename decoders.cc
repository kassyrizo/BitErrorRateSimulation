#include <armadillo>
#include <time.h>
#include <sys/time.h>
#include <stdexcept>
#include <stdlib.h>
#include <iomanip>
#include "decoders.h"


////////////////////////////////
// Puncturing Class Functions //
////////////////////////////////

CodePuncturer::CodePuncturer(vector<int> bits, int codeWordSize)
{
	sizeOfCodeWord = codeWordSize;
	totalNumberOfPuncturedBits = bits.size();
	bitsThatCanBePunctured = new int[totalNumberOfPuncturedBits]();

	for(int i = 0; i < totalNumberOfPuncturedBits; i++)
		bitsThatCanBePunctured[i] = bits[i];

	return;
}

CodePuncturer::~CodePuncturer()
{
	delete [] bitsThatCanBePunctured;
}

void CodePuncturer::Puncture(int* &codeWord, int numberOfBitsToPuncture)
{
	if(numberOfBitsToPuncture > totalNumberOfPuncturedBits)
		throw std::invalid_argument("cannot puncture that many bits!");

	int* newPuncturedCodeWord = new int[sizeOfCodeWord - numberOfBitsToPuncture]();
	int* bitsToPuncture = new int[numberOfBitsToPuncture]();

	for(int i = 0; i < numberOfBitsToPuncture; i++)
		bitsToPuncture[i] = bitsThatCanBePunctured[i];

	sortBits(bitsToPuncture, numberOfBitsToPuncture);

	int* puncturePtr = bitsToPuncture;
	int  bitsPunctured = 0;
	for(int i = sizeOfCodeWord - 1; i >= 0; i--)
	{
		if((*puncturePtr) != i)
			newPuncturedCodeWord[i - numberOfBitsToPuncture + bitsPunctured] = codeWord[i];
		else
		{
			bitsPunctured++;
			if(bitsPunctured < numberOfBitsToPuncture)
				puncturePtr++;
		}
	}

	delete [] bitsToPuncture;
	delete [] codeWord;
	codeWord = newPuncturedCodeWord;

	return;
}

void CodePuncturer::UnPuncture(double* &puncturedWord, int numOfBitsPunctured)
{
	if(numOfBitsPunctured > totalNumberOfPuncturedBits)
			throw std::invalid_argument("cannot unpuncture that many bits!");

	double* newUnpuncturedCodeWord = new double[sizeOfCodeWord]();
	int* bitsToUnpuncture = new int[numOfBitsPunctured]();

	for(int i = 0; i < numOfBitsPunctured; i++)
		bitsToUnpuncture[i] = bitsThatCanBePunctured[i];

	sortBits(bitsToUnpuncture, numOfBitsPunctured);

	int* puncturePtr = (&bitsToUnpuncture[numOfBitsPunctured - 1]);
	int  bitsPutBack = 0;
	for(int i = 0; i < sizeOfCodeWord; i++)
	{
		if( (*puncturePtr) == i)
		{
			newUnpuncturedCodeWord[i] = 0;
			bitsPutBack++;
			if(bitsPutBack < numOfBitsPunctured)
			puncturePtr--;
		}
		else
			newUnpuncturedCodeWord[i] = puncturedWord[i - bitsPutBack];
	}

	delete [] bitsToUnpuncture;
	delete [] puncturedWord;
	puncturedWord = newUnpuncturedCodeWord;

	return;
}

void CodePuncturer::ChangePuncturedBits(int numBitsToPuncture, int n, int k)
{
	// This is for random puncturing
	// Right now, it's configured for randomly puncturing ONLY parity bits
	totalNumberOfPuncturedBits = numBitsToPuncture;
	delete [] bitsThatCanBePunctured;
	bitsThatCanBePunctured = new int[numBitsToPuncture]();
	int number = 0;
	bool newBit = true;

	for(int i = 0; i < numBitsToPuncture; i++)
	{
		newBit = true;
		number = std::rand() % (n-k); // this is the line that makes it puncture only parity bits. change this so to (n) instead of (n-k) if all bits punctured

		//cout << "\n" << number << endl;

		for(int j = 0; j < i; j++)
			if(bitsThatCanBePunctured[j] == number)
				newBit = false;

		if(newBit)
		{
			bitsThatCanBePunctured[i] = number;
			newBit = false;
		}
		else
			--i;
	}

	/*cout << "Random punctured bits: ";
	for(int i = 0; i < numBitsToPuncture; i++)
		cout <<  bitsThatCanBePunctured[i] << " ";
	cout << endl;*/
 }

void CodePuncturer::sortBits(int* bitsToSort, int numBitsToSort)
{
	// sort the bits in ascending order
	int temp;
	for(int i = 0; i < numBitsToSort; i++)
	{
		for(int j = 0; j < numBitsToSort - 1; j++)
		{
			if(bitsToSort[j] < bitsToSort[j+1])
			{
				temp = bitsToSort[j+1];
				bitsToSort[j+1] = bitsToSort[j];
				bitsToSort[j] = temp;
			}
		}
	}

	return;
}

////////////////////////////////////
/// Block Encoder Class Functions //
////////////////////////////////////

BlockEncoder::BlockEncoder(umat* Gptr)
{
	rows = (*Gptr).n_rows;
	cols = (*Gptr).n_cols;
	colSupport = new int*[cols]();
	colOnes = new int[cols]();

	// Allocate the memory for the columns in G
	for(int i = 0; i < cols; i++)
	{
		uvec ones = find( (*Gptr).col(i) );
		colSupport[i] = new int[ones.size()]();
		colOnes[i] = ones.size();
		for(int j = 0; j < (int)ones.size(); j++)
			colSupport[i][j] = ones[j];
	}
}

BlockEncoder::~BlockEncoder()
{
	for(int i = 0; i < cols; i++)
		delete [] colSupport[i];

	delete [] colSupport;
	delete [] colOnes;
}

void BlockEncoder::Encode(umat u, int* &codeWord)
{
	// The encoded word is basically the xor (mod 2 addition) of
	// the bits in u with row support in each column in G
	for(int i = 0; i < cols - rows; i++)
	{
		codeWord[i] = 0; // reset the bit value

		for(int j = 0; j < colOnes[i]; j++)
			codeWord[i] ^= u[ colSupport[i][j] ];
	}

	// Concat u on the end, since this is a systematic code
	for(int i = cols - rows; i < cols; i++)
		codeWord[i] = u[ i - (cols - rows) ];

	return;
}

///////////////////////////////////////////
// Sum Product Algorithm Class Functions //
///////////////////////////////////////////

SumProductDecoder::SumProductDecoder(umat* Hptr, int maxNumIterations)
{
	// Max number of iterations for the sum product algorithm
	totalIterations = maxNumIterations;
	numIterations = 0;

	// Parity matrix dimensions
	rows = (*Hptr).n_rows;
	cols = (*Hptr).n_cols;

	// Create 2D dynamic array for rowSupport and colSupport
	// Dynamic memory allocation for row support and column support of H matrix
	rowSupport = new int*[rows](); // each pointer in rowSupport points to an array that holds the column index of the ones in the H matrix for that particular row
	colSupport = new int*[cols](); // each pointer in colSupport points to an array that holds the row index of the ones in the H matrix for that particular column

	rowOnes = new int[rows](); // hold the number of ones in each row
	colOnes = new int[cols](); // hold the number of ones in each column

	// Allocate the memory for the rows in H
	for(int i = 0; i < rows; i++)
	{
		uvec ones = find((*Hptr).row(i));
		rowSupport[i] = new int[ones.size()]();
		rowOnes[i] = ones.size();

		for(int j = 0; j < (int)ones.size(); j++)
			rowSupport[i][j] = ones[j];
	}

	// Allocate the memory for the columns in  H
	for(int i = 0; i < cols; i++)
	{
		uvec ones = find((*Hptr).col(i));
		colSupport[i] = new  int[ones.size()]();
		colOnes[i] = ones.size();

		for(int j = 0; j < (int)ones.size(); j++)
			colSupport[i][j] = ones[j];
	}

	hardReceivedVector = new int[cols]();
	syndrome = new int[rows]();

	r = new double[cols]();
	eps_tot = new double[cols]();
	Y = new double*[rows]();
	Z = new double*[rows]();
	E = new double*[rows]();
	epsilon = new double*[rows]();

	//Memory allocation and initialization for SPA variables
	for(int i = 0; i < cols; i++)
	{
		r[i]   = 0;
		eps_tot[i] = 0;
	}

	for(int i = 0; i < rows; i++)
	{
		Y[i] = new double[rowOnes[i]]();
		Z[i] = new double[rowOnes[i]]();
		E[i] = new double[rowOnes[i]]();
		epsilon[i] = new double[rowOnes[i]]();

		for(int j = 0; j < rowOnes[i]; j++)
		{
			Y[i][j] = 0;
			Z[i][j] = 0;
			E[i][j] = 0;
			epsilon[i][j] = 0;
		}
	}
}

SumProductDecoder::~SumProductDecoder()
{
	// Free up the memory from dynamic allocations
	for(int i = 0; i < rows; i++)
		delete [] rowSupport[i];
	delete [] rowSupport;

	for(int i = 0; i < cols; i++)
		delete [] colSupport[i];
	delete [] colSupport;

	delete [] rowOnes;
	delete [] colOnes;

	// Delete memory allocations
	for(int i = 0; i < rows; i++)
	{
		delete [] Y[i];
		delete [] Z[i];
		delete [] E[i];
		delete [] epsilon[i];
	}

	delete [] hardReceivedVector;
	delete [] syndrome;

	delete [] Y;
	delete [] Z;
	delete [] E;
	delete [] epsilon;
	delete [] r;
	delete [] eps_tot;
}

void SumProductDecoder::Decode(double* receivedVector, double noiseVariance, int* &estimatedCodeWord)
{
	pSoftReceivedVector = receivedVector; // set the soft Received vector to point at the received vector

	/////////////////////////////////////
	// Begin the Sum Product Algorithm //
	/////////////////////////////////////

	//////////////////////////////////////////
	// Initialize the Sum Product Algorithm //
	//////////////////////////////////////////
	numIterations = 0;

	for(int i = 0; i < cols; i++)
	{
		r[i] = 2 * pSoftReceivedVector[i] / noiseVariance;
		eps_tot[i] = 0;
		estimatedCodeWord[i]  = r[i] > 0 ? 1 : -1;
	}

	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < rowOnes[i]; j++)
		{
			Y[i][j] = r[rowSupport[i][j]];
			Z[i][j] = Y[i][j];
			E[i][j] = 0;
			epsilon[i][j] = 0;
		}
	}

	while(numIterations < totalIterations)
	{
		// Check the syndrome to see if it is all 0's, if it is, don't do the SPA
		for(int i = 0; i < cols; i++)
			hardReceivedVector[i] = -1*estimatedCodeWord[i] > 0;

		if(calculateSyndrome())
			break;

		// For each row, calculate the extrinsic information provided by the OTHER bits with support in the row, not including the bit itself
		double eps;
		for(int i = 0; i < rows; i++)
		{
			for(int j = 0; j < rowOnes[i]; j++)
			{
				eps = 1;
				for(int l = 0; l < rowOnes[i]; l++)
				{
					if(rowSupport[i][j] != rowSupport[i][l])
						eps *= tanh( Z[i][l]/2 );
				}

				eps = log( (1+eps)/(1-eps) );
				epsilon[i][j] = eps;
			}
		}

		// Calculate the Extrinsic Matrix
		double sum = 0;
		for(int i = 0; i < cols; i++)
		{
			for(int j = 0; j < colOnes[i]; j++)
			{
				sum = 0;
				for(int l = 0; l < colOnes[i]; l++)
				{
					if(colSupport[i][j] != colSupport[i][l])
					{
						for(int t = 0; t < rowOnes[ colSupport[i][l] ]; t++)
							if (rowSupport[ colSupport[i][l] ][t] == i)
							{
								sum += epsilon[ colSupport[i][l] ][t];
								break;
							}
					}
				}

				for(int t = 0; t < rowOnes[ colSupport[i][j] ]; t++)
					if(rowSupport[ colSupport[i][j]][t] == i)
					{
						E[ colSupport[i][j] ][t] = sum;
						break;
					}
			}
		}

		// Calculate epsilon_total
		for(int i = 0; i < cols; i++)
		{
			for(int t = 0; t < rowOnes[ colSupport[i][0] ]; t++)
				if(rowSupport[ colSupport[i][0] ][t] == i)
				{
					eps_tot[i] = E[ colSupport[i][0] ][t] + epsilon[ colSupport[i][0] ][t];
					break;
				}
		}

		////////////////////////////////////////////////
		// Step 2: Form the next Z, r, and z matrices //
		////////////////////////////////////////////////

		for(int i = 0; i < rows; i++)
			for(int j = 0; j < rowOnes[i]; j++)
				Z[i][j] = Y[i][j]+E[i][j];

		for(int i = 0; i < cols; i++)
		{
			r[i] += eps_tot[i];
			estimatedCodeWord[i] = r[i] > 0 ? 1 : -1;
		}

		numIterations++;
	}

	return;
}

bool SumProductDecoder::calculateSyndrome()
{
	bool syndromeCheck = true;

	// The syndrome is basically the xor (mod 2 addition) of
	// the bits codeword with support in each row in H
	for(int i = 0; i < rows; i++)
	{
		syndrome[i] = 0; // reset the bit value

		for(int j = 0; j < rowOnes[i]; j++)
			syndrome[i] ^= hardReceivedVector[ rowSupport[i][j] ];

		if(syndrome[i] == 1)
			syndromeCheck = false;
	}

	return syndromeCheck;
}

////////////////////////////////////////////
// Gallager's Algorithm A Class Functions //
////////////////////////////////////////////

GallagerADecoder::GallagerADecoder(umat* Hptr, int maxNumIterations)
{
	// Max number of iterations for the sum product algorithm
	totalIterations = maxNumIterations;
	numIterations = 0;

	// Parity matrix dimensions
	rows = (*Hptr).n_rows;
	cols = (*Hptr).n_cols;

	// Create 2D dynamic array for rowSupport and colSupport
	// Dynamic memory allocation for row support and column support of H matrix
	rowSupport = new int*[rows](); // each pointer in rowSupport points to an array that holds the column index of the ones in the H matrix for that particular row
	colSupport = new int*[cols](); // each pointer in colSupport points to an array that holds the row index of the ones in the H matrix for that particular column

	rowOnes = new int[rows](); // hold the number of ones in each row
	colOnes = new int[cols](); // hold the number of ones in each column

	// Allocate the memory for the rows in H
	for(int i = 0; i < rows; i++)
	{
		uvec ones = find((*Hptr).row(i));
		rowSupport[i] = new int[ones.size()]();
		rowOnes[i] = ones.size();

		for(int j = 0; j < (int)ones.size(); j++)
			rowSupport[i][j] = ones[j];
	}

	// Allocate the memory for the columns in  H
	for(int i = 0; i < cols; i++)
	{
		uvec ones = find((*Hptr).col(i));
		colSupport[i] = new  int[ones.size()]();
		colOnes[i] = ones.size();

		for(int j = 0; j < (int)ones.size(); j++)
			colSupport[i][j] = ones[j];
	}

	initHardDecision = new int[cols]();
	hardDecisionVector = new int[cols]();
	erasedBits = new int[cols]();
	syndrome = new int[rows]();

	varNodeVals  = new int*[rows]();
	parityChecks = new int*[rows]();

	for(int i = 0; i < rows; i++)
	{
		varNodeVals[i]  = new int[rowOnes[i]]();
		parityChecks[i] = new int[rowOnes[i]]();
	}

	peelingVector = new double[cols]();

}

GallagerADecoder::~GallagerADecoder()
{
	// Free up the memory from dynamic allocations
	for(int i = 0; i < rows; i++)
		delete [] rowSupport[i];
	delete [] rowSupport;

	for(int i = 0; i < cols; i++)
		delete [] colSupport[i];

	delete [] colSupport;

	delete [] rowOnes;
	delete [] colOnes;

	for(int i = 0; i < rows; i++)
	{
		delete [] parityChecks[i];
		delete [] varNodeVals[i];
	}

	delete [] initHardDecision;
	delete [] hardDecisionVector;
	delete [] erasedBits;
	delete [] syndrome;
	delete [] parityChecks;
	delete [] varNodeVals;
	delete [] peelingVector;
}

void GallagerADecoder::Decode(double* receivedVectorPtr, int* &estimatedCodeWord)
{
	for(int i = 0; i < cols; i++)
	{
		initHardDecision[i] = receivedVectorPtr[i] < 0;
		hardDecisionVector[i] = receivedVectorPtr[i] < 0;
	}

	/*cout << "hrdDecision: ";
	for(int i = 0; i < cols; i++)
		cout << hardDecisionVector[i] << " ";
	cout << endl;*/

	// Initialize the algorithm using the hard decision guess
	for(int i = 0; i < rows; i++)
		for(int j = 0; j < rowOnes[i]; j++)
			varNodeVals[i][j] = initHardDecision[rowSupport[i][j]];

	/*cout << "varNodeVals: ";
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < rowOnes[i]; j++)
			cout << varNodeVals[i][j] << " ";
		cout << "\n             ";
	}
	cout << endl;*/

	for(numIterations = 0; numIterations < totalIterations; numIterations++)
	{
		if(calculateSyndrome())
			break;

		// calculate the parity check equations (i.e. extrinsic information)
		for(int i = 0; i < rows; i++)
		{
			for(int j = 0; j < rowOnes[i]; j++)
			{
				parityChecks[i][j] = 0;
				for(int l = 0; l < rowOnes[i]; l++)
				{
					if(rowSupport[i][j] != rowSupport[i][l])
						parityChecks[i][j] ^= varNodeVals[i][l];
				}
			}
		}

		/*cout << "parityChecks: ";
		for(int i = 0; i < rows; i++)
		{
			for(int j = 0; j < rowOnes[i]; j++)
				cout << parityChecks[i][j] << " ";
			cout << "\n              ";
		}
		cout << endl;*/

		// Check the check node information for each variable node
		int value = 0;

		for(int i = 0; i < cols; i++) // Iterate through the columns (i.e. variable nodes)
		{
			for(int j = 0; j < colOnes[i]; j++) // Iterate through the one values in the columns
			{
				value = 100; // set to an initial value, as an indication that this hasn't been set yet.

				for(int l = 0; l < colOnes[i]; l++)
				{
					if(j != l) // Make sure we're only sending information from all OTHER check nodes
					{
						int columnIndex = 0;
						// populate a new vector with all the values, then check if all values are equal
						for(int t = 0; t < rowOnes[ colSupport[i][l] ]; t++) // initially set value to paritycheck[first 1 in column][column we're looking at]
						{
							if (rowSupport[ colSupport[i][l] ][t] == i)
							{
								columnIndex = t;
								break;
							}
						}

						if(value == 100)
							value = parityChecks[ colSupport[i][l] ][columnIndex]; // set value if it hasn't been set before

						else if(value != parityChecks[ colSupport[i][l] ][columnIndex])
						{
							value = initHardDecision[i];
							break;
						}
					}
				}

				// Update the variable node value
				for(int t = 0; t < rowOnes[ colSupport[i][j] ]; t++)
					if(rowSupport[ colSupport[i][j]][t] == i)
					{
						varNodeVals[ colSupport[i][j] ][t] = value;
						break;
					}
			}
		}

		/*cout << "varNodeVals: ";
		for(int i = 0; i < rows; i++)
		{
			for(int j = 0; j < rowOnes[i]; j++)
				cout << varNodeVals[i][j] << " ";
			cout << "\n             ";
		}
		cout << endl;*/

		// hard decision here

		int numOnes = 0;
		int numZeros = 0;

		for(int i = 0; i < cols; i++)
		{
			numOnes = 0;
			numZeros = 0;

			for(int j = 0; j < colOnes[i]; j++)
			{
				for(int t = 0; t < rowOnes[ colSupport[i][j] ]; t++)
				{
					if( (rowSupport[ colSupport[i][j] ][t]) == i )
					{	if(varNodeVals[ colSupport[i][j] ][t] == 0)
							numZeros++;
						else
							numOnes++;

						break;
					}
				}
			}

			if( (colOnes[i] % 2) == 0 ) // if the degree of the variable node is even, also take into account the initial channel information
			{
				if(initHardDecision[i] == 0)
					numZeros++;
				else
					numOnes++;
			}

			hardDecisionVector[i] = numZeros > numOnes? 0 : 1;
		}

		/*cout << std::setw(15) << "hardDecision: ";
		for(int i = 0; i < cols; i++)
			cout << std::setw(3) << hardDecisionVector[i] << " ";
		cout << endl;*/
	}

	for(int i = 0; i < cols; i++)
		estimatedCodeWord[i] = -2*hardDecisionVector[i]+1;

	return;
}

void GallagerADecoder::DecodeWithErasures(double* receivedVectorPtr, int* &estimatedCodeWord, int k)
{
	// In the erasure case, the message alphabet is {1, 0, -1}
	// First, use the Peeling Decoder to fill in the missing (erased) bits
	numErasedBits = 0;
	bool erased = true;

	// Fill in the hard decision vector with erasures
	for(int i = 0; i < cols; i++)
	{
		if(receivedVectorPtr[i] > 0)
			peelingVector[i] = 1.0;
		else if(receivedVectorPtr[i] < 0)
			peelingVector[i] = -1.0;
		else if(receivedVectorPtr[i] == 0) // check if the bit is erased
			peelingVector[i] = 0.0;
	}

	while(erased)
	{
		// Fill in the erased bits vector first (will only be in parity bits, this code is in systematic form)
		numErasedBits = 0;

		for(int i = 0; i < cols - k; i++)
		{
			if(peelingVector[i] == 0.0) // check if the bit is erased
			{
				numErasedBits++;
				int* temp = new int[numErasedBits];

				for(int j = 0; j < numErasedBits-1; j++)
					temp[j] = erasedBits[j];
				temp[numErasedBits-1] = i;

				delete [] erasedBits;
				erasedBits = temp;
			}
		}

		/*cout << std::setw(15) << "hardDecision: ";
		for(int i = 0; i < cols; i++)
			cout << std::setw(3) << peelingVector[i] << " ";
		cout << endl;

		cout << "Erased Bits: ";
		for(int i = 0; i < numErasedBits; i++)
			cout << erasedBits[i] << " ";
		cout << endl;*/

		// Create a matrix to store the bit guesses
		peelingGuesses = new int*[numErasedBits]();

		for(int i = 0; i < numErasedBits; i++)
			peelingGuesses[i] = new int[ 3 ](); // change this from hard coding, basically increase the positive, negative or erasure counters

		for(int i = 0; i < numErasedBits; i++) // for each erased bit
		{
			for(int j = 0; j < colOnes[ erasedBits[i] ]; j++) // find the index of the parity check equations for each erased bit
			{
				int product = 1;

				for(int k = 0; k < rowOnes[ colSupport[ erasedBits[i] ][j] ]; k++) // find the ones in the parity check equation
				{
					if(rowSupport[ colSupport[ erasedBits[i] ][j] ][k] != erasedBits[i])
						product *= peelingVector[ rowSupport[ colSupport[ erasedBits[i] ][j] ][k] ];

					if(product == 0)
						break;
				}

				if(product > 0)
					peelingGuesses[i][POSITIVE]++;
				else if(product < 0)
					peelingGuesses[i][NEGATIVE]++;
				else
					peelingGuesses[i][ERASURE]++;
			}

		}

		/*cout << "peeling: \n";
		for(int i = 0; i < numErasedBits; i++)
		{
			for(int j = 0; j < colOnes[ erasedBits[i] ]; j++)
				cout << peelingGuesses[i][j] << " ";
			cout << endl;
		}
		cout << endl;*/

		// Make the guesses here based on the peeling matrix
		bitsPutBack = 0;

		for(int i = 0; i < numErasedBits; i++)
		{
			positives = peelingGuesses[i][POSITIVE];
			negatives = peelingGuesses[i][NEGATIVE];
			erasures  = peelingGuesses[i][ERASURE];

			// This decision criteria can probably be tweaked and played with
			// If the number of positives is greater than the number

			if(positives > negatives/*&& !negatives*/)
			{
				peelingVector[ erasedBits[i] ] = 1.0;
				bitsPutBack++;
			}
			else if(negatives > positives/*&& !positives*/)
			{
				peelingVector[ erasedBits[i] ] = -1.0;
				bitsPutBack++;
			}
		}

		// If all of the bits have been erased, stop the loop, or if nothing was updated, end the loop
		if( (numErasedBits - bitsPutBack) == 0 || (bitsPutBack == 0) )
			erased = false;

		for(int i = 0; i < numErasedBits; i++)
			delete [] peelingGuesses[i];
		delete [] peelingGuesses;
	}

	/*cout << std::setw(15) << "hardDecision: ";
	for(int i = 0; i < cols; i++)
		cout << std::setw(3) << peelingVector[i] << " ";
	cout << endl;*/

	// Set any bits that have not been peeled to zero
	if( bitsPutBack == 0)
	{
		// Since no bits were put back, the erasedBits haven't changed. Fill them with 1's
		for(int i = 0; i < numErasedBits; i++)
			peelingVector[ erasedBits[i] ] = 1.0;

		/*cout << std::setw(15) << "hardDecision: ";
		for(int i = 0; i < cols; i++)
			cout << std::setw(3) << peelingVector[i] << " ";
		cout << endl;*/
	}


	Decode(peelingVector, estimatedCodeWord);

	return;
}

bool GallagerADecoder::calculateSyndrome()
{
	bool syndromeCheck = true;

	// The syndrome is basically the xor (mod 2 addition) of
	// the bits codeword with support in each row in H
	for(int i = 0; i < rows; i++)
	{
		syndrome[i] = 0; // reset the bit value

		for(int j = 0; j < rowOnes[i]; j++)
			syndrome[i] ^= hardDecisionVector[ rowSupport[i][j] ];

		if(syndrome[i] == 1)
			syndromeCheck = false;
	}

	return syndromeCheck;
}

/////////////////////////////////////////////////
// Weighted Bit Flip Algorithm Class Functions //
/////////////////////////////////////////////////

WeightedBitFlipDecoder::WeightedBitFlipDecoder(umat* Hptr, int maxNumIterations)
{
	// Max number of iterations for the sum product algorithm
	totalIterations = maxNumIterations;

	// Parity matrix dimensions
	rows = (*Hptr).n_rows;
	cols = (*Hptr).n_cols;

	// Create 2D dynamic array for rowSupport and colSupport
	// Dynamic memory allocation for row support and column support of H matrix
	rowSupport = new int*[rows](); // each pointer in rowSupport points to an array that holds the column index of the ones in the H matrix for that particular row
	colSupport = new int*[cols](); // each pointer in colSupport points to an array that holds the row index of the ones in the H matrix for that particular column

	rowOnes = new int[rows](); // hold the number of ones in each row
	colOnes = new int[cols](); // hold the number of ones in each column

	// Allocate the memory for the rows in H
	for(int i = 0; i < rows; i++)
	{
		uvec ones = find((*Hptr).row(i));
		rowSupport[i] = new int[ones.size()]();
		rowOnes[i] = ones.size();

		for(int j = 0; j < (int)ones.size(); j++)
			rowSupport[i][j] = ones[j];
	}

	// Allocate the memory for the columns in  H
	for(int i = 0; i < cols; i++)
	{
		uvec ones = find((*Hptr).col(i));
		colSupport[i] = new int[ones.size()]();
		colOnes[i] = ones.size();

		for(int j = 0; j < (int)ones.size(); j++)
			colSupport[i][j] = ones[j];
	}

	syndromeReliability = new double[rows]();
	bitReliability = new double[cols]();


	hardReceivedVector = new int[cols]();
	syndrome = new int[rows]();
}

WeightedBitFlipDecoder::~WeightedBitFlipDecoder()
{
	// Free up the memory from dynamic allocations
	for(int i = 0; i < rows; i++)
		delete [] rowSupport[i];
	delete [] rowSupport;

	for(int i = 0; i < cols; i++)
		delete [] colSupport[i];
	delete [] colSupport;

	delete [] rowOnes;
	delete [] colOnes;

	delete [] syndromeReliability;
	delete [] bitReliability;

	delete [] hardReceivedVector;
	delete [] syndrome;
}

void WeightedBitFlipDecoder::Decode(double* receivedVectorPtr, int* &estimatedCodeWord )
{
	pSoftReceivedVector = receivedVectorPtr;

	for(int i = 0; i < cols; i++)
		hardReceivedVector[i] = pSoftReceivedVector[i] < 0;

    // Initialize the algorithm
    for(int i = 0; i < rows; i++)
    	syndromeReliability[i] = 0;

    for(int i = 0; i < cols; i++)
    	bitReliability[i] = 0;

    /* Calculate the reliability measures of the syndrome components
     *
     * First, find the support of each row (position of any ones in the row).
     * Then, find the value of the received vector with support in the row that
     * has the minimum magnitude
     */
    for(int row = 0; row < rows; row++)
    {
    	for(int col = 0; col < rowOnes[row]; col++)
    	{
    		if( col == 0 || syndromeReliability[row] > abs( pSoftReceivedVector[rowSupport[row][col]] ))
    			syndromeReliability[row] = abs(pSoftReceivedVector[rowSupport[row][col]]);
    	}
    }

    for(int iteration = 0; iteration < totalIterations; iteration++)
    {
        if(calculateSyndrome())
            break;

        /* Calculate the reliability measure for each bit
         *
         * First, find the support of each column (position of any ones in the row)
         * Then, calculate the reliability component from each parity equation on
         * each bit. Finally, sum all reliability components together.
         */
        for(int i = 0; i < cols; i++)
        	bitReliability[i] = 0;

        for(int col = 0; col < cols; col++)
            for(int row = 0; row < colOnes[col]; row++)
                bitReliability[col] += (2.0 * (double)(syndrome[colSupport[col][row]]) - 1.0) * syndromeReliability[colSupport[col][row]];

        // Locate the value where the bit reliability measure is the largest
        index = 0;

        for(int i = 0; i < cols; i++)
        {
        	if(i == 0 || bitReliability[i] > bitReliability[index])
        		index = i;
        }

        // Toggle that bit
        hardReceivedVector[index] = 1 - hardReceivedVector[index];
    }

    for(int i = 0; i < cols; i++)
    	estimatedCodeWord[i] = -2*hardReceivedVector[i]+1;
}

bool WeightedBitFlipDecoder::calculateSyndrome()
{
	bool syndromeCheck = true;

	// The syndrome is basically the xor (mod 2 addition) of
	// the bits codeword with support in each row in H
	for(int i = 0; i < rows; i++)
	{
		syndrome[i] = 0; // reset the bit value

		for(int j = 0; j < rowOnes[i]; j++)
			syndrome[i] ^= hardReceivedVector[ rowSupport[i][j] ];

		if(syndrome[i] == 1)
			syndromeCheck = false;
	}

	return syndromeCheck;
}

/////////////////////////////////////////////////////////////
// Gradient Descent Bit Flipping Algorithm Class Functions //
/////////////////////////////////////////////////////////////

GDBFDecoder::GDBFDecoder(umat* Hptr, int maxNumIterations)
{
	// Max number of iterations for the sum product algorithm
	totalIterations = maxNumIterations;

	// Parity matrix dimensions
	rows = (*Hptr).n_rows;
	cols = (*Hptr).n_cols;

	// Create 2D dynamic array for rowSupport and colSupport
	// Dynamic memory allocation for row support and column support of H matrix
	rowSupport = new int*[rows](); // each pointer in rowSupport points to an array that holds the column index of the ones in the H matrix for that particular row
	colSupport = new int*[cols](); // each pointer in colSupport points to an array that holds the row index of the ones in the H matrix for that particular column

	rowOnes = new int[rows](); // hold the number of ones in each row
	colOnes = new int[cols](); // hold the number of ones in each column

	// Allocate the memory for the rows in H
	for(int i = 0; i < rows; i++)
	{
		uvec ones = find((*Hptr).row(i));
		rowSupport[i] = new int[ones.size()]();
		rowOnes[i] = ones.size();

		for(int j = 0; j < (int)ones.size(); j++)
			rowSupport[i][j] = ones[j];
	}

	// Allocate the memory for the columns in  H
	for(int i = 0; i < cols; i++)
	{
		uvec ones = find((*Hptr).col(i));
		colSupport[i] = new int[ones.size()]();
		colOnes[i] = ones.size();

		for(int j = 0; j < (int)ones.size(); j++)
			colSupport[i][j] = ones[j];
	}

	decisionVector = new int[cols]();
	syndromeComponents = new int[rows]();
	inverseFunctions = new double[cols]();

	hardDecisionVector = new int[cols]();
	syndrome = new int[rows]();
}

GDBFDecoder::~GDBFDecoder()
{
	// Free up the memory from dynamic allocations
	for(int i = 0; i < rows; i++)
		delete [] rowSupport[i];
	delete [] rowSupport;

	for(int i = 0; i < cols; i++)
		delete [] colSupport[i];
	delete [] colSupport;

	delete [] rowOnes;
	delete [] colOnes;

	delete [] decisionVector;
	delete [] syndromeComponents;
	delete [] inverseFunctions;

	delete [] hardDecisionVector;
	delete [] syndrome;
}

void GDBFDecoder::Decode(double* receivedVectorPtr, int* &estimatedCodeWord )
{
	pSoftReceivedVector = receivedVectorPtr;

	for(int i = 0; i < cols; i++)
		hardDecisionVector[i] = pSoftReceivedVector[i] < 0;

	/*cout << "hardDecision: ";
	for(int i = 0; i < cols; i++)
		cout << hardDecisionVector[i] << " ";
	cout << endl;*/

    // Initialize the algorithm
    for(int i = 0; i < cols; i++)
    	decisionVector[i] = pSoftReceivedVector[i] > 0 ? 1 : -1;

    /*cout << "decisionVector: ";
	for(int i = 0; i < cols; i++)
		cout << decisionVector[i] << " ";
	cout << endl;*/

    for(int iteration = 0; iteration < totalIterations; iteration++)
    {
        if(calculateSyndrome())
            break;

        // compute the syndrome components
        for(int row = 0; row < rows; row++)
        {
        	syndromeComponents[row] = 1;
        	for(int col = 0; col < rowOnes[row]; col++)
        		syndromeComponents[row] *= decisionVector[rowSupport[row][col]];
        }

        /*cout << "syndromeComponents: ";
		for(int i = 0; i < rows; i++)
			cout << syndromeComponents[i] << " ";
		cout << endl;*/

        // compute the inverse functions
        for(int col = 0; col < cols; col++)
        {
        	inverseFunctions[col] = pSoftReceivedVector[col]*decisionVector[col];

        	for(int row = 0; row < colOnes[col]; row++)
        		inverseFunctions[col] += syndromeComponents[colSupport[col][row]];
        }

        /*cout << "inverseFunctions: ";
		for(int i = 0; i < cols; i++)
			cout << inverseFunctions[i] << " ";
		cout << endl;*/

        // single bit flip version
        // flip the bit for which the inverse function is the smallest
        // Locate the value where the inverse functions measure is the smallest
        index = 0;

        for(int i = 0; i < cols; i++)
        {
        	if(i == 0 || inverseFunctions[i] < inverseFunctions[index])
        		index = i;
        }

        //cout << "index: " << index << endl;

        // Toggle that bit
        decisionVector[index] *= -1;

        /*cout << "decisionVector: ";
		for(int i = 0; i < cols; i++)
			cout << decisionVector[i] << " ";
		cout << endl;*/

        for(int i = 0; i < cols; i++)
        	hardDecisionVector[i] = decisionVector[i] < 0;

        /*cout << "hardDecisionVector: ";
		for(int i = 0; i < cols; i++)
			cout << hardDecisionVector[i] << " ";
		cout << endl;*/
    }

    for(int i = 0; i < cols; i++)
    	estimatedCodeWord[i] = -2*hardDecisionVector[i]+1;

    /*cout << "estimatedCodeWord: ";
	for(int i = 0; i < cols; i++)
		cout << estimatedCodeWord[i] << " ";
	cout << endl;*/
}

bool GDBFDecoder::calculateSyndrome()
{
	bool syndromeCheck = true;

	// The syndrome is basically the xor (mod 2 addition) of
	// the bits codeword with support in each row in H
	for(int i = 0; i < rows; i++)
	{
		syndrome[i] = 0; // reset the bit value

		for(int j = 0; j < rowOnes[i]; j++)
			syndrome[i] ^= hardDecisionVector[ rowSupport[i][j] ];

		if(syndrome[i] == 1)
			syndromeCheck = false;
	}

	return syndromeCheck;
}

////////////////////////////////////////
// Bit Flip Algorithm Class Functions //
////////////////////////////////////////

BitFlipDecoder::BitFlipDecoder(umat* Hptr, int maxNumIterations)
{
	// Max number of iterations for the sum product algorithm
	totalIterations = maxNumIterations;

	// Parity matrix dimensions
	rows = (*Hptr).n_rows;
	cols = (*Hptr).n_cols;

	// Create 2D dynamic array for rowSupport and colSupport
	// Dynamic memory allocation for row support and column support of H matrix
	rowSupport = new int*[rows](); // each pointer in rowSupport points to an array that holds the column index of the ones in the H matrix for that particular row
	colSupport = new int*[cols](); // each pointer in colSupport points to an array that holds the row index of the ones in the H matrix for that particular column

	rowOnes = new int[rows](); // hold the number of ones in each row
	colOnes = new int[cols](); // hold the number of ones in each column

	// Allocate the memory for the rows in H
	for(int i = 0; i < rows; i++)
	{
		uvec ones = find((*Hptr).row(i));
		rowSupport[i] = new int[ones.size()]();
		rowOnes[i] = ones.size();

		for(int j = 0; j < (int)ones.size(); j++)
			rowSupport[i][j] = ones[j];
	}

	// Allocate the memory for the columns in  H
	for(int i = 0; i < cols; i++)
	{
		uvec ones = find((*Hptr).col(i));
		colSupport[i] = new  int[ones.size()]();
		colOnes[i] = ones.size();

		for(int j = 0; j < (int)ones.size(); j++)
			colSupport[i][j] = ones[j];
	}

	failedChecks = new int[cols]();

	hardReceivedVector = new int[cols]();
	syndrome = new int[rows]();

	syndromeCheck = false;
}

BitFlipDecoder::~BitFlipDecoder()
{
	// Free up the memory from dynamic allocations
	for(int i = 0; i < rows; i++)
		delete [] rowSupport[i];
	delete [] rowSupport;

	for(int i = 0; i < cols; i++)
		delete [] colSupport[i];
	delete [] colSupport;

	delete [] rowOnes;
	delete [] colOnes;

	delete [] failedChecks;

	delete [] hardReceivedVector;
	delete [] syndrome;
}

void BitFlipDecoder::Decode(double* receivedVectorPtr, int* &estimatedCodeWord)
{
	pSoftReceivedVector = receivedVectorPtr;

	for(int i = 0; i < cols; i++)
		hardReceivedVector[i] = pSoftReceivedVector[i] < 0;

	for(int iterations = 0; iterations < totalIterations; iterations++)
	{
		syndromeCheck = true;

		// Calculate the syndrome
		// The syndrome is basically the xor (mod 2 addition) of
		// the bits codeword with support in each row in H
		for(int i = 0; i < rows; i++)
		{
			syndrome[i] = 0; // reset the bit value

			for(int j = 0; j < rowOnes[i]; j++)
				syndrome[i] ^= hardReceivedVector[ rowSupport[i][j] ];

			if(syndrome[i] == 1)
				syndromeCheck = false;
		}

		if(syndromeCheck)
			break;

		maxFailedChecks = 0;

		// If the syndrome is non-zero, find the  number of failed check sums for each bit
		for(int i = 0; i < cols; i++)
		{
			failedChecks[i] = 0;

			for(int j = 0; j < colOnes[i]; j++)
				failedChecks[i] += syndrome[ colSupport[i][j] ];

			if(failedChecks[i] > maxFailedChecks)
				maxFailedChecks = failedChecks[i];
		}

		// Flip the bits with the most number of failed checks
		for(int i = 0; i < cols; i++)
		{
			if(failedChecks[i] == maxFailedChecks)
				hardReceivedVector[i] = 1 - hardReceivedVector[i];
		}
	}

	for(int i = 0; i < cols; i++)
		estimatedCodeWord[i] = -2*hardReceivedVector[i]+1;

	return;
}

double get_wall_time(){
	struct timeval time;

	if(gettimeofday(&time,NULL)){
		return 0;
	}
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

double get_cpu_time(){
	return (double)clock() / CLOCKS_PER_SEC;
}
