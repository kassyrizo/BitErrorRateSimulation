#include <armadillo>
#include <cmath>
#include <algorithm>
#include "Puncture.h"


vector< vector<int> > total_combs;

vector <int> PunctureCode(umat H)
{
    cout << "Begin Puncturing..." << endl;
    umat newH = H;

    // Variables for Algorithm One
    vector <int> INDEX;
    vector <int> REMAIN;
    vector <int> NEIBOR;
    vector<int> J;

    // Variables for Algorithm Two
    vector<int> S;
    vector<int> I;
    vector<int> B;
    vector< vector<int> > Si;

    // Return value, the bits to be punctured
    vector <int> puncturedBits;
    vector <int> cantBePunctured;


    /*************************************************/
    /*************************************************/
    /************** BEGIN ALGORITHM ONE **************/
    /*************************************************/
    /*************************************************/

    #if defined(DEBUG)
        cout << "\n************ BEGIN ALGORITHM 1 ************" << endl;
    #endif
    do
    {
        /************************************************************************/
        /********* Step 1: REMAIN = {variable nodes in H}, NEIBOR = {0} *********/
        /************************************************************************/
    	NEIBOR.clear();
        REMAIN.clear();

        // Set REMAIN to be the indices of the remaining variable nodes in H
        for(int row = 0; row < (int)newH.n_rows; row++)
        {
            uvec remaining = find(newH.row(row)); // Find the indices of all remaining variable nodes in this row (i.e. 1's in the row)

            for(int j = 0; j < (int)remaining.size(); j++)
            {
                if(!(find(REMAIN.begin(), REMAIN.end(), remaining.at(j)) != REMAIN.end()))
                    REMAIN.push_back(remaining.at(j)); // If those variable nodes are not already in REMAIN, add them
            }
        }
        sort(REMAIN.begin(), REMAIN.end());

        #if defined(DEBUG)
            cout << "\n\nStep 1: Initialization" << endl;
            cout << "INDEX: ";
            vector_print(INDEX);
            cout << "REMAIN: ";
            vector_print(REMAIN);
            cout << "NEIBOR: ";
            vector_print(NEIBOR);
        #endif
            
        do
        {
            /***************************************************************************/
            /********* Step 3: Find vertex, v, uniformly at random from REMAIN *********/
            /***************************************************************************/
            mat r(1,1);
            r.randu();
            int v = floor(r(0)*(REMAIN.size()));
            v = REMAIN[v];

            #if defined(DEBUG)
                cout << "\n\nStep 3: Find vertex v uniformly at random from REMAIN" << endl;
                cout << "v: " << v << endl;
            #endif

            /***********************************************************************************/
            /********* Step 4: Set INDEX = INDEX + {v} and NEIBOR = NEIBOR + Ng(Ng(v)) *********/
            /***********************************************************************************/
            INDEX.push_back(v);

            // Set NEIBOR = NEIBOR + NG(NG(v))
            uvec rowIndices = find(newH.col(v)); // Find the rows (i.e. check equations) that involve node v

            // For each check node (row) that was found, add all of it's neighbours (indices of 1 positions) to NEIBOR
            for(int i = 0; i < (int)rowIndices.size(); i++)
            {
                // Find the indices of all neighbours (i.e. the positions of the 1's in the current row)
                uvec neighbors = find(newH.row(rowIndices.at(i)));

                for(int j = 0; j < (int)neighbors.size(); j++)
                {
                    if(!(find(NEIBOR.begin(), NEIBOR.end(), neighbors.at(j)) != NEIBOR.end()))
                        NEIBOR.push_back(neighbors.at(j)); // If those neighbors are not already in NEIBOR, add them
                }
            }

            sort(NEIBOR.begin(), NEIBOR.end()); // sort for easier debugging

            #if defined(DEBUG)
                cout << "\n\nStep 4: INDEX = INDEX+v, NEIBOR = NEIBOR+Ng(Ng(v))" << endl;

                cout << "INDEX: ";
                vector_print(INDEX);

                cout << "NEIBOR: ";
                vector_print(NEIBOR);
            #endif

            /************************************************/
            /********* Step 5: REMAIN=REMAIN-NEIBOR *********/
            /************************************************/

            vector<int> difference(REMAIN.size()); // holds the set differencec between REMAIN and NEIBOR
            vector<int>::iterator difference_it; // Iterator for the set difference operation

            sort(REMAIN.begin(), REMAIN.end()); // set_difference requires vectors involved to be sorted. NEIBOR was sorted earlier

            difference_it = set_difference(REMAIN.begin(), REMAIN.end(), NEIBOR.begin(), NEIBOR.end(), difference.begin()); // get the set difference
            difference.resize(difference_it-difference.begin());

            REMAIN = difference; // REMAIN is now everything that was in REMAIN that was not in NEIBOR (REMAIN = REMAIN\NEIBOR}

            #if defined(DEBUG)
                cout << "\n\n Step 5: REMAIN=REMAIN-NEIBOR" << endl;
                cout << "REMAIN: ";
                vector_print(REMAIN);
            #endif

        /*****************************************************************************/
        /********* Step 2: If REMAIN = {0}, go to Step 3, else, go to Step 6 *********/
        /*****************************************************************************/
        }while(!REMAIN.empty());

        /**********************************************************************/
        /********* Step 6:  Reinitialize, NEIBOR = {1,2....n} - INDEX *********/
        /**********************************************************************/

        // Set NEIBOR to be the indices of the variable nodes that are not in INDEX
        for(int row = 0; row < (int)newH.n_rows; row++)
        {
            uvec remaining = find(newH.row(row)); // Find the indices of all remaining variable nodes in this row

            for(int j = 0; j < (int)remaining.size(); j++)
            {
                if(!(find(NEIBOR.begin(), NEIBOR.end(), remaining.at(j)) != NEIBOR.end()))
                    NEIBOR.push_back(remaining.at(j)); // If the variable node is not already in NEIBOR, add it
            }
        }

        vector<int> tempIndex = INDEX; // holds a temporary copy of the INDEX (I didn't want to sort the actual INDEX, since I'm not sure if order matters with puncturing at this time)
        vector<int> difference(NEIBOR.size()); // holds the set difference between tempINDEX and NEIBOR
        vector<int>::iterator difference_it; // Iterator for the set difference operation

        sort(NEIBOR.begin(), NEIBOR.end()); // set_difference requires sorted vectors to work properly
        sort(tempIndex.begin(), tempIndex.end());

        difference_it = set_difference(NEIBOR.begin(), NEIBOR.end(), tempIndex.begin(), tempIndex.end(), difference.begin()); // get the set difference
        difference.resize(difference_it-difference.begin());

        NEIBOR = difference; //NEIBOR is now everything that was in NEIBOR that was not in INDEX

        #if defined(DEBUG)
            cout << "\n\nStep 6: Reinitialize, NEIBOR = {1,2....n} - INDEX" << endl;
            cout << "NEIBOR: ";
            vector_print(NEIBOR);
        #endif

        /**********************************************************************************************/
        /********* Step 7:  Set J to be set of check nodes involving only bit nodes in NEIBOR *********/
        /**********************************************************************************************/

        J.clear();
        for(int row = 0; row < (int)newH.n_rows; row++)
        {
            uvec ones = find(newH.row(row)); // indices of 1's in row
            vector<int> intersection(newH.n_cols); // holds the set intersection between ones and tempINDEX (i.e. check for any ones involved in INDEX)
            vector<int>::iterator intersection_it; // Iterator for the set intersection operation

            sort(ones.begin(), ones.end()); // set_intersection requires sorted vectors to work properly

            intersection_it = std::set_intersection(ones.begin(), ones.end(), tempIndex.begin(), tempIndex.end(), intersection.begin()); // get the set intersection
            intersection.resize(intersection_it-intersection.begin());

            // If none of the bit nodes in INDEX are present in the row, add it to J
            if(intersection.empty())
                J.push_back(row);
        }

        #if defined(DEBUG)
            cout << "\n\nStep 7: Set J to be set of check nodes involving only bit nodes in NEIBOR" << endl;
            cout << "J: ";
            vector_print(J);
        #endif

        /*************************************************************************************************/
        /********* Step 8/Step9: Set G' = G[Ng(J)+J], H' is parity matrix that corresponds to G' *********/
        /*************************************************************************************************/

        umat tempH; // matrix for the temporary H matrix
        tempH.set_size(J.size(), newH.n_cols); // the new H will have J rows, and n columns

        for(int index = 0; index < (int)J.size(); index++)
            tempH.row(index) = newH.row(J[index]);

        newH.resize(tempH.n_rows, tempH.n_cols);
        newH = tempH; // set H'

        #if defined(DEBUG)
            cout << "\n\nStep 8/Step 9: Set G' = G[Ng(J)+J], H' is parity matrix that corresponds to G'" << endl;
            cout << "H': " << endl;
            cout << newH;
        #endif

        S.push_back(INDEX.size()-1);//Hold the last index of the set, this is for algorithm 2

    /*****************************************************************************/
    /********* Step 10: If H' = {0}, halt algorithm 1, else go to Step 1 *********/
    /*****************************************************************************/
    }while(!newH.empty()); // END ALGORITHM ONE

    /*************************************************/
    /*************************************************/
    /************** BEGIN ALGORITHM TWO **************/
    /*************************************************/
    /*************************************************/
    //cout << "INDEX: ";
    //vector_print(INDEX);

    #if defined(DEBUG)
        cout << "\n************ TERMINATE ALGORITHM 1 ************" << endl;
        cout << "\n\n************ BEGIN ALGORITHM 2 ************" << endl;
    #endif

    printf( " Number of punctured bits:        ");
    /********************************************************************************/
    /********* Step 1: Generate INDEX and Si, i = 1,...,L using Algorithm 1 *********/
    /********************************************************************************/

    // Generate Si. The vector S holds the indices of the last element of each set in Si
    for(int index = 0; index < (int)S.size(); index++)
    {
        vector<int> newS; // new vector to add to Si
        vector<int>::iterator fill_it; // iterator for filling vector<int>

        if(index == 0)
            fill_it = INDEX.begin();
        else
            fill_it = INDEX.begin()+S[index-1] + 1;

        newS.assign(fill_it, INDEX.begin() + S[index] + 1);

        Si.push_back(newS); // add new vector to list
    }

    #if defined(DEBUG)
        cout << "\n\nStep 1: Generate INDEX and Si, i = 1,...,L using Algorithm 1" << endl;
        cout << "INDEX: ";
        vector_print(INDEX);
        cout << "S: ";
        vector_print(S);
        cout << "Si: \n";
        vector_print(Si);
    #endif


    vector< vector<int> >::iterator Si_iterator = Si.begin(); // iterator for determining when to terminate the dowhile loop

    do
    {
        /****************************************************************************************************************************/
        /********* Step 2: Set I to be the set of all rows of H that involve bit nodes from exactly one set Sj, 0 <= j <= L *********/
        /****************************************************************************************************************************/

        // Find the rows in H that involve only bits from the current set, Sj
        for(int row = 0; row < (int)H.n_rows; row++)
        {
            uvec ones = find(H.row(row)); // find the positions of the ones in the current row
            vector<int> intersection(H.n_cols); // holds the set intersection of the position of the current ones and the values in Sj
            vector<int>::iterator intersection_it;

            sort(ones.begin(), ones.end()); // set_intersection requires sorted vectors to work properly
            sort(Si_iterator->begin(), Si_iterator->end());

            intersection_it = set_intersection(ones.begin(), ones.end(), Si_iterator->begin(), Si_iterator->end(), intersection.begin()); // get the intersection
            intersection.resize(intersection_it-intersection.begin());

            if(!intersection.empty())
            {
                I.push_back(row);
            }

        }




        #if defined(DEBUG)
            cout << "\n\nStep 2: Set I to be the set of all rows of H that involve bit nodes from exactly one set, Sj" << endl;
            cout << "I: ";
            vector_print(I);
        #endif

        do
        {
            /**************************************************************************************************************/
            /********* Step 3: Set B to be the set of all those bit nodes that are not in INDEX that are involved *********/
            /********* in the least number of equations represented in the rows I of H and such that no stopping  *********/
            /********* set is contained in INDEX + {v}                                                            *********/
            /**************************************************************************************************************/

            vector<int> counter;
            B.clear();

            // B is the set of all nodes not in INDEX involved in rows, I
            for(int v = 0; v < (int)H.n_cols; v++)
            {
                // B is not in INDEX
                if(!(find(INDEX.begin(), INDEX.end(), v) != INDEX.end()))
                {
                    // B is in I
                    for(int row = 0; row < (int)I.size(); row++)
                    {
                        if(H(I[row],v) == 1)
                        {
                            if(!(find(B.begin(), B.end(), v) != B.end())) // if v hasn't been added, ad it
                            {
                                B.push_back(v);
                                counter.push_back(1);
                            }
                            else
                                counter.back() = counter.back()+1;
                        }
                    }
                }
            }

            #if defined(DEBUG)
                cout << "\n\nStep 3, Part 1: Set B to be the set of all nodes not in INDEX involved in rows I" << endl;
                cout << "B: ";
                vector_print(B);
                cout << "counter: ";
                vector_print(counter);
            #endif

            // Find minimum number of times that variables in B show up in check equations I
            int minimum = *std::min_element(counter.begin(), counter.end());

            // Erase all elements in B that do not show up the minimum amount of time in I
            for(int j = 0; j < (int)counter.size(); j++)
            {
                if(counter[j] != minimum)
                {
                    B.erase(B.begin()+j);
                    counter.erase(counter.begin()+j);
                    --j;
                }
            }

            // Erase all elements in B that have already been found to create stopping sets
            for(int j = 0; j < (int)cantBePunctured.size(); j++)
            {
                if(find(B.begin(), B.end(), cantBePunctured[j]) != B.begin())
                    B.erase(remove(B.begin(), B.end(), cantBePunctured[j]), B.end());
            }

            #if defined(DEBUG)
                cout << "\n\nStep 3, Part 2: Set B to be the bit nodes involved in the least number of equations" << endl;
                cout << "minimum: " << minimum << endl;
                cout << "B: ";
                vector_print(B);
            #endif

            //Check each value of B, make sure they don't form a stopping set with the bits in INDEX
            vector<int> tempINDEX = INDEX;

            for(int k = 0; k < (int)B.size(); k++)
            {
                tempINDEX.push_back(B[k]); // tempINDEX = INDEX

                #if defined(DEBUG)
                    cout << "\n\nStep 3, Part 3: Find stopping Sets, B = " << B[k] << endl;
                #endif

                // Use the stopping set algorithm found in a paper to find the stopping sets
                // First step is initialization, set p_prime = puncturing status of all variable nodes
                // Puncturing status is 0 if bit is punctured, 1 if bit is unpunctured

                umat p_prime = ones<umat>(1, H.n_cols);
                for(int index = 0; index < (int)INDEX.size(); index++)
                    p_prime[INDEX[index]] = 0;

                // Then, set p'(v*) = 0, in this case, v* = B[k], since we are trying to puncture this bit
                p_prime[B[k]] = 0;

                // Set the maximum number of iterations to the size of (INDEX) + B[k]
                int max_iterations = INDEX.size()+1;

                #if defined(DEBUG)
                    cout << "\nStopping Set Check Algorithm, Step 1: Initialization, set p'(v) = p(v) and set p'(v*) = 0" << endl;
                    cout << "p_prime: " << p_prime << endl;
                    cout << "maximum number of iterations: " << max_iterations << endl;
                #endif

                // Run the algorithm
                for(int iteration = max_iterations; iteration > 0; iteration--)
                {
                    //umat tempp_prime = p_prime;

                    // for each punctured bit
                    uvec puncBits = find(p_prime == 0);

                    for(int index = 0; index < (int)puncBits.size(); index++)
                    {
                        // find all connected check nodes
                        uvec connectedCheckNodes = find(H.col(puncBits[index]));

                        // check the neighbouring variable nodes to each connected check node
                        for(int chkNode = 0; chkNode < (int)connectedCheckNodes.size(); chkNode++)
                        {
                            // find the neighbouring check nodes that aren't the node in question
                            uvec varNodes = find(H.row(connectedCheckNodes[chkNode]));
                            vector<int> connectedVarNodes;
                            connectedVarNodes.assign(varNodes.begin(), varNodes.end());
                            connectedVarNodes.erase(remove(connectedVarNodes.begin(), connectedVarNodes.end(), puncBits[index]), connectedVarNodes.end());

                            int canBePunctured = 1;

                            for(int node = 0; node < (int)connectedVarNodes.size(); node++)
                            {
                                if(p_prime.at(connectedVarNodes.at(node)) == 0)
                                    canBePunctured = 0;
                            }

                            if(canBePunctured)
                            {
                                p_prime[puncBits[index]] = 1;
                                break;
                            }
                        }
                    }

                    // If there are no punctured bits left, stop running the algorithm
                    uvec anyZeros = find(p_prime == 0);
                    if(anyZeros.empty())
                        break;

                    /*// If no bit was updated, stop running the algorithm
                    umat equals = (tempp_prime == p_prime);
                    if(all(vectorise(equals)))
                        break;*/
                }

                uvec anyZeros = find(p_prime == 0);

                // if there are any zeros, we can't puncture this bit and a stopping set was found.
                if(!anyZeros.empty())
                {
                    #if defined(DEBUG)
                        cout << "A stopping set was found! " << p_prime << endl;
                    #endif
                        B.erase(B.begin()+k);
                        cantBePunctured.push_back(B[k]);
                        --k;
                }
            }

            #if defined(DEBUG)
                cout << "Final B: ";
                vector_print(B);
            #endif

            if(B.empty())
                break;

            /********************************************************************************************/
            /********* Step 5: Select b from B randomly, and set INDEX = INDEX + {b}, B = B\{b} *********/
            /********************************************************************************************/

            mat r(1,1);
            r.randu();
            int b = floor(r(0)*(B.size()));
            b = B[b];

            INDEX.push_back(b);

            B.erase(remove(B.begin(), B.end(), b), B.end());

            #if defined(DEBUG)
                cout << "Step 5: Select b from B randomly, and set INDEX = INDEX + {b}, B = B\{b}" << endl;
                cout << "b: " << b << endl;
                cout << "INDEX: ";
                vector_print(INDEX);
                cout << "B: ";
                vector_print(B);
            #endif


            printf("\b\b\b\b\b\b\b%7ld", INDEX.size());
            fflush(0);

        /******************************************************************************************/
        /********* Step 4: If B = {0}, try again with new Si, otherwise go back to Step 3 *********/
        /******************************************************************************************/
        }while(!B.empty());

        I.clear();
        Si_iterator++;

    }while(Si_iterator != Si.end());

    #if defined(DEBUG)
        cout << "\n\nStep 6: Output INDEX and halt" << endl;
    #endif
    cout << "\nFinished Puncturing!" << endl;
    cout << "Total Number of Punctured Bits: " << INDEX.size() << endl;
    puncturedBits = INDEX;

    return puncturedBits;
}

void vector_print(vector<int> vectorToPrint)
{
    vector<int>::iterator it;

    for(it=vectorToPrint.begin(); it!=vectorToPrint.end();++it)
        cout << *it << " ";
    cout << endl;
}

void vector_print(vector< vector<int> > vectorToPrint)
{
    vector< vector<int> >::iterator vector_it;
    vector<int>::iterator int_it;

    for(vector_it = vectorToPrint.begin(); vector_it != vectorToPrint.end(); ++vector_it)
    {
        for(int_it = vector_it->begin(); int_it != vector_it->end(); ++int_it)
            cout << *int_it << " ";
        cout << endl;
    }
}

