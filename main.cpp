#include <iostream>
#include <vector>
#include <math.h>


/// Compute n!
int factorial(int n)
{
    if (n <= 1)
        return 1;

    return n*factorial(n-1);
}


/// Compute b^e when these are integers (avoids floating point errors).
/// Note: returns 0^0 = 1 and returns -1 (error) if b or e negative.
int pow(int b, int e)
{
    if (e==0) return 1;
    if (b==0) return 0;
    if (e<0) return -1;
    if (b<0) return -1;
    int result = 1;
    for (int i = 0; i < e; i++)
        result *= b;
    return result;
}


/// Replace 0 with 1, or 1 with 0.
/// (convert the int to a bool, then switch it using the ! operator, then convert back to int)
int complement(int a)
{
    return int ( ! bool (a) );
}


/// This needs to be defined before the functions below.
struct corner{
        int vertex;
        int firstSymbol[2];
        int secondSymbol[2];
    };


/// Return true if the corners have the the same data, false if not.
bool is_equal(std::vector<corner>::iterator a, corner* b)
{
    if (a->vertex != b->vertex) return false;
    if (a->firstSymbol[0] != b->firstSymbol[0]) return false;
    if (a->firstSymbol[1] != b->firstSymbol[1]) return false;
    if (a->secondSymbol[0] != b->secondSymbol[0]) return false;
    if (a->secondSymbol[1] != b->secondSymbol[1]) return false;

    return true;
}

bool is_equal(corner* a, corner* b)
{
    if (a->vertex != b->vertex) return false;
    if (a->firstSymbol[0] != b->firstSymbol[0]) return false;
    if (a->firstSymbol[1] != b->firstSymbol[1]) return false;
    if (a->secondSymbol[0] != b->secondSymbol[0]) return false;
    if (a->secondSymbol[1] != b->secondSymbol[1]) return false;

    return true;
}


/// Print a corner's contents.
void print_corner(corner* a)
{
    std::cout << "Vertex: " << a->vertex << ", " << "First symbol: " << a->firstSymbol[0]
        << " " << a->firstSymbol[1] << ", " << "Second symbol: " << a->secondSymbol[0]
        << " " << a->secondSymbol[1] << std::endl;
}

void print_corner(std::vector<corner>::iterator a)
{
    std::cout << "Vertex: " << a->vertex << ", " << "First symbol: " << a->firstSymbol[0]
        << " " << a->firstSymbol[1] << ", " << "Second symbol: " << a->secondSymbol[0]
        << " " << a->secondSymbol[1] << std::endl;
}


/// Stores all strings of <length> numbers, where each number is between 0 and <greatest>.
/// On first call, should have <place> = 0.
void StoreAllStrings(std::vector<std::vector<int>> &storage, int length, int greatest, int place )
{
    int counter = 0;
    int entry = 0;
    int numParts = pow(greatest+1, place+1);
    int partSize = pow(greatest+1, length - (place+1));

    for(int part = 0; part < numParts; part++)
    {
        for (counter = part*partSize; counter < (part+1)*partSize; counter++)
            storage[counter].push_back(entry);

        entry = (entry+1)%(greatest+1);
    }

    place++;
    if (place < length)
        StoreAllStrings(storage, length, greatest, place);
}



/// Instead of storing ALL strings, this will simply get the next string. It returns TRUE if a
/// new string was constructed and FALSE otherwise (i.e. it returns FALSE if the last string it
/// would construct is passed to it as an argument).
bool GetNextString(int *pString, int length, int greatest)
{
    /// Check if any of the digits can be increased. If so, increment the rightmost such digit.
    /// Then reset all digits to the right of that to be 0.
    /// Then return true (signifying that a new string was generated).
    for (int place = length-1; place > -1; place--)
        if (pString[place] < greatest)
        {
            pString[place] = pString[place] + 1;  /// Increment rightmost nonmaxed digit.
            if (place < length-1)                 /// Set all digits to the right of that to zero.
                for (int i = place+1; i < length; i++)
                    pString[i] = 0;

            return true;
        }

    return false;  /// If you make it to this line, then all the digits are set to the greatest
                   /// possible value. Hence, the last string was fed to the function.
}





int main()
{
    ///******************************************************
    ///     Recording the information about the group
    ///     and Cayley graph
    ///******************************************************


    /// Format for graphTable:
    /// 1. group elements are assigned a number 0, 1, ... , GROUPORDER-1.
    /// 2. row i, column 2j lists the product (element i)(generator j)
    /// 3. row i, column 2j+1 lists the product (element i)(generator j)^-1
    /// In other words, row i lists the product of element i with the zeroth generator, the inverse
    /// of the zeroth generator, the first generator, the inverse of the first generator, and so on.
    /// This contains all the information for the Cayley graph.

    /// Currently only a few graphs are stored, corresponding to standard 2-generator presentations
    /// of a few groups. You can choose which graph you want by changing GRAPHCHOICE. Make sure
    /// to also change GROUPORDER accordingly.

    const int GRAPHCHOICE = 4;  ///     1 = Z_2 X Z_2  Order 4
                                ///     2 = Z_2 X Z_3  Order 6
                                ///     3 = D_6  Order 6
                                ///     4 = Q_8  Order 8
                                ///     5 = Z_3 X Z_3  Order 9
                                ///     6 = A_4    Order 12 (note: this took over 24 hours to run
                                ///                          on my laptop.)

    int tempOrder;
    switch(GRAPHCHOICE) {
        case 1:
            tempOrder = 4;
            break;
        case 2:
            tempOrder = 6;
            break;
        case 3:
            tempOrder = 6;
            break;
        case 4:
            tempOrder = 8;
            break;
        case 5:
            tempOrder = 9;
            break;
        case 6:
            tempOrder = 12;
            break;
        default:
            std::cout << "Invalid graph choice; ending program." << std::endl;
            return 1;
            break;
    }
    const int GROUPORDER = tempOrder;
    const int NUMGENERATORS = 2;

    /// Depending on the number of systems that must be generated, it may be impossible to store them.
    /// E.g. GRAPHCHOICE = 6 gives a memory error. Therefore, you can choose whether to store
    /// all the systems or simply to generate the system, record the genus, and then delete
    /// the system.
    /// OPERATIONMODE = 0 will have the program store all systems.
    /// OPERATIONMODE = 1 will have the program store only the number of systems of each genus.

    const int OPERATIONMODE = 1;



    int graphTable [GROUPORDER][2*NUMGENERATORS];

    if (GRAPHCHOICE == 1)
    {
        int temp [][2*NUMGENERATORS] =
        {
            {1, 1, 2, 2},
            {0, 0, 3, 3},
            {3, 3, 0, 0},
            {2, 2, 1, 1}
        };
        for(int i = 0; i < GROUPORDER; i++)
            for(int j = 0; j < 2*NUMGENERATORS; j++)
                graphTable[i][j] = temp[i][j];
    }
    else if (GRAPHCHOICE == 2)
    {
        int temp [][2*NUMGENERATORS] =
        {
            {1, 1, 2, 3},
            {0, 0, 4, 5},
            {4, 4, 3, 0},
            {5, 5, 0, 2},
            {2, 2, 5, 1},
            {3, 3, 1, 4}
        };
        for(int i = 0; i < GROUPORDER; i++)
            for(int j = 0; j < 2*NUMGENERATORS; j++)
                graphTable[i][j] = temp[i][j];
    }
    else if (GRAPHCHOICE == 3)
    {
        int temp [][2*NUMGENERATORS] =
        {
            {3, 3, 1, 2},
            {5, 5, 2, 0},
            {4, 4, 0, 1},
            {0, 0, 4, 5},
            {2, 2, 5, 3},
            {1, 1, 3, 4}
        };
        for(int i = 0; i < GROUPORDER; i++)
            for(int j = 0; j < 2*NUMGENERATORS; j++)
                graphTable[i][j] = temp[i][j];
    }
    else if (GRAPHCHOICE == 4)
    {
        int temp [][2*NUMGENERATORS] =
        {
            {1,3,4,6},
            {2,0,7,5},
            {3,1,6,4},
            {0,2,5,7},
            {5,7,2,0},
            {6,4,1,3},
            {7,5,0,2},
            {4,6,3,1}
        };
        for(int i = 0; i < GROUPORDER; i++)
            for(int j = 0; j < 2*NUMGENERATORS; j++)
                graphTable[i][j] = temp[i][j];
    }
    else if (GRAPHCHOICE == 5)
    {
        int temp [][2*NUMGENERATORS] =
        {
            {1,2,3,6},
            {2,0,4,7},
            {0,1,5,8},
            {4,5,6,0},
            {5,3,7,1},
            {3,4,8,2},
            {7,8,0,3},
            {8,6,1,4},
            {6,7,2,5}
        };
        for(int i = 0; i < GROUPORDER; i++)
            for(int j = 0; j < 2*NUMGENERATORS; j++)
                graphTable[i][j] = temp[i][j];
    }
    else if (GRAPHCHOICE == 6)
    {
        int temp [][2*NUMGENERATORS] =
        {
            {1,2,3,4},
            {2,0,5,6},
            {0,1,7,8},
            {8,11,4,0},
            {9,5,0,3},
            {4,9,6,1},
            {10,7,1,5},
            {6,10,8,2},
            {11,3,2,7},
            {5,4,11,10},
            {7,6,9,11},
            {3,8,10,9}
        };
        for(int i = 0; i < GROUPORDER; i++)
            for(int j = 0; j < 2*NUMGENERATORS; j++)
                graphTable[i][j] = temp[i][j];
    }
    else
    {
        std::cout << "Nope.\n";
            return 1;
    }



    ///******************************************************
    ///     Vertex configurations
    ///******************************************************

    /// Throughout, let n = GROUPORDER and m = NUMGENERATORS.

    /// A system of vertex configurations is a map from G to the order lists
    /// of the symbols x_0^+, x_0^-, ... , x_(m-1)^+, x_(m-1)^- for the generators x_i.
    /// x_i^+ corresponds to (i, 0) and x_i^- corresponds to (i, 1).
    /// Then a configuration for element k is an array
    /// {k, k_1, 0 or 1, k_2, 0 or 1, ... , k_(2m), 0 or 1}
    /// where each {k_j, 0 or 1} represents the symbol x_(k_j)^(+ or -).
    /// A system of vertex configurations is an array of these - one for each k,
    /// 1 <= k <= n. We can delete the first coordinate k, and just use
    /// the index in the array as the group element.

    /// Sample system of vertex configurations:
    int system [4][8] =
    {
        {1, 0, 1, 1, 2, 0, 2, 1},
        {2, 0, 1, 0, 2, 1, 1, 1},
        {1, 1, 1, 0, 2, 0, 2, 1},
        {1, 0, 1, 1, 2, 0, 2, 1}
    };

    /// There are (2*m - 1)! possible vertex configurations, and
    /// ((2*m - 1)!)^n possible systems of vertex configurations. Thus even for Z_2 X Z_2,
    /// with n=4 and m=2, there are 6^4 = 1296 possible systems.


    ///********************************************************
    ///     Generating all systems of vertex configurations
    ///********************************************************


    /// 1.
    /// Find the (2*m - 1)! possible vertex configurations and store them in possibleConfigurations.

    /// possibleConfigurations will store the (2*m - 1)! vertex configurations, which are
    /// implemented as arrays of 4*m integers.

    /// (Right now, this is specific to the case of 2 generators only.)
    int possibleConfigurations [6][8] = {
        {0,0,0,1,1,0,1,1},
        {0,0,0,1,1,1,1,0},
        {0,0,1,0,0,1,1,1},
        {0,0,1,0,1,1,0,1},
        {0,0,1,1,0,1,1,0},
        {0,0,1,1,1,0,0,1}
    };

    /// 2.
    /// List all possible systems of vertex configurations. These are maps from G to
    /// possibleConfigurations. We can implement that as an array of n integers, where each
    /// integer is in the range [ 0, (2*NUMGENERATORS - 1)! ). That is, such an array A corresponds to
    /// such a map f according to f(i) <-> A[i].
    /// There are ((2*NUMGENERATORS - 1)!)^n of these.

    /// NOTE: path splits depending on choice of OPERATIONMODE. If OPMODE is 0, then we store
    /// all strings. If OPMODE is 1, then we generate them, find and record the genus, and
    /// delete them, one by one.


    int numChoices = factorial(2*NUMGENERATORS - 1);
    int numArrays = pow(numChoices, GROUPORDER);
    int* pSystemGenus = new int [numArrays];   /// Points to a int array that holds the system genera.
                                               /// For OPMODE == 0, this stores genus of each system.

    unsigned long int* pGenusCounts = new unsigned long [101];   /// Stores number of systems of any genus up to 100.
                                         /// Used when OPMODE == 1, instead of storing genus of each
                                         /// system.
    for (int i = 0; i < 100; i++) pGenusCounts[i] = 0;


    if (OPERATIONMODE == 0)  {  /// Store all strings mode
    std::vector< std::vector<int> > possibleSystems(numArrays);
    StoreAllStrings(possibleSystems, GROUPORDER, numChoices-1, 0);




    ///**********************************************************
    ///     Finding the surface corresponding to each system
    ///**********************************************************



    /// A corner is a subsequence {..., y, z, ...} of one of the configurations. We also
    /// allow {z, ... , y} to be a corner, since these are defined up to cyclic permutation.
    /// y, z stand for symbols like x^+, so y, z are pairs {i, 0 or 1}, {j, 0 or 1}.

    /// Corners are implemented as structures, with
    /// int vertex = vertex of the corner, an element of the group
    /// int firstSymbol[2], secondSymbol[2] = the two symbols stored as
    ///     arrays {i, 0 or 1} and {j, 0 or 1}. i, j, indicate the ith, jth generators.

    /// The struct must actually be declared above, before main(), in order to define the is_equal
    /// function. But here is the definition, commented out:
    ///
    ///        struct corner{
    ///            int vertex;
    ///            int firstSymbol[2];
    ///            int secondSymbol[2];
    ///        };



    int currentSystem [GROUPORDER];
    corner definingCorner, currentCorner;


    /// Loop over the systems in possibleSystems. For each system, determine the number of circuits
    /// from the system, which determines the genus of the surface. Store the genus in systemGenus.
    for (int i = 0; i < numArrays; i++)
    {

        /// This loop should 1. get the next system from possibleSystems, 2. fill a vector corners
        /// with all the corners in the system, 3. apply the circuit procedure to count the number
        /// of circuits in the system, 4. compute the genus of the resulting surface and store it
        /// in systemGenus.




        /// 1. Get the next system.
        for (int j = 0; j < GROUPORDER; j++)
            currentSystem[j] = possibleSystems[i][j];

        /// 2. Fill a vector with all corners in the system.
        std::vector<corner> corners;

        /// Loop over vertices.
        for (int j = 0; j < GROUPORDER; j++)
        {
            currentCorner.vertex = j;

            /// currentSystem[j] is the index in possibleConfigurations of the configuration for
            /// vertex j, so possibleConfigurations[currentSystem[j]] is the configuration
            /// for vertex j.
            currentCorner.firstSymbol[0] = possibleConfigurations[currentSystem[j]][4*NUMGENERATORS-2];
            currentCorner.firstSymbol[1] = possibleConfigurations[currentSystem[j]][4*NUMGENERATORS-1];
            currentCorner.secondSymbol[0] = possibleConfigurations[currentSystem[j]][0];
            currentCorner.secondSymbol[1] = possibleConfigurations[currentSystem[j]][1];

            corners.push_back(currentCorner);

            /// Loop over corners for vertex j.
            for (int k = 0; k < 2*NUMGENERATORS-1; k+=2)
            {
                currentCorner.firstSymbol[0] = possibleConfigurations[currentSystem[j]][k];
                currentCorner.firstSymbol[1] = possibleConfigurations[currentSystem[j]][k+1];
                currentCorner.secondSymbol[0] = possibleConfigurations[currentSystem[j]][k+2];
                currentCorner.secondSymbol[1] = possibleConfigurations[currentSystem[j]][k+3];

                corners.push_back(currentCorner);
            }
        } /// corners now contains all corners in the system.


        /// numCircuits is the number of circuits derived from the system (used to compute genus).
        int numCircuits = 0;

        /// 3. Count circuits using circuit procedure.
        /// Loop until all corners have been used in at least one circuit (i.e. corners is empty).
        while (!corners.empty())
        {

            /// This loop should 1. get a defining corner from corners, 2. inductively delete all
            /// corners in the circuit generated by defining corner from corners, and
            /// 3. increase numCircuits.



            /// 1. Get definingCorner (accessing/erasing elements is most efficient at end of vector).
            currentCorner = definingCorner = corners.back();
            corners.pop_back();



            /// 2. Inductively delete corners in the circuit generated by definingCorner.
            bool stopProcedure = false;
            while (!stopProcedure)
            {

                /// This loop should 1. set currentCorner.vertex to be the next vertex in the circuit,
                /// 2. set currentCorner.firstSymbol to be the inverse of currentCorner.secondSymbol,
                /// 3. search the configuration of currentCorner.vertex for the inverse of
                /// currentCorner.secondSymbol, then 4. set currentCorner.secondSymbol to be the next
                /// symbol in the configuration. Then 5. if this corner is the defining corner then
                /// set bool to stop loop, otherwise remove the corner from corners.


                /// 1. Set currentCorner.vertex.
                /// If secondSymbol = {j, d} then the next symbol is the product of
                /// currentCorner.vertex and the jth generator (if d == 0, or its inverse if d == 1).
                if (currentCorner.secondSymbol[1] == 0)
                    currentCorner.vertex =
                    graphTable[currentCorner.vertex][2*currentCorner.secondSymbol[0]];
                    else
                        currentCorner.vertex =
                        graphTable[currentCorner.vertex][2*currentCorner.secondSymbol[0] + 1];


                /// 2. Set currentCorner.firstSymbol to the inverse of currentCorner.secondSymbol.
                currentCorner.firstSymbol[0] = currentCorner.secondSymbol[0];
                currentCorner.firstSymbol[1] = complement(currentCorner.secondSymbol[1]);


                /// 3. Search configuration of currentCorner.vertex for inverse of second symbol,
                /// which is now stored as the first symbol.
                /// Note: the configuration of currentCorner. in currentSystem is the array
                /// possibleConfigurations[currentSystem[currentCorner.vertex]].
                int a; /// will store the index (in the configuration) of the first symbol.
                for (a = 0; a < 4*NUMGENERATORS; a += 2)
                    if (    possibleConfigurations[currentSystem[currentCorner.vertex]][a]
                            == currentCorner.firstSymbol[0]
                        &&  possibleConfigurations[currentSystem[currentCorner.vertex]][a+1]
                            == currentCorner.firstSymbol[1] )
                            break;


                /// 4. Set currentCorner.secondSymbol to be the next symbol in the configuration.
                /// If a is the index of the last symbol in the array, then the "next" symbol is
                /// the first symbol in the array (we consider these to be "cyclic lists").
                if (a == 4*NUMGENERATORS - 2)
                {
                    currentCorner.secondSymbol[0] =
                        possibleConfigurations[currentSystem[currentCorner.vertex]][0];
                    currentCorner.secondSymbol[1] =
                        possibleConfigurations[currentSystem[currentCorner.vertex]][1];
                }
                else
                {
                    currentCorner.secondSymbol[0] =
                        possibleConfigurations[currentSystem[currentCorner.vertex]][a+2];
                    currentCorner.secondSymbol[1] =
                        possibleConfigurations[currentSystem[currentCorner.vertex]][a+3];
                }


                /// 5. If currentCorner == definingCorner set bool to stop loop.
                /// Otherwise search for currentCorner in corners and erase it.
                if (is_equal(&currentCorner, &definingCorner))
                    stopProcedure = true;
                else
                {
                    /// Search for the corner in corners then erase it.
                    for (std::vector<corner>::iterator p = corners.begin(); p != corners.end(); p++)
                        if (is_equal(p, &currentCorner))
                        {
                            corners.erase(p);
                            break;
                        }
                }


            } /// End of circuit procedure.



            /// 3. Increase numCircuits.
            numCircuits++;


        } /// End of computation of numCircuits.


                /// This module can be used to print the system if it has a given number of
                /// circuits, in this case 6.
                    /*
                    if (numCircuits == 6)
                    {
                        std::cout << "\nSYSTEM:\n";
                        for(int h = 0; h < 4; h++)
                        {
                            std::cout << "Vertex " << h << ": ";
                            for(int y = 0; y < 8; y++)
                                std::cout << possibleConfigurations[currentSystem[h]][y] << ", ";
                            std::cout << std::endl;
                        }
                    }
                    */


        /// 4.
        /// numCircuits = number of faces, n = number of vertices,
        /// m*n = number of edges. So the Euler characteristic is n - m*n + numCircuits = 2-2g
        /// where g is the genus, so the genus is (n - m*n + numCircuits - 2)/(-2).

        pSystemGenus[i] = (GROUPORDER - NUMGENERATORS*GROUPORDER + numCircuits - 2)/(-2);


    } /// End of computation of the genus of each surface.

    } /// End of computation for OPERATIONMODE == 0 (store all strings mode).




    /// If OPERATIONMODE == 1, the computation is almost the same. The difference is that instead
    /// of getting each system from an array containing all the systems, we just generate the next
    /// system until we've found the genus of each one. This is done by the GetNextString function,
    /// which will return FALSE when it is given the last system.

    else if (OPERATIONMODE == 1) {   /// String by string mode.

    int currentSystem [GROUPORDER];  /// currentSystem will hold each system; after its genus is
                                     /// is computed, it will be replaced by the next system.

    for (int i = 0; i < GROUPORDER; i++)  /// Initialize to string of all zeros.
        currentSystem[i] = 0;

    int * pCurrentSystem = currentSystem;
    corner definingCorner, currentCorner;



    /// Loop generates the next system, then computes and records its genus.
    /// If the last system was generated, GetNextString() will return false and exit the loop.
    int systemCounter = 0;  /// This is used to store the genus as pSystemGenus[systemCounter].
    bool isNewString = true;
    while (isNewString)
    {
        std::vector<corner> corners;

        /// Get all corners.
        for (int j = 0; j < GROUPORDER; j++)
        {
            currentCorner.vertex = j;

            currentCorner.firstSymbol[0] = possibleConfigurations[currentSystem[j]][4*NUMGENERATORS-2];
            currentCorner.firstSymbol[1] = possibleConfigurations[currentSystem[j]][4*NUMGENERATORS-1];
            currentCorner.secondSymbol[0] = possibleConfigurations[currentSystem[j]][0];
            currentCorner.secondSymbol[1] = possibleConfigurations[currentSystem[j]][1];

            corners.push_back(currentCorner);

            /// Loop over corners for vertex j.
            for (int k = 0; k < 2*NUMGENERATORS-1; k+=2)
            {
                currentCorner.firstSymbol[0] = possibleConfigurations[currentSystem[j]][k];
                currentCorner.firstSymbol[1] = possibleConfigurations[currentSystem[j]][k+1];
                currentCorner.secondSymbol[0] = possibleConfigurations[currentSystem[j]][k+2];
                currentCorner.secondSymbol[1] = possibleConfigurations[currentSystem[j]][k+3];

                corners.push_back(currentCorner);
            }
        } /// corners now contains all corners in the system.

        int numCircuits = 0;

        /// Count circuits.
        /// Loop until all corners have been used in at least one circuit (i.e. corners is empty).
        while (!corners.empty())
        {
            /// 1. Get definingCorner.
            currentCorner = definingCorner = corners.back();
            corners.pop_back();

            /// 2. Inductively delete corners in the circuit generated by definingCorner.
            bool stopProcedure = false;
            while (!stopProcedure)
            {
                if (currentCorner.secondSymbol[1] == 0)
                    currentCorner.vertex =
                    graphTable[currentCorner.vertex][2*currentCorner.secondSymbol[0]];
                    else
                        currentCorner.vertex =
                        graphTable[currentCorner.vertex][2*currentCorner.secondSymbol[0] + 1];

                currentCorner.firstSymbol[0] = currentCorner.secondSymbol[0];
                currentCorner.firstSymbol[1] = complement(currentCorner.secondSymbol[1]);

                int a; /// will store the index (in the configuration) of the first symbol.
                for (a = 0; a < 4*NUMGENERATORS; a += 2)
                    if (    possibleConfigurations[currentSystem[currentCorner.vertex]][a]
                            == currentCorner.firstSymbol[0]
                        &&  possibleConfigurations[currentSystem[currentCorner.vertex]][a+1]
                            == currentCorner.firstSymbol[1] )
                            break;

                if (a == 4*NUMGENERATORS - 2)
                {
                    currentCorner.secondSymbol[0] =
                        possibleConfigurations[currentSystem[currentCorner.vertex]][0];
                    currentCorner.secondSymbol[1] =
                        possibleConfigurations[currentSystem[currentCorner.vertex]][1];
                }
                else
                {
                    currentCorner.secondSymbol[0] =
                        possibleConfigurations[currentSystem[currentCorner.vertex]][a+2];
                    currentCorner.secondSymbol[1] =
                        possibleConfigurations[currentSystem[currentCorner.vertex]][a+3];
                }

                if (is_equal(&currentCorner, &definingCorner))
                    stopProcedure = true;
                else
                {
                    for (std::vector<corner>::iterator p = corners.begin(); p != corners.end(); p++)
                        if (is_equal(p, &currentCorner))
                        {
                            corners.erase(p);
                            break;
                        }
                }
            } /// End of circuit procedure.

            numCircuits++;

        } /// End of computation of numCircuits.


                /// This module can be used to print the system if it has a given number of
                /// circuits, in this case 6.
                    /*
                    if (numCircuits == 6)
                    {
                        std::cout << "\nSYSTEM:\n";
                        for(int h = 0; h < 4; h++)
                        {
                            std::cout << "Vertex " << h << ": ";
                            for(int y = 0; y < 8; y++)
                                std::cout << possibleConfigurations[currentSystem[h]][y] << ", ";
                            std::cout << std::endl;
                        }
                    }
                    */

        /// Compute and record the system's genus.
        int gen = (GROUPORDER - NUMGENERATORS*GROUPORDER + numCircuits - 2)/(-2);
        pGenusCounts[gen]++;

        /// Get the next system; if return is false, set isNewString = false to terminate loop.
        isNewString = GetNextString(pCurrentSystem, GROUPORDER, numChoices-1);



    } /// End of computation of the genus of each surface.

    } /// End of computation for OPERATIONMODE == 1 (string by string mode)


    else  ///OPERATIONMODE must be 0 or 1, otherwise give an error message.
    {
        std::cout << "Bad OPERATIONMODE.\n";
        return 1;
    }



    ///********************************************
    ///     Presenting the results
    ///********************************************


    /// Count the number of surfaces of each genus.
    /// Then print the counts.


    if (OPERATIONMODE == 0){

    int maxGenus, minGenus;
    maxGenus = minGenus = 0;

    for (int i = 0; i < numArrays; i++)
    {
        if (pSystemGenus[i] > maxGenus)
            maxGenus = pSystemGenus[i];
        if (pSystemGenus[i] < minGenus)
            minGenus = pSystemGenus[i];
    }

    int surfaceGenusCounts[maxGenus - minGenus + 1];

    for (int i = 0; i < maxGenus - minGenus + 1; i++)
        surfaceGenusCounts[i] = 0;

    for (int i = 0; i < numArrays; i++)
        for (int j = 0; j < maxGenus - minGenus + 1; j++)
            if (pSystemGenus[i] == minGenus + j)
                surfaceGenusCounts[j]++;

    for (int i = 0; i < maxGenus - minGenus + 1; i++)
        std::cout << "Genus " << minGenus + i << ": " << surfaceGenusCounts[i]
        << "  -----  " << double(surfaceGenusCounts[i])/double(numArrays) << std::endl;

    std::cout << "\nNumber of flip-inequivalent decompositions:\n";
    for (int i = 0; i < maxGenus - minGenus + 1; i++)
        std::cout << "Genus " << minGenus + i << ": " << surfaceGenusCounts[i]/2
        << "  -----  " << (double(surfaceGenusCounts[i])/double(numArrays)) << std::endl;


    } /// End of OPMODE == 0


    else if (OPERATIONMODE == 1){

        for (int i = 0; i < 100; i++)
            if (pGenusCounts[i] > 0)
                std::cout << "Genus " << i << ": " << pGenusCounts[i]
                << "  -----  " << double(pGenusCounts[i])/double(numArrays) << std::endl;

        std::cout << "\nNumber of flip-inequivalent decompositions:\n";
        for (int i = 0; i < 100; i++)
            if (pGenusCounts[i] > 0)
                std::cout << "Genus " << i << ": " << pGenusCounts[i]/2
                << "  -----  " << (double(pGenusCounts[i])/double(numArrays)) << std::endl;


    } /// End of OPMODE == 1

    else {std::cout << "Bad OPERATIONMODE.\n"; return 1;}

}
