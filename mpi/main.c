### DEMET YAYLA

#include <iostream>
#include<fstream>
#include<set>
#include<iterator>
#include<map>
#include<queue>
#include <metis.h>
#include <time.h>
#include <mpi.h>

#define epsilon 0.000001
#define alpha 0.2
#define linecount 16741171 //number of lines in graph.txt
#define FILENAME "grapht.txt"

using namespace std;

/*
rth row and jth column (r and j starting from 0)
a. increase row_begin(r+1) by 1
b. sum from 0 to r-1. values[sum+1] is the new node's position
c. col_indices[sum+1] is its equivalent column reference

3 r 4 c
in-degree --> row sum
oudegree --> column sum

If there is an edge from k to l, the kth column lth row is 1

*/
class CSR{
    public:
        vector<int> row_begin;
        CSR(int nodeNum);
        double* values;
        vector<int> col_indices;
};

CSR::CSR(int nodeNum){
    this->row_begin = vector<int>(nodeNum+1);
    this->col_indices = vector<int>();
}

//creates a CSR object representing the matrice given in graph.txt. The returned matrix is a directed graph's representation
void createCSR(CSR* csr, int numOfNodes, set<int>* adjacencyList,  vector<int>* outdegrees){
    //the names for the three fields of CSR matrix format are preserved as in hw description
    printf("--> Creating CSR format\n");
    csr = new CSR(numOfNodes);
    csr->row_begin[0] = 0; //sets the initial element of row_begin as 0
    int priorLen = 0; //while adding to csr, keep the prior len of the previous row
    int priorNonzero = 0; //holds the index of the former non-zero element
    for(int i = 0; i < numOfNodes; i++){
        if(adjacencyList[i].size() == 0){ //if a row has zero adjacencies its corresponding row_begin entry is equal to zero
            csr->row_begin[i] = 0; 
        }else{
            csr->row_begin[i] = priorNonzero + priorLen; //finding the current row's (i'th) row_begin value
            priorNonzero = priorNonzero + priorLen;
            priorLen = adjacencyList[i].size();
            //adjacencyList holds the id for adjacent nodes thus the indexes of the non-zero entries in a row. Thus appending it to col_indices directly forms a correct col_indices below
            std::copy(adjacencyList[i].begin(), adjacencyList[i].end(), std::back_inserter(csr->col_indices));
        }
    }

    csr->row_begin[numOfNodes] = priorNonzero + priorLen; //number of non-zero entries in a matrix

    int filledSize = csr->col_indices.size(); //size of csr->col_indices is equal to size of csr->values which is filledSize
    csr->values = new double[filledSize];
    for(int i = 0; i < filledSize; i++){
        if(outdegrees->at(csr->col_indices[i])==0){
            continue; //to eliminate arithmetic division by zero error
        }
        csr->values[i] = 1. / outdegrees->at(csr->col_indices[i]); //line provides that the sum of values for a column should give 1
    }

    printf("--> Created CSR format\n");
}


//the output value for the function is idx_t* partitioning given in arguments. it is an array for which i'th element corresponds to the processor index - 1 the i'th node is assigned to.
int partitionWithMETIS(int worldSize, idx_t *xadj, idx_t* adjncy, int numOfNodes, idx_t* partitioning) {
    // Undirected Graph Matrix for METIS

    //readMetisArrays(xadj, adjncy);

    idx_t nvtxs = numOfNodes; //number of vertices in the graph
    idx_t ncon = 1; //number of balancing constraints. It should be at least 1
    idx_t nparts = worldSize; //number of parts to partition the graph

    idx_t metisProgramState;
    idx_t  adjwgt = 1;	

    return METIS_PartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, NULL, NULL, &adjwgt, &nparts, NULL, NULL, NULL, &metisProgramState, partitioning);
}


//forms the adjacency list of undirected and directed versions of the graph represented in text file. It also keeps track of the outdegrees of each node later to normalize the values for each column
void formAdjacancyMatricesAndDictionary(int numOfNodes, set<int>* adjacencyList, set<int>* undirectedAdjacencyList, map<string, int>* nodeIndexes, vector<int>* outdegrees){
    printf("--> Start forming metadata\n");
    ifstream input;
	input.open(FILENAME);

    set<string> setOfAddresses; //collects a set of hash strings corresponding to nodes given in file
    
    string adr;
    int count = (int)linecount * 2; //number of lines
	while (count-- != 0) {
        input >> adr; //hash string for each node in file
        setOfAddresses.insert(adr);
    }

    input.close();

    cout << "adress size: " << numOfNodes << endl; //numOfNodes is actually equal to the size of setOfAddresses. I use a precalculated static value for it throughout the program

    set<string>::iterator itr = setOfAddresses.begin();
    
    //below create an index map. Each unique hash string is mapped to an integer value representation
    for (int i = 0; i < numOfNodes ; i++) {
        nodeIndexes->insert({*itr, i});
        itr++;
    }

    input.open(FILENAME);


    string from;
    string to;
    count = (int)linecount ;
    while (count-- != 0) {
        input >> from;
        input >> to;

        //finding the corresponding integer values for both nodes below
        int fromIndex = nodeIndexes->find(from)->second; 
        int toIndex = nodeIndexes->find(to)->second;

        outdegrees->at(fromIndex) += 1; //outdegree is incremented for a node (column)

        adjacencyList[toIndex].insert(fromIndex); //inserted to directed adjacency list
        //below inserted to undirected adjacency list. METIS partitioning requires this
        if(toIndex != fromIndex){
            undirectedAdjacencyList[toIndex].insert(fromIndex);
            undirectedAdjacencyList[fromIndex].insert(toIndex);
        }
    }
    input.close();

    printf("--> End forming metadata\n");
}


//method is for checking if the interest rates vector for nodes converged to a desired value. Return value from this being true implies the convergence to desired interval has occured
bool complyWithConstraint(vector<double> r_prev, vector<double>* r_init){
    double sum = 0; 
    //below sums differences between r^(i+1) - r^(i)
    for(int i = 0; i < r_init->size(); i++){
        double signedDiff = r_init->at(i) - r_prev.at(i); 
        sum += signedDiff < 0 ? signedDiff * -1 : signedDiff;
    }
    return sum <= epsilon;
}



//
int main(int arg, char* argv[])
{
    MPI_Init(&arg, &argv); //initialize MPI

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); //get the processor rank
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); //get the number of processors

    int timeElapsed = 0;

    if(world_rank == 0){ //if we are in processor zero (master processor)

        CSR* csr; //define csr for the graph represented in file
        map<string, int>* nodeIndexes = new map<string, int>(); //map of hash string of a node and its corresponding integer ID

        const int numOfNodes = 1850066; //number of distinct nodes in given internet graph. Hardcoded, calculated with a function call
        set<idx_t> adjacencyList[numOfNodes]; //adjacency list for the directed internet graph
        set<idx_t> undirectedAdjacencyList[numOfNodes]; //adjacency list for the undirected internet graph
        
        //below to be used when calling partitionWithMETIS function
        idx_t *xadj = new idx_t[numOfNodes];
        vector<idx_t> adjncy;
        
        vector<int>* outdegrees = new vector<int>(numOfNodes); //beware can have elements with zero value
        idx_t *partitioning = new idx_t[numOfNodes]; //partitioning array to be filled by METIS partitioning function call

        formAdjacancyMatricesAndDictionary(numOfNodes, adjacencyList, undirectedAdjacencyList,nodeIndexes ,outdegrees);

        printf("--> Creating xadj and adjncy arrays\n");

        int count = 0;
        xadj[0] = 0;
        for(int i = 0; i < numOfNodes; i++){
            if(i != 0)
                xadj[i] = xadj[i-1] + count;

            count += undirectedAdjacencyList[i].size();
            std::copy(undirectedAdjacencyList[i].begin(), undirectedAdjacencyList[i].end(), std::back_inserter(adjncy));
        }

        printf("--> Done creating xadj and adjncy arrays\n");

        idx_t* prt = new idx_t[adjncy.size()];	
        for(int a = 0; a < adjncy.size(); a++)	
            prt[a] = adjncy.at(a);	

        partitionWithMETIS(world_size - 1, xadj, prt, numOfNodes, partitioning); //partitions of nodes are created to be distributed among processors

        createCSR(csr, numOfNodes, adjacencyList, outdegrees); //csr format is created with the adjacency list and outdegrees found from formAdjacancyMatricesAndDictionary call

        vector<double>* r_init = new vector<double>(numOfNodes, 1.0); //the vector holding rankings of each node

        vector<int>* processorIndexes = new vector<int>[world_size]; //is a hashmap for processors pointing to the vector of the nodes assigned to each corresponding processor. i'th processorIndexes vector
                                                                     //corresponds to i'th processor's assigned nodes to be processed. zeroth vector will be size zero: master node, won't be assigned any.

        //below we collect the nodes that are assigned to the processors to the processor's bucket
        for(int i = 0; i < numOfNodes; i++)
            processorIndexes[partitioning[i] + 1].push_back(i); //all processor indexes assigned to the node by metis partitioner are incremented by one to exlude zero'th index processor

        timeElapsed = MPI_Wtime();
        while(true){
            vector<double> r_prev; //copy of the current r_init. r_init will change below
        
            for(int i = 1; i < world_size; i++){
                
                int numOfAdjs = processorIndexes[i].size(); //send the len of processorIndexes[i] which signifies num of nodes assigned to it
                MPI_Send( &numOfAdjs , 1 ,MPI_INT, i, 3141, MPI_COMM_WORLD); //msg sent to slaves with tag -1 means this is the num of nodes assigned to the slave processor
                for(int j = 0; j < processorIndexes[i].size(); j++){  //i being the rank of the processor, j being the node's index in assigned nodes vector of processor array
                    int nodeIndex = processorIndexes[i].at(j); //node index to be assigned to the processor i
                    int getArrayLength = sizeof(csr->values) / sizeof(int); //number of entries in internet graph matrix
                    int lenOfAdj = nodeIndex == numOfNodes - 1 ? getArrayLength - csr->row_begin.at(nodeIndex) : csr->row_begin.at(nodeIndex + 1) - csr->row_begin.at(nodeIndex); 
                    //lenOfAdj above is the length of the partition of values array of CSR matrix to be sent to processor for the corresponding node j assigned to processor i

                    vector<int> v(csr->values[csr->row_begin.at(nodeIndex)], csr->values[csr->row_begin.at(nodeIndex)] + lenOfAdj); //v is the values subarray corresponding to node j to be sent to processor i
                    v.push_back(processorIndexes[i].at(j)); //here, v's last element is the id of the node j being assigned to processor i (done to reduce number of MPI_Sends)
                    vector<int>::const_iterator first = csr->col_indices.begin() + csr->row_begin.at(nodeIndex);
                    vector<int>::const_iterator last = first + lenOfAdj;
                    vector<int> k(first, last); //sub col_indices vector corresponding to the sub values array v above

                    int size = v.size(); //inconsistent size of vs throughout j s. thus we are sending the size for them before actually sending v and k for j.
                    
                    MPI_Send( &size, 1, MPI_INT, i, 0, MPI_COMM_WORLD); //send size of sub array of values and col_indices corresponding to node j
                    MPI_Send( &v , size , MPI_INT, i, 5, MPI_COMM_WORLD); //send sub array values corresponding to node j
                    MPI_Send( &k , size - 1, MPI_INT, i, 8, MPI_COMM_WORLD); //send sub vector col_indices corresponding to node j
                }

                MPI_Send(&numOfNodes, 1, MPI_INT, i , 3, MPI_COMM_WORLD); //sends size of rank vector
                MPI_Send(r_init, numOfNodes, MPI_DOUBLE, i, 4, MPI_COMM_WORLD); //send the rank vector itself
            }

            MPI_Status status[world_size-1];
            MPI_Request requests[world_size-1];
            for (int i = 1; i < world_size; i++) {
                int recvData; //dummy receive data sent with signal
                MPI_Irecv(&recvData, 1, MPI_INT, i, 478, MPI_COMM_WORLD, &requests[i-1]); //recieve termination signals from processors async
            }

             MPI_Waitall(world_size -1, requests, status); //wait for all async signals signaling the termination of processes in slave processors
            
            std::copy(r_init->begin(), r_init->end(), std::back_inserter(r_prev)); //copy r_init values to r_prev for future calculations of convergence
            if(complyWithConstraint(r_prev, r_init)) break; //checking if convergence has occured
        }

        timeElapsed = MPI_Wtime() - timeElapsed;

        cout << "Time elapsed for the parallelized operations: " << timeElapsed << endl << endl; //printing out the time elapsed

        //below, first five nodes with highest rankings is printed to the screen with their corresponding ranking values
        map<double, int, greater<double>> mp;
        for(int k = 0; k < numOfNodes; k++)
            mp[r_init->at(k)] = k;
        
        int c = 0;
        for (map<double, int>::iterator it = mp.begin(); c != 5; it++, c++)
            cout << it->first << " => " << it->second << endl;

    }

    //MPI_Barrier(MPI_COMM_WORLD);

    //START OF PARALLEL WORK DONE
    if(world_rank != 0){ 
        printf("--> starting to collect data being sent to processor %d", world_rank);
        int numOfNodesAssigned; //number of nodes to be processed
        MPI_Recv(&numOfNodesAssigned, 1, MPI_INT, 0, 3141, MPI_COMM_WORLD ,MPI_STATUS_IGNORE); //gets the number of nodes it will be receiving
        vector<int> adjacencies[numOfNodesAssigned]; //value sub-arrays map for each assigned node
        vector<int> indexOfNodeAssigned[numOfNodesAssigned]; //col_indices sub-vectors map for each assigned node
        for(int i = 0 ; i < numOfNodesAssigned; i++){
            int lenOfAdj;
            MPI_Recv(&lenOfAdj, 1, MPI_INT, 0, 0, MPI_COMM_WORLD ,MPI_STATUS_IGNORE); //receive size of value and col_indices sub-arrays
            MPI_Recv(&(adjacencies[i]), lenOfAdj, MPI_INT, 0, 5, MPI_COMM_WORLD ,MPI_STATUS_IGNORE); //receive the values sub-array for corresponding node
            MPI_Recv(&indexOfNodeAssigned[i], lenOfAdj - 1, MPI_INT, 0, 8, MPI_COMM_WORLD ,MPI_STATUS_IGNORE); //receive the col_indices sub-array for corresponding node
        }

        int lenOfR; //length of the ranks vector of nodes in internet graph
        MPI_Recv(&lenOfR, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //receive the value of lenOfR

        vector<double>* r_i_prev; //ranks vector of nodes in internet graph
        MPI_Recv(r_i_prev, lenOfR, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //receive the ranks vector of nodes in internet graph

        vector<double> r_prev_copy;
        std::copy(r_i_prev->begin(), r_i_prev->end(), std::back_inserter(r_prev_copy));

        //below, matrix multiplication is done for each rank vector index corresponding to the assigned nodes
        for(int i = 0; i < numOfNodesAssigned; i++){
            int rowToBeMultiplied = adjacencies[i].back(); //the index of the node the recalculation is being done. We sent it by appending it to values sub-array
            adjacencies[i].pop_back();
            
            int newRValue = 0;
            //below, the multiplication happens for the node with index rowToBeMultiplied
            for(int j = 0; j < adjacencies[i].size() ; j++){
                int correspondingRIndex = indexOfNodeAssigned[i].at(j);
                int valueFromCSR = adjacencies[i].at(j);
                double toBeReplaced = alpha * valueFromCSR * r_i_prev->at(correspondingRIndex) + (1 - alpha);
                newRValue += toBeReplaced;
            }

            r_i_prev->at(rowToBeMultiplied) = newRValue;
        }

        int gibberish = 9;//dummy data of termination signal below
        MPI_Send((void *) &gibberish, 1, MPI_INT, 0, 478, MPI_COMM_WORLD); //termination signal sent from slaves to master processor

        printf("--> die arbeit hat fertig bleiben für die process %d", world_rank);
    }
    
    //END OF PARALLEL WORK DONE

    //MPI_Barrier(MPI_COMM_WORLD);
    
    return 0;
}
