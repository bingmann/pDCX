/************************************************************************************
 *      Program:    pDC3
 *      Author:     Fabian Kulla
 *      Mail:       @
 *      Description: parallel suffix array construction
 *
 *
 *      Different measurements and output for different DEBUGLEVEL
 *      Debug0: no timemeasuring
 *      Debug1: only total runtime
 *      Debug2: timing local sorting
 *      Debug3: timing communication may include synchronisation
 *      Debug4: timing communication without synchronisation-time (MPI_Barrier used)
 *      Debug5: Lastverteilung
 *      Debug6: not really used
 *
**************************************************************************************/

#define DEBUGLEVEL 1
#define SORT 0

static const unsigned int X	= 7;
static const unsigned int D	= 3;
static const unsigned int DC[3] = { 0, 1, 3 };

static const unsigned int DCH[X] = { 3, 6, 6, 5, 6, 5, 4 };	// additional chars in tuple
static const unsigned int DCD[X] = { 0, 0, 1, 0, 3, 2, 1 };	// depth to sort chars before using first rank

static const unsigned int inDC[X] = { 1, 1, 0, 1, 0, 0, 0 };

static const bool debug = true;

#include <mpi.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <fcntl.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <limits.h>
#include <inttypes.h>

#include "tuple.h"

#include "util.h"

#include <cassert>
#include <cmath>
#include <vector>
#include <limits>

#include "yuta-sais-lite.h"
#include "sachecker.h"

typedef unsigned int	uint;

#define MSGTAG 42

typedef TupleA<X>	TupleX;


#define DBG_ARRAY(dbg,text,X)  do {						\
    if (dbg)									\
    {										\
	std::cout << text << " line " << __LINE__ << " : i - " #X "[i]\n"; \
 	for (unsigned int i = 0; i < X.size(); ++i)				\
	    std::cout << i << " : " << X[i] << "\n";				\
    }										\
} while(0)

#define DBG_ARRAY2(dbg,text,X,Xsize)  do {					\
    if (dbg)									\
    {										\
	std::cout << text << " line " << __LINE__ << " : i - " #X "[i]\n"; \
 	for (unsigned int i = 0; i < (unsigned int)(Xsize); ++i)		\
	    std::cout << i << " : " << X[i] << "\n";				\
    }										\
} while(0)

uint* sortS0S1S2( uint* inbuffer, Pair* P, uint localSize, uint* n, uint* imod3, uint* salen, uint& half);
inline void getTuple(uint* inbuffer, Pair* recvBufPair, Quintuple* S0, Quadruple* S1, Quintuple* S2, uint localSize, uint* n,  uint* imod3, uint& half);

inline int namelex( TupleX* in, Pair* P, uint arraylen );

void writesa( uint* salen, uint* suffixarray, char* filenameOut );
void writeTimes( uint filelength );

int myproc, nprocs;
const int ROOT = 0;

/**
 *
 * @param in
 * @param out       :   out-array.length > in-array.length => some fields are not used
 * @param length    :   length of in-array
 * @param n12all    :   splitter
 * @param half      :   1/2 length of out-array
 */
void sortPair( Pair* in, Pair* out, uint length, uint* n12all, uint half )
{
    if ( myproc == 0 )
	for ( uint i = 0 ; i < length; i++ ) {
	    if ( in[ i ].index >= n12all[ nprocs - 1 ] ) {
		//          assert( (in[ i ].index - n12all[ nprocs -1 ] + half) < 2*half);
		//          std::cout << "in[i].index -n12all[nprocs-1] +half " << in[ i ].index - n12all[ nprocs-1 ] + half << std::endl;
		out[ in[ i ].index - n12all[ nprocs - 1 ] + half ] = in[ i ];
	    } else {
		//          std::cout << "in[i].index " << in[ i ].index <<" half "<<2*half<< std::endl;
		//          if(in[ i ].index > 2*half) assert(0);
		//          if ((int)in[ i ].index < 0 ) assert(0);
		out[ in[ i ].index ] = in[ i ];
	    }
	}
    else
	for ( uint i = 0 ; i < length; i++ ) {
	    if ( in[ i ].index >= n12all[ nprocs - 1 ] ) {
		//          assert( (in[ i ].index - n12all[myproc+ nprocs -1 ] + half )< 2*half);
		//          std::cout << "in[i].index -n12all[myproc-1+nprocs] +half " << in[ i ].index - n12all[ myproc -1 + nprocs ] + half << std::endl;
		out[ in[ i ].index - n12all[ myproc - 1 + nprocs ] + half ] = in[ i ];
	    } else {
		//          std::cout << "in[i].index -n12all[myproc-1] " << in[ i ].index - n12all[ myproc-1 ] << std::endl;
		out[ in[ i ].index - n12all[ myproc - 1 ] ] = in[ i ];
	    }
	}
}

MPI_Datatype MPI_TUPLEX;
MPI_Datatype MPI_PAIR;
MPI_Datatype MPI_TUPLE_SX;


//MPI_Datatype MPI_QUADRUPLE;

MPI_Status status;

int samplefactor;

/**used for timing*/
Debug1(
    double dc3StartTime;  // complete time
    double dc3FinishTime=0;
)

Debug2(
    double starttime;
    double sortingTime = 0;     // time for stl-sort
    double mergeTime = 0;       // time for my mergesort
    double alltoallvTime = 0;   // time for MPI_Alltoallv
    double permuteTime=0;       // time for permuting Pairs
    double mpiComStartTime;
    double mpiComTime=0;            // time for communication without MPI_Alltoallv
    double sampleStartTime;
    double sampleTime=0;        // time used for sampling (includes some mpiComTime-time)

    double sortS0S1S2TimeStart=0;
    double sortS0S1S2TimeEnd= 0; // time for Step 3 of DC3 (includes other times)

    const int maxRekursion=35;  // no chance that there will be more than 35 rekursions (worst case 4G * a)
    double rekursionStart[maxRekursion];
    double rekursionEnd[maxRekursion];
    int countRek=0;

    uint counter=0;
    double name2tuple[maxRekursion];
)

Debug6(double stlSortTime;)
/** end timing */

void mpi_init_datatypes()
{
    //Zu den Tupeln passenden MPI-Datentyp
    MPI_Datatype mpi_tupel[1] = { MPI_UNSIGNED };

    int mpi_blocklen[1] = { X+1 };
    MPI_Aint mpi_displ[1] = { 0 };
    MPI_Type_struct( 1, mpi_blocklen, mpi_displ, mpi_tupel, &MPI_TUPLEX );
    MPI_Type_commit( &MPI_TUPLEX );

    mpi_blocklen[ 0 ] = 2;
    MPI_Type_struct( 1, mpi_blocklen, mpi_displ, mpi_tupel, &MPI_PAIR );
    MPI_Type_commit( &MPI_PAIR );

    mpi_blocklen[ 0 ] = (X-1) + D + 1;
    MPI_Type_struct( 1, mpi_blocklen, mpi_displ, mpi_tupel, &MPI_TUPLE_SX );
    MPI_Type_commit( &MPI_TUPLE_SX );

#if 0
    //Zu den Tupeln passenden MPI-Datentyp
    MPI_Datatype mpi_tupel[ 2 ] = {MPI_UNSIGNED, MPI_UNSIGNED};

    int mpi_blocklen[ 2 ] = {4, 1};
    MPI_Aint mpi_displ[ 2 ] = {0, 4 * sizeof( uint ) };
    MPI_Type_struct( 2, mpi_blocklen, mpi_displ, mpi_tupel, &MPI_QUINTUPLE );
    MPI_Type_commit( &MPI_QUINTUPLE );

    mpi_blocklen[ 0 ] = 3;
    mpi_displ[ 1 ] = 3 * sizeof( uint );
    MPI_Type_struct( 2, mpi_blocklen, mpi_displ, mpi_tupel, &MPI_QUADRUPLE );
    MPI_Type_commit( &MPI_QUADRUPLE );

#endif
}

template < typename Comparator >
class LoserTree
{
private:

    /// the tree of size n-1
    std::vector<int>	m_tree;

    /// the comparator object of this tree
    const Comparator&	m_less;

public:

    LoserTree(unsigned int size, const Comparator& less)
	: m_less(less)
    {
	// initialize loser tree by performing comparisons

	unsigned int treesize = (1 << (int)(log2(size - 1) + 2)) - 1;
	m_tree.resize(treesize, -1);

	// fill in lowest level: ascending numbers until each sequence
	// finishes.
	int levelput = m_tree.size() / 2;
	for (unsigned int i = 0; i < size; ++i)
	    m_tree[levelput + i] = i;

	int levelsize = levelput + 1;

	// construct higher levels iteratively from bottom up
	while ( levelsize > 1 )
	{
	    levelsize = levelsize / 2;
	    int levelget = levelput;
	    levelput /= 2;

	    for (int i = 0; i < levelsize; ++i)
	    {
		if ( m_tree[levelget + 2*i + 1] < 0 )
		    m_tree[levelput + i] = m_tree[levelget + 2*i];
		else if (m_less( m_tree[levelget + 2*i], m_tree[levelget + 2*i + 1] ))
		    m_tree[levelput + i] = m_tree[levelget + 2*i];
		else
		    m_tree[levelput + i] = m_tree[levelget + 2*i + 1];
	    }
	}
    }

    int top() const
    {
	return m_tree[0];
    }

    void replay(bool done)
    {
	int top = m_tree[0];

	int p = (m_tree.size() / 2) + top;
	
	if (done) top = -1;	// mark sequence as done

	while( p > 0 )
	{
	    m_tree[p] = top;

	    p -= (p+1) % 2;	// round down to left node position
    
	    if (m_tree[p] < 0)
		top = m_tree[p+1];
	    else if (m_tree[p+1] < 0)
		top = m_tree[p];
	    else if ( m_less( m_tree[p], m_tree[p+1] ) )
		top = m_tree[p];
	    else
		top = m_tree[p+1];

	    p /= 2;
	}

	m_tree[p] = top;
    }

    void print() const
    {
	int levelsize = 1;
	int j = 0;

	for (unsigned int i = 0; i < m_tree.size(); ++i)
	{
	    if (i >= j + levelsize) {
		std::cout << "\n";
		j = i; levelsize *= 2;
	    }
	    std::cout << m_tree[i] << " ";
	}
	std::cout << "\n";
    }
};

static const int cmpDepth[X][X] = 
{
    {  0,  0,  1,  0,  3,  3,  1 },
    { -1,  0,  6,  0,  6,  2,  2 },
    { -1, -1,  1,  5,  6,  5,  1 },
    { -1, -1, -1,  0,  4,  5,  4 },
    { -1, -1, -1, -1,  3,  3,  4 },
    { -1, -1, -1, -1, -1,  2,  2 },
    { -1, -1, -1, -1, -1, -1,  1 },
};

static const int cmpRanks[X][X][2] = 
{
    { {  0, 0 }, {  0, 0 }, {  1, 0 }, {  0, 0 }, {  2, 0 }, {  2, 1 }, {  1, 0 } },
    { { -1,-1 }, {  0, 0 }, {  2, 2 }, {  0, 0 }, {  2, 2 }, {  1, 0 }, {  1, 1 } },
    { { -1,-1 }, { -1,-1 }, {  0, 0 }, {  1, 2 }, {  2, 2 }, {  1, 2 }, {  0, 0 } },
    { { -1,-1 }, { -1,-1 }, { -1,-1 }, {  0, 0 }, {  1, 1 }, {  2, 2 }, {  1, 2 } },
    { { -1,-1 }, { -1,-1 }, { -1,-1 }, { -1,-1 }, {  0, 0 }, {  0, 1 }, {  1, 2 } },
    { { -1,-1 }, { -1,-1 }, { -1,-1 }, { -1,-1 }, { -1,-1 }, {  0, 0 }, {  0, 1 } },
    { { -1,-1 }, { -1,-1 }, { -1,-1 }, { -1,-1 }, { -1,-1 }, { -1,-1 }, {  0, 0 } },
};

bool cmpTupleSXCompare(const TupleSX& t1, const TupleSX& t2)
{
    unsigned int v1 = t1.index % X, v2 = t2.index % X;

    if (v1 <= v2)
    {
	int depth = cmpDepth[v1][v2];
	assert(depth >= 0);

	std::cout << "cmp " << v1 << "(" << t1.index << ") ? " << v2 << "(" << t2.index << ") - depth " << depth << "\n";

	for (int d = 0; d < depth; ++d)
	{
	    if (t1.chars[d] == t2.chars[d]) continue;
	    return (t1.chars[d] < t2.chars[d]);
	}

	const int* rs = cmpRanks[v1][v2];

	std::cout << "break tie using ranks " << rs[0] << " - " << rs[1] << " = " << t1.ranks[ rs[0] ] << " - " << t2.ranks[ rs[1] ] << "\n";

	return (t1.ranks[ rs[0] ] < t2.ranks[ rs[1] ]);
    }
    else // (v1 > v2)
    {
	int depth = cmpDepth[v2][v1];
	assert(depth >= 0);

	std::cout << "cmp " << v1 << "(" << t1.index << ") ? " << v2 << "(" << t2.index << ") - depth " << depth << "\n";

	for (int d = 0; d < depth; ++d)
	{
	    if (t1.chars[d] == t2.chars[d]) continue;
	    return (t1.chars[d] < t2.chars[d]);
	}

	const int* rs = cmpRanks[v2][v1];

	std::cout << "break tie using ranks " << rs[1] << " - " << rs[0] << " = " << t1.ranks[ rs[1] ] << " - " << t2.ranks[ rs[0] ] << "\n";

	return (t1.ranks[ rs[1] ] < t2.ranks[ rs[0] ]);
    }
}

template <int X>
struct TupleSXMerge
{
    const std::vector<TupleSX>*	 m_S;
    
    std::vector<unsigned int>	 m_ptr;

    TupleSXMerge(const std::vector<TupleSX>* S)
	: m_S(S), m_ptr(X, 0)
    {
    }

    inline bool operator()(int v1, int v2) const
    {
	assert( v1 < v2 );

	int depth = cmpDepth[v1][v2];
	assert(depth >= 0);

	const TupleSX& t1 = m_S[v1][ m_ptr[v1] ];
	const TupleSX& t2 = m_S[v2][ m_ptr[v2] ];
	
	assert( t1.index % X == v1 );
	assert( t2.index % X == v2 );

	std::cout << "cmp " << v1 << "(" << t1.index << ") ? " << v2 << "(" << t2.index << ") - depth " << depth << "\n";

	for (int d = 0; d < depth; ++d)
	{
	    if (t1.chars[d] == t2.chars[d]) continue;
	    return (t1.chars[d] < t2.chars[d]);
	}

	const int* rs = cmpRanks[v1][v2];

	{ // self-check matrices

	    int d;
	    uint p1 = v1, p2 = v2;
	    int r1 = 0, r2 = 0;

	    for (d = 0; d < X; ++d)
	    {
		if (inDC[ p1 % X ] && inDC[ p2 % X ])
		    break;

		if (inDC[ p1 % X ]) r1++;
		if (inDC[ p2 % X ]) r2++;

		assert( t1.chars[d] == t2.chars[d] );
		p1++, p2++;
	    }

	    assert( d == depth );
	    assert( r1 == rs[0] );
	    assert( r2 == rs[1] );
	}
	
	std::cout << "break tie using ranks " << rs[0] << " - " << rs[1] << " = " << t1.ranks[ rs[0] ] << " - " << t2.ranks[ rs[1] ] << "\n";

	return (t1.ranks[ rs[0] ] < t2.ranks[ rs[1] ]);
    }
};

/**
 *
 * @param inbuffer
 * @param filelength
 * @param salen
 * @return suffixarray
 */
uint* dc3( uint* input, uint globalSize, uint* salen, uint localSizeGeneral )
{
    if ( myproc == ROOT )
	std::cout << "******************** DCX ********************" << std::endl;


    //uint localSize = ( globalSize + nprocs - 1 ) / nprocs ;

    uint samplesize = (uint)sqrt(2.0 * localSizeGeneral / 3.0 / nprocs) * samplefactor;

// Debug5(  if ( myproc == ROOT ) std::cout << "samplefactor=" << k <<" |Tuple|="<<2*(localSizeGeneral/3)<<" in %="<<100*k/(2*localSizeGeneral/3.0)<<std::endl;)
//  int k=(uint) ceil(samplefactor*(sqrt(2.0*localSizeGeneral/3.0)/100));

    if( samplesize >= 2*(localSizeGeneral/3)) samplesize = 2*(localSizeGeneral/3)-1;

// if ( myproc == ROOT ) std::cout << "samplefactor=" << k <<" |Tuple|="<<2*(localSizeGeneral/3)<<" in %="<<100*k/(2*localSizeGeneral/3.0)<<std::endl;

    unsigned int globalMultipleOfX = (globalSize + X - 1) / X;	// rounded up number of global multiples of X

    unsigned int localOffset = myproc * localSizeGeneral;
    unsigned int localSize = (myproc != nprocs-1) ? localSizeGeneral : globalSize - localOffset;
    unsigned int localSizeReal = (myproc != nprocs-1) ? localSizeGeneral + (X-1) : globalSize - localOffset;

    printf("localSize = %d\n", localSize);

    // **********************************************************************
    // * calculate build DC-tuple array and sort locally

    const unsigned int M = localSize / X;		// multiples of complete X chars in local area size
    const unsigned int M2 = (localSize + X - 1) / X;	// number of incomplete X chars in local area size

    std::cout << "localSize = " << localSize << " - M = " << M << " - M2 = " << M2 << "\n";

    std::vector<TupleX> R (D * M2);			// create D*M2 tuples which might include up to D-1 dummies

    uint j = 0 ;
    for (uint i = 0; i < localSize; i += X)
    {
	for (uint d = 0; d < D; ++d)
	{
	    R[j].index = localOffset + i + DC[d];

	    for (uint x = i + DC[d], y = 0; y < X; ++x, ++y)
		R[j].chars[y] = (x < localSizeReal) ? input[x] : 0;

	    ++j;
	}
    }

    assert( j == D * M2 );

    // **********************************************************************
    // {{{ Sample sort of array R
    {	
	std::sort(R.begin(), R.end());		// sort locally

	DBG_ARRAY(1, "Sample suffixes", R);

	// **********************************************************************
	// * select equidistance samples and redistribute sorted DC-tuples

	// select samples
	TupleX* samplebuf = new TupleX[ samplesize ];

	double dist = ( double ) R.size() / samplesize;
	for ( uint i = 0; i < samplesize; i++ )
	    samplebuf[i] = R[ int( i * dist ) ];

	// root proc collects all samples
	TupleX* samplebufall;
	if (myproc == ROOT) samplebufall = new TupleX[ nprocs * samplesize ];

	MPI_Gather( samplebuf, samplesize, MPI_TUPLEX, samplebufall, samplesize, MPI_TUPLEX, ROOT, MPI_COMM_WORLD );

	delete[] samplebuf;

	// root proc sorts samples as splitters

	TupleX* splitterbuf = new TupleX[ nprocs ];

	if ( myproc == ROOT )
	{
	    std::sort( samplebufall, samplebufall + samplesize * nprocs );

	    DBG_ARRAY2(1, "Sample splitters", samplebufall, nprocs * samplesize);

	    for ( int i = 0; i < nprocs; i++ )		// pick splitters
		splitterbuf[i] = samplebufall[ i * samplesize ];

	    DBG_ARRAY2(1, "Selected splitters", splitterbuf, nprocs);

	    delete[] samplebufall;
	}

	// distribute splitters

	MPI_Bcast( splitterbuf, nprocs, MPI_TUPLEX, ROOT, MPI_COMM_WORLD );

	// find nearest splitters in locally sorted tuple list

	uint* splitterpos = new uint[ nprocs + 1 ];

	splitterpos[0] = 0;
	for ( int i = 1; i < nprocs; i++ )
	{
	    //splitterpos[i] = findPos( R.data(), splitterbuf[i].data(), R.size(), cmpSplitterGreaterS12 );

	    std::vector<TupleX>::const_iterator it = std::lower_bound(R.begin(), R.end(), splitterbuf[i]);

	    splitterpos[i] = it - R.begin();
	}
	splitterpos[ nprocs ] = R.size();

	DBG_ARRAY2(1, "Splitters positions", splitterpos, nprocs+1);

	delete [] splitterbuf;

	// boardcast number of element in each division

	int* sendcnt = new int[ nprocs ];
	int* recvcnt = new int[ nprocs ];

	for ( int i = 0; i < nprocs; i++ )
	{
	    sendcnt[i] = splitterpos[i+1] - splitterpos[i] ;
	    assert( sendcnt[i] >= 0 );
	}

	delete[] splitterpos;

	MPI_Alltoall( sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT , MPI_COMM_WORLD );

	// calculate number of received tuples

	int* sendoff = new int[ nprocs + 1 ];
	int* recvoff = new int[ nprocs + 1 ];

	sendoff[ 0 ] = recvoff[ 0 ] = 0;
	for ( int i = 1; i < nprocs + 1; i++ )
	{
	    sendoff[i] = sendoff[i-1] + sendcnt[i-1];
	    recvoff[i] = recvoff[i-1] + recvcnt[i-1];
	}

	std::vector<TupleX> recvbuf ( recvoff[ nprocs ] );

	MPI_Alltoallv( R.data(), sendcnt, sendoff, MPI_TUPLEX, recvbuf.data(), recvcnt, recvoff, MPI_TUPLEX, MPI_COMM_WORLD );

	std::swap(R, recvbuf);

	delete [] sendcnt;
	delete [] recvcnt;
	delete [] sendoff;

	if ( R.size() )
	{
	    /** merge */
	    TupleX * helparray4 = new TupleX[ R.size() ];

	    // binary merge sort!
	    mergesort( R.data(), helparray4, 0, nprocs, recvoff, cmpTupleXLeq<X> );

	    // really no need for all this sorting: just name using loser tree!

	    delete[] helparray4;
	}
	else {
	    std::cout << myproc << " hat nichts bekommen" << std::endl;
	}
    }
    // }}}

    // recvbuf contain DC-tuples in unsorted order

    DBG_ARRAY(1, "Sorted sample suffixes", R);

    /** Lexicographical naming   */
    std::vector<Pair> P ( R.size() );

    int recursion;

    {
	// **********************************************************************
	// *** naming with local names

	TupleX temp;	// get last tuple from previous process as basis for name comparison (cyclicly)
	MPI_Sendrecv( &(R[R.size()-1]), 1, MPI_TUPLEX, (myproc + 1) % nprocs, MSGTAG,
		      &temp, 1, MPI_TUPLEX, (myproc - 1 + nprocs) % nprocs, MSGTAG,
		      MPI_COMM_WORLD, &status );

	uint name = 0;
	for ( uint i = 0; i < R.size(); i++ ) {
	    if ( !( R[i] == temp ) ) {
		name++;
		std::cout << "Giving name " << name << " to " << R[i] << "\n";
		temp = R[i];
	    }
	    P[i].name = name;
	    P[i].index = R[i].index;
	}
	R.clear();		// Why?: because it is easier to recreate the tuples later on

	DBG_ARRAY(1, "Local Names", P);

	// **********************************************************************
	// *** renaming with global names: calculate using prefix sum
	
	uint namesglob;
	MPI_Scan( &name, &namesglob, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );

	// update local names
	for ( uint i = 0; i < P.size(); i++ )
	    P[i].name += (namesglob - name);

	DBG_ARRAY(1, "Global Names", P);

	// determine whether recursion is necessary: last proc checks its highest name
	if (myproc == nprocs - 1)
	{
	    std::cout << "last name: " << P.back().name << " -? " << D * globalMultipleOfX << "\n";
	    recursion = (P.back().name != D * globalMultipleOfX);
	    std::cout << "recursion: " << recursion << "\n";
	}

	MPI_Bcast( &recursion, 1, MPI_INT, nprocs - 1, MPI_COMM_WORLD );
    }

    if ( recursion )
    {
	if ( myproc == ROOT )
	    std::cout << "---------------------   RECURSION ---------------- " << localSize << std::endl;

	uint namesGlobalSize = D * globalMultipleOfX;
	uint namesLocalSize = ( namesGlobalSize + nprocs - 1 ) / nprocs;	// rounded up 

	if (namesGlobalSize > 1000)
	{
	    // **********************************************************************
	    // {{{ Sample sort of array P by (i mod X, i div X)

	    std::sort(P.begin(), P.end(), Pair::cmpIndexModDiv<X>);		// sort locally

	    DBG_ARRAY(1, "Global Names sorted cmpModDiv", P);

	    uint* splitterpos = new uint[nprocs+1];
	    int* sendcnt = new int[nprocs];
	    int* sendoff = new int[nprocs+1];
	    int* recvcnt = new int[nprocs];
	    int* recvoff = new int[nprocs+1];

	    // use equidistance splitters from 0..namesGlobalSize (because indexes are known in advance)
	    splitterpos[ 0 ] = 0;
	    Pair ptemp;
	    ptemp.name=0;
	    for ( int i = 1; i < nprocs; i++ ) {
		ptemp.index = i * namesLocalSize;	// TODO: check range (maybe index doesnt start at 0)?

		unsigned int x = ptemp.index;

		unsigned int divM = ptemp.index / globalMultipleOfX;
		ptemp.index = DC[divM] + X * (ptemp.index - divM * globalMultipleOfX);
		std::cout << "splitter: " << ptemp.index << " = " << x << " - " << divM << "\n";

		std::vector<Pair>::const_iterator it = std::lower_bound(P.begin(), P.end(), ptemp, Pair::cmpIndexModDiv<X>);
		splitterpos[i] = it - P.begin();
	    }
	    splitterpos[ nprocs ] = P.size();

	    for ( int i = 0; i < nprocs; i++ )
	    {
		sendcnt[ i ] = splitterpos[ i + 1 ] - splitterpos[ i ];
		assert( sendcnt[ i ] >= 0 );
	    }

	    MPI_Alltoall( sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT , MPI_COMM_WORLD );
	
	    sendoff[0] = recvoff[0] = 0;
	    for ( int i = 1; i <= nprocs; i++ ) {
		sendoff[i] = sendoff[i - 1] + sendcnt[i - 1];
		recvoff[i] = recvoff[i - 1] + recvcnt[i - 1];
	    }
	    recvoff[nprocs] = recvoff[nprocs - 1] + recvcnt[nprocs - 1];

	    std::vector<Pair> recvBufPair ( recvoff[ nprocs ] );
	    unsigned int recvBufPairSize = recvoff[ nprocs ];

	    MPI_Alltoallv( P.data(), sendcnt, sendoff, MPI_PAIR, recvBufPair.data(), recvcnt, recvoff, MPI_PAIR, MPI_COMM_WORLD );

	    //delete[] P;
	    P.clear();

	    Pair* helparray2 = new Pair[ recvBufPairSize ];
	    mergesort( recvBufPair.data(), helparray2, 0, nprocs, recvoff, Pair::cmpIndexModDiv<X> );
	    delete[] helparray2;

	    // TODO: merge and reduce at once

	    uint* namearray = new uint[ recvBufPairSize ];
	    for (unsigned int i = 0; i < recvBufPairSize; ++i)
		namearray[i] = recvBufPair[i].index;
	
	    DBG_ARRAY2(1, "Pairs P (globally sorted by indexModDiv)", recvBufPair.data(), recvBufPairSize);

	    // }}} end Sample sort of array P

	    uint rSAsize = 0;
	    uint* rSA = dc3( namearray, namesGlobalSize, &rSAsize, namesLocalSize );

	    delete [] namearray;

	    DBG_ARRAY2(1, "Recursive SA", rSA, rSAsize);
	}
	else
	{
	    int Psize = P.size();

	    int* recvcnt = new int[nprocs];
	    int* recvoff = new int[nprocs+1];

	    MPI_Gather( &Psize, 1, MPI_INT, recvcnt, 1, MPI_INT, ROOT, MPI_COMM_WORLD );

	    if (myproc == ROOT)
	    {
		recvoff[0] = 0;
		for ( int i = 1; i <= nprocs; i++ ) {
		    recvoff[i] = recvoff[i-1] + recvcnt[i-1];
		}
		assert( recvoff[nprocs] == namesGlobalSize );
	    }

	    std::vector<Pair> Pall;

	    if (myproc == ROOT)
		Pall.resize( namesGlobalSize );

	    MPI_Gatherv(P.data(), P.size(), MPI_PAIR,
			Pall.data(), recvcnt, recvoff, MPI_PAIR, ROOT, MPI_COMM_WORLD);

	    if (myproc == ROOT)
	    {
		uint maxname = Pall.back().name;

		std::sort(Pall.begin(), Pall.end(), Pair::cmpIndexModDiv<X>);		// sort locally

		DBG_ARRAY(1, "Global Names sorted cmpModDiv", Pall);

		uint* namearray = new uint[ Pall.size() ];
		for (unsigned int i = 0; i < Pall.size(); ++i)
		    namearray[i] = Pall[i].name;

		DBG_ARRAY2(1, "Recursion input", namearray, Pall.size());

		int* rSA = new int [ Pall.size() ];

		yuta_sais_lite::saisxx< uint*, int*, int >( namearray, rSA, Pall.size(), maxname+1 );

		delete [] namearray;

		DBG_ARRAY2(1, "Recursive SA", rSA, Pall.size());

		// generate rank array - same as updating pair array with correct names
		for (unsigned int i = 0; i < Pall.size(); ++i)
		{
		    Pall[ rSA[i] ].name = i+1;
		}
		
		DBG_ARRAY(1, "Fixed Global Names sorted cmpModDiv", Pall);

		std::swap(P, Pall);
	    }

	    // **********************************************************************
	    // *** redistribute pairs partitioned by localSize
	    
	    if (myproc == ROOT)
	    {
		std::sort(P.begin(), P.end());			// sort locally by index

		DBG_ARRAY(1, "Fixed Global Names sorted index", P);
	    }


	    if (myproc == ROOT)
	    {
		uint* splitterpos = new uint[nprocs+1];

		// use equidistance splitters from 0..globalSize (because indexes are fixed)
		splitterpos[ 0 ] = 0;
		Pair ptemp; ptemp.name=0;
		for ( int i = 1; i < nprocs; i++ ) {
		    ptemp.index = i * localSizeGeneral;
	    
		    std::vector<Pair>::const_iterator it = std::lower_bound(P.begin(), P.end(), ptemp);
		    splitterpos[i] = it - P.begin();
		}
		splitterpos[ nprocs ] = P.size();

		recvcnt[0] = splitterpos[1];
		recvoff[0] = 0;
		for ( int i = 1; i < nprocs; i++ )
		{
		    recvcnt[i] = splitterpos[i+1] - splitterpos[i];
		    recvoff[i] = recvoff[i-1] + recvcnt[i-1];
		    assert( recvcnt[i] >= 0 );
		}
		recvoff[nprocs] = recvoff[nprocs-1] + recvcnt[nprocs-1];

		DBG_ARRAY2(1, "recvcnt", recvcnt, nprocs);
		DBG_ARRAY2(1, "recvoff", recvoff, nprocs+1);
	    }

	    MPI_Bcast( recvcnt, nprocs, MPI_INT, ROOT, MPI_COMM_WORLD );

	    std::vector<Pair> recvBufPair ( recvcnt[ myproc ] + D );
	    unsigned int recvBufPairSize = recvcnt[ myproc ];

	    MPI_Scatterv( P.data(), recvcnt, recvoff, MPI_PAIR, recvBufPair.data(), recvBufPairSize, MPI_PAIR, ROOT, MPI_COMM_WORLD );

	    //delete[] P;
	    P.clear();

	    DBG_ARRAY2(1, "Pairs P (globally sorted by index)", recvBufPair.data(), recvBufPairSize);

	    // **********************************************************************
	    // *** every PE needs D additional sample suffix ranks to complete the final tuple

	    Pair temp[ D ];

	    MPI_Sendrecv( recvBufPair.data(), D, MPI_PAIR, ( myproc - 1 + nprocs ) % nprocs, MSGTAG,
			  &temp, D, MPI_PAIR, ( myproc + 1 ) % nprocs, MSGTAG,
			  MPI_COMM_WORLD, &status );

	    if ( myproc == nprocs - 1 )	// last processor gets sentinel tuples with maximum ranks
	    {
		for (unsigned int i = 0; i < D; ++i)
		{
		    recvBufPair[ recvBufPairSize + i ].name = INT_MAX;
		    recvBufPair[ recvBufPairSize + i ].index = recvBufPair[ recvBufPairSize - D ].index + X + DC[i];
		}
	    }
	    else	// other processors get D following tuples with indexes from the DC
	    {
		for (unsigned int i = 0; i < D; ++i)
		{
		    recvBufPair[ recvBufPairSize + i ].name = temp[i].name;
		    recvBufPair[ recvBufPairSize + i ].index = recvBufPair[ recvBufPairSize - D ].index + X + DC[i];
		    assert( recvBufPair[ recvBufPairSize + i ].index == temp[i].index );
		}
	    }

	    //std::cout<<"CPU ["<<myproc<<"] nach rekursion"<<std::endl;

	    //std::sort( recvBufPair.begin(), recvBufPair.end(), cmpIndexModDiv );

	    std::swap(recvBufPair, P);

	    DBG_ARRAY(1, "Pairs P (globally sorted by index + extra tuples)", P);

	    delete[] recvcnt;
	    delete[] recvoff;

	} // end use sequential suffix sorter

	if ( myproc == ROOT )
	    std::cout << "---------------------   END  RECURSION  --------------- " << std::endl;
    }
    else {
	if ( myproc == ROOT )
	    std::cout << "---------------------   keine  Recursion---------------- " << localSize << std::endl;

	// **********************************************************************
	// *** sample sort pairs P by index

	std::sort( P.begin(), P.end() );

	DBG_ARRAY(1, "Pairs P (sorted by index)", P);

	uint* splitterpos = new uint[nprocs+1];
	int* sendcnt = new int[nprocs];
	int* sendoff = new int[nprocs+1];
	int* recvcnt = new int[nprocs];
	int* recvoff = new int[nprocs+1];

	// use equidistance splitters from 0..globalSize (because names are unique)
	splitterpos[ 0 ] = 0;
	Pair ptemp;
	ptemp.name=0;
	for ( int i = 1; i < nprocs; i++ ) {
	    ptemp.index = i * localSizeGeneral;	// TODO: check range (maybe index doesnt start at 0)?
	    
	    std::vector<Pair>::const_iterator it = std::lower_bound(P.begin(), P.end(), ptemp);
	    splitterpos[i] = it - P.begin();
	}
	splitterpos[ nprocs ] = P.size();

	for ( int i = 0; i < nprocs; i++ )
	{
	    sendcnt[ i ] = splitterpos[ i + 1 ] - splitterpos[ i ];
	    assert( sendcnt[ i ] >= 0 );
	}

	MPI_Alltoall( sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT , MPI_COMM_WORLD );

	sendoff[ 0 ] = recvoff[ 0 ] = 0;
	for ( int i = 1; i <= nprocs; i++ ) {
	    sendoff[ i ] = sendoff[ i - 1 ] + sendcnt[ i - 1 ];
	    recvoff[ i ] = recvoff[ i - 1 ] + recvcnt[ i - 1 ];
	}
	recvoff[ nprocs ] = recvoff[ nprocs - 1 ] + recvcnt[ nprocs - 1 ];

	std::vector<Pair> recvBufPair ( recvoff[ nprocs ] + D );
	unsigned int recvBufPairSize = recvoff[ nprocs ];

	MPI_Alltoallv( P.data(), sendcnt, sendoff, MPI_PAIR, recvBufPair.data(), recvcnt, recvoff, MPI_PAIR, MPI_COMM_WORLD );

	//delete[] P;
	P.clear();

	Pair* helparray2 = new Pair[ recvBufPairSize ];
	mergesort( recvBufPair.data(), helparray2, 0, nprocs, recvoff, cmpLessIndex );
	delete[] helparray2;

	
	DBG_ARRAY2(1, "Pairs P (globally sorted by index)", recvBufPair.data(), recvBufPairSize);

	/* unknown reason?
	if ( myproc == nprocs - 1 ) {
	    if ( recvoff[ nprocs ] == 1 ) {       				// kann nur bei letztem PE auftreten?
		if ( recvBufPair[ recvoff[ nprocs ] - 1 ].index % 3 == 1 )
		    recvBufPair[ recvoff[ nprocs ] ].index = recvBufPair[ recvoff[ nprocs ] - 1 ].index + 1;
		else
		    recvBufPair[ recvoff[ nprocs ] ].index = recvBufPair[ recvoff[ nprocs ] - 1 ].index + 2;
	    }
	}
	*/

	// **********************************************************************
	// *** every PE needs D additional sample suffix ranks to complete the final tuple

	Pair temp[ D ];

	MPI_Sendrecv( recvBufPair.data(), D, MPI_PAIR, ( myproc - 1 + nprocs ) % nprocs, MSGTAG,
		      &temp, D, MPI_PAIR, ( myproc + 1 ) % nprocs, MSGTAG,
		      MPI_COMM_WORLD, &status );

	if ( myproc == nprocs - 1 )	// last processor gets sentinel tuples with maximum ranks
	{
	    for (unsigned int i = 0; i < D; ++i)
	    {
		recvBufPair[ recvBufPairSize + i ].name = INT_MAX;
		recvBufPair[ recvBufPairSize + i ].index = recvBufPair[ recvBufPairSize - D ].index + X + DC[i];
	    }
	}
	else	// other processors get D following tuples with indexes from the DC
	{
	    for (unsigned int i = 0; i < D; ++i)
	    {
		recvBufPair[ recvBufPairSize + i ].name = temp[i].name;
		recvBufPair[ recvBufPairSize + i ].index = recvBufPair[ recvBufPairSize - D ].index + X + DC[i];
		assert( recvBufPair[ recvBufPairSize + i ].index == temp[i].index );
	    }
	}

        //std::cout<<"CPU ["<<myproc<<"] nach rekursion"<<std::endl;

	//std::sort( recvBufPair.begin(), recvBufPair.end(), cmpIndexModDiv );

	std::swap(recvBufPair, P);

	DBG_ARRAY(1, "Pairs P (globally sorted by index + extra tuples)", P);

	delete[] sendcnt;
	delete[] sendoff;
	delete[] recvcnt;
	delete[] recvoff;
    }

    // P contains pairs of index and global rank: sorted by index and partitioned by localSize

    // **********************************************************************
    // *** Generate tuple arrays of samples and non-samples

    std::vector<TupleSX> S [X];

    for (unsigned int k = 0; k < X; ++k)
	S[k].resize(M2);

    unsigned int dp = 0;	// running index into P incremented when crossing DC-indexes

    for (unsigned int i = 0; i < M2; ++i)
    {
	for (unsigned int k = 0; k < X; ++k)
	{
	    S[k][i].index = myproc * localSizeGeneral + i * X + k;

	    for (unsigned int c = 0; c < X-1; ++c)
		S[k][i].chars[c] = (i*X + k + c < localSizeReal) ? input[ i*X + k + c ] : 0;

	    for (unsigned int d = 0; d < D; ++d)
		S[k][i].ranks[d] = P[ dp + d ].name;

	    if (DC[dp % D] == k) ++dp;
	}
    }

    for (unsigned int k = 0; k < X; ++k)
    {
	//DBG_ARRAY(1, "S" << k, S[k]);
    }
    
    // **********************************************************************
    // *** Sample sort tuple arrays
    {
	for (unsigned int k = 0; k < X; ++k)
	{
	    std::sort(S[k].begin(), S[k].end(), cmpTupleSXdepth<6>);	// TODO: sort less
	}

	// select equidistant samples
    
	TupleSX* samplebuf = new TupleSX[ X * samplesize ];

	double dist = ( double ) M2 / samplesize;
	for ( uint i = 0, j = 0; i < samplesize; i++ )
	{
	    for (unsigned int k = 0; k < X; ++k)
	    {
		samplebuf[j++] = S[k][ int( i * dist ) ];
	    }
	}
    
	// root proc collects all samples
	TupleSX* samplebufall;
	if (myproc == ROOT) samplebufall = new TupleSX[ nprocs * X * samplesize ];

	MPI_Gather( samplebuf, X * samplesize, MPI_TUPLE_SX, samplebufall, X * samplesize, MPI_TUPLE_SX, ROOT, MPI_COMM_WORLD );

	delete[] samplebuf;
	
	// root proc sorts samples as splitters

	TupleSX* splitterbuf = new TupleSX[ nprocs ];

	if ( myproc == ROOT )
	{
	    std::sort( samplebufall, samplebufall + nprocs * X * samplesize, cmpTupleSXCompare );

	    DBG_ARRAY2(1, "Sample splitters", samplebufall, nprocs * X * samplesize);

	    for ( int i = 0; i < nprocs; i++ )		// pick splitters
		splitterbuf[i] = samplebufall[ i * X * samplesize ];

	    DBG_ARRAY2(1, "Selected splitters", splitterbuf, nprocs);

	    delete[] samplebufall;
	}

	// distribute splitters

	MPI_Bcast( splitterbuf, nprocs, MPI_TUPLE_SX, ROOT, MPI_COMM_WORLD );

	// find nearest splitters in each of the locally sorted tuple list

	uint** splitterpos = new uint*[X];

	for (unsigned int k = 0; k < X; ++k)
	{
	    splitterpos[k] = new uint[ nprocs + 1 ];

	    splitterpos[k][0] = 0;
	    for ( int i = 1; i < nprocs; i++ )
	    {
		std::vector<TupleSX>::const_iterator it = std::lower_bound(S[k].begin(), S[k].end(), splitterbuf[i], cmpTupleSXCompare );

		splitterpos[k][i] = it - S[k].begin();
	    }
	    splitterpos[k][ nprocs ] = M2;
	}

	for (unsigned int k = 0; k < X; ++k)
	{
	    DBG_ARRAY2(1, "Splitters S." << k, splitterpos[k], nprocs+1);
	}

	delete [] splitterbuf;
	
	// boardcast number of element in each division

	int** sendcnt = new int*[X];
	int** recvcnt = new int*[X];

	for (unsigned int k = 0; k < X; ++k)
	{
	    sendcnt[k] = new int[ nprocs ];
	    recvcnt[k] = new int[ nprocs ];

	    for ( int i = 0; i < nprocs; i++ )
	    {
		sendcnt[k][i] = splitterpos[k][i+1] - splitterpos[k][i];

		assert( sendcnt[k][i] >= 0 );
	    }

	    delete [] splitterpos[k];

	    MPI_Alltoall( sendcnt[k], 1, MPI_INT, recvcnt[k], 1, MPI_INT , MPI_COMM_WORLD );
	}

	delete [] splitterpos;

	// calculate number of received tuples

	int** sendoff = new int*[X];
	int** recvoff = new int*[X];
	unsigned int totalsize = 0;

	for (unsigned int k = 0; k < X; ++k)
	{
	    sendoff[k] = new int[ nprocs + 1 ];
	    recvoff[k] = new int[ nprocs + 1 ];

	    sendoff[k][ 0 ] = recvoff[k][ 0 ] = 0;
	    for ( int i = 1; i < nprocs + 1; i++ )
	    {
		sendoff[k][i] = sendoff[k][i-1] + sendcnt[k][i-1];
		recvoff[k][i] = recvoff[k][i-1] + recvcnt[k][i-1];

		assert(sendoff[k][i] >= 0);
		assert(recvoff[k][i] >= 0);
	    }

	    std::vector<TupleSX> recvbuf ( recvoff[k][ nprocs ] );

	    MPI_Alltoallv( S[k].data(), sendcnt[k], sendoff[k], MPI_TUPLE_SX, recvbuf.data(), recvcnt[k], recvoff[k], MPI_TUPLE_SX, MPI_COMM_WORLD );

	    std::swap(S[k], recvbuf);

	    delete [] sendcnt[k];
	    delete [] recvcnt[k];
	    delete [] sendoff[k];

	    totalsize += S[k].size();
	}

	// merge received array parts

	for (unsigned int k = 0; k < X; ++k)
	{
	    if ( S[k].size() )
	    {
		/** merge */
		TupleSX * helparray4 = new TupleSX[ S[k].size() ];

		// binary merge sort!
		mergesort( S[k].data(), helparray4, 0, nprocs, recvoff[k], cmpTupleSXdepth<6> );
	    }
	}

	for (unsigned int k = 0; k < X; ++k)
	{
	    DBG_ARRAY(1, "After samplesort S" << k, S[k]);
	}

	std::vector<TupleSX> suffixarray (totalsize);
	uint* suffixarray2 = new uint[ totalsize ];
	int j = 0;

	TupleSXMerge<X> tuplecmp(S);
	LoserTree< TupleSXMerge<X> > LT(X, tuplecmp);

	int top;

	while( (top = LT.top()) >= 0 )
	{
	    suffixarray[j] = S[top][ tuplecmp.m_ptr[top] ];
	    suffixarray2[j] = suffixarray[j].index;

	    std::cout << "Winning tuple: " << suffixarray[j] << "\n";

	    if (suffixarray[j].index < globalSize) ++j;

	    tuplecmp.m_ptr[top]++;

	    LT.replay( tuplecmp.m_ptr[top] >= S[top].size() );
	}

	DBG_ARRAY2(1, "Suffixarray merged", suffixarray, j);

	*salen = j;
	return suffixarray2;
    }
}

void writesa( uint* salen, uint* suffixarray, char* filenameOut )
{
    /** Suffixarry in Datei schreiben */
    MPI_Request request;

    uint* saLenAll = new uint[ nprocs ];
    MPI_Gather( salen, 1, MPI_UNSIGNED, saLenAll, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

    if ( myproc == ROOT )
    {
	int fd = open(filenameOut, O_WRONLY | O_CREAT | O_TRUNC , S_IRUSR | S_IWUSR);

	if (fd < 0) { std::cout << "Fehler bei der Ausgabe" << std::endl; return; }

	uint maxLen = 0;
	for ( int i = 0;i < nprocs;i++ ) if ( maxLen < saLenAll[ i ] )   maxLen = saLenAll[ i ];

	write(fd, &suffixarray[1], saLenAll[0] * sizeof(uint));

	std::cout << "Von 0 Daten geschrieben" << std::endl;
	delete[] suffixarray;

	uint* outbuffer = new uint[ maxLen ];

	for (int k = 1; k < nprocs; k++)
	{
	    MPI_Recv( outbuffer, saLenAll[k], MPI_UNSIGNED, k, MSGTAG, MPI_COMM_WORLD, &status );
	    std::cout << "salen von " << k << " ist " << saLenAll[ k ] << std::endl;

	    write(fd, outbuffer, saLenAll[k] * sizeof(uint));

//          Debug6(for (uint i=0 ; i<saLenAll[k] ; i++)
//                      std::cout<< outbuffer[i]<<" ";
//                  std::cout << std::endl;)

	    std::cout << "Von " << k << " Daten geschrieben" << std::endl;
	}
	delete[] saLenAll;
	delete[] outbuffer;

	if (close(fd) != 0) std::cout << "Fehler beim Schließen der Datei";

    }
    else {
	MPI_Isend( suffixarray, *salen , MPI_UNSIGNED, 0, MSGTAG , MPI_COMM_WORLD , &request );
	MPI_Wait( &request, &status );
	delete[] suffixarray;
	std::cout << myproc << " ist fertig" << std::endl;
    }
}

void getTuple(uint* inbuffer, Pair* recvBufPair, Quintuple* S0, Quadruple* S1, Quintuple* S2, uint localSize, uint* n,  uint* imod3, uint& half)
{
    /** Tupel erstellen
	 *    S0[j]={inbuffer[i], inbuffer[i+1], P[i*2/3], P[i*2/3+1],i};              // mit i mod 3 == 0
	 *    S1[j]={P[i*2/3], inbuffer[i+1], P[i*2/3+1],i};                          // mit i mod 3 == 1
	 *    S2[j]={P[i*2/3+1], inbuffer[i+2], inbuffer[i+3], P[i*2/3+2],i};        // mit i mod 3 == 2
	 */
    uint pairOffset[ 6 ];
    pairOffset[ 0 ] = imod3[ 0 ] / 2;
    pairOffset[ 1 ] = half + ( imod3[ 0 ] + 1 ) / 2;   //imod3[0]/2+(imod3[2]?0:1);
    pairOffset[ 2 ] = 0;
    pairOffset[ 3 ] = half + imod3[ 1 ] / 2; //(imod3[2]?0:1);
    pairOffset[ 4 ] = ( imod3[ 2 ] + 1 ) / 2; //(imod3[2]?1:0);
    pairOffset[ 5 ] = half ;
    uint displacement = myproc * localSize;
    uint temp;
    uint n012 = n[0]+n[1]+n[2];
    uint index[3]={0,0,0};
    while (index[0] + index[1] +index[2] < n012 ){
	if    (index[0] < n[0]){
			temp=imod3[ 0 ] + 3 * index[0];
			S0[ index[0] ].name[ 0 ] = inbuffer[ temp ];
			S0[ index[0] ].name[ 1 ] = inbuffer[ temp + 1 ];
			S0[ index[0] ].name[ 2 ] = recvBufPair[ index[0] + pairOffset[ 0 ] ].name;
			S0[ index[0] ].name[ 3 ] = recvBufPair[ index[0] + pairOffset[ 1 ] ].name;
			S0[ index[0] ].index = temp + displacement;
			index[0]++;
	}
	if    (index[1] < n[1]){
			temp=imod3[ 1 ] + 3 * index[1];
			S1[ index[1] ].name[ 0 ] = recvBufPair[ index[1] + pairOffset[ 2 ] ].name;
			S1[ index[1] ].name[ 1 ] = inbuffer[ temp ];
			S1[ index[1] ].name[ 2 ] = recvBufPair[ index[1] + pairOffset[ 3 ] ].name;
			S1[ index[1] ].index = temp + displacement;
			index[1]++;
	}
	if    (index[2] < n[2]){
			temp = 3 * index[2] + imod3[ 2 ];
			S2[ index[2] ].name[ 0 ] = recvBufPair[ index[2] + pairOffset[ 5 ] ].name;
			S2[ index[2] ].name[ 1 ] = inbuffer[ temp ];
			S2[ index[2] ].name[ 2 ] = inbuffer[ temp + 1 ];
			S2[ index[2] ].name[ 3 ] = recvBufPair[ index[2] + pairOffset[ 4 ] ].name;
			S2[ index[2] ].index = temp + displacement;
			index[2]++;
	}
    }
    delete[] inbuffer;
    delete[] recvBufPair;
}
#if 0
uint* sortS0S1S2( uint* inbuffer, Pair* recvBufPair, uint localSize, uint* n, uint* imod3, uint* salen, uint& half )
{
    Debug2( sortS0S1S2TimeStart = MPI_Wtime() );
//  uint k=(uint) ceil(samplefactor*(sqrt((double)localSize)/300));
    uint k = (uint)(sqrt((double) localSize / nprocs) *samplefactor/3.0);
    if(k>=localSize/3){std::cout<<"k="<<k<<" zu groß, jetzt k="<<localSize/3-1<<std::endl; k=localSize/3-1;}

//  uint k1 = (uint)sqrt( localSize / (3.0 * nprocs) *samplefactor);
//   if ( myproc == ROOT ) std::cout<<" k "<<3*k<<" % "<<3.0*k/localSize*100.0<<" Tuple " <<localSize<<std::endl;
    //    std::cout << myproc << " 0:" << imod3[ 0 ] << " 1:" << imod3[ 1 ] << " 2:" << imod3[ 2 ] << " localSize " <<  localSize << std::endl;
    Quintuple* S0 = new Quintuple[ n[ 0 ] ];
    Quadruple* S1 = new Quadruple[ n[ 1 ] ];
    Quintuple* S2 = new Quintuple[ n[ 2 ] ];
    getTuple( inbuffer,  recvBufPair,  S0,  S1,  S2, localSize,  n,  imod3,  half);

    /**Samplesort für S1 */
    Debug2( starttime = MPI_Wtime() );
    Switch0(
	std::sort( S0, S0 + n[ 0 ], cmpNameRank ); //Übergebe Funktion
	std::sort( S1, S1 + n[ 1 ], cmpFirstNameS1 ); //Übergebe Funktion
	std::sort( S2, S2 + n[ 2 ], cmpFirstNameS2 ); //Übergebe Funktion
	);
    Switch1(
	std::stable_sort( S0, S0 + n[ 0 ], cmpNameRank ); //Übergebe Funktion
	std::stable_sort( S1, S1 + n[ 1 ], cmpFirstNameS1 ); //Übergebe Funktion
	std::stable_sort( S2, S2 + n[ 2 ], cmpFirstNameS2 ); //Übergebe Funktion
	);
    Debug2( sortingTime += MPI_Wtime() - starttime );
    Debug6(stlSortTime=MPI_Wtime() - starttime;
//  MPI_Reduce(&MergeSA12, &stlSortTime, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
//  MPI_Reduce(&MergeSA12, &stlSortTime, 1, MPI_DOUBLE, MPI_MIN, ROOT, MPI_COMM_WORLD);
	   std::cout<<myproc<<"   "<<stlSortTime<<" sec "<< n[0]+n[1]+n[2]<<" Tuple sort S0, S1, S2 "<<std::endl;
	);

    /** Sampling */
    Quintuple* samplebufS0 = new Quintuple[ k ];
    Quintuple* samplebufS1 = new Quintuple[ k ];
    Quintuple* samplebufS2 = new Quintuple[ k ];
    double d[3] = {( double ) n[ 0 ] / k,( double ) n[ 1] / k,( double ) n[ 2 ] / k};

    for ( uint i = 0; i < k; i++ ) {
	samplebufS0[ i ] = S0[ int( i * d[0] ) ];
	samplebufS0[i].index = samplebufS0[i].index/3;

	samplebufS1[ i ].name[0] = S1[ int( i * d[1] ) ].name[0];
	samplebufS1[ i ].name[1] = S1[ int( i * d[1] ) ].name[1];
	samplebufS1[ i ].name[2] = S1[ int( i * d[1] ) ].name[2];
	samplebufS1[i].index = (S1[i].index / 3) + MOD0;
	//if(samplebufS1[i].index>MOD0+(localSize+1)/3)std::cout<<"Fehler samplebufS1[i].index>MOD0+n[i] "<< samplebufS1[i].index<<" "<< S1[i].index/3 <<" "<<MOD0<<" "<<(localSize+1)/3 <<std::endl;

	samplebufS2[ i ] = S2[ int( i * d[2] ) ];
	samplebufS2[i].index = samplebufS2[i].index/3+MOD1;
    }

    /**Proc 0 sammelt Samples ein*/
    Quintuple* samplerecvbufS0=NULL;
    Quintuple* samplerecvbufS1=NULL;
    Quintuple* samplerecvbufS2=NULL;

    if ( myproc == ROOT ) {
	samplerecvbufS0 = new Quintuple[ k * nprocs ];
	samplerecvbufS1 = new Quintuple[ k * nprocs +1];
	samplerecvbufS2 = new Quintuple[ k * nprocs +1];
    }

    Debug2(sampleStartTime=MPI_Wtime() );
    Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

    MPI_Gather( samplebufS0, k, MPI_QUINTUPLE, samplerecvbufS0, k, MPI_QUINTUPLE, ROOT, MPI_COMM_WORLD );
    MPI_Gather( samplebufS1, k, MPI_QUINTUPLE, samplerecvbufS1, k, MPI_QUINTUPLE, ROOT, MPI_COMM_WORLD );
    MPI_Gather( samplebufS2, k, MPI_QUINTUPLE, samplerecvbufS2, k, MPI_QUINTUPLE, ROOT, MPI_COMM_WORLD );

    Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime);

    delete[] samplebufS0;
    delete[] samplebufS1;
    delete[] samplebufS2;

    Quintuple* pivbuf= new Quintuple[ nprocs ];

    /**Proc 0 ermittelt Splitter*/
    if (myproc==ROOT)
    {
	Switch0(
	    std::sort( samplerecvbufS0, samplerecvbufS0 + k * nprocs, cmpNameRank );
	    std::sort( samplerecvbufS1, samplerecvbufS1 + k * nprocs, cmpFirstNameS2 );
	    std::sort( samplerecvbufS2, samplerecvbufS2 + k * nprocs, cmpFirstNameS2 );
	    );
	Switch1(
	    std::stable_sort( samplerecvbufS0, samplerecvbufS0 + k * nprocs, cmpNameRank );
	    std::stable_sort( samplerecvbufS1, samplerecvbufS1 + k * nprocs, cmpFirstNameS2 );
	    std::stable_sort( samplerecvbufS2, samplerecvbufS2 + k * nprocs, cmpFirstNameS2 );
	    );
	samplerecvbufS1[k*nprocs].name[0]=MAX_INT;
	samplerecvbufS1[k*nprocs].name[1]=MAX_INT;
	samplerecvbufS1[k*nprocs].name[2]=MAX_INT;
	samplerecvbufS1[k*nprocs].name[3]=MAX_INT;
	samplerecvbufS2[k*nprocs].name[0]=MAX_INT_SAMPLE;
	samplerecvbufS2[k*nprocs].name[1]=MAX_INT_SAMPLE;
	samplerecvbufS2[k*nprocs].name[2]=MAX_INT_SAMPLE;
	samplerecvbufS2[k*nprocs].name[3]=MAX_INT_SAMPLE;
	uint k3=3*k;
	Quintuple* Samples= new Quintuple[k3*nprocs+2];
//std::cout<<myproc<<" 1006"<<std::endl;
	merge2(samplerecvbufS0,samplerecvbufS1,samplerecvbufS2,Samples, k*nprocs);
//      std::cout<<"samples"<<std::endl;
//      for(int i =0;i<k3*nprocs+2;i++) Samples[i].print();
//std::cout<<myproc<<" 1008"<<std::endl;
	delete[] samplerecvbufS0;
	delete[] samplerecvbufS1;
	delete[] samplerecvbufS2;
//      std::cout<<"pivbuf"<<std::endl;
	for ( int i = 0; i < nprocs; i++ ) {
	    pivbuf[ i ] = Samples[ i * k3 ];
//          pivbuf[ i ].print();
	}
	delete[] Samples;
    }

    /**Splitter verteilen*/
    Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );
    MPI_Bcast( pivbuf, nprocs, MPI_QUINTUPLE, ROOT, MPI_COMM_WORLD );
// cerr <<myproc<<"----- 1031"<<std::endl;
    Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime);
    Debug2( sampleTime+=MPI_Wtime()-sampleStartTime );
//std::cout<<myproc<<" 1021"<<std::endl;

    /**Einordnen der lokal sortierten Elemente anhand von Splitter*/
    uint* splitterpos0= new uint[nprocs+1];
    uint* splitterpos1= new uint[nprocs+1];
    uint* splitterpos2= new uint[nprocs+1];
    splitterpos0[0]=0;
    splitterpos1[0]=0;
    splitterpos2[0]=0;
    MPI_Barrier(MPI_COMM_WORLD);
// cerr <<myproc<<"------- 1041"<<std::endl;

    for ( int i = 1; i < nprocs; i++ ){
	if(pivbuf[ i ].index<MOD0){
	    //             if (myproc==0){std::cout<<"0      -- ";pivbuf[i] .print();}
	    pivbuf[i].index=pivbuf[i].index*3;
	    splitterpos0[ i ] = binsearch( S0, pivbuf[ i ], n[ 0 ], cmpS0S0 );
	    splitterpos1[ i ] = binsearch( S1, pivbuf[ i ], n[ 1 ], cmpS0S1 );
	    splitterpos2[ i ] = binsearch( S2, pivbuf[ i ], n[ 2 ], cmpS0S2 );
	}
	else if(pivbuf[ i ].index<MOD1){
//              if (myproc==0){std::cout<<"1      -- ";pivbuf[i] .print();}
	    pivbuf[ i ].index=(pivbuf[i].index-MOD0)*3+1;
	    splitterpos0[ i ] = binsearch( S0, pivbuf[ i ], n[ 0 ], cmpS1S0 );
	    splitterpos1[ i ] = binsearch( S1, pivbuf[ i ], n[ 1 ], cmpS1S1 );
	    splitterpos2[ i ] = binsearch( S2, pivbuf[ i ], n[ 2 ], cmpS1S2 );
	}
	else{
//              if (myproc==0){std::cout<<"2      -- ";pivbuf[i] .print();}
	    pivbuf[i].index=(pivbuf[i].index-MOD1)*3+2;
	    splitterpos0[ i ] = binsearch( S0, pivbuf[ i ], n[ 0 ], cmpS2S0 );
	    splitterpos1[ i ] = binsearch( S1, pivbuf[ i ], n[ 1 ], cmpS2S1 );
	    splitterpos2[ i ] = binsearch( S2, pivbuf[ i ], n[ 2 ], cmpS2S2 );
	}
    }

    delete[] pivbuf;
    splitterpos0[ nprocs ] = n[ 0 ];
    splitterpos1[ nprocs ] = n[ 1 ];
    splitterpos2[ nprocs ] = n[ 2 ];

    int* sendcnt0= new int[nprocs];
    int* sendcnt1= new int[nprocs];
    int* sendcnt2= new int[nprocs];
    int* recvcnt0= new int[nprocs];
    int* recvcnt1= new int[nprocs];
    int* recvcnt2= new int[nprocs];


    for ( int i = 0; i < nprocs; i++ ) {
	sendcnt0[ i ] = splitterpos0[ i + 1 ] > splitterpos0[ i ] ? splitterpos0[ i + 1 ] - splitterpos0[ i ] : 0;
	sendcnt1[ i ] = splitterpos1[ i + 1 ] > splitterpos1[ i ] ? splitterpos1[ i + 1 ] - splitterpos1[ i ] : 0;
	sendcnt2[ i ] = splitterpos2[ i + 1 ] > splitterpos2[ i ] ? splitterpos2[ i + 1 ] - splitterpos2[ i ] : 0;
    }

    Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

    MPI_Alltoall( sendcnt0, 1, MPI_INT, recvcnt0, 1, MPI_INT , MPI_COMM_WORLD );
    MPI_Alltoall( sendcnt1, 1, MPI_INT, recvcnt1, 1, MPI_INT , MPI_COMM_WORLD );
    MPI_Alltoall( sendcnt2, 1, MPI_INT, recvcnt2, 1, MPI_INT , MPI_COMM_WORLD );

    Debug3( mpiComTime += MPI_Wtime()-mpiComStartTime );

    int* sendoff0= new int[nprocs+1];
    int* sendoff1= new int[nprocs+1];
    int* sendoff2= new int[nprocs+1];
    int* recvoff0= new int[nprocs+1];
    int* recvoff1= new int[nprocs+1];
    int* recvoff2= new int[nprocs+1];

    sendoff0[ 0 ] = recvoff0[ 0 ] = 0;
    sendoff1[ 0 ] = recvoff1[ 0 ] = 0;
    sendoff2[ 0 ] = recvoff2[ 0 ] = 0;

    for ( int i = 1; i < nprocs + 1; i++ ) {
	sendoff0[ i ] = sendoff0[ i - 1 ] + sendcnt0[ i - 1 ];
	sendoff1[ i ] = sendoff1[ i - 1 ] + sendcnt1[ i - 1 ];
	sendoff2[ i ] = sendoff2[ i - 1 ] + sendcnt2[ i - 1 ];
	recvoff0[ i ] = recvoff0[ i - 1 ] + recvcnt0[ i - 1 ];
	recvoff1[ i ] = recvoff1[ i - 1 ] + recvcnt1[ i - 1 ];
	recvoff2[ i ] = recvoff2[ i - 1 ] + recvcnt2[ i - 1 ];
    }

    n[ 0 ] = recvoff0[ nprocs ];

    Quintuple* recvbufS0 = new Quintuple[ n[0] ];
    Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)     starttime = MPI_Wtime() );
    MPI_Alltoallv( S0, sendcnt0, sendoff0, MPI_QUINTUPLE, recvbufS0, recvcnt0, recvoff0, MPI_QUINTUPLE, MPI_COMM_WORLD );
    Debug3( alltoallvTime += MPI_Wtime() - starttime );
    delete[] S0;

    /** merge */
    Quintuple* helparray5 = new Quintuple[ n[0] ];
    Debug2(starttime = MPI_Wtime() );

    mergesort( recvbufS0, helparray5, 0, nprocs, recvoff0, cmpSplitterLeqS0 );

    Debug2( mergeTime += MPI_Wtime() - starttime );

    delete[] helparray5;

    /** Samplesort S0 end */

    /** ein 'MAX'-Element am Ende mehr */
    n[ 1 ] = recvoff1[ nprocs ];
    Quadruple* recvbufS1 = new Quadruple[ n[1] + 1 ];
    Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)     starttime = MPI_Wtime() );
    MPI_Alltoallv( S1, sendcnt1, sendoff1, MPI_QUADRUPLE, recvbufS1, recvcnt1, recvoff1, MPI_QUADRUPLE, MPI_COMM_WORLD );
    Debug3( alltoallvTime += MPI_Wtime() - starttime );
    delete[] S1;

    /** merge */
    Quadruple* helparray4 = new Quadruple[ n[1] ];
    Debug2( starttime = MPI_Wtime() );
    mergesort( recvbufS1, helparray4, 0, nprocs, recvoff1, cmpFirstNameS1 );
    Debug2( mergeTime += MPI_Wtime() - starttime );
    delete[] helparray4;

    recvbufS1[ n[1] ].name[ 0 ] = MAX_INT_SAMPLE;
    recvbufS1[ n[1] ].name[ 1 ] = MAX_INT_SAMPLE;
    recvbufS1[ n[1] ].name[ 2 ] = MAX_INT_SAMPLE;

    n[ 2 ] = recvoff2[ nprocs ];
    Quintuple* recvbufS2 = new Quintuple[ n[ 2 ] + 1 ];

    Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)     starttime = MPI_Wtime() );
    MPI_Alltoallv( S2, sendcnt2 , sendoff2 , MPI_QUINTUPLE, recvbufS2, recvcnt2, recvoff2, MPI_QUINTUPLE, MPI_COMM_WORLD );
    Debug3( alltoallvTime += MPI_Wtime() - starttime );
    delete[] S2;

    /** merge */
    helparray5 = new Quintuple[ n[ 2 ] ];
    Debug2( starttime = MPI_Wtime() );
    mergesort( recvbufS2, helparray5, 0, nprocs, recvoff2, cmpFirstNameS2 );
    Debug2( mergeTime += MPI_Wtime() - starttime );
    delete[] helparray5;

    recvbufS2[ n[ 2 ] ].name[ 0 ] = MAX_INT_SAMPLE;
    recvbufS2[ n[ 2 ] ].name[ 1 ] = MAX_INT_SAMPLE;
    recvbufS2[ n[ 2 ] ].name[ 2 ] = MAX_INT_SAMPLE;
    recvbufS2[ n[ 2 ] ].name[ 3 ] = MAX_INT_SAMPLE;

    /** Samplesort S1, S2 end */
    *salen = n[ 0 ] + n[ 1 ] + n[ 2 ];
    uint* suffixarray = new uint[ *salen ];

    delete[] sendcnt0;
    delete[] sendoff0;
    delete[] recvcnt0;
    delete[] recvoff0;
    delete[] splitterpos0;
    delete[] sendcnt1;
    delete[] sendoff1;
    delete[] recvcnt1;
    delete[] recvoff1;
    delete[] splitterpos1;
    delete[] sendcnt2;
    delete[] sendoff2;
    delete[] recvcnt2;
    delete[] recvoff2;
    delete[] splitterpos2;

    Debug2( starttime = MPI_Wtime() );
    merge( recvbufS0, recvbufS1, recvbufS2, suffixarray, n );
    Debug2( mergeTime += MPI_Wtime() - starttime );

    Debug5( //Lastverteilung
	double sortSAll= MPI_Wtime() - starttime;
	double maxsortAll;
	double minsortAll;
	uint max;
	uint min;
//  Debug6(std::cout << myproc << " SortS0S1S2 Alltoallv: before "<< localSize <<" after " << *salen << std::endl; )
	MPI_Reduce(salen, &max, 1, MPI_UNSIGNED, MPI_MAX, ROOT, MPI_COMM_WORLD);
	MPI_Reduce(salen, &min, 1, MPI_UNSIGNED, MPI_MIN, ROOT, MPI_COMM_WORLD);
	MPI_Reduce(&sortSAll, &maxsortAll, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
	MPI_Reduce(&sortSAll, &minsortAll, 1, MPI_DOUBLE, MPI_MIN, ROOT, MPI_COMM_WORLD);

	uint space=11;
/**
   min -- max -- max-min -- %-Diff -- %-Tuple -- merge-min  -- merge-max
*/
	if (myproc==ROOT) {std::cout<< setw( space )<<min<< setw( space )<<max<< setw( space )<<max-min<< setw( space );
//  std::cout.setf( ios::scientific, ios::floatfield );
	    std::cout << (double)(max-min)/min*100<< setw( space )<<(300.0*k)/localSize<< setw( space )<<minsortAll<< setw( space )<< maxsortAll<<setw (space)<<(maxsortAll-minsortAll)/minsortAll*100<<" Alltoallv all Tuple"<<std::endl;}
	);

    delete[] recvbufS0;
    delete[] recvbufS1;
    delete[] recvbufS2;

//  Debug6( if ( myproc == 0 ) cerr << " ende sortS0S1S2 " << std::endl; )
    Debug2( sortS0S1S2TimeEnd+= (MPI_Wtime()-sortS0S1S2TimeStart) );
// std::cout<<myproc<<" hat "<<n[ 0 ] + n[ 1 ] + n[ 2 ] <<" zu sortieren in "<<(MPI_Wtime()-sortS0S1S2TimeStart) <<std::endl;
    return suffixarray;
}
#endif

/**
 * main
 * @param argc
 * @param argv
 * @return
 */

int main( int argc, char **argv )
{
    MPI_Init( &argc, &argv );

    MPI_Comm_rank( MPI_COMM_WORLD, &myproc );
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

    if ( nprocs <= 1 )
    {
	std::cerr << "Error: requires more than one MPI processor (use -np 2)." << std::endl;
	return -1;
    }

    mpi_init_datatypes();

    // **********************************************************************
    // * Read input file size

    uint filelength;
    const char* filename = NULL;

    if ( myproc == ROOT )
    {
	if ( argc < 2 ) {
	    std::cout << "No input file! Call using mpdrun -np 4 ./pDCX <input-file> [output-file]" << std::endl;
	    return 0;
	}

	filename = argv[1];
	std::ifstream infile( filename );
	if (!infile.good()) {
	    perror("Cannot read input file");
	    return -1;
	}

	// determine file size
	infile.seekg( 0, std::ios::end );
	filelength = infile.tellg();

	//if (argc > 3 ) if (atoll(argv[3])< filelength && atoll(argv[3])>0) filelength= atoll(argv[3]);
    }

    MPI_Bcast( &filelength, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD );	// synchronized broadcast

    samplefactor = nprocs;	// TODO

    // **********************************************************************
    // * Calculate local input size and send data portion

    uint localSize = (uint)ceil( (double)filelength / nprocs );	// divide by processors rounding up
    localSize += X - localSize % X;				// round up to nearest multiple of X

    assert( localSize * nprocs >= filelength );

    if ( myproc == ROOT )
    {
	std::cout << "Total input size = " << filelength << " bytes. localSize = " << localSize << std::endl;
    }

    // Verteilung der Daten nach inbuffer<- enthält lokal benötigte Zeichen +2 wegen Überschneidung und zum Erstellen der S0-S1 Tupel

    // Send data portions to _other_ n-1 processors -- including an overlap of +(X-1) characters for the final tuples

    int overlap = X-1;
    char* inputChars = new char[ localSize + overlap ];

    if ( myproc == ROOT )
    {
	std::ifstream infile( filename, std::ios::binary );

	// Einlesen der Elemente für die Prozessoren 1 bis n-2
	for ( int p = 1; p < nprocs; p++ )
	{
	    infile.seekg( p * localSize, std::ios::beg );

	    uint readsize = (p != nprocs-1) ? localSize + overlap : filelength - (p * localSize);

	    std::cout << "Read portion " << p << " from pos " << p * localSize << " of length " << readsize << std::endl;

	    infile.read( inputChars, readsize );
	    MPI_Send( inputChars, readsize, MPI_CHAR, p, MSGTAG, MPI_COMM_WORLD );
	}

	if (!infile.good()) {
	    perror("Error reading file.");
	    return -1;
	}

	// Einlesen der Elemente für Prozessor ROOT

	std::cout << "Read portion 0 from pos 0 of length " << localSize + overlap << std::endl;

	infile.seekg( 0, std::ios::beg );
	infile.read( inputChars, localSize + overlap );

	if (!infile.good()) {
	    perror("Error reading file.");
	    return -1;
	}

	Debug5(
	    std::cout << "         min   -- max     -- Diff -- %-Diff -- %-Tuple -- merge-min--merge-max-- %-Diff alltoall-min alltoall-max"<<std::endl
	    );
    }
    else // not ROOT
    {
	MPI_Recv( inputChars, localSize + overlap, MPI_CHAR, ROOT, MSGTAG, MPI_COMM_WORLD, &status );
    }

    // convert input characters into uint data

    uint* input = new uint[ localSize + overlap ];

    for ( unsigned int i = 0; i < localSize + overlap; i++ ) {
	input[i] = inputChars[i];
    }

    delete[] inputChars;

    MPI_Barrier( MPI_COMM_WORLD );

    // **********************************************************************
    // * Construct suffix array recursively
    
    uint    salen;
    uint*   suffixarray = dc3( input, filelength, &salen, localSize );

    if ( argc > 2 ) {
	//writesa( &salen, suffixarray, argv[ 2 ] );
    }

    // **********************************************************************
    // * Collect and check suffix array

    {
	int* recvcnt = new int[nprocs];
	int* recvoff = new int[nprocs+1];

	MPI_Gather( &salen, 1, MPI_INT, recvcnt, 1, MPI_INT, ROOT, MPI_COMM_WORLD );

	if (myproc == ROOT)
	{
	    recvoff[0] = 0;
	    for ( int i = 1; i <= nprocs; i++ ) {
		recvoff[i] = recvoff[i-1] + recvcnt[i-1];
	    }
	    assert( recvoff[nprocs] == filelength );
	}

	std::vector<uint> SA;

	if (myproc == ROOT)
	    SA.resize( filelength );

	MPI_Gatherv(suffixarray, salen, MPI_INT,
		    SA.data(), recvcnt, recvoff, MPI_INT, ROOT, MPI_COMM_WORLD);

	if (myproc == ROOT)
	{
	    DBG_ARRAY(1, "Suffixarray collected", SA);
	    
	    std::ifstream infile( filename );
	    if (!infile.good()) {
		perror("Cannot read input file");
		return -1;
	    }

	    std::vector<char> string(filelength);

	    infile.read((char*)string.data(), filelength);
	    if (!infile.good()) {
		perror("Cannot read input file");
		return -1;
	    }

	    std::cout << "result suffix array: \n";

	    for (unsigned int i = 0; i < SA.size(); ++i)
	    {
		std::cout << i << " : " << SA[i] << " : ";

		for (unsigned int j = 0; SA[i]+j < filelength; ++j)
		{
		    std::cout << string[SA[i]+j] << " ";
		}

		std::cout << "\n";
	    }

	    assert( sachecker::sa_checker(string, SA) );
	}
    }

    MPI_Finalize();
    return 0;
}
