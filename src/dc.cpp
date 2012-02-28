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

static const bool debug = true;

#include <mpi.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <fcntl.h>
#include <sys/stat.h>
#include <stdlib.h>

#include "tuple.h"

#include "util.h"

#include <cassert>
#include <cmath>
#include <vector>

typedef unsigned int	uint;

#define MAX_INT ~1U
#define MAX_INT_SAMPLE ~0U
//#define MOD2_INDEX ~1UL // ~1UL % 3 = 4294967294 % 3 = 2  max mod 2 Index for UNSIGNED INTEGER
#define MIN_INT 0UL
#define MSGTAG 42
#define MOD0 ~0U/3
#define MOD1 2*(~0U/3)


uint* dc3( uint* inbuffer, uint filelength, uint* salen );
uint* sortS0S1S2( uint* inbuffer, Pair* P, uint localArraylen, uint* n, uint* imod3, uint* salen, uint& half);
inline void getTuple(uint* inbuffer, Pair* recvBufPair, Quintuple* S0, Quadruple* S1, Quintuple* S2, uint localArraylen, uint* n,  uint* imod3, uint& half);

inline int namelex( TupleA<X>* in, Pair* P, uint arraylen );

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

/**
 *
 * @param in
 * @param out
 * @param length
 * @param n12all
 * @param names
 */
void sortPairIndex( Pair* in, Pair* out, uint length, uint* n12all ,uint names)
{
	if (myproc==0)
		for ( uint i = 0 ; i < length; i++ )
			out[( in[ i ].index % 3 - 1 ) * ( names / 2 ) + in[ i ].index / 3 ]=in[i];
	else
		for ( uint i = 0 ; i < length; i++ ) {
			out[( in[ i ].index % 3 - 1 ) * ( names / 2 ) + in[ i ].index / 3 - n12all[myproc-1]]=in[i];
		}

}

MPI_Datatype MPI_TUPLEX;

MPI_Datatype MPI_QUADRUPLE;
MPI_Datatype MPI_PAIR;

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
    MPI_Datatype mpi_tupel[2] = { MPI_UNSIGNED, MPI_UNSIGNED };

    int mpi_blocklen[2] = { X, 1 };
    MPI_Aint mpi_displ[2] = { 0, X * sizeof( uint ) };
    MPI_Type_struct( 2, mpi_blocklen, mpi_displ, mpi_tupel, &MPI_TUPLEX );
    MPI_Type_commit( &MPI_TUPLEX );

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

    mpi_blocklen[ 0 ] = 1;
    mpi_displ[ 1 ] = sizeof( uint );
    MPI_Type_struct( 2, mpi_blocklen, mpi_displ, mpi_tupel, &MPI_PAIR );
    MPI_Type_commit( &MPI_PAIR );
#endif
}


/**
 *
 * @param inbuffer
 * @param filelength
 * @param salen
 * @return suffixarray
 */
uint* dc3( uint* input, uint N, uint* salen, uint localArraylen )
{
    // Debug6( std::cout << "dc3 gestartet" << std::endl; )

    /** Tupel einlesen*/

    // Anzahl mod 1, mod 2, mod 3 Tupel berechnen. Ist für jedes PE +-1 Element gleich
    //uint localArraylen = ( filelength + nprocs - 1 ) / nprocs ;

    uint samplesize = (uint)sqrt(2.0 * localArraylen / 3.0 / nprocs) * samplefactor;

// Debug5(  if ( myproc == ROOT ) std::cout << "samplefactor=" << k <<" |Tuple|="<<2*(localArraylen/3)<<" in %="<<100*k/(2*localArraylen/3.0)<<std::endl;)
//  int k=(uint) ceil(samplefactor*(sqrt(2.0*localArraylen/3.0)/100));

    if( samplesize >= 2*(localArraylen/3)) samplesize = 2*(localArraylen/3)-1;

// if ( myproc == ROOT ) std::cout << "samplefactor=" << k <<" |Tuple|="<<2*(localArraylen/3)<<" in %="<<100*k/(2*localArraylen/3.0)<<std::endl;

    unsigned int myAoffset = myproc * localArraylen;
    unsigned int myArraylen = (myproc != nprocs-1) ? localArraylen : N - myAoffset;

    printf("myArraylen = %d\n", myArraylen);

    // **********************************************************************
    // * calculate build DC-tuple array and sort locally

    const unsigned int M = myArraylen / X;		// multiples of complete X chars in local area size
    const unsigned int M2 = (X * M == myArraylen ? M : M+1); // number of incomplete X chars in local area size
    unsigned int Rsize = D * M2;

    std::cout << "myAlen = " << myArraylen << " - M = " << M << " - M2 = " << M2 << "\n";

    TupleA<X>* R = new TupleA<X>[ Rsize ];		// create D*M2 tuples which might include up to M-1 dummies

    uint j = 0 ;
    for (uint i = 0; i < myArraylen; i += X)
    {
	for (uint d = 0; d < D; ++d)
	{
	    R[j].index = myAoffset + i + DC[d];

	    for (uint x = i + DC[d], y = 0; y < X; ++x, ++y)
		R[j].chars[y] = (x < myArraylen) ? input[x] : 0;

	    ++j;
	}
    }

    assert( j == D * M2 );

    std::sort(&R[0], &R[Rsize]);	// sort samples locally

    if (debug && 1)
    {
	// print out sorted suffixes
	std::cout << "Sorted sample X-tuples: i - R[i] - T[R[i]] - pos in R2\n";
 
	for (unsigned int i = 0; i < Rsize; ++i)
	{
	    std::cout << i << " : " << R[i].index << " : " << R[i].str() << "\n";
	}
    }

    // **********************************************************************
    // * select equidistance samples and redistribute sorted DC-tuples

    // select samples
    TupleA<X>* samplebuf = new TupleA<X>[ samplesize ];

    double dist = ( double ) Rsize / samplesize;
    for ( uint i = 0; i < samplesize; i++ )
	samplebuf[i] = R[ int( i * dist ) ];

    // root proc collects all samples
    TupleA<X>* samplebufall;
    if (myproc == ROOT) samplebufall = new TupleA<X>[ nprocs * samplesize ];

    MPI_Gather( samplebuf, samplesize, MPI_TUPLEX, samplebufall, samplesize, MPI_TUPLEX, ROOT, MPI_COMM_WORLD );

    delete[] samplebuf;

    // root proc sorts samples as splitters

    TupleA<X>* splitterbuf = new TupleA<X>[ nprocs ];

    if ( myproc == ROOT )
    {
	std::sort( samplebufall, samplebufall + samplesize * nprocs );

	if (debug && 1)
	{
	    // print out sorted suffixes
	    std::cout << "Available Splitters: i - samplebufall[i]\n";
 
	    for (unsigned int i = 0; i < nprocs * samplesize; ++i)
	    {
		std::cout << i << " : " << samplebufall[i].str() << "\n";
	    }
	}

	for ( int i = 0; i < nprocs; i++ )		// pick splitters
	    splitterbuf[i] = samplebufall[ i * samplesize ];

	if (debug && 1)
	{
	    // print out sorted suffixes
	    std::cout << "Selected Splitters: i - splitterbuf[i]\n";
 
	    for (int i = 0; i < nprocs; ++i)
	    {
		std::cout << i << " : " << splitterbuf[i].str() << "\n";
	    }
	}

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

	TupleA<X>* it = std::lower_bound(&R[0], &R[Rsize], splitterbuf[i]);

	splitterpos[i] = it - &R[0];
    }
    splitterpos[ nprocs ] = Rsize;

    if (debug && 1)
    {
	// print out sorted suffixes
	std::cout << "Splitter Positions: i - splitterpos[i]\n";
 
	for (int i = 0; i < nprocs+1; ++i)
	{
	    std::cout << i << " : " << splitterpos[i] << "\n";
	}
    }

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

    TupleA<X>* recvbuf;
    if ( recvoff[ nprocs ] )
	recvbuf = new TupleA<X>[ recvoff[ nprocs ] ];
    else
	recvbuf = NULL;

    MPI_Alltoallv( R, sendcnt, sendoff, MPI_TUPLEX, recvbuf, recvcnt, recvoff, MPI_TUPLEX, MPI_COMM_WORLD );

    std::swap(R, recvbuf);
    Rsize = recvoff[ nprocs ];

    delete [] recvbuf;
    delete [] sendcnt;
    delete [] recvcnt;
    delete [] sendoff;

    // recvbuf contain DC-tuples in unsorted order

    if (debug && 1)
    {
	// print out sorted suffixes
	std::cout << "Redistributed DC-tuples: i - R[i]\n";
 
	for (unsigned int i = 0; i < Rsize; ++i)
	{
	    std::cout << i << " : " << R[i].str() << "\n";
	}
    }

    Pair* P;
    int recursion = false;

    if ( Rsize )
    {
	/** merge */
	TupleA<X> * helparray4 = new TupleA<X>[ Rsize ];

	// binary merge sort!
	mergesort( R, helparray4, 0, nprocs, recvoff, cmpTupleXLeq<X> );

	// really no need for all this sorting: just name using loser tree!

#if 0
	{
	    Debug5(  //Lastverteilung
		double MergeSA12= MPI_Wtime() - starttime;
		double maxMergeSA12;
		double minMergeSA12;
		double maxAlltoallSA12;
		double minAlltoallSA12;
		uint max;
		uint min;
		Debug6(std::cout << myproc << " first MPI_Alltoallv: before "<< n12<<" after " << recvoff[ nprocs ]<<" mergetime "<<MergeSA12 << std::endl;
		       for(int i=0;i<nprocs;i++)std::cout<< recvcnt[i]<<" r "; std::cout<<std::endl; )
		MPI_Reduce(&recvoff[nprocs], &max, 1, MPI_UNSIGNED, MPI_MAX, ROOT, MPI_COMM_WORLD);
		MPI_Reduce(&recvoff[nprocs], &min, 1, MPI_UNSIGNED, MPI_MIN, ROOT, MPI_COMM_WORLD);
		MPI_Reduce(&AlltoallSA12, &maxAlltoallSA12, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
		MPI_Reduce(&AlltoallSA12, &minAlltoallSA12, 1, MPI_DOUBLE, MPI_MIN, ROOT, MPI_COMM_WORLD);
		MPI_Reduce(&MergeSA12, &maxMergeSA12, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
		MPI_Reduce(&MergeSA12, &minMergeSA12, 1, MPI_DOUBLE, MPI_MIN, ROOT, MPI_COMM_WORLD);

		uint space=11;
/**
   min -- max -- max-min -- %-Diff -- %-Tuple -- merge-min -- merge-max -- %-Diff -- alltoall-min  -- alltoall-max
*/

		if (myproc==ROOT) {std::cout<< setw( space )<<min<< setw( space )<<max<< setw( space )<<max-min<< setw( space );
		    std::cout.precision(4);
		    std::cout << (double)(max-min)/min*100<< setw( space )<<100*k/(2*localArraylen/3.0)<< setw( space )<< minMergeSA12<< setw( space )<<maxMergeSA12<<setw (space)<<(maxMergeSA12-minMergeSA12)/minMergeSA12*100 << setw( space )<<minAlltoallSA12<< setw( space )<< maxAlltoallSA12<<"  SA12"<<std::endl;}
		);
	}
#endif
	delete[] helparray4;

	/** Lexicographical naming   */
	P = new Pair[ Rsize ];

	recursion = namelex( R, P, Rsize ) ;
    }
    else {
	std::cout << myproc << " hat nichts bekommen" << std::endl;
	P = NULL;
    }

    MPI_Bcast( &recursion, 1, MPI_INT, nprocs - 1, MPI_COMM_WORLD );

    /** recursion ! */
    if ( recursion ) {
#if 0
//      if(myproc==ROOT) std::cout<<"Rekursion"<<std::endl;
	uint names = (uint)(2*((double)N /3.0)) + 1;
	uint newlocalArraylen = ( names + nprocs - 1 ) / nprocs;
	uint arraylen=recvoff[nprocs];
	uint* tmparray= new uint[arraylen];
	Pair* sendarray= new Pair[arraylen];
	uint* t = new uint[nprocs];

	int* sendcnt= new int[nprocs];
	int* sendoff= new int[nprocs+1];
	int* recvcnt= new int[nprocs];

	for (int i=0 ; i<nprocs ; i++) {
	    t[i] = (i+1) * newlocalArraylen;
	    sendcnt[i] = 0;
	}
	uint half_names = names / 2;

	for (uint i = 0 ; i < arraylen; i++) {
	    if (P[i].index%3==2) {
		tmparray[i] = findPosSATest(t, half_names + P[ i ].index / 3 ,nprocs, sendcnt);
	    }
	    else {
		tmparray[i] = findPosSATest(t, P[ i ].index / 3 ,nprocs, sendcnt);
	    }
	}
	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

	MPI_Alltoall( sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT , MPI_COMM_WORLD );

	Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime );

	uint* index= new uint[nprocs];

	sendoff[0] = recvoff[0] = 0;

	for (int i = 0 ; i < nprocs ; i++) {
	    index[i] = 0;
	    sendoff[i+1] = sendcnt[i] + sendoff[i];
	    recvoff[i+1] = recvcnt[i] + recvoff[i];
	}

	for (uint i = 0; i < arraylen; i++)
	    sendarray[sendoff[tmparray[i]] + index[tmparray[i]]++]=P[i];

	delete[] index;
	delete[] P;
	delete[] tmparray;

	Pair* in = new Pair[recvoff[nprocs]];
	Pair* recvBufPair= new Pair[recvoff[nprocs]+2];

	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    starttime=MPI_Wtime() );

	MPI_Alltoallv( sendarray, sendcnt, sendoff, MPI_PAIR, in, recvcnt, recvoff, MPI_PAIR, MPI_COMM_WORLD );

	Debug3( alltoallvTime+= MPI_Wtime()-starttime );

	delete[] sendarray;

	Debug2( starttime=MPI_Wtime() );

	sortPairIndex(in, recvBufPair,recvoff[nprocs], t, names);

	Debug2( permuteTime+=MPI_Wtime()-starttime );

	delete[] t;
	delete[] in;

	/** Sehr unwahrscheinlicher Fall: die ersten zwei Elemente werden    an das PE[nprocs-2] gegeben jedoch hat
	    das letzte PE nur ein Element.  */
	if ( myproc == nprocs - 1 && recvoff[ nprocs ] < 2 ) {
	    Debug0( std::cout << "Sehr dummer unwahrscheinlicher Fall eingetreten" << std::endl; );
	    recvBufPair[ recvoff[ nprocs ] ].name = MIN_INT;
	    recvBufPair[ recvoff[ nprocs ] ].index = recvBufPair[ recvoff[ nprocs ] - 1 ].index + 3;
	}

	/** der recvBufPair des letzten PEs enthält am Ende ungültige Zeichen */
	uint* newInBuffer = new uint[ newlocalArraylen + 2 ];
	for ( int i = 0;i < recvoff[nprocs]; i++ )
	    newInBuffer[ i ] = recvBufPair[ i ].name;


	/**Jedes PE muss seinem linken Nachbarn noch die ersten zwei Elemente schicken. Die letzte PE bekommt diese von PE 0 und mus diese danach mit 0 und MAX_INT überschreiben*/
	Pair temp[ 2 ];
	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

	MPI_Sendrecv( recvBufPair, 2, MPI_PAIR, ( myproc - 1 + nprocs ) % nprocs, MSGTAG, temp, 2, MPI_PAIR, ( myproc + 1 ) % nprocs , MSGTAG, MPI_COMM_WORLD, &status );

	Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime );

	delete[] recvBufPair;

	if ( myproc == nprocs - 1 ) { //letzte PE hat naturgemäß zu wenig Elemente
	    newInBuffer[ recvoff[nprocs] ] = 0;
	    for ( uint i = recvoff[nprocs];i < newlocalArraylen+2; i++ )
		newInBuffer[ i ] = MAX_INT;
	}
	else if ((uint)recvoff[nprocs]!=newlocalArraylen) { //vorletzte PE kann zu wenig Elemente haben?!? Führt vorher schon zu Fehler. s. Rekursion start.
	    newInBuffer[ recvoff[nprocs] ] = 0;
	    for ( uint i = recvoff[nprocs];i < newlocalArraylen+2; i++ )
		newInBuffer[ i ] = MAX_INT;
	}
	else { //alle anderen PEs
	    newInBuffer[ newlocalArraylen ] = temp[ 0 ].name;
	    newInBuffer[ newlocalArraylen + 1 ] = temp[ 1 ].name;
	}


	uint* sa12rec;
	uint sa12len = 0;

	// cerr<<myproc<<"alles ok"<<std::endl;

	Debug2( rekursionStart[countRek++]= MPI_Wtime() );

	sa12rec = dc3( newInBuffer, names, &sa12len, localArraylen );

// cerr<<myproc<<"alles ok"<<std::endl;

	Debug2( rekursionEnd[--countRek]=   MPI_Wtime() - rekursionStart[countRek] );

	uint* n1all = new uint[ nprocs ];
	uint* n2all = new uint[ nprocs ];

	uint salenall;

	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

	MPI_Allgather( &( n[ 1 ] ), 1, MPI_UNSIGNED, n1all, 1, MPI_UNSIGNED, MPI_COMM_WORLD );
	MPI_Allgather( &( n[ 2 ] ), 1, MPI_UNSIGNED, n2all, 1, MPI_UNSIGNED, MPI_COMM_WORLD );
	MPI_Scan(&sa12len, &salenall, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );

	Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime );

	uint* n12all = new uint[ 2 * nprocs ];

	n12all[ 0 ] = n1all[ 0 ];
	for ( int i = 1; i < nprocs ;i++ )
	    n12all[ i ] = n1all[ i ] + n12all[ i - 1 ];
	for ( int i = nprocs ; i < 2*nprocs; i++ )
	    n12all[ i ] = n2all[ i - nprocs ] + n12all[ i - 1 ];
	n12all[ 2 * nprocs - 1 ] ++; // die letzte 'Grenze wird um eins erhöht, damit das Element dem letzten PE gehört

	delete[] n1all;
	delete[] n2all;

	for ( int i = 0 ; i < nprocs ; i++ )  sendcnt[ i ] = 0;
	P = new Pair[ sa12len ];
	salenall=salenall-sa12len+1;

	for ( uint i = 0; i < sa12len; i++ ) {
	    P[ i ].name = i + salenall;
	    P[ i ].index = sa12rec[ i ];
	    sa12rec[ i ] = findPosSA( n12all, sa12rec[ i ], nprocs, sendcnt );

//          sendcnt[sa12rec[i]]++;
	    Debug0( if( sa12rec[ i ] >= (uint)nprocs )
			std::cout << myproc << " Hier stimmt was nicht " << sa12rec[ i ] << " sollte kleiner sein als #PEs" << std::endl;
		);
	};
//      delete[] salenAll;
	sendarray = new Pair[ sa12len ];
	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

	MPI_Alltoall( sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT, MPI_COMM_WORLD );

	Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime );

	sendoff[ 0 ] = 0;
	recvoff[ 0 ] = 0;

	index = new uint[nprocs];

	for ( int i = 0 ; i < nprocs; i++ ) {
	    index[ i ] = 0;
	    sendoff[ i + 1 ] = sendcnt[ i ] + sendoff[ i ];
	    recvoff[ i + 1 ] = recvoff[ i ] + recvcnt[ i ];
	}

	for ( uint i = 0; i < sa12len; i++ )
	    sendarray[ sendoff[ sa12rec[ i ] ] + index[ sa12rec[ i ] ] ++ ] = P[ i ];

	delete[] index;
	delete[] sa12rec;
	delete[] P;

	Pair* recvarray = new Pair[ recvoff[ nprocs ] + 2 ];

	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)     starttime = MPI_Wtime() );

	MPI_Alltoallv( sendarray, sendcnt, sendoff, MPI_PAIR, recvarray, recvcnt, recvoff, MPI_PAIR, MPI_COMM_WORLD );

	Debug3( alltoallvTime += MPI_Wtime() - starttime );

	delete[] sendarray;
	uint half = ( localArraylen + 2 ) / 3 + 1;
	P = new Pair[ 2 * half ];

	Debug2( starttime = MPI_Wtime() );

	sortPair( recvarray, P, recvoff[ nprocs ], n12all, half );

	Debug2( permuteTime += MPI_Wtime() - starttime );

	delete[] recvarray;
	delete[] n12all;

	Pair* sendTmp = new Pair[ 2 ];
	Pair* recvTmp = new Pair[ 2 ];

	sendTmp[ 0 ] = P[ 0 ];
	sendTmp[ 1 ] = P[ half ];

	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

	MPI_Sendrecv( sendTmp, 2, MPI_PAIR, ( myproc - 1 + nprocs ) % nprocs, MSGTAG, recvTmp, 2, MPI_PAIR, ( myproc + 1 ) % nprocs , MSGTAG, MPI_COMM_WORLD, &status );

	Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime );

	delete[] sendTmp;

	P[ n[ 1 ] ] = recvTmp[ 0 ];
	P[ half + n[ 2 ] ] = recvTmp[ 1 ];

	delete[] recvTmp;
	delete[] sendcnt; //Hier richtig?
	delete[] sendoff;
	delete[] recvcnt;
	delete[] recvoff;

	uint* satmp = sortS0S1S2( input, P, localArraylen, n, imod3, salen, half );
	return satmp;
	/** end recursion */
#endif
    }
    else {
	if ( myproc == ROOT )
	    std::cout << "---------------------   keine  Recursion---------------- " << localArraylen << std::endl;

	/** Samplesort P nach Index */
	std::sort( P, P + Rsize );

/** testing */
	uint* splitterpos = new uint[nprocs+1];
	int* sendcnt = new int[nprocs];
	int* sendoff = new int[nprocs+1];
	int* recvcnt = new int[nprocs];
//      int* recvoff = new int[nprocs+1];

/** end *****/
	splitterpos[ 0 ] = 0;
	Pair ptemp;
	ptemp.name=0;
	for ( int i = 0; i < nprocs - 1; i++ ) {
	    ptemp.index = (i+1) * localArraylen;
	    splitterpos[ i + 1 ] = findPos( P, ptemp , recvoff [ nprocs ], cmpGreaterIndex );
	}
	splitterpos[ nprocs ] = recvoff[ nprocs ];

	for ( int i = 0; i < nprocs; i++ )
	    sendcnt[ i ] = splitterpos[ i + 1 ] > splitterpos[ i ] ? splitterpos[ i + 1 ] - splitterpos[ i ] : 0;

	MPI_Alltoall( sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT , MPI_COMM_WORLD );

	sendoff[ 0 ] = recvoff[ 0 ] = 0;
	for ( int i = 1; i <= nprocs; i++ ) {
	    sendoff[ i ] = sendoff[ i - 1 ] + sendcnt[ i - 1 ];
	    recvoff[ i ] = recvoff[ i - 1 ] + recvcnt[ i - 1 ];
	}
	recvoff[ nprocs ] = recvoff[ nprocs - 1 ] + recvcnt[ nprocs - 1 ];
	Pair* recvBufPair = new Pair[ recvoff[ nprocs ] + 2 ];

	MPI_Alltoallv( P, sendcnt, sendoff, MPI_PAIR, recvBufPair, recvcnt, recvoff, MPI_PAIR, MPI_COMM_WORLD );

	delete[] P;

	Pair* helparray2 = new Pair[ recvoff[ nprocs ] ];
	mergesort( recvBufPair, helparray2, 0, nprocs, recvoff, cmpLessIndex );
	delete[] helparray2;

	if ( myproc == nprocs - 1 ) {
	    if ( recvoff[ nprocs ] == 1 ) {       //kann nur bei letztem PE auftreten?
		if ( recvBufPair[ recvoff[ nprocs ] - 1 ].index % 3 == 1 )
		    recvBufPair[ recvoff[ nprocs ] ].index = recvBufPair[ recvoff[ nprocs ] - 1 ].index + 1;
		else
		    recvBufPair[ recvoff[ nprocs ] ].index = recvBufPair[ recvoff[ nprocs ] - 1 ].index + 2;
	    }
	}

	/**jedes PE erhält zwei Elemente am Ende mehr*/
	Pair temp[ 2 ];
	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

	MPI_Sendrecv( recvBufPair, 2, MPI_PAIR, ( myproc - 1 + nprocs ) % nprocs, MSGTAG, &temp, 2, MPI_PAIR, ( myproc + 1 ) % nprocs , MSGTAG, MPI_COMM_WORLD, &status );

	Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime);

	if ( myproc == nprocs - 1 ) {
	    recvBufPair[ recvoff[ nprocs ] ].name = MAX_INT;
	    recvBufPair[ recvoff[ nprocs ] + 1 ].name = MAX_INT;
	    if ( recvBufPair[ recvoff[ nprocs ] - 1 ].index % 3 == 1 ) {
		recvBufPair[ recvoff[ nprocs ] ].index = recvBufPair[ recvoff[ nprocs ] - 1 ].index + 1;
		recvBufPair[ recvoff[ nprocs ] + 1 ].index = recvBufPair[ recvoff[ nprocs ] - 1 ].index + 3;
	    } else {
		recvBufPair[ recvoff[ nprocs ] ].index = recvBufPair[ recvoff[ nprocs ] - 1 ].index + 2;
		recvBufPair[ recvoff[ nprocs ] + 1 ].index = recvBufPair[ recvoff[ nprocs ] - 1 ].index + 3;
	    }
	}
	else {
	    recvBufPair[ recvoff[ nprocs ] ] = temp[ 0 ];
	    recvBufPair[ recvoff[ nprocs ] + 1 ] = temp[ 1 ];
	    if ( myproc == nprocs - 2 ) { //wenn das letzte PE nur mod 2 Tupel hat,...
		if ( recvBufPair[ recvoff[ nprocs ] ].index == recvBufPair[ recvoff[ nprocs ] + 1 ].index )
		    recvBufPair[ recvoff[ nprocs ] ].index = recvBufPair[ recvoff[ nprocs ] ].index - 1;
	    }
	}

//      std::cout<<"CPU ["<<myproc<<"] nach rekursion"<<std::endl;
	Debug2( starttime = MPI_Wtime() );

	std::sort( recvBufPair, recvBufPair + recvoff[ nprocs ] + 2, cmpIndexModDiv );

	Debug2( sortingTime += MPI_Wtime() - starttime );
	Debug6(stlSortTime=MPI_Wtime() - starttime;
//  MPI_Reduce(&MergeSA12, &stlSortTime, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
//  MPI_Reduce(&MergeSA12, &stlSortTime, 1, MPI_DOUBLE, MPI_MIN, ROOT, MPI_COMM_WORLD);
	       std::cout<<myproc<<"   "<<stlSortTime<<" sec "<<recvoff[ nprocs ] + 2<<" Pair -- sort Pair %/ rekursion "<<std::endl;
	    );

#if 0

	uint half= n[1]+1;

	delete[] sendcnt; //Hier richtig?
	delete[] sendoff;
	delete[] recvcnt;
	delete[] recvoff;

	uint* satmp = sortS0S1S2( input, recvBufPair, localArraylen, n, imod3, salen, half );
	return satmp;
#endif

	/**Samplesort P end*/
    }
}

/**
 * Global naming of all Quadruple: same tuple => same name
 * @param in Quadruple*
 * @param P Pair* out-param
 * @param arraylen number of Quadruples
 */
inline int namelex( TupleA<X>* in, Pair* P, uint arraylen )
{
    /** naming with local names */
    TupleA<X> temp;

    MPI_Sendrecv( &( in[ arraylen - 1 ] ), 1, MPI_TUPLEX, ( myproc + 1) % nprocs, MSGTAG,
		  &temp, 1, MPI_TUPLEX, ( myproc - 1 +nprocs) % nprocs , MSGTAG, MPI_COMM_WORLD, &status );

    uint name = 0;
    for ( uint i = 0; i < arraylen; i++ ) {
	if ( !( in[ i ] == temp ) ) {
	    name++;
	    temp = in[ i ];
	}
	P[ i ].name = name;
	P[ i ].index = in[ i ].index;
    }
    delete[] in;

    /** renaming with global names */
    uint namesglob;

    MPI_Scan( &name, &namesglob, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );

    for ( uint i = 0; i < arraylen; i++ )
	P[i].name += ( namesglob - name );

    // todo: determine recursion?

    return false;
}

void writesa( uint* salen, uint* suffixarray, char* filenameOut )
{
    /** Suffixarry in Datei schreiben */
    MPI_Request request;

    uint* saLenAll=new uint[ nprocs ];
    MPI_Gather( salen, 1, MPI_UNSIGNED, saLenAll, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

    if ( myproc == ROOT ) {
	int fd=open(filenameOut, O_WRONLY | O_CREAT | O_TRUNC , S_IRUSR | S_IWUSR);
	if(fd==-1) {std::cout << "Fehler bei der Ausgabe"<<std::endl; return ;}
	FILE* fp = fdopen(fd,"w");
	if (!fp)std::cout <<"fehler beim öffnen"<< std::endl;
	uint maxLen = 0;
	for ( int i = 0;i < nprocs;i++ ) if ( maxLen < saLenAll[ i ] )   maxLen = saLenAll[ i ];

	fwrite(&suffixarray[1],sizeof(uint),saLenAll[0]-1,fp); //abschneiden der 'null' am anfang
	std::cout << "Von 0 Daten geschrieben" << std::endl;
	delete[] suffixarray;
	uint* outbuffer = new uint[ maxLen ];
	for ( int k = 1;k < nprocs;k++ ) {
	    MPI_Recv( outbuffer, saLenAll[ k ], MPI_UNSIGNED, k, MSGTAG, MPI_COMM_WORLD, &status );
			std::cout << "salen von " << k << " ist " << saLenAll[ k ] << std::endl;
	    fwrite(outbuffer,sizeof(uint),saLenAll[k],fp);
//          Debug6(for (uint i=0 ; i<saLenAll[k] ; i++)
//                      std::cout<< outbuffer[i]<<" ";
//                  std::cout << std::endl;)

	    std::cout << "Von " << k << " Daten geschrieben" << std::endl;
	}
	delete[] saLenAll;
	delete[] outbuffer;
	if(fclose(fp)!=0) std::cout <<"Fehler beim Schließen der Datei";

    }
	else {
	MPI_Isend( suffixarray, *salen , MPI_UNSIGNED, 0, MSGTAG , MPI_COMM_WORLD , &request );
	MPI_Wait( &request, &status );
	delete[] suffixarray;
	std::cout << myproc << " ist fertig" << std::endl;
    }
}

Debug2(
void writeTimes( uint N )
{
    ofstream outfile( "timing", ios::app );
    assert( outfile.good() );
	std::cout.precision(4);
/*  if ( myproc == ROOT ) {
	outfile << nprocs << " CPUs und Dateigroesse ist " << N << " Bytes" << std::endl;
	ofstream sampleTime("sampling",ios::app);
	sampleTime<< nprocs << " CPUs und Dateigroesse ist " << N << " Bytes" << std::endl;
	int i=0;
	sampleTime<<"Sampling 1: ";
	while(samplingTime1[i])
	sampleTime<<samplingTime1[i++]<<" ";
	sampleTime<<std::endl<<"Sampling 2: ";
	i=0;
	while(samplingTime2[i])
	sampleTime<<samplingTime2[i++]<<" ";
	sampleTime<<std::endl;
    }*/
	Debug2(
		if ( myproc == ROOT ) {
			outfile << nprocs << " CPUs und Dateigroesse ist " << N << " Bytes" << std::endl;
			ofstream rekursion("rekursion",ios::app);
			rekursion<< nprocs << " CPUs und Dateigroesse ist " << N << " Bytes" << std::endl;
			int i=0;
			rekursion<<"Rekursionszeit: ";
			while(rekursionEnd[i++]>0 && i<maxRekursion)        rekursion<<rekursionEnd[i-1]-rekursionEnd[i]<<" ";
			i=0;
			rekursion<<std::endl<<"Rekursionszeit kummuliert: ";
			while(rekursionEnd[i]>0 && i<maxRekursion)  rekursion<<rekursionEnd[i++]<<" ";
			rekursion<<std::endl<<"Verhältnis Namen zu Tupel: ";
			i=0;
			while(name2tuple[i]>0 && i<maxRekursion)    rekursion<<name2tuple[i++]<<" ";
			rekursion<<std::endl;
			rekursion.close();
		}
		);

	outfile << setw( 4 ) << myproc;
    outfile << setprecision( 5 );
    outfile.setf( ios::fmtflags( 0 ), ios::floatfield );
    int space = 13;
    outfile << setw( space ) << dc3FinishTime;
    outfile << setw( space ) << sortingTime ;
    outfile << setw( space ) << mergeTime;
    outfile << setw( space ) << permuteTime;
    outfile << setw( space ) << alltoallvTime;
    outfile << setw( space ) << mpiComTime;
    outfile << setw( space ) << dc3FinishTime - ( permuteTime + sortingTime + mergeTime + alltoallvTime + mpiComTime );
    outfile << setw( space ) << sortS0S1S2TimeEnd;
    outfile.setf( ios::scientific, ios::floatfield );
    outfile << setw( space ) << dc3FinishTime / N;
    outfile << std::endl;
    outfile.close();
    space = 11;
    std::cout << setw( 4 ) << myproc;
    std::cout << setw( space ) << dc3FinishTime;
    std::cout << setw( space ) << sortingTime;
    std::cout << setw( space ) << mergeTime;
    std::cout << setw( space ) << permuteTime;
    std::cout << setw( space ) << alltoallvTime;
    std::cout << setw( space ) << mpiComTime;
    std::cout << setw( space ) << dc3FinishTime - ( sortingTime + mergeTime + permuteTime + alltoallvTime + mpiComTime);
    std::cout.setf( ios::scientific, ios::floatfield );
    std::cout << setw( space ) << dc3FinishTime / N;
    std::cout << std::endl;
}
)

void getTuple(uint* inbuffer, Pair* recvBufPair, Quintuple* S0, Quadruple* S1, Quintuple* S2, uint localArraylen, uint* n,  uint* imod3, uint& half)
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
    uint displacement = myproc * localArraylen;
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
uint* sortS0S1S2( uint* inbuffer, Pair* recvBufPair, uint localArraylen, uint* n, uint* imod3, uint* salen, uint& half )
{
    Debug2( sortS0S1S2TimeStart = MPI_Wtime() );
//  uint k=(uint) ceil(samplefactor*(sqrt((double)localArraylen)/300));
    uint k = (uint)(sqrt((double) localArraylen / nprocs) *samplefactor/3.0);
    if(k>=localArraylen/3){std::cout<<"k="<<k<<" zu groß, jetzt k="<<localArraylen/3-1<<std::endl; k=localArraylen/3-1;}

//  uint k1 = (uint)sqrt( localArraylen / (3.0 * nprocs) *samplefactor);
//   if ( myproc == ROOT ) std::cout<<" k "<<3*k<<" % "<<3.0*k/localArraylen*100.0<<" Tuple " <<localArraylen<<std::endl;
    //    std::cout << myproc << " 0:" << imod3[ 0 ] << " 1:" << imod3[ 1 ] << " 2:" << imod3[ 2 ] << " localArraylen " <<  localArraylen << std::endl;
    Quintuple* S0 = new Quintuple[ n[ 0 ] ];
    Quadruple* S1 = new Quadruple[ n[ 1 ] ];
    Quintuple* S2 = new Quintuple[ n[ 2 ] ];
    getTuple( inbuffer,  recvBufPair,  S0,  S1,  S2, localArraylen,  n,  imod3,  half);

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
	//if(samplebufS1[i].index>MOD0+(localArraylen+1)/3)std::cout<<"Fehler samplebufS1[i].index>MOD0+n[i] "<< samplebufS1[i].index<<" "<< S1[i].index/3 <<" "<<MOD0<<" "<<(localArraylen+1)/3 <<std::endl;
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
//  Debug6(std::cout << myproc << " SortS0S1S2 Alltoallv: before "<< localArraylen <<" after " << *salen << std::endl; )
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
	    std::cout << (double)(max-min)/min*100<< setw( space )<<(300.0*k)/localArraylen<< setw( space )<<minsortAll<< setw( space )<< maxsortAll<<setw (space)<<(maxsortAll-minsortAll)/minsortAll*100<<" Alltoallv all Tuple"<<std::endl;}
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
    else if (nprocs % 2 != 0)
    {
	std::cerr << "Error: requires an even number of  MPI processors." << std::endl;
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

    uint localArraylen = (uint)ceil( (double)filelength / nprocs );	// divide by processors rounding up
    localArraylen += X - localArraylen % X;				// round up to nearest multiple of X

    assert( localArraylen * nprocs >= filelength );

    if ( myproc == ROOT )
    {
	std::cout << "Total input size = " << filelength << " bytes. localArraylen = " << localArraylen << std::endl;
    }

    // Verteilung der Daten nach inbuffer<- enthält lokal benötigte Zeichen +2 wegen Überschneidung und zum Erstellen der S0-S1 Tupel

    // Send data portions to _other_ n-1 processors -- including an overlap of +(X-1) characters for the final tuples

    int overlap = X-1;
    char* inputChars = new char[ localArraylen + overlap ];

    if ( myproc == ROOT )
    {
	std::ifstream infile( filename, std::ios::binary );

	// Einlesen der Elemente für die Prozessoren 1 bis n-2
	for ( int p = 1; p < nprocs; p++ )
	{
	    infile.seekg( p * localArraylen, std::ios::beg );

	    uint readsize = (p != nprocs-1) ? localArraylen + overlap : filelength - (p * localArraylen);

	    std::cout << "Read portion " << p << " from pos " << p * localArraylen << " of length " << readsize << std::endl;

	    infile.read( inputChars, readsize );
	    MPI_Send( inputChars, readsize, MPI_CHAR, p, MSGTAG, MPI_COMM_WORLD );
	}

	if (!infile.good()) {
	    perror("Error reading file.");
	    return -1;
	}

	// Einlesen der Elemente für Prozessor ROOT

	std::cout << "Read portion 0 from pos 0 of length " << localArraylen + overlap << std::endl;

	infile.seekg( 0, std::ios::beg );
	infile.read( inputChars, localArraylen + overlap );

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
	MPI_Recv( inputChars, localArraylen + overlap, MPI_CHAR, ROOT, MSGTAG, MPI_COMM_WORLD, &status );
    }

    // convert input characters into uint data

    uint* input = new uint[ localArraylen + overlap ];

    for ( unsigned int i = 0; i < localArraylen + overlap; i++ ) {
	input[i] = inputChars[i];
    }

    delete[] inputChars;

    MPI_Barrier( MPI_COMM_WORLD );

    // **********************************************************************
    // * Construct suffix array recursively
    
    Debug2( rekursionStart[countRek++]= MPI_Wtime() );
    Debug1( dc3StartTime = MPI_Wtime() );

    uint    salen;
    uint*   suffixarray = dc3( input, filelength, &salen, localArraylen );

    Debug1( dc3FinishTime = MPI_Wtime() - dc3StartTime );
    Debug2( rekursionEnd[--countRek]=   MPI_Wtime() - rekursionStart[countRek] );

    Debug1( if ( myproc == ROOT ) std::cout << "Zeit für Erstellen von Suffixarry = "<<dc3FinishTime<<std::endl );

    Debug2(
	if (myproc==ROOT) //Zeitmessung
	    MPI_Recv( name2tuple, 35, MPI_DOUBLE, nprocs - 1, MSGTAG, MPI_COMM_WORLD, &status );
	if (myproc==nprocs-1)
	    MPI_Send( name2tuple, 35, MPI_DOUBLE, 0, MSGTAG , MPI_COMM_WORLD );
	int token = 0;
	if ( myproc == ROOT ) {
	    std::cout << " CPU    Gesamtzeit    Sort    Merge     Permute   MPI_Alltoallv  MPI_Com   Rest   Zeit/Zeichen " << std::endl;
	    token = 1;
	}
	if ( !token )
	    MPI_Recv( &token, 1, MPI_INT, myproc - 1, MSGTAG, MPI_COMM_WORLD, &status );
	writeTimes( filelength );
	if ( myproc != nprocs - 1 )
	    MPI_Send( &token, 1, MPI_INT, ( myproc + 1 ) % nprocs, MSGTAG , MPI_COMM_WORLD );
	);

    if ( argc > 3 ) writesa( &salen, suffixarray, argv[ 3 ] );

    MPI_Finalize();
    return 0;
}
