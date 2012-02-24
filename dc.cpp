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

typedef unsigned int	uint;

#define MAX_INT ~1U
#define MAX_INT_SAMPLE ~0U
//#define MOD2_INDEX ~1UL // ~1UL % 3 = 4294967294 % 3 = 2  max mod 2 Index for UNSIGNED INTEGER
#define MIN_INT 0UL
#define MSGTAG 42
#define MOD0 ~0U/3
#define MOD1 2*(~0U/3)


uint* dc3( uint* inbuffer, uint filelength, uint* salen );
uint* sortS0S1S2( uint* inbuffer, Pair* P, uint locArraylen, uint* n, uint* imod3, uint* salen, uint& half);
inline void getTuple(uint* inbuffer, Pair* recvBufPair, Quintuple* S0, Quadruple* S1, Quintuple* S2, uint locArraylen, uint* n,  uint* imod3, uint& half);

inline void namelex( Quadruple* in, Pair* P, uint arraylen );

void writesa( uint* salen, uint* suffixarray, char* filenameOut );
void writeTimes( uint filelength );

int myproc, size;
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
	    if ( in[ i ].index >= n12all[ size - 1 ] ) {
		//          assert( (in[ i ].index - n12all[ size -1 ] + half) < 2*half);
		//          cout << "in[i].index -n12all[size-1] +half " << in[ i ].index - n12all[ size-1 ] + half << endl;
		out[ in[ i ].index - n12all[ size - 1 ] + half ] = in[ i ];
	    } else {
		//          cout << "in[i].index " << in[ i ].index <<" half "<<2*half<< endl;
		//          if(in[ i ].index > 2*half) assert(0);
		//          if ((int)in[ i ].index < 0 ) assert(0);
		out[ in[ i ].index ] = in[ i ];
	    }
	}
    else
	for ( uint i = 0 ; i < length; i++ ) {
	    if ( in[ i ].index >= n12all[ size - 1 ] ) {
		//          assert( (in[ i ].index - n12all[myproc+ size -1 ] + half )< 2*half);
		//          cout << "in[i].index -n12all[myproc-1+size] +half " << in[ i ].index - n12all[ myproc -1 + size ] + half << endl;
		out[ in[ i ].index - n12all[ myproc - 1 + size ] + half ] = in[ i ];
	    } else {
		//          cout << "in[i].index -n12all[myproc-1] " << in[ i ].index - n12all[ myproc-1 ] << endl;
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
MPI_Datatype MPI_QUINTUPLE;
MPI_Datatype MPI_QUADRUPLE;
MPI_Datatype MPI_PAIR;

MPI_Status status;

int samplefactor;

/**used for testing*/
int ausgabeCPU;
/** end testing */

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

void mpi_get_init()
{
    MPI_Comm_rank( MPI_COMM_WORLD, &myproc );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    assert( size > 1 );
    assert( ( size & 0x1 ) == 0 );


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
}


/**
 *
 * @param inbuffer
 * @param filelength
 * @param salen
 * @return suffixarray
 */
uint* dc3( uint* inbuffer, uint filelength, uint* salen )
{
	// Debug6( cout << "dc3 gestartet" << endl; )

    /** Tupel einlesen*/

    //Anzahl mod 1, mod 2, mod 3 Tupel berechnen. Ist für jedes PE +-1 Element gleich
    uint locArraylen = ( filelength + size - 1 ) / size ;
    uint k = (uint)sqrt(2.0 * locArraylen / 3.0 / size) * samplefactor;

// Debug5(  if ( myproc == ROOT ) cout << "samplefactor=" << k <<" |Tuple|="<<2*(locArraylen/3)<<" in %="<<100*k/(2*locArraylen/3.0)<<endl;)
//  int k=(uint) ceil(samplefactor*(sqrt(2.0*locArraylen/3.0)/100));

    if( k >= 2*(locArraylen/3)) k = 2*(locArraylen/3)-1;

// if ( myproc == ROOT ) cout << "samplefactor=" << k <<" |Tuple|="<<2*(locArraylen/3)<<" in %="<<100*k/(2*locArraylen/3.0)<<endl;

    uint n[ 3 ];
    n[ 0 ] = n[ 1 ] = n[ 2 ] = ( locArraylen ) / 3;
    uint imod3[ 3 ] = {0, 1, 2};
    uint displacement = myproc * locArraylen;

    /* Neun verschiedene Möglichkeiten, da Arraygröße mod 3 und Startwert mod 3 kombiniert werden muss
     * für Startwert mod 3 = 0 ist Arraygröße mod 3 egal.
     * für Startwert mod 3 = 1 ist ein Tupel mehr vorhanden
     * für Startwert mod 3 = 2 sind zwei Tupel mehr vorhanden
     */
    switch ( locArraylen % 3 ) {
	case 0:
		break; //nothing to do
	case 1:
		switch ( displacement % 3 ) {
		case 0:
			n[ 0 ] ++;
			break;
		case 1:
			n[ 1 ] ++;
			imod3[ 0 ] = 2;
			imod3[ 1 ] = 0;
			imod3[ 2 ] = 1;
			break;
		case 2:
			n[ 2 ] ++;
			imod3[ 0 ] = 1;
			imod3[ 1 ] = 2;
			imod3[ 2 ] = 0;
			break;
		}
		break;
	case 2:
		switch ( displacement % 3 ) {
		case 0:
			n[ 0 ] ++;
			n[ 1 ] ++;
			break;
		case 1:
			n[ 1 ] ++;
			n[ 2 ] ++;
			imod3[ 0 ] = 2;
			imod3[ 1 ] = 0;
			imod3[ 2 ] = 1;
			break;
		case 2:
			n[ 0 ] ++;
			n[ 2 ] ++;
			imod3[ 0 ] = 1;
			imod3[ 1 ] = 2;
			imod3[ 2 ] = 0;
			break;
		}
		break;
    }
    if ( myproc == size - 1 ) {
	n[ 0 ] = ( filelength + 2 ) / 3 - ( locArraylen * myproc + 2 ) / 3;
	n[ 1 ] = ( filelength + 1 ) / 3 - ( locArraylen * myproc + 1 ) / 3;
	n[ 2 ] = ( filelength + 0 ) / 3 - ( locArraylen * myproc + 0 ) / 3;
    }

    uint n12 = n[ 1 ] + n[ 2 ];
    // Debug6( if ( myproc == ausgabeCPU ) cout << myproc << " erstellt Quadruple - n[0] " << n[ 0 ] << " n[1] " << n[ 1 ] << " n[2] " << n[ 2 ] << endl; )
    Quadruple* sa12 = new Quadruple[ n12 ];

    {
	uint j1 = 0,j2=n[1];
	/** zwei recht ähnliche Schleifen, mit dem Unterschied, dass eine davon 1 mal länger laufen kann */
	uint x1 = imod3[1], x2 =imod3[2];
	while (j1<n[1] || j2<n12){
	    if(j1<n[1]){// i mod 3 = 1
		sa12[ j1 ].index = x1 + displacement;
		sa12[ j1 ].name[ 0 ] = inbuffer[ x1++ ];
		sa12[ j1 ].name[ 1 ] = inbuffer[ x1++ ];
		sa12[ j1++ ].name[ 2 ] = inbuffer[ x1++ ];
	    }
	    if(j2<n12){// i mod 3 = 2
		sa12[ j2 ].index = x2 + displacement;
		sa12[ j2 ].name[ 0 ] = inbuffer[ x2++ ];
		sa12[ j2 ].name[ 1 ] = inbuffer[ x2++ ];
		sa12[ j2++ ].name[ 2 ] = inbuffer[ x2++ ];
	    }
	}
    }

    /** Samplesort Quadruple start */
	Debug2( starttime = MPI_Wtime() );

	Switch0( std::sort( sa12, sa12 + n12 ) );
	Switch1( std::stable_sort( sa12, sa12 + n12 ) );

	Debug2( sortingTime += MPI_Wtime() - starttime );
	Debug6(	stlSortTime=MPI_Wtime() - starttime;
//  MPI_Reduce(&MergeSA12, &stlSortTime, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
//  MPI_Reduce(&MergeSA12, &stlSortTime, 1, MPI_DOUBLE, MPI_MIN, ROOT, MPI_COMM_WORLD);
			cout<<myproc<<"   "<<stlSortTime<<" sec "<<n12<<" Tuple sort SA12 "<<endl;
		);

	/** Sampling */
	Quadruple* samplebuf = new Quadruple[ k ];
    double dist = ( double ) n12 / k;
    for ( uint i = 0; i < k; i++ )  samplebuf[ i ] = sa12[ int( i * dist ) ];


    /**Proc 0 sammelt Samples ein*/
	Debug2( sampleStartTime = MPI_Wtime() );
	Quadruple* samplerecvbuf=NULL;
    if ( myproc == ROOT )   samplerecvbuf = new Quadruple[ k * size ];
	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );
	MPI_Gather( samplebuf, k, MPI_QUADRUPLE, samplerecvbuf, k, MPI_QUADRUPLE, 0, MPI_COMM_WORLD );
	Debug3( mpiComTime += MPI_Wtime()-mpiComStartTime );
	/**Proc 0 ermittelt Splitter*/
	delete[] samplebuf;

    Quadruple* pivbuf= new Quadruple[ size ];
    if ( myproc == ROOT ) {
	Switch0( sort( samplerecvbuf, samplerecvbuf + k * size ) );
		Switch1( stable_sort( samplerecvbuf, samplerecvbuf + k * size ) );
		for ( int i = 0; i < size; i++ ) pivbuf[ i ] = samplerecvbuf[ i * k ];
	delete[] samplerecvbuf;
    }

    /**Splitter verteilen*/
	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );
	MPI_Bcast( pivbuf, size, MPI_QUADRUPLE, ROOT, MPI_COMM_WORLD );
	Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime );

	Debug2( sampleTime+=MPI_Wtime()-sampleStartTime );

	/** Einordnen der lokal sortierten Elemente anhand von Splitter */
	uint* pivpos= new uint[ size + 1 ];

    pivpos[0] = 0;
    for ( int i = 1; i < size; i++ )
	pivpos[i] = findPos( sa12, pivbuf[i], n12, cmpSplitterGreaterS12 );
    pivpos[ size ] = n12;
    delete[] pivbuf;

    /** Ermitteln der Anzahl der zu versendenden bzw. empfangenden Elemente */
    int* sendcnt;
    int* recvcnt;
    sendcnt = new int[ size ];
    recvcnt = new int[ size ];
    for ( int i = 0; i < size; i++ )
	sendcnt[ i ] = pivpos[ i + 1 ] > pivpos[ i ] ? pivpos[ i + 1 ] - pivpos[ i ] : 0;
    delete[] pivpos;

	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD));    mpiComStartTime = MPI_Wtime() );

	MPI_Alltoall( sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT , MPI_COMM_WORLD );

	Debug3( mpiComTime += MPI_Wtime()-mpiComStartTime );

	int* sendoff;
    int* recvoff;
    sendoff = new int[ size + 1 ];
    recvoff = new int[ size + 1 ];

    sendoff[ 0 ] = recvoff[ 0 ] = 0;
    for ( int i = 1; i < size + 1; i++ ) {
	sendoff[ i ] = sendoff[ i - 1 ] + sendcnt[ i - 1 ];
	recvoff[ i ] = recvoff[ i - 1 ] + recvcnt[ i - 1 ];
    }

    Quadruple* recvbuf;
    if ( recvoff[ size ] )
	recvbuf = new Quadruple[ recvoff[ size ] ];
    else
	recvbuf = NULL;

	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD)); starttime = MPI_Wtime());

	MPI_Alltoallv( sa12, sendcnt, sendoff, MPI_QUADRUPLE, recvbuf, recvcnt, recvoff, MPI_QUADRUPLE, MPI_COMM_WORLD );

	Debug3( alltoallvTime += MPI_Wtime() - starttime );

	delete[] sa12;
    delete[] sendcnt;
    delete[] recvcnt;
    delete[] sendoff;

	Debug5(  //Lastverteilung
		double AlltoallSA12 = MPI_Wtime() - starttime
		);

	Pair* P;

    if ( recvoff[ size ] ) {
		/** merge */
	Quadruple * helparray4 = new Quadruple[ recvoff[ size ] ];
		Debug2( starttime = MPI_Wtime() );

		mergesort( recvbuf, helparray4, 0, size, recvoff, cmpSplitterLeqS12 );

		Debug2( mergeTime += MPI_Wtime() - starttime );

		{
			Debug5(  //Lastverteilung
				double MergeSA12= MPI_Wtime() - starttime;
				double maxMergeSA12;
				double minMergeSA12;
				double maxAlltoallSA12;
				double minAlltoallSA12;
				uint max;
				uint min;
				Debug6(cout << myproc << " first MPI_Alltoallv: before "<< n12<<" after " << recvoff[ size ]<<" mergetime "<<MergeSA12 << endl;
					   for(int i=0;i<size;i++)cout<< recvcnt[i]<<" r "; cout<<endl; )
				MPI_Reduce(&recvoff[size], &max, 1, MPI_UNSIGNED, MPI_MAX, ROOT, MPI_COMM_WORLD);
				MPI_Reduce(&recvoff[size], &min, 1, MPI_UNSIGNED, MPI_MIN, ROOT, MPI_COMM_WORLD);
				MPI_Reduce(&AlltoallSA12, &maxAlltoallSA12, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
				MPI_Reduce(&AlltoallSA12, &minAlltoallSA12, 1, MPI_DOUBLE, MPI_MIN, ROOT, MPI_COMM_WORLD);
				MPI_Reduce(&MergeSA12, &maxMergeSA12, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
				MPI_Reduce(&MergeSA12, &minMergeSA12, 1, MPI_DOUBLE, MPI_MIN, ROOT, MPI_COMM_WORLD);

				uint space=11;
/**
   min -- max -- max-min -- %-Diff -- %-Tuple -- merge-min -- merge-max -- %-Diff -- alltoall-min  -- alltoall-max
*/

				if (myproc==ROOT) {cout<< setw( space )<<min<< setw( space )<<max<< setw( space )<<max-min<< setw( space );
					cout.precision(4);
					cout << (double)(max-min)/min*100<< setw( space )<<100*k/(2*locArraylen/3.0)<< setw( space )<< minMergeSA12<< setw( space )<<maxMergeSA12<<setw (space)<<(maxMergeSA12-minMergeSA12)/minMergeSA12*100 << setw( space )<<minAlltoallSA12<< setw( space )<< maxAlltoallSA12<<"  SA12"<<endl;}
				);
		}
	delete[] helparray4;

		/** Lexicographical naming   */
	P = new Pair[ recvoff[ size ] ];
	namelex( recvbuf, P, recvoff[ size ] ) ;

    }
	else {
	cout << myproc << " hat nichts bekommen" << endl;
	P = NULL;
    }


    /** recursion ? */
    int recursion = false;

    Debug2(
		if (myproc==size-1) {
			name2tuple[counter++] = 1 - (double)((uint)(2*((double)filelength /3.0))-P[ recvoff[ size ] - 1 ].name )/((uint)(2*((double)filelength /3.0))); });

	if ( myproc == size - 1 )       //Was passiert wenn das letzte PE keine Elemente hat?
		if ( P[ recvoff[ size ] - 1 ].name < (uint)(2*((double)filelength /3.0)) )
			recursion = true;

	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

	MPI_Bcast( &recursion, 1, MPI_INT, size - 1, MPI_COMM_WORLD );

	Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime );

	/** recursion ! */
	if ( recursion ) {
//      if(myproc==ROOT) cout<<"Rekursion"<<endl;
		uint names = (uint)(2*((double)filelength /3.0)) + 1;
		uint newlocArraylen = ( names + size - 1 ) / size;
		uint arraylen=recvoff[size];
		uint* tmparray= new uint[arraylen];
		Pair* sendarray= new Pair[arraylen];
		uint* t = new uint[size];

		int* sendcnt= new int[size];
		int* sendoff= new int[size+1];
		int* recvcnt= new int[size];

		for (int i=0 ; i<size ; i++) {
			t[i] = (i+1) * newlocArraylen;
			sendcnt[i] = 0;
		}
		uint half_names = names / 2;

		for (uint i = 0 ; i < arraylen; i++) {
			if (P[i].index%3==2) {
				tmparray[i] = findPosSATest(t, half_names + P[ i ].index / 3 ,size, sendcnt);
			}
			else {
				tmparray[i] = findPosSATest(t, P[ i ].index / 3 ,size, sendcnt);
			}
		}
		Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

		MPI_Alltoall( sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT , MPI_COMM_WORLD );

		Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime );

		uint* index= new uint[size];

		sendoff[0] = recvoff[0] = 0;

		for (int i = 0 ; i < size ; i++) {
			index[i] = 0;
			sendoff[i+1] = sendcnt[i] + sendoff[i];
			recvoff[i+1] = recvcnt[i] + recvoff[i];
		}

		for (uint i = 0; i < arraylen; i++)
			sendarray[sendoff[tmparray[i]] + index[tmparray[i]]++]=P[i];

		delete[] index;
		delete[] P;
		delete[] tmparray;

		Pair* in = new Pair[recvoff[size]];
		Pair* recvBufPair= new Pair[recvoff[size]+2];

		Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    starttime=MPI_Wtime() );

		MPI_Alltoallv( sendarray, sendcnt, sendoff, MPI_PAIR, in, recvcnt, recvoff, MPI_PAIR, MPI_COMM_WORLD );

		Debug3( alltoallvTime+= MPI_Wtime()-starttime );

		delete[] sendarray;

		Debug2( starttime=MPI_Wtime() );

		sortPairIndex(in, recvBufPair,recvoff[size], t, names);

		Debug2( permuteTime+=MPI_Wtime()-starttime );

		delete[] t;
		delete[] in;

		/** Sehr unwahrscheinlicher Fall: die ersten zwei Elemente werden    an das PE[size-2] gegeben jedoch hat
		    das letzte PE nur ein Element.  */
		if ( myproc == size - 1 && recvoff[ size ] < 2 ) {
			Debug0( cout << "Sehr dummer unwahrscheinlicher Fall eingetreten" << endl; );
			recvBufPair[ recvoff[ size ] ].name = MIN_INT;
			recvBufPair[ recvoff[ size ] ].index = recvBufPair[ recvoff[ size ] - 1 ].index + 3;
		}

		/** der recvBufPair des letzten PEs enthält am Ende ungültige Zeichen */
		uint* newInBuffer = new uint[ newlocArraylen + 2 ];
		for ( int i = 0;i < recvoff[size]; i++ )
			newInBuffer[ i ] = recvBufPair[ i ].name;


		/**Jedes PE muss seinem linken Nachbarn noch die ersten zwei Elemente schicken. Die letzte PE bekommt diese von PE 0 und mus diese danach mit 0 und MAX_INT überschreiben*/
		Pair temp[ 2 ];
		Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

		MPI_Sendrecv( recvBufPair, 2, MPI_PAIR, ( myproc - 1 + size ) % size, MSGTAG, temp, 2, MPI_PAIR, ( myproc + 1 ) % size , MSGTAG, MPI_COMM_WORLD, &status );

		Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime );

		delete[] recvBufPair;

		if ( myproc == size - 1 ) { //letzte PE hat naturgemäß zu wenig Elemente
			newInBuffer[ recvoff[size] ] = 0;
			for ( uint i = recvoff[size];i < newlocArraylen+2; i++ )
				newInBuffer[ i ] = MAX_INT;
		}
		else if ((uint)recvoff[size]!=newlocArraylen) { //vorletzte PE kann zu wenig Elemente haben?!? Führt vorher schon zu Fehler. s. Rekursion start.
			newInBuffer[ recvoff[size] ] = 0;
			for ( uint i = recvoff[size];i < newlocArraylen+2; i++ )
				newInBuffer[ i ] = MAX_INT;
		}
		else { //alle anderen PEs
			newInBuffer[ newlocArraylen ] = temp[ 0 ].name;
			newInBuffer[ newlocArraylen + 1 ] = temp[ 1 ].name;
		}


		uint* sa12rec;
		uint sa12len = 0;

	// cerr<<myproc<<"alles ok"<<endl;

		Debug2( rekursionStart[countRek++]= MPI_Wtime() );

		sa12rec = dc3( newInBuffer, names, &sa12len );

// cerr<<myproc<<"alles ok"<<endl;

		Debug2( rekursionEnd[--countRek]=   MPI_Wtime() - rekursionStart[countRek] );

		uint* n1all = new uint[ size ];
		uint* n2all = new uint[ size ];

		uint salenall;

		Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

		MPI_Allgather( &( n[ 1 ] ), 1, MPI_UNSIGNED, n1all, 1, MPI_UNSIGNED, MPI_COMM_WORLD );
		MPI_Allgather( &( n[ 2 ] ), 1, MPI_UNSIGNED, n2all, 1, MPI_UNSIGNED, MPI_COMM_WORLD );
		MPI_Scan(&sa12len, &salenall, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );

		Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime );

		uint* n12all = new uint[ 2 * size ];

		n12all[ 0 ] = n1all[ 0 ];
		for ( int i = 1; i < size ;i++ )
			n12all[ i ] = n1all[ i ] + n12all[ i - 1 ];
		for ( int i = size ; i < 2*size; i++ )
			n12all[ i ] = n2all[ i - size ] + n12all[ i - 1 ];
		n12all[ 2 * size - 1 ] ++; // die letzte 'Grenze wird um eins erhöht, damit das Element dem letzten PE gehört

		delete[] n1all;
		delete[] n2all;

		for ( int i = 0 ; i < size ; i++ )  sendcnt[ i ] = 0;
		P = new Pair[ sa12len ];
		salenall=salenall-sa12len+1;

		for ( uint i = 0; i < sa12len; i++ ) {
			P[ i ].name = i + salenall;
			// if (myproc==ausgabeCPU) cout <<i + 1 + salenAll[myproc+1] - sa12len<<" saAll "<< salenAll[myproc] <<" sa12len "<< sa12len<< endl;
			P[ i ].index = sa12rec[ i ];
			sa12rec[ i ] = findPosSA( n12all, sa12rec[ i ], size, sendcnt );

//          sendcnt[sa12rec[i]]++;
			Debug0( if( sa12rec[ i ] >= (uint)size )
						cout << myproc << " Hier stimmt was nicht " << sa12rec[ i ] << " sollte kleiner sein als #PEs" << endl;
				);
		};
//      delete[] salenAll;
		sendarray = new Pair[ sa12len ];
		Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

		MPI_Alltoall( sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT, MPI_COMM_WORLD );

		Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime );

		sendoff[ 0 ] = 0;
		recvoff[ 0 ] = 0;

		index = new uint[size];

		for ( int i = 0 ; i < size; i++ ) {
			index[ i ] = 0;
			sendoff[ i + 1 ] = sendcnt[ i ] + sendoff[ i ];
			recvoff[ i + 1 ] = recvoff[ i ] + recvcnt[ i ];
		}

		for ( uint i = 0; i < sa12len; i++ )
			sendarray[ sendoff[ sa12rec[ i ] ] + index[ sa12rec[ i ] ] ++ ] = P[ i ];

		delete[] index;
		delete[] sa12rec;
		delete[] P;

		Pair* recvarray = new Pair[ recvoff[ size ] + 2 ];

		Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)     starttime = MPI_Wtime() );

		MPI_Alltoallv( sendarray, sendcnt, sendoff, MPI_PAIR, recvarray, recvcnt, recvoff, MPI_PAIR, MPI_COMM_WORLD );

		Debug3( alltoallvTime += MPI_Wtime() - starttime );

		delete[] sendarray;
		uint half = ( locArraylen + 2 ) / 3 + 1;
		P = new Pair[ 2 * half ];

		Debug2( starttime = MPI_Wtime() );

		sortPair( recvarray, P, recvoff[ size ], n12all, half );

		Debug2( permuteTime += MPI_Wtime() - starttime );

		delete[] recvarray;
		delete[] n12all;

		Pair* sendTmp = new Pair[ 2 ];
		Pair* recvTmp = new Pair[ 2 ];

		sendTmp[ 0 ] = P[ 0 ];
		sendTmp[ 1 ] = P[ half ];

		Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

		MPI_Sendrecv( sendTmp, 2, MPI_PAIR, ( myproc - 1 + size ) % size, MSGTAG, recvTmp, 2, MPI_PAIR, ( myproc + 1 ) % size , MSGTAG, MPI_COMM_WORLD, &status );

		Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime );

		delete[] sendTmp;

		P[ n[ 1 ] ] = recvTmp[ 0 ];
		P[ half + n[ 2 ] ] = recvTmp[ 1 ];

		delete[] recvTmp;
		delete[] sendcnt; //Hier richtig?
		delete[] sendoff;
		delete[] recvcnt;
		delete[] recvoff;

		uint* satmp = sortS0S1S2( inbuffer, P, locArraylen, n, imod3, salen, half );
		return satmp;
		/** end recursion */
	}
	else {
//      Debug0(if ( myproc == ROOT )
//          cout << "---------------------   keine  Recursion---------------- " << locArraylen << endl;)

		/** Samplesort P nach Index */
		Debug2( starttime = MPI_Wtime() );

		Switch0(sort( P, P + recvoff[ size ] ) ); // cmpLessIndex
		Switch1(stable_sort( P, P + recvoff[ size ] ) ); // cmpLessIndex

		Debug2( sortingTime += MPI_Wtime() - starttime );
		Debug6(stlSortTime=MPI_Wtime() - starttime;
//  MPI_Reduce(&MergeSA12, &stlSortTime, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
//  MPI_Reduce(&MergeSA12, &stlSortTime, 1, MPI_DOUBLE, MPI_MIN, ROOT, MPI_COMM_WORLD);
			   cout<<myproc<<"   "<<stlSortTime<<" sec "<< recvoff[size]<<" Pair -- sort Pair Index "<<endl;
			);

/** testing */
		uint* pivpos= new uint[size+1];
		int* sendcnt= new int[size];
		int* sendoff= new int[size+1];
		int* recvcnt= new int[size];
//      int* recvoff= new int[size+1];

/** end *****/
		pivpos[ 0 ] = 0;
		Pair ptemp;
		ptemp.name=0;
		for ( int i = 0; i < size - 1; i++ ) {
			ptemp.index = (i+1) * locArraylen;
			pivpos[ i + 1 ] = findPos( P, ptemp , recvoff [ size ], cmpGreaterIndex );
		}
		pivpos[ size ] = recvoff[ size ];

		for ( int i = 0; i < size; i++ )
			sendcnt[ i ] = pivpos[ i + 1 ] > pivpos[ i ] ? pivpos[ i + 1 ] - pivpos[ i ] : 0;

		Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime());

		MPI_Alltoall( sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT , MPI_COMM_WORLD );

		Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime);

		sendoff[ 0 ] = recvoff[ 0 ] = 0;
		for ( int i = 1; i <= size; i++ ) {
			sendoff[ i ] = sendoff[ i - 1 ] + sendcnt[ i - 1 ];
			recvoff[ i ] = recvoff[ i - 1 ] + recvcnt[ i - 1 ];
		}
		recvoff[ size ] = recvoff[ size - 1 ] + recvcnt[ size - 1 ];
		Pair* recvBufPair = new Pair[ recvoff[ size ] + 2 ];

		Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)     starttime = MPI_Wtime() );

		MPI_Alltoallv( P, sendcnt, sendoff, MPI_PAIR, recvBufPair, recvcnt, recvoff, MPI_PAIR, MPI_COMM_WORLD );

		Debug3( alltoallvTime += MPI_Wtime() - starttime );

		delete[] P;

		Pair* helparray2 = new Pair[ recvoff[ size ] ];
		Debug2(starttime=MPI_Wtime());

		mergesort( recvBufPair, helparray2, 0, size, recvoff, cmpLessIndex );

		Debug2(mergeTime += MPI_Wtime() - starttime);

		delete[] helparray2;

		if ( myproc == size - 1 ) {
			if ( recvoff[ size ] == 1 ) {       //kann nur bei letztem PE auftreten?
				if ( recvBufPair[ recvoff[ size ] - 1 ].index % 3 == 1 )
					recvBufPair[ recvoff[ size ] ].index = recvBufPair[ recvoff[ size ] - 1 ].index + 1;
				else
					recvBufPair[ recvoff[ size ] ].index = recvBufPair[ recvoff[ size ] - 1 ].index + 2;
			}
		}

		/**jedes PE erhält zwei Elemente am Ende mehr*/
		Pair temp[ 2 ];
		Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

		MPI_Sendrecv( recvBufPair, 2, MPI_PAIR, ( myproc - 1 + size ) % size, MSGTAG, &temp, 2, MPI_PAIR, ( myproc + 1 ) % size , MSGTAG, MPI_COMM_WORLD, &status );

		Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime);

		if ( myproc == size - 1 ) {
			recvBufPair[ recvoff[ size ] ].name = MAX_INT;
			recvBufPair[ recvoff[ size ] + 1 ].name = MAX_INT;
			if ( recvBufPair[ recvoff[ size ] - 1 ].index % 3 == 1 ) {
				recvBufPair[ recvoff[ size ] ].index = recvBufPair[ recvoff[ size ] - 1 ].index + 1;
				recvBufPair[ recvoff[ size ] + 1 ].index = recvBufPair[ recvoff[ size ] - 1 ].index + 3;
			} else {
				recvBufPair[ recvoff[ size ] ].index = recvBufPair[ recvoff[ size ] - 1 ].index + 2;
				recvBufPair[ recvoff[ size ] + 1 ].index = recvBufPair[ recvoff[ size ] - 1 ].index + 3;
			}
		}
		else {
			recvBufPair[ recvoff[ size ] ] = temp[ 0 ];
			recvBufPair[ recvoff[ size ] + 1 ] = temp[ 1 ];
			if ( myproc == size - 2 ) { //wenn das letzte PE nur mod 2 Tupel hat,...
				if ( recvBufPair[ recvoff[ size ] ].index == recvBufPair[ recvoff[ size ] + 1 ].index )
					recvBufPair[ recvoff[ size ] ].index = recvBufPair[ recvoff[ size ] ].index - 1;
			}
		}

//      cout<<"CPU ["<<myproc<<"] nach rekursion"<<endl;
		Debug2( starttime = MPI_Wtime() );

		Switch0(sort( recvBufPair, recvBufPair + recvoff[ size ] + 2, cmpIndexModDiv ) );
		Switch1(stable_sort( recvBufPair, recvBufPair + recvoff[ size ] + 2, cmpIndexModDiv ) );

		Debug2( sortingTime += MPI_Wtime() - starttime );
		Debug6(stlSortTime=MPI_Wtime() - starttime;
//  MPI_Reduce(&MergeSA12, &stlSortTime, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
//  MPI_Reduce(&MergeSA12, &stlSortTime, 1, MPI_DOUBLE, MPI_MIN, ROOT, MPI_COMM_WORLD);
			   cout<<myproc<<"   "<<stlSortTime<<" sec "<<recvoff[ size ] + 2<<" Pair -- sort Pair %/ rekursion "<<endl;
			);

		uint half= n[1]+1;

		delete[] sendcnt; //Hier richtig?
		delete[] sendoff;
		delete[] recvcnt;
		delete[] recvoff;

		uint* satmp = sortS0S1S2( inbuffer, recvBufPair, locArraylen, n, imod3, salen, half );

		return satmp;

		/**Samplesort P end*/
	}
}


/**
 * Global naming of all Quadruple: same tuple => same name
 * @param in Quadruple*
 * @param P Pair* out-param
 * @param arraylen number of Quadruples
 */
inline void namelex( Quadruple* in, Pair* P, uint arraylen )
{
    /** naming with local names */
    Quadruple temp;
	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD));    mpiComStartTime = MPI_Wtime() );

	MPI_Sendrecv( &( in[ arraylen - 1 ] ), 1, MPI_QUADRUPLE, ( myproc + 1) % size, MSGTAG, &temp, 1, MPI_QUADRUPLE, ( myproc - 1 +size) % size , MSGTAG, MPI_COMM_WORLD, &status );

	Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime);

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
	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

	MPI_Scan( &name, &namesglob, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );

	Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime );

	for ( uint i = 0; i < arraylen; i++ )
		P[ i ].name += ( namesglob - name );
}

void writesa( uint* salen, uint* suffixarray, char* filenameOut )
{
    /** Suffixarry in Datei schreiben */
    MPI_Request request;

    uint* saLenAll=new uint[ size ];
    MPI_Gather( salen, 1, MPI_UNSIGNED, saLenAll, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

    if ( myproc == ROOT ) {
	int fd=open(filenameOut, O_WRONLY | O_CREAT | O_TRUNC , S_IRUSR | S_IWUSR);
	if(fd==-1) {cout << "Fehler bei der Ausgabe"<<endl; return ;}
	FILE* fp = fdopen(fd,"w");
	if (!fp)cout <<"fehler beim öffnen"<< endl;
	uint maxLen = 0;
	for ( int i = 0;i < size;i++ ) if ( maxLen < saLenAll[ i ] )   maxLen = saLenAll[ i ];

	fwrite(&suffixarray[1],sizeof(uint),saLenAll[0]-1,fp); //abschneiden der 'null' am anfang
	cout << "Von 0 Daten geschrieben" << endl;
	delete[] suffixarray;
	uint* outbuffer = new uint[ maxLen ];
	for ( int k = 1;k < size;k++ ) {
	    MPI_Recv( outbuffer, saLenAll[ k ], MPI_UNSIGNED, k, MSGTAG, MPI_COMM_WORLD, &status );
			cout << "salen von " << k << " ist " << saLenAll[ k ] << endl;
	    fwrite(outbuffer,sizeof(uint),saLenAll[k],fp);
//          Debug6(for (uint i=0 ; i<saLenAll[k] ; i++)
//                      cout<< outbuffer[i]<<" ";
//                  cout << endl;)

	    cout << "Von " << k << " Daten geschrieben" << endl;
	}
	delete[] saLenAll;
	delete[] outbuffer;
	if(fclose(fp)!=0) cout <<"Fehler beim Schließen der Datei";

    }
	else {
	MPI_Isend( suffixarray, *salen , MPI_UNSIGNED, 0, MSGTAG , MPI_COMM_WORLD , &request );
	MPI_Wait( &request, &status );
	delete[] suffixarray;
	cout << myproc << " ist fertig" << endl;
    }
}

Debug2(
void writeTimes( uint filelength )
{
    ofstream outfile( "timing", ios::app );
    assert( outfile.good() );
	cout.precision(4);
/*  if ( myproc == ROOT ) {
	outfile << size << " CPUs und Dateigroesse ist " << filelength << " Bytes" << endl;
	ofstream sampleTime("sampling",ios::app);
	sampleTime<< size << " CPUs und Dateigroesse ist " << filelength << " Bytes" << endl;
	int i=0;
	sampleTime<<"Sampling 1: ";
	while(samplingTime1[i])
	sampleTime<<samplingTime1[i++]<<" ";
	sampleTime<<endl<<"Sampling 2: ";
	i=0;
	while(samplingTime2[i])
	sampleTime<<samplingTime2[i++]<<" ";
	sampleTime<<endl;
    }*/
	Debug2(
		if ( myproc == ROOT ) {
			outfile << size << " CPUs und Dateigroesse ist " << filelength << " Bytes" << endl;
			ofstream rekursion("rekursion",ios::app);
			rekursion<< size << " CPUs und Dateigroesse ist " << filelength << " Bytes" << endl;
			int i=0;
			rekursion<<"Rekursionszeit: ";
			while(rekursionEnd[i++]>0 && i<maxRekursion)        rekursion<<rekursionEnd[i-1]-rekursionEnd[i]<<" ";
			i=0;
			rekursion<<endl<<"Rekursionszeit kummuliert: ";
			while(rekursionEnd[i]>0 && i<maxRekursion)  rekursion<<rekursionEnd[i++]<<" ";
			rekursion<<endl<<"Verhältnis Namen zu Tupel: ";
			i=0;
			while(name2tuple[i]>0 && i<maxRekursion)    rekursion<<name2tuple[i++]<<" ";
			rekursion<<endl;
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
    outfile << setw( space ) << dc3FinishTime / filelength;
    outfile << endl;
    outfile.close();
    space = 11;
    cout << setw( 4 ) << myproc;
    cout << setw( space ) << dc3FinishTime;
    cout << setw( space ) << sortingTime;
    cout << setw( space ) << mergeTime;
    cout << setw( space ) << permuteTime;
    cout << setw( space ) << alltoallvTime;
    cout << setw( space ) << mpiComTime;
    cout << setw( space ) << dc3FinishTime - ( sortingTime + mergeTime + permuteTime + alltoallvTime + mpiComTime);
    cout.setf( ios::scientific, ios::floatfield );
    cout << setw( space ) << dc3FinishTime / filelength;
    cout << endl;
}
)

void getTuple(uint* inbuffer, Pair* recvBufPair, Quintuple* S0, Quadruple* S1, Quintuple* S2, uint locArraylen, uint* n,  uint* imod3, uint& half)
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
    uint displacement = myproc * locArraylen;
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

uint* sortS0S1S2( uint* inbuffer, Pair* recvBufPair, uint locArraylen, uint* n, uint* imod3, uint* salen, uint& half )
{
	Debug2( sortS0S1S2TimeStart = MPI_Wtime() );
//     if ( myproc == ausgabeCPU )cout << "------------------sortS0S1S2-------------------" << endl;
//  uint k=(uint) ceil(samplefactor*(sqrt((double)locArraylen)/300));
	uint k = (uint)(sqrt((double) locArraylen / size) *samplefactor/3.0);
    if(k>=locArraylen/3){cout<<"k="<<k<<" zu groß, jetzt k="<<locArraylen/3-1<<endl; k=locArraylen/3-1;}

//  uint k1 = (uint)sqrt( locArraylen / (3.0 * size) *samplefactor);
//   if ( myproc == ROOT ) cout<<" k "<<3*k<<" % "<<3.0*k/locArraylen*100.0<<" Tuple " <<locArraylen<<endl;
    //    cout << myproc << " 0:" << imod3[ 0 ] << " 1:" << imod3[ 1 ] << " 2:" << imod3[ 2 ] << " locArraylen " <<  locArraylen << endl;
    Quintuple* S0 = new Quintuple[ n[ 0 ] ];
    Quadruple* S1 = new Quadruple[ n[ 1 ] ];
    Quintuple* S2 = new Quintuple[ n[ 2 ] ];
    getTuple( inbuffer,  recvBufPair,  S0,  S1,  S2, locArraylen,  n,  imod3,  half);

    /**Samplesort für S1 */
	Debug2( starttime = MPI_Wtime() );
	Switch0(
		sort( S0, S0 + n[ 0 ], cmpNameRank ); //Übergebe Funktion
		sort( S1, S1 + n[ 1 ], cmpFirstNameS1 ); //Übergebe Funktion
		sort( S2, S2 + n[ 2 ], cmpFirstNameS2 ); //Übergebe Funktion
		);
	Switch1(
		stable_sort( S0, S0 + n[ 0 ], cmpNameRank ); //Übergebe Funktion
		stable_sort( S1, S1 + n[ 1 ], cmpFirstNameS1 ); //Übergebe Funktion
		stable_sort( S2, S2 + n[ 2 ], cmpFirstNameS2 ); //Übergebe Funktion
		);
	Debug2( sortingTime += MPI_Wtime() - starttime );
	Debug6(stlSortTime=MPI_Wtime() - starttime;
//  MPI_Reduce(&MergeSA12, &stlSortTime, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
//  MPI_Reduce(&MergeSA12, &stlSortTime, 1, MPI_DOUBLE, MPI_MIN, ROOT, MPI_COMM_WORLD);
		   cout<<myproc<<"   "<<stlSortTime<<" sec "<< n[0]+n[1]+n[2]<<" Tuple sort S0, S1, S2 "<<endl;
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
	//if(samplebufS1[i].index>MOD0+(locArraylen+1)/3)cout<<"Fehler samplebufS1[i].index>MOD0+n[i] "<< samplebufS1[i].index<<" "<< S1[i].index/3 <<" "<<MOD0<<" "<<(locArraylen+1)/3 <<endl;
	samplebufS2[ i ] = S2[ int( i * d[2] ) ];
	samplebufS2[i].index = samplebufS2[i].index/3+MOD1;
    }

    /**Proc 0 sammelt Samples ein*/
    Quintuple* samplerecvbufS0=NULL;
    Quintuple* samplerecvbufS1=NULL;
    Quintuple* samplerecvbufS2=NULL;

    if ( myproc == ROOT ) {
	samplerecvbufS0 = new Quintuple[ k * size ];
	samplerecvbufS1 = new Quintuple[ k * size +1];
	samplerecvbufS2 = new Quintuple[ k * size +1];
    }

	Debug2(sampleStartTime=MPI_Wtime() );
	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

	MPI_Gather( samplebufS0, k, MPI_QUINTUPLE, samplerecvbufS0, k, MPI_QUINTUPLE, 0, MPI_COMM_WORLD );
    MPI_Gather( samplebufS1, k, MPI_QUINTUPLE, samplerecvbufS1, k, MPI_QUINTUPLE, 0, MPI_COMM_WORLD );
    MPI_Gather( samplebufS2, k, MPI_QUINTUPLE, samplerecvbufS2, k, MPI_QUINTUPLE, 0, MPI_COMM_WORLD );

	Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime);

	delete[] samplebufS0;
    delete[] samplebufS1;
    delete[] samplebufS2;

    Quintuple* pivbuf= new Quintuple[ size ];

    /**Proc 0 ermittelt Splitter*/
    if (myproc==ROOT) {
		Switch0(
			sort( samplerecvbufS0, samplerecvbufS0 + k * size, cmpNameRank );
			sort( samplerecvbufS1, samplerecvbufS1 + k * size, cmpFirstNameS2 );
			sort( samplerecvbufS2, samplerecvbufS2 + k * size, cmpFirstNameS2 );
			);
		Switch1(
			stable_sort( samplerecvbufS0, samplerecvbufS0 + k * size, cmpNameRank );
			stable_sort( samplerecvbufS1, samplerecvbufS1 + k * size, cmpFirstNameS2 );
			stable_sort( samplerecvbufS2, samplerecvbufS2 + k * size, cmpFirstNameS2 );
			);
		samplerecvbufS1[k*size].name[0]=MAX_INT;
	samplerecvbufS1[k*size].name[1]=MAX_INT;
	samplerecvbufS1[k*size].name[2]=MAX_INT;
	samplerecvbufS1[k*size].name[3]=MAX_INT;
	samplerecvbufS2[k*size].name[0]=MAX_INT_SAMPLE;
	samplerecvbufS2[k*size].name[1]=MAX_INT_SAMPLE;
	samplerecvbufS2[k*size].name[2]=MAX_INT_SAMPLE;
	samplerecvbufS2[k*size].name[3]=MAX_INT_SAMPLE;
	uint k3=3*k;
	Quintuple* Samples= new Quintuple[k3*size+2];
//cout<<myproc<<" 1006"<<endl;
	merge2(samplerecvbufS0,samplerecvbufS1,samplerecvbufS2,Samples, k*size);
//      cout<<"samples"<<endl;
//      for(int i =0;i<k3*size+2;i++) Samples[i].print();
//cout<<myproc<<" 1008"<<endl;
	delete[] samplerecvbufS0;
	delete[] samplerecvbufS1;
	delete[] samplerecvbufS2;
//      cout<<"pivbuf"<<endl;
	for ( int i = 0; i < size; i++ ) {
	    pivbuf[ i ] = Samples[ i * k3 ];
//          pivbuf[ i ].print();
	}
	delete[] Samples;
    }

    /**Splitter verteilen*/
	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );
	MPI_Bcast( pivbuf, size, MPI_QUINTUPLE, ROOT, MPI_COMM_WORLD );
// cerr <<myproc<<"----- 1031"<<endl;
	Debug3(mpiComTime += MPI_Wtime()-mpiComStartTime);
	Debug2( sampleTime+=MPI_Wtime()-sampleStartTime );
//cout<<myproc<<" 1021"<<endl;

	/**Einordnen der lokal sortierten Elemente anhand von Splitter*/
	uint* pivpos0= new uint[size+1];
    uint* pivpos1= new uint[size+1];
    uint* pivpos2= new uint[size+1];
    pivpos0[0]=0;
    pivpos1[0]=0;
    pivpos2[0]=0;
	MPI_Barrier(MPI_COMM_WORLD);
// cerr <<myproc<<"------- 1041"<<endl;

    for ( int i = 1; i < size; i++ ){
		if(pivbuf[ i ].index<MOD0){
			//             if (myproc==0){cout<<"0      -- ";pivbuf[i] .print();}
			pivbuf[i].index=pivbuf[i].index*3;
			pivpos0[ i ] = binsearch( S0, pivbuf[ i ], n[ 0 ], cmpS0S0 );
			pivpos1[ i ] = binsearch( S1, pivbuf[ i ], n[ 1 ], cmpS0S1 );
			pivpos2[ i ] = binsearch( S2, pivbuf[ i ], n[ 2 ], cmpS0S2 );
		}
		else if(pivbuf[ i ].index<MOD1){
//              if (myproc==0){cout<<"1      -- ";pivbuf[i] .print();}
			pivbuf[ i ].index=(pivbuf[i].index-MOD0)*3+1;
			pivpos0[ i ] = binsearch( S0, pivbuf[ i ], n[ 0 ], cmpS1S0 );
			pivpos1[ i ] = binsearch( S1, pivbuf[ i ], n[ 1 ], cmpS1S1 );
			pivpos2[ i ] = binsearch( S2, pivbuf[ i ], n[ 2 ], cmpS1S2 );
		}
		else{
//              if (myproc==0){cout<<"2      -- ";pivbuf[i] .print();}
			pivbuf[i].index=(pivbuf[i].index-MOD1)*3+2;
			pivpos0[ i ] = binsearch( S0, pivbuf[ i ], n[ 0 ], cmpS2S0 );
			pivpos1[ i ] = binsearch( S1, pivbuf[ i ], n[ 1 ], cmpS2S1 );
			pivpos2[ i ] = binsearch( S2, pivbuf[ i ], n[ 2 ], cmpS2S2 );
		}
    }

    delete[] pivbuf;
    pivpos0[ size ] = n[ 0 ];
    pivpos1[ size ] = n[ 1 ];
    pivpos2[ size ] = n[ 2 ];

    int* sendcnt0= new int[size];
    int* sendcnt1= new int[size];
    int* sendcnt2= new int[size];
    int* recvcnt0= new int[size];
    int* recvcnt1= new int[size];
    int* recvcnt2= new int[size];


    for ( int i = 0; i < size; i++ ) {
	sendcnt0[ i ] = pivpos0[ i + 1 ] > pivpos0[ i ] ? pivpos0[ i + 1 ] - pivpos0[ i ] : 0;
	sendcnt1[ i ] = pivpos1[ i + 1 ] > pivpos1[ i ] ? pivpos1[ i + 1 ] - pivpos1[ i ] : 0;
	sendcnt2[ i ] = pivpos2[ i + 1 ] > pivpos2[ i ] ? pivpos2[ i + 1 ] - pivpos2[ i ] : 0;
    }

	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)    mpiComStartTime = MPI_Wtime() );

	MPI_Alltoall( sendcnt0, 1, MPI_INT, recvcnt0, 1, MPI_INT , MPI_COMM_WORLD );
	MPI_Alltoall( sendcnt1, 1, MPI_INT, recvcnt1, 1, MPI_INT , MPI_COMM_WORLD );
	MPI_Alltoall( sendcnt2, 1, MPI_INT, recvcnt2, 1, MPI_INT , MPI_COMM_WORLD );

	Debug3( mpiComTime += MPI_Wtime()-mpiComStartTime );

	int* sendoff0= new int[size+1];
    int* sendoff1= new int[size+1];
    int* sendoff2= new int[size+1];
    int* recvoff0= new int[size+1];
    int* recvoff1= new int[size+1];
    int* recvoff2= new int[size+1];

    sendoff0[ 0 ] = recvoff0[ 0 ] = 0;
    sendoff1[ 0 ] = recvoff1[ 0 ] = 0;
    sendoff2[ 0 ] = recvoff2[ 0 ] = 0;

    for ( int i = 1; i < size + 1; i++ ) {
	sendoff0[ i ] = sendoff0[ i - 1 ] + sendcnt0[ i - 1 ];
	sendoff1[ i ] = sendoff1[ i - 1 ] + sendcnt1[ i - 1 ];
	sendoff2[ i ] = sendoff2[ i - 1 ] + sendcnt2[ i - 1 ];
	recvoff0[ i ] = recvoff0[ i - 1 ] + recvcnt0[ i - 1 ];
	recvoff1[ i ] = recvoff1[ i - 1 ] + recvcnt1[ i - 1 ];
	recvoff2[ i ] = recvoff2[ i - 1 ] + recvcnt2[ i - 1 ];
    }

    n[ 0 ] = recvoff0[ size ];

    Quintuple* recvbufS0 = new Quintuple[ n[0] ];
	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)     starttime = MPI_Wtime() );
	MPI_Alltoallv( S0, sendcnt0, sendoff0, MPI_QUINTUPLE, recvbufS0, recvcnt0, recvoff0, MPI_QUINTUPLE, MPI_COMM_WORLD );
	Debug3( alltoallvTime += MPI_Wtime() - starttime );
	delete[] S0;

    /** merge */
    Quintuple* helparray5 = new Quintuple[ n[0] ];
	Debug2(starttime = MPI_Wtime() );

	mergesort( recvbufS0, helparray5, 0, size, recvoff0, cmpSplitterLeqS0 );

	Debug2( mergeTime += MPI_Wtime() - starttime );

	delete[] helparray5;

    /** Samplesort S0 end */

    /** ein 'MAX'-Element am Ende mehr */
    n[ 1 ] = recvoff1[ size ];
    Quadruple* recvbufS1 = new Quadruple[ n[1] + 1 ];
	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)     starttime = MPI_Wtime() );
	MPI_Alltoallv( S1, sendcnt1, sendoff1, MPI_QUADRUPLE, recvbufS1, recvcnt1, recvoff1, MPI_QUADRUPLE, MPI_COMM_WORLD );
	Debug3( alltoallvTime += MPI_Wtime() - starttime );
	delete[] S1;

    /** merge */
    Quadruple* helparray4 = new Quadruple[ n[1] ];
	Debug2( starttime = MPI_Wtime() );
	mergesort( recvbufS1, helparray4, 0, size, recvoff1, cmpFirstNameS1 );
	Debug2( mergeTime += MPI_Wtime() - starttime );
	delete[] helparray4;

    recvbufS1[ n[1] ].name[ 0 ] = MAX_INT_SAMPLE;
    recvbufS1[ n[1] ].name[ 1 ] = MAX_INT_SAMPLE;
    recvbufS1[ n[1] ].name[ 2 ] = MAX_INT_SAMPLE;

    n[ 2 ] = recvoff2[ size ];
    Quintuple* recvbufS2 = new Quintuple[ n[ 2 ] + 1 ];

	Debug3( Debug4(MPI_Barrier(MPI_COMM_WORLD);)     starttime = MPI_Wtime() );
	MPI_Alltoallv( S2, sendcnt2 , sendoff2 , MPI_QUINTUPLE, recvbufS2, recvcnt2, recvoff2, MPI_QUINTUPLE, MPI_COMM_WORLD );
	Debug3( alltoallvTime += MPI_Wtime() - starttime );
	delete[] S2;

    /** merge */
    helparray5 = new Quintuple[ n[ 2 ] ];
	Debug2( starttime = MPI_Wtime() );
	mergesort( recvbufS2, helparray5, 0, size, recvoff2, cmpFirstNameS2 );
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
	delete[] pivpos0;
	delete[] sendcnt1;
	delete[] sendoff1;
	delete[] recvcnt1;
	delete[] recvoff1;
	delete[] pivpos1;
	delete[] sendcnt2;
	delete[] sendoff2;
	delete[] recvcnt2;
	delete[] recvoff2;
	delete[] pivpos2;

	Debug2( starttime = MPI_Wtime() );
	merge( recvbufS0, recvbufS1, recvbufS2, suffixarray, n );
	Debug2( mergeTime += MPI_Wtime() - starttime );

	Debug5( //Lastverteilung
		double sortSAll= MPI_Wtime() - starttime;
		double maxsortAll;
		double minsortAll;
		uint max;
		uint min;
//  Debug6(cout << myproc << " SortS0S1S2 Alltoallv: before "<< locArraylen <<" after " << *salen << endl; )
		MPI_Reduce(salen, &max, 1, MPI_UNSIGNED, MPI_MAX, ROOT, MPI_COMM_WORLD);
		MPI_Reduce(salen, &min, 1, MPI_UNSIGNED, MPI_MIN, ROOT, MPI_COMM_WORLD);
		MPI_Reduce(&sortSAll, &maxsortAll, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
		MPI_Reduce(&sortSAll, &minsortAll, 1, MPI_DOUBLE, MPI_MIN, ROOT, MPI_COMM_WORLD);

		uint space=11;
/**
   min -- max -- max-min -- %-Diff -- %-Tuple -- merge-min  -- merge-max
*/
		if (myproc==ROOT) {cout<< setw( space )<<min<< setw( space )<<max<< setw( space )<<max-min<< setw( space );
//  cout.setf( ios::scientific, ios::floatfield );
			cout << (double)(max-min)/min*100<< setw( space )<<(300.0*k)/locArraylen<< setw( space )<<minsortAll<< setw( space )<< maxsortAll<<setw (space)<<(maxsortAll-minsortAll)/minsortAll*100<<" Alltoallv all Tuple"<<endl;}
		);

	delete[] recvbufS0;
    delete[] recvbufS1;
    delete[] recvbufS2;

//  Debug6( if ( myproc == 0 ) cerr << " ende sortS0S1S2 " << endl; )
	Debug2( sortS0S1S2TimeEnd+= (MPI_Wtime()-sortS0S1S2TimeStart) );
// cout<<myproc<<" hat "<<n[ 0 ] + n[ 1 ] + n[ 2 ] <<" zu sortieren in "<<(MPI_Wtime()-sortS0S1S2TimeStart) <<endl;
	return suffixarray;
}

/**
 * main
 * @param argc
 * @param argv
 * @return
 */

int main( int argc, char **argv )
{
    MPI_Init( &argc, &argv );

    mpi_get_init();

    uint filelength;

    char* filename=NULL;
    if ( myproc == ROOT ) {
	// Eingabe Datei
	if ( argc < 2 ) {
	    cout << "Keine Eingabedatei! Aufruf z.B. mpdrun -np 4 ./dc INPUT SAMPLESIZE INPUTSIZE" << endl;
	    return 0;
	}
	//  cout << argv[ 1 ] << endl;
	filename = argv[ 1 ];
	ifstream infile( filename );
	assert( infile );

	// Groesse der Datei
	infile.seekg( 0L, ios::end );
	filelength = ( uint ) infile.tellg() + 1 ; // + 1 wegen Abschlusszeichen
	infile.seekg( 0L, ios::beg );
	infile.close();
	if (argc > 3 ) if (atoll(argv[3])< filelength && atoll(argv[3])>0) filelength= atoll(argv[3]);
    }

    MPI_Bcast( &filelength, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD );

    if ( argc > 2 )
	samplefactor = atoi( argv[ 2 ] ) ;
    else
	samplefactor = size;

    if ( argc > 5 )
	ausgabeCPU = atoi( argv[ 5 ] );

    uint locArraylen = (uint)ceil( (double)filelength / size);

    if ( myproc == ROOT )
	cout << "Groesse der Datei ist " << filelength << " bytes " << (uint)(2*((double)filelength /3.0))<< " locArraylen " << locArraylen << endl;

    // Verteilung der Daten nach inbuffer<- enthält lokal benötigte Zeichen +2 wegen Überschneidung und zum Erstellen der S0-S1 Tupel

    int overlap = 2;
    char* inbuffer = new char[ locArraylen + overlap ];

    if ( myproc == ROOT ) {
	ifstream infile( filename, ios::binary );
	infile.seekg( locArraylen + overlap, ios::beg );
	//Einlesen der Elemente für die Prozessoren 1 bis n-2
	for ( int i = 1; i < size - 1; i++ ) {
	    infile.seekg( -overlap, ios::cur );
//          cout << endl << " Einlesen von Daten ab Position " << i * locArraylen << " ";
	    infile.read( inbuffer, locArraylen + overlap );
	    MPI_Send( inbuffer, locArraylen + overlap, MPI_CHAR, i, MSGTAG , MPI_COMM_WORLD );
	}
	//Einlesen der Elemente für Prozessor n-1
	infile.seekg( -overlap, ios::cur );
	uint sonderfallPn = filelength - 1 - locArraylen * ( size - 1 ) ;
//      cout << endl <<-overlap<<" "<<sonderfallPn<< " Einlesen von Daten ab Position " << ( size - 1 ) * locArraylen << " ";
	infile.read( inbuffer, sonderfallPn );
	MPI_Send( inbuffer, locArraylen + overlap , MPI_CHAR, size - 1, MSGTAG , MPI_COMM_WORLD );

	if (!infile.good()) {cout<<"Error reading File!"<<endl; assert(0);}
	//Einlesen der Elemente für Prozessor ROOT
	infile.seekg( 0L, ios::beg );
//      cout << endl << " Einlesen von Daten ab Position " << 0L << endl;
	infile.read( inbuffer, locArraylen + overlap );
	infile.close();
		Debug5(
			cout << "         min   -- max     -- Diff -- %-Diff -- %-Tuple -- merge-min--merge-max-- %-Diff alltoall-min alltoall-max"<<endl
			);
	}
	else { // not ROOT
	MPI_Recv( inbuffer, locArraylen + overlap, MPI_CHAR, ROOT, MSGTAG, MPI_COMM_WORLD, &status );
	}

    uint* integerBuf = new uint[ locArraylen + overlap ];
    read( inbuffer, locArraylen + overlap, integerBuf );

    if ( myproc == size - 1 )   // Abschlusszeichen setzen
	integerBuf[ filelength - 1 - locArraylen * ( size - 1 ) ] = 0;

    delete[] inbuffer;

    uint salen;
    MPI_Barrier( MPI_COMM_WORLD );

	Debug2( rekursionStart[countRek++]= MPI_Wtime() );
	Debug1( dc3StartTime = MPI_Wtime() );

	uint*   suffixarray = dc3( integerBuf, filelength, &salen );

	Debug1( dc3FinishTime = MPI_Wtime() - dc3StartTime );
	Debug2( rekursionEnd[--countRek]=   MPI_Wtime() - rekursionStart[countRek] );

	Debug1( if ( myproc == ROOT ) cout << "Zeit für Erstellen von Suffixarry = "<<dc3FinishTime<<endl );

	Debug2(
		if (myproc==ROOT) //Zeitmessung
			MPI_Recv( name2tuple, 35, MPI_DOUBLE, size - 1, MSGTAG, MPI_COMM_WORLD, &status );
		if (myproc==size-1)
			MPI_Send( name2tuple, 35, MPI_DOUBLE, 0, MSGTAG , MPI_COMM_WORLD );
		int token = 0;
		if ( myproc == ROOT ) {
			cout << " CPU    Gesamtzeit    Sort    Merge     Permute   MPI_Alltoallv  MPI_Com   Rest   Zeit/Zeichen " << endl;
			token = 1;
		}
		if ( !token )
			MPI_Recv( &token, 1, MPI_INT, myproc - 1, MSGTAG, MPI_COMM_WORLD, &status );
		writeTimes( filelength );
		if ( myproc != size - 1 )
			MPI_Send( &token, 1, MPI_INT, ( myproc + 1 ) % size, MSGTAG , MPI_COMM_WORLD );
		);
	if ( argc > 3 ) writesa( &salen, suffixarray, argv[ 3 ] );

    MPI_Finalize();
    return 0;
}
