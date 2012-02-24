#include <iostream>
#include <fstream>
#include <cassert>

#include <mpi.h>
#include "tuple.h"
#define ROOT 0
#define MSGTAG 42

using namespace std;
typedef unsigned int uint;
typedef unsigned long long ull;

int myproc,size;
ull filelength;


/** timing */
double starttimeSA, endtimeSA, starttimeMAIN, endtimeMAIN;


MPI_Datatype MPI_TRIPLE;
MPI_Datatype MPI_PAIR;
MPI_Status status;

void mpi_get_init() {
	MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	MPI_Datatype mpi_tupel[ 2 ] = {MPI_UNSIGNED, MPI_UNSIGNED};

	int mpi_blocklen[ 2 ] = {2, 1};
	MPI_Aint mpi_displ[ 2 ] = {0, 2 * sizeof( uint ) };
	MPI_Type_struct( 2, mpi_blocklen, mpi_displ, mpi_tupel, &MPI_TRIPLE );
	MPI_Type_commit( &MPI_TRIPLE );
	
	mpi_blocklen[ 0 ] =  1;
	mpi_displ[ 1 ] =   sizeof( uint ) ;
	MPI_Type_struct( 2, mpi_blocklen, mpi_displ, mpi_tupel, &MPI_PAIR );
	MPI_Type_commit( &MPI_PAIR );
}

void print_sendcount(int* sendcnt){
	int token = 0;
	if ( myproc == ROOT ) 	token = 1;
	if ( !token ) 
		MPI_Recv( &token, 1, MPI_INT, myproc - 1, MSGTAG, MPI_COMM_WORLD, &status );
	
	cout<<"if(myproc=="<<myproc<<"){  ";
	for (int i=0 ; i<size ; i++)
		cout<<"sendcnt["<<i<<"]="<<sendcnt[i]<<"; ";
	cout<<"}"<<endl;
	if ( myproc != size - 1 )
		MPI_Send( &token, 1, MPI_INT, ( myproc + 1 ) % size, MSGTAG , MPI_COMM_WORLD );
}

template <class T1,class T2> void initArray(T1* array, uint n, T2 x){
		for (uint i=0 ; i<n ; i++)	array[i]=(T1)x;
}

inline void pairsort( Pair* in, Pair* out, uint length, uint Slen) {
	uint tmp=myproc*Slen;
	if (myproc==size-1)	length = (uint)filelength - myproc*Slen;
	for ( uint i = 0 ; i < length; i++ ) 
		out[( in[ i ].index-tmp) ]=in[i];
}

inline void triplesort( Triple* in, Triple* out, uint length, uint Slen) {
	uint tmp=myproc*Slen;
	if (myproc==size-1)	length = (uint)filelength - myproc*Slen;
	for ( uint i = 0 ; i < length; i++ ) 
		out[( in[ i ].name[0]-tmp) ]=in[i];
}

namespace check{
/**
 * liefert entweder das gesuchte Element oder das nächst kleinere. wird verwendet für die aufteilung der Elemente
 * auf Prozessoren
 * @param array 
 * @param element gesuchte Element
 * @param right Arraygröße
 * @return Position
 */
uint binsearch(uint* array, uint element, int right) {
	int left=0;
	uint  middle;
	while (left < right){
		middle= (left+right) / 2;
		if (array[middle] < element) left= middle+1;
		else 	right= middle; 
	}
	if (array[left] == element) return left;
	return left-1;
}

}
uint isPermutation(Pair* P, uint Slen){
	uint index=0;
	uint start=P[index].index;
	for (uint i=0 ; i<Slen ; i++){
		if(++start!=P[++index].index) return 0;
// 		cout<<start<<" "<<index<<" "<<P[index].index<<endl;
	}
	return 1;

}

int sa_check(uint* SA, char* S, uint SAlen, uint Slen ){
//	mpi_get_init();
	uint globSlen=(uint)((filelength+size-1)/size);
	starttimeSA=MPI_Wtime();	
	uint* splitter= new uint[size+1];
	for (int i=0 ; i<=size ; i++)	splitter[i]=i*globSlen;
	splitter[size]++;
	uint suffixe;
	//Suffixe nicht gleich verteilt
/* 	MPI_Scan( &SAlen, &suffixe, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );
 	suffixe-=SAlen;*/
	//Suffixe neu eingelesen und damit gleich verteilt
	suffixe=myproc*globSlen;

	Pair* P= new Pair[SAlen];
	uint* receiver= new uint[SAlen];
	int sendcnt[size],recvcnt[size];
 	initArray(sendcnt,size,0);

	for (uint i=0 ; i<SAlen ; i++){
		P[i].name=i+suffixe;
		P[i].index=SA[i];
//   		P[i].print();		
		receiver[i]=check::binsearch(splitter, SA[i], size);
// 		cout<<receiver[i]<<" ";
		sendcnt[receiver[i]]++;
	}

/** ********* testing */
	
//	print_sendcount(sendcnt);

/** ******* testing end ***/

delete[] SA;
// cout << endl;
	MPI_Alltoall( sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT , MPI_COMM_WORLD );
	int sendoff[size+1],recvoff[size+1];
	sendoff[ 0 ] = recvoff[ 0 ] = 0;
	for ( int i = 1; i < size + 1; i++ ) {
		sendoff[ i ] = sendoff[ i - 1 ] + sendcnt[ i - 1 ];
		recvoff[ i ] = recvoff[ i - 1 ] + recvcnt[ i - 1 ];
	}

// 	if(myproc<2){
// 		cout<<endl<<"Sendcount:";	for (int i=0 ; i<size ; i++) cout<<" sc "<<sendcnt[i];
// 		cout<<endl<<"Recvcount:";	for (int i=0 ; i<size ; i++) cout<<" rc "<<recvcnt[i];
// 		cout<<endl<<"Sendoff:";		for (int i=0 ; i<=size ; i++) cout<<" so "<<sendoff[i];
// 		cout<<endl<<"Recvoff:";		for (int i=0 ; i<=size ; i++) cout<<" ro "<<recvoff[i];
// 	}

	Pair* Psend= new Pair[SAlen];
	uint* x= new uint[size];
 	initArray(x,size,0);
	for (uint i=0 ; i<SAlen ; i++){
		Psend[x[receiver[i]] + sendoff[receiver[i]]]=P[i];
		x[receiver[i]]++;
	}
	delete[] P;
	delete[] receiver;

	Pair* Ptmp= new Pair[recvoff[ size ]];
	MPI_Alltoallv( Psend, sendcnt, sendoff, MPI_PAIR, Ptmp, recvcnt, recvoff, MPI_PAIR, MPI_COMM_WORLD );
	delete[] Psend;
	uint error=0;
	if((uint)recvoff[size]!=Slen){
		cout<<myproc<<" recvoff[size]= "<<recvoff[size]<<"!="<<Slen<<"=Slen"<<endl;
		error=1;
	}
	MPI_Allgather( &error, 1, MPI_UNSIGNED, &error, 1, MPI_UNSIGNED, MPI_COMM_WORLD );
		//Bedingung für gleich verteilt und Permutation
	if(error) return 1;

	P=new Pair[recvoff[size]+1];
	pairsort(Ptmp, P, recvoff[size],globSlen);
	delete[] Ptmp;
	MPI_Sendrecv( &P[0], 1, MPI_PAIR, ( myproc - 1 + size ) % size, MSGTAG, &P[recvoff[size]], 1, MPI_PAIR, ( myproc + 1 ) % size , MSGTAG, MPI_COMM_WORLD, &status );

//  	for (uint i=0 ; i<Slen+1 ; i++)	P[i].print(); cout << endl;

	if(myproc==size-1) P[Slen].index=P[Slen-1].index+1; //für Permutation benötigt, damit alle Elemente überprüft werden können
	if(!isPermutation(P,Slen)){
		cout<<endl<<myproc<<" keine Permutation"<<endl;
		error=1;
	}
	MPI_Allgather( &error, 1, MPI_UNSIGNED, &error, 1, MPI_UNSIGNED, MPI_COMM_WORLD );
	if(error) return 1;

	receiver = new uint[Slen];
	Triple* I= new Triple[Slen];

	initArray(sendcnt,size,0);
	for (uint i=0 ; i<Slen; i++){
		I[i].name[0]=P[i].name;
		I[i].name[1]=P[i+1].name;
		I[i].index=(uint)S[i];
//  		I[i].print();
		receiver[i]= check::binsearch(splitter, I[i].name[0], size);	
		sendcnt[receiver[i]]++;
	}

/** ********* testing */
	
//	print_sendcount(sendcnt);

/** ******* testing end ***/

delete[] P;
delete[] S;
//   	cout<<endl<<myproc<<"I erstellt"<<endl;
// 	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Alltoall( sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT , MPI_COMM_WORLD );
	sendoff[ 0 ] = recvoff[ 0 ] = 0;
	for ( int i = 1; i < size + 1; i++ ) {
		sendoff[ i ] = sendoff[ i - 1 ] + sendcnt[ i - 1 ];
		recvoff[ i ] = recvoff[ i - 1 ] + recvcnt[ i - 1 ];
	}
/*	if(myproc<2){
		cout<<endl<<"Sendcount:";	for (int i=0 ; i<size ; i++) cout<<" sc "<<sendcnt[i];
		cout<<endl<<"Recvcount:";	for (int i=0 ; i<size ; i++) cout<<" rc "<<recvcnt[i];
		cout<<endl<<"Sendoff:";		for (int i=0 ; i<=size ; i++) cout<<" so "<<sendoff[i];
		cout<<endl<<"Recvoff:";		for (int i=0 ; i<=size ; i++) cout<<" ro "<<recvoff[i];
	}*/
	Triple* Isend=new Triple[Slen];
 	initArray(x,size,0);
	for (uint i=0 ; i<SAlen ; i++){
		Isend[x[receiver[i]] + sendoff[receiver[i]]]=I[i];
		x[receiver[i]]++;
	}
	delete[] I;
	delete[] receiver;
	
	Triple* Itmp= new Triple[recvoff[ size ]];
	MPI_Alltoallv( Isend, sendcnt, sendoff, MPI_TRIPLE, Itmp, recvcnt, recvoff, MPI_TRIPLE, MPI_COMM_WORLD );
	delete[] Isend;

	I = new Triple[Slen+1];
	triplesort(Itmp, I, recvoff[size],globSlen);
	MPI_Sendrecv( I, 1, MPI_TRIPLE, ( myproc - 1 + size ) % size, MSGTAG, &(I[Slen]), 1, MPI_TRIPLE, ( myproc + 1 ) % size , MSGTAG, MPI_COMM_WORLD, &status );


	if(myproc==size-1)	Slen--; //letztes Element hat keinen Nachfolger 
	for (uint i=0 ; i<Slen ; i++){
		if(I[i].index==I[i+1].index){
			if(I[i].name[1]==I[i+1].name[1]) { cout<<myproc <<" Fehler bei "<<i<<" von "<<Slen<<endl; error++;}
		}else{
			if(I[i].index>I[i+1].index){ cout<<myproc <<" Fehler bei "<<i<<endl;
				error++;}
		}
	}
	endtimeSA=MPI_Wtime()-starttimeSA;
	return error;
}

int readFile(char* S, uint Slen, char* filename, uint filelength){
	MPI_Status status;
	
	if ( myproc == ROOT ) {
		ifstream infile( filename, ios::binary );
		char* Stmp=new char[Slen];
		infile.read( S, Slen );
		//Einlesen der Elemente für die Prozessoren 1 bis n-2
		for ( int i = 1; i < size - 1; i++ ) {
			infile.read(Stmp, Slen);
			MPI_Send( Stmp, Slen, MPI_CHAR, i, MSGTAG , MPI_COMM_WORLD );
		}
		//Einlesen der Elemente für Prozessor n-1
		uint PEmax = filelength - Slen * ( size - 1 ) ;
		infile.read(Stmp, PEmax );
		MPI_Send( Stmp, PEmax, MPI_CHAR, size-1, MSGTAG , MPI_COMM_WORLD );
		delete[] Stmp;
		infile.close();
	} else   //not ROOT
		MPI_Recv( S, Slen, MPI_CHAR, ROOT, MSGTAG, MPI_COMM_WORLD, &status );
	return 0;
}

int readSA(uint* S, uint Slen, char* filename, uint filelength){
	MPI_Status status;
	if ( myproc == ROOT ) {
		ifstream infile( filename, ios::binary );
		uint* Stmp=new uint[Slen];
		infile.read( (char*)S, 4*Slen );
		//Einlesen der Elemente für die Prozessoren 1 bis n-2
		for ( int i = 1; i < size - 1; i++ ) {
			infile.read((char*)Stmp, 4*Slen);
			MPI_Send( Stmp, Slen, MPI_UNSIGNED, i, MSGTAG , MPI_COMM_WORLD );
		}
		uint PEmax = filelength - Slen * ( size - 1 ) ;
		infile.read((char*)Stmp, 4*PEmax );
		MPI_Send( Stmp, PEmax, MPI_UNSIGNED, size-1, MSGTAG , MPI_COMM_WORLD );
		delete[] Stmp;
		infile.close();
	} else   //not ROOT
		MPI_Recv( S, Slen, MPI_UNSIGNED, ROOT, MSGTAG, MPI_COMM_WORLD, &status );
	return 0;
}

ull get_filelength(char* filename){
	ull filelength;
	ifstream infile( filename );
	if(!infile)cout<<"kann Eingabedatei nicht öffnen"<<endl;
	assert( infile );
	infile.seekg( 0L, ios::end );
	filelength = ( ull ) infile.tellg() ; 
	infile.seekg( 0L, ios::beg );
	infile.close();
	return filelength;
}

template <class T>
	void printfile(T* S,uint len){
		cout<<myproc<<" Ausgabe Zeichenkette start: ";
		for (uint i=0 ; i<len ; i++){
			cout<<(uint)S[i]<<" ";
		}
		cout<<" ende" <<endl;
}

	void printfile(char* S,uint len){
		cout<<myproc<<" Ausgabe Zeichenkette start: ";
		for (uint i=0 ; i<len ; i++){
			cout<<S[i]<<" ";
		}
		cout<<" ende" <<endl;
}


int main( int argc, char **argv ) {

	MPI_Init( &argc, &argv );
	mpi_get_init();
	starttimeMAIN=MPI_Wtime();
		uint error=0;
		if ( myproc == ROOT && argc<3){
				cout<<"sa_check  SUFFIXARRAY   STRING"<<endl;
				error=1;
		}
	
		if(myproc==ROOT){
			filelength=get_filelength(argv[2]);
			if(4*filelength!=get_filelength(argv[1])){
				cout	<<"SA entspricht nicht S oder Sonderzeichen $ ist in SA enthalten"<<endl
						<<"Dateigröße SA="<< get_filelength(argv[1]) <<" != "<<filelength <<" 4* Dateigröße S"<<endl;
				error=1;
			}
		}
		MPI_Bcast( &error, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD );

		if(error) MPI_Finalize();

		MPI_Bcast( &filelength, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD );
	
		uint Slen;
		if(myproc==size-1) Slen =(uint)( filelength  - (( filelength + size - 1 ) / size ) * ( size - 1 ) );
		else 	Slen=(uint)(( filelength + size - 1 ) / size );
//		if (myproc==0)cout<<"Slen "<<Slen;

		char* S=new char[Slen];
		uint* SA=new uint[Slen];
		readSA(SA, Slen, argv[1], (uint)filelength);
		readFile(S, Slen, argv[2], (uint)filelength);
		if(argc>3 && argv[3][0]=='p'){
			printfile(S,Slen);
			printfile(SA,Slen);
		}
		cout<<myproc<<" hat "<<Slen<<" von S[" <<filelength<<"] bekommen"<<endl;
		uint result =	sa_check(SA, S, Slen, Slen );
		if(result)	{
			cout<<myproc<<" hat "<< result<<"  FEHLER  im SA entdeckt"<<endl;
		}else{
			endtimeMAIN=MPI_Wtime()-starttimeMAIN;
			if (myproc==ROOT) 
				cout	<<"Gesamtzeit inklusive Einlesen Dateien: " <<endtimeMAIN<<"s (pro Zeichen " <<endtimeMAIN/filelength<<" )" <<endl
				<<"Zeit für Überprüfung: " <<endtimeSA<<"s (pro Zeichen "<<endtimeSA/filelength<<" )" <<endl;
		}
		MPI_Finalize();
		return 0;
}


