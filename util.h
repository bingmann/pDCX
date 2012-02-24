// // // // this files contains all the application independent little
// functions and macros used for the optimizer.
// In particular Peters debug macros and Dags stuff
// from dbasic.h cdefs, random,...

//////////////// stuff originally from debug.h ///////////////////////////////
// (c) 1997 Peter Sanders
// some little utilities for debugging adapted
// to the paros conventions


#ifndef UTIL
#define UTIL
typedef unsigned char uchar;
typedef unsigned int uint;
using namespace std;

// default debug level. will be overidden e.g. if debug.h is included
#ifndef DEBUGLEVEL
#define DEBUGLEVEL 3
#endif

#if DEBUGLEVEL >= 0
#define Debug0(A) A
#else
#define Debug0(A) 
#endif
#if DEBUGLEVEL >= 1
#define Debug1(A) A
#else
#define Debug1(A) 
#endif
#if DEBUGLEVEL >= 2
#define Debug2(A) A
#else
#define Debug2(A) 
#endif
#if DEBUGLEVEL >= 3
#define Debug3(A) A
#else
#define Debug3(A) 
#endif
#if DEBUGLEVEL >= 4
#define Debug4(A) A
#else
#define Debug4(A) 
#endif
#if DEBUGLEVEL >= 5
#define Debug5(A) A
#else
#define Debug5(A) 
#endif
#if DEBUGLEVEL >= 6
#define Debug6(A) A
#else
#define Debug6(A) 
#endif

#define Assert(c) if(!(c))\
  {cout << "\nAssertion violation " << __FILE__ << ":" << __LINE__ << endl;}
#define Assert0(C) Debug0(Assert(C))
#define Assert1(C) Debug1(Assert(C))
#define Assert2(C) Debug2(Assert(C))
#define Assert3(C) Debug3(Assert(C))
#define Assert4(C) Debug4(Assert(C))
#define Assert5(C) Debug5(Assert(C))

#define Error(s) {cout << "\nError:" << s << " " << __FILE__ << ":" << __LINE__ << endl;}

#if SORT == 0
#define Switch0(A) A
#else
#define Switch0(A) 
#endif
#if SORT == 1
#define Switch1(A) A
#else
#define Switch1(A) 
#endif


uchar* int2uchar(int number){
	uchar* temp = new uchar[4];
	for(int i=0;i<4;i++)
		temp[i]=(uchar) ((number >> 8*i) & 0xFF);

/*	cout<<endl<<"Zahl "<<number<<" = ";
	for(int i=0;i<4;i++)
		cout<< (int)temp[i]<<" ";
	cout<<endl; */

	return temp;
}

uint uchar2int(uchar* string){
	uint temp=0;
	for(int i=3;i>0;i--){
		temp+=(uint)string[i];
		temp<<=8;
	}
	temp+=(uint)string[0];
	
/*	cout<<endl<<"String ";
	for(int i=0;i<4;i++)
		cout<< (int)string[i]<<" ";
	cout<<" = " <<temp<<endl;*/

	return temp;
}

void int2uchar(int number, uchar* temp){
	for(int i=0;i<4;i++)
		temp[i]=(uchar) ((number >> 8*i) & 0xFF);
}

void uchar2int(uchar* string, uint temp){
	for(int i=3;i>0;i--){
		temp+=(uint)string[i];
		temp<<=8;
	}
	temp+=(uint)string[0];
}

double pow( int a, int b ) {
	double temp = 1;
	if(b>=0){
		for ( int i = 0; i < b ; i++ ) 
			temp *= a;
		return temp;
	}else{
		for ( int i = 0; i > b ; i-- ) 
			temp /= a;
		return temp;
	}
}

void printBits(uint t){
	for (int i=31 ; i>=0 ; i--){
		if ( t & (uint)pow(2,i)) cout<<"1";
		else cout<<"0";
		if (i%4==0) cout<<".";
	}

}

/**
 * changes character array into unsigned integer array
 * @param inbuffer character array
 * @param arraylen length of inbuffer
 * @return uint array
 */
uint* read( char* inbuffer, int arraylen ){
	uint * temp = new uint[ arraylen ];
	for ( int i = 0; i < arraylen; i++ ) {
		temp[ i ] = (unsigned char)inbuffer[ i ];
	}
	return temp;
}

void read( char* inbuffer, int arraylen, uint* temp){
	for ( int i = 0; i < arraylen; i++ ) 
		temp[ i ] = (unsigned char)inbuffer[ i ] + 1;
	
}

#endif
