#include <iostream>
#include <iomanip>
#include "tuple.h"
#include <cassert>

typedef unsigned char uchar;

#define MOD0 ~0U/3
#define MOD1 2*(~0U/3)
#define SPACE 12 

/** Quadruple */
//bool Quadruple::less ( const Quadruple& a ) const { return * this < a; }
bool Quadruple::operator< ( const Quadruple& a ) const {
    if ( name[ 0 ] < a.name[ 0 ] )
        return true;
    else if ( name[ 0 ] == a.name[ 0 ] ) {
        if ( name[ 1 ] < a.name[ 1 ] )
            return true;
        else if ( name[ 1 ] == a.name[ 1 ] ) {
            if ( name[ 2 ] < a.name[ 2 ] )
                return true;
            else if ( name[ 2 ] == a.name[ 2 ] ) {
                if ( index < a.index )
                    return true;
			}
		}
	}
    return false;
}
bool Quadruple::operator> ( const Quadruple& a ) const {
    if ( name[ 0 ] > a.name[ 0 ] )
        return true;
    else if ( name[ 0 ] == a.name[ 0 ] ) {
        if ( name[ 1 ] > a.name[ 1 ] )
            return true;
        else if ( name[ 1 ] == a.name[ 1 ] ) {
            if ( name[ 2 ] > a.name[ 2 ] )
                return true;
            else if ( name[ 2 ] == a.name[ 2 ] ) {
                if ( index > a.index )
                    return true;
			}
		}
	}
    return false;
}

// bool Quadruple::equal ( const Quadruple& a ) const { return index == a.index; }
bool Quadruple::operator== ( const Quadruple& a ) const {
    if ( name[ 0 ] != a.name[ 0 ] )
        return false;
    else if ( name[ 1 ] != a.name[ 1 ] )
        return false;
    else if ( name[ 2 ] != a.name[ 2 ] )
        return false;
    return true;
}

bool Quadruple::operator<= ( const Quadruple& a ) const {
    if ( * this > a )
        return false;
    return true;
}


void Quadruple::print() const {
    for ( int i = 0;i < 3;i++ )
        std::cout << name[ i ] << " ";
    std::cout << index << std::endl;
}
void Quadruple::printchar() const {
    std::cout << ( char ) name[ 0 ] << ( char ) name[ 1 ] << ( char ) name[ 2 ] << " " << index << std::endl;
}
void Quadruple::printcharS1() const {
    std::cout << name[ 0 ] << " " << ( char ) name[ 1 ] << " " << name[ 2 ] << " " << index << std::endl;
}


/** Quintuple */
void Quintuple::print() const {
    for ( int i = 0;i < 4;i++ )
        std::cout << name[ i ] << " ";
    std::cout << index << std::endl;
}

void Quintuple::printcharS0() const {
    std::cout << ( char ) name[ 0 ] << " " << ( char ) name[ 1 ] << " " << name[ 2 ] << " " << name[ 3 ] << " " << index << std::endl;
}

void Quintuple::printcharS2() const {
    std::cout << name[ 0 ] << " " << ( char ) name[ 1 ] << " " << ( char ) name[ 2 ] << " " << name[ 3 ] << " " << index << std::endl;
}

void Quintuple::printS() const {
	if(index<MOD0)    printS0();
	else if(index<MOD1)   printS1();
	else printS2();
        
}

void Quintuple::printS0() const {
    std::cout << std::setw( SPACE )<< name[ 0 ] << std::setw( SPACE ) << name[ 1 ] << std::setw( SPACE ) << name[ 2 ] << std::setw( SPACE ) << index <<" 0"<< std::endl;
}
void Quintuple::printS1() const {
    std::cout << std::setw( SPACE )<< name[ 1 ] << std::setw( SPACE )<<" "<< std::setw( SPACE ) << name[ 0 ] << std::setw( SPACE ) << index <<" 1"<< std::endl;
}
void Quintuple::printS2() const {
    std::cout << std::setw( SPACE )<< name[ 1 ] << std::setw( SPACE )<< name[ 2 ] << std::setw( SPACE )<< name[ 0 ] << std::setw( SPACE )<< index <<" 2"<< std::endl;
}


bool Quintuple::less ( const Quintuple& a ) const {
    return * this < a;
}
bool Quintuple::operator< ( const Quintuple& a ) const {
    return ( name[ 0 ] < a.name[ 0 ] ||
             ( ( name[ 0 ] == a.name[ 0 ] ) && ( name[ 2 ] < a.name[ 2 ] ) ) );
}

bool Quintuple::equal ( const Quintuple& a ) const {
    return index == a.index;
}
bool Quintuple::operator== ( const Quintuple& a ) const {
    return index == a.index;
    //  return (( name[ 0 ] == a.name[ 0 ] ) && ( name[ 2 ] == a.name[ 2 ] ));
}
bool Quintuple::operator<= ( const Quintuple& a ) const {
    if ( *this < a )
        return true;
    if ( *this == a )
        return true;
    return false;
}




/**
 * merges S0, S1, S2 -Tuple
 * @param S0 : mod 3 = 0-Tuple
 * @param S1 : mod 3 = 1-Tuple
 * @param S2 : mod 3 = 2-Tuple
 * @param SA : suffixarray
 * @param n  : array: n0 = |mod 3 = 0-Tuple| ; n1 = |mod 3 = 1-Tuple| ; n2 = |mod 3 = 2-Tuple|  
 */
void merge( Quintuple* S0, Quadruple* S1, Quintuple* S2, uint* SA, uint* n )
{
    uint x0 = 0, x1 = 0, x2 = 0;
    uint nsa = 0;
    //  std::cout<<"x0:"<<x0<<" n0:"<<n[0]<<std::endl;
    while ( x0 != n[ 0 ] ) {
        /*  S0[x0].print();
            S1[x1].print();
            S2[x2].print();
        */
        if ( cmpS0LessS1( S0[ x0 ], S1[ x1 ] ) )                                                                // S0 < S1
            if ( cmpS0LessS2( S0[ x0 ], S2[ x2 ] ) )
                SA[ nsa++ ] = S0[ x0++ ].index;                 // S0 < S2
            else
                SA[ nsa++ ] = S2[ x2++ ].index;                                                         // S2 < S0 < S1
        else if ( S1[ x1 ].name[ 0 ] < S2[ x2 ].name[ 0 ] )
            SA[ nsa++ ] = S1[ x1++ ].index; // S1 < S0 und S1 < S2
        else
            SA[ nsa++ ] = S2[ x2++ ].index;                                                             // S2 < S1 < S0
        //      std::cout<<"SA["<<nsa-1<<"]="<<SA[nsa-1 ]<<std::endl;
        assert( x1 <= n[ 1 ] && x2 <= n[ 2 ] );
    }
    //  std::cout<<"x1:"<<x1<<" n1:"<<n[1]<<std::endl;
    while ( x1 != n[ 1 ] ) { //sind noch S1-Tuple übrig
        if ( S1[ x1 ].name[ 0 ] < S2[ x2 ].name[ 0 ] )
            SA[ nsa++ ] = S1[ x1++ ].index;
        else
            SA[ nsa++ ] = S2[ x2++ ].index;
        //      std::cout<<"SA["<<nsa-1<<"]="<<SA[nsa-1 ]<<std::endl;
        assert( x2 <= n[ 2 ] );
    }
    //  std::cout<<"x2:"<<x2<<" n2:"<<n[2]<<std::endl;
    while ( x2 != n[ 2 ] ) { //sind noch S2-Tuple übrig
        SA[ nsa++ ] = S2[ x2++ ].index;
        //      std::cout<<"SA["<<nsa-1<<"]="<<SA[nsa-1 ]<<std::endl;
    }
}


void findPosPair( Pair* P, uint splitter, uint& start, uint size, uint names ) {
    /*  find position of key (or next larger) in sorted array */
    uint left = start;
    uint end = size;
    size -= start;
    while ( size > 0 ) {
        start = ( left + size / 2 >= end ) ? end - 1 : left + size / 2;
        if ( splitter > ( P[ start ].index % 3 - 1 ) * ( names / 2 ) + P[ start ].index / 3 )
            left = start + 1;
        size /= 2;
    }

    while ( start < end && splitter > ( P[ start ].index % 3 - 1 ) * ( names / 2 ) + P[ start ].index / 3 )
        ++start;
    //  cerr<<"start "<<start<< " size "<<end<<std::endl;
    //start=(start>end)?end:start;
}


uint findPosSA( uint* array, uint splitter, uint size, int* sendcnt ) {
    uint left = 0, end = size, done = 0 ;
    size *= 2;
    uint middle = ( left + size / 2 ); //>=2*end? 2*end-1 :(left + size/2);
    while ( !done && size > 0 ) {
        if ( splitter == array[ middle ] ) {
            done = 1;
            middle++;
        } else {
            if ( splitter > array[ middle ] )
                left = middle + 1;
            size /= 2;
            middle = ( left + size / 2 ) >= 2 * end ? 2 * end - 1 : ( left + size / 2 );
        }
    }
    while ( !done && middle < 2 * end && ( splitter > array[ middle ] ) )
        ++middle;
    middle = ( middle ) % end;
    sendcnt[ middle ] ++;
    return middle;
}

uint findPosSATest( uint* array, uint splitter, uint size, int* sendcnt ) {
    uint left = 0, end = size, done = 0 ;
    uint middle = ( left + size / 2 ); //>=end? end-1 :(left + size/2);
    //  std::cout<<"splitter "<<splitter<<" array[middle] "<<array[middle]<<std::endl;
    while ( !done && size > 0 ) {
        if ( splitter == array[ middle ] ) {
            done = 1;
            middle++;
        } else {
            //          std::cout<<" array[middle] "<<array[middle]<<" middle "<<middle<<std::endl;
            if ( splitter > array[ middle ] )
                left = middle + 1;
            size /= 2;
            middle = ( left + size / 2 ) >= end ? end - 1 : ( left + size / 2 );
        }
    }
    while ( !done && middle < end && ( splitter > array[ middle ] ) )
        ++middle;
    //  std::cout<<" array[middle] "<<array[middle]<<" middle nach while "<<middle<<std::endl;

    middle = ( middle ) % end;
    //  std::cout<<" array[middle] "<<array[middle]<<" middle%end "<<middle<<std::endl;

    sendcnt[ middle ] ++;
    return middle;
}

uint findPosSATest2( uint* array, uint splitter, uint size, int* sendcnt ) {
    uint left = 0, end = size, done = 0 ;
    size *= 2;
    uint middle = ( left + size / 2 ); //>=2*end? 2*end-1 :(left + size/2);
    std::cout << "splitter " << splitter << " array[middle] " << array[ middle ] << std::endl;
    while ( !done && size > 0 ) {
        if ( splitter == array[ middle ] ) {
            done = 1;
            middle++;
        } else {
            std::cout << " array[middle] " << array[ middle ] << " middle " << middle << std::endl;
            if ( splitter > array[ middle ] )
                left = middle + 1;
            size /= 2;
            middle = ( left + size / 2 ) >= 2 * end ? 2 * end - 1 : ( left + size / 2 );
        }
    }
    while ( !done && middle < 2 * end && ( splitter > array[ middle ] ) )
        ++middle;
    std::cout << " array[middle] " << array[ middle ] << " middle nach while " << middle << std::endl;

    middle = ( middle ) % end;
    std::cout << " array[middle] " << array[ middle ] << " middle%end " << middle << std::endl;

    sendcnt[ middle ] ++;
    return middle;
}


/** find Position */

bool cmpGreaterUint( const uint& a, const uint& b ) {
    return a > b ;
}

/**
 * Splitter > S12
 */
bool cmpSplitterGreaterS12( const Quadruple& splitter, const Quadruple& S12 ) {
    return splitter > S12 ;
}
/**
 * Splitter > S12
 */
bool cmpSplitterGreaterS0( const Quintuple& splitter, const Quintuple& S0 ) {
    return ( splitter.name[ 0 ] > S0.name[ 0 ] ) ||
           ( ( splitter.name[ 0 ] == S0.name[ 0 ] ) && ( splitter.name[ 2 ] > S0.name[ 2 ] ) );
}
/**
 * S0 > S1 
 */
bool cmpS0GreaterS1( const Quintuple& S0, const Quadruple& S1 ) {
    return ( S0.name[ 0 ] > S1.name[ 1 ] ) ||
           ( ( S0.name[ 0 ] == S1.name[ 1 ] ) && ( S0.name[ 2 ] > S1.name[ 2 ] ) );
}
/**
 * S0 > S2 
 */
bool cmpS0GreaterS2( const Quintuple& S0, const Quintuple& S2 ) {
    return ( S0.name[ 0 ] > S2.name[ 1 ] ) ||
           ( ( S0.name[ 0 ] == S2.name[ 1 ] ) &&
             ( ( S0.name[ 1 ] > S2.name[ 2 ] ) || ( ( S0.name[ 1 ] == S2.name[ 2 ] ) && ( S0.name[ 3 ] > S2.name[ 3 ] ) ) ) );
}

/** mergesort */
bool cmpSplitterLeqS12( const Quadruple& splitter, const Quadruple& S12 ) {
    return ( splitter == S12 || splitter < S12 );
}
bool cmpSplitterLeqS0( const Quintuple& splitter, const Quintuple& S0 ) {
    return ( splitter.name[ 0 ] < S0.name[ 0 ] ) ||
           ( ( splitter.name[ 0 ] == S0.name[ 0 ] ) && ( splitter.name[ 2 ] <= S0.name[ 2 ] ) );
}
/**
 * S0 < S1 
 */
bool cmpS0LessS1( const Quintuple& S0, const Quadruple& S1 ) {
    return ( S0.name[ 0 ] < S1.name[ 1 ] ) ||
           ( ( S0.name[ 0 ] == S1.name[ 1 ] ) && ( S0.name[ 2 ] < S1.name[ 2 ] ) );
}
bool cmpS0LessS1( const Quintuple& S0, const Quintuple& S1 ) {
    return ( S0.name[ 0 ] < S1.name[ 1 ] ) ||
           ( ( S0.name[ 0 ] == S1.name[ 1 ] ) && ( S0.name[ 2 ] < S1.name[ 2 ] ) );
}

/**
 * S0 < S2 
 */
bool cmpS0LessS2( const Quintuple& S0, const Quintuple& S2 ) {
    return ( S0.name[ 0 ] < S2.name[ 1 ] ) ||
           ( ( S0.name[ 0 ] == S2.name[ 1 ] ) &&
             ( ( S0.name[ 1 ] < S2.name[ 2 ] ) || ( ( S0.name[ 1 ] == S2.name[ 2 ] ) && ( S0.name[ 3 ] < S2.name[ 3 ] ) ) ) );
}

/** sort Vergleichsfunktionen für Quadruple und Quintuple */
bool cmpNameRank ( const Quintuple &a, const Quintuple &b ) {
    return ( a.name[ 0 ] < b.name[ 0 ] ) ||
           ( ( a.name[ 0 ] == b.name[ 0 ] ) && ( a.name[ 2 ] < b.name[ 2 ] ) );
}
/** sort & mergesort */
bool cmpFirstNameS1( const Quadruple& a, const Quadruple& b ) {
    return ( a.name[ 0 ] < b.name[ 0 ] );

}
bool cmpFirstNameS2( const Quintuple& a, const Quintuple& b ) {
    return ( a.name[ 0 ] < b.name[ 0 ] );
}

  /** ************************ **/
 /** ************************ **/
/** used for Testing 12.11.05 */

bool cmpSplitterS0( const Quadruple& a, const Quintuple& b ) {
    return ( a.name[ 1 ] > b.name[ 0 ] ) ||
            (( a.name[ 1 ] == b.name[ 0 ] ) && ( a.name[ 2 ] > b.name[ 2 ] ));
}

bool cmpSplitterS1( const Quadruple& a, const Quadruple& b ) {
    return ( a.name[ 0 ] > b.name[ 0 ] );

}
bool cmpSplitterS2( const Quadruple& a, const Quintuple& b ) {
    return ( a.name[ 0 ] > b.name[ 0 ] );
}



 /** ************************ **/


/**
 * 
 * @param S zwei sortierte Teilarrays( mod0 und mod1/mod2)
 * @param n0Len 
 * @param size 
 */
void merge( Quintuple* S, Quintuple* helparray, uint n0Len, uint size )
{
    uint x0 = 0, x1 = n0Len, index = 0;
    while ( (x0 != n0Len) && (x1 != size) ) {

        if ( S[ x1 ].index < MOD1 )   //mod 1 Tuple
            if ( cmpS0LessS1( S[ x0 ], S[ x1 ] ) )
                helparray[ index++ ] = S[ x0++ ];
            else
                helparray[ index++ ] = S[ x1++ ];
        else                                    //mod 2 Tuple
            if ( cmpS0LessS2( S[ x0 ], S[ x1 ] ) )
                helparray[ index++ ] = S[ x0++ ];
            else
                helparray[ index++ ] = S[ x1++ ];
//          std::cout<<" x0 "<<x0<<" x1 "<<x1<<" n0Len-"<<n0Len<<std::endl;

    }
    while ( x0 < n0Len) {
        helparray[ index++ ] = S[ x0++ ];
//          std::cout<<" x0 "<<x0<<" x1 "<<x1<<" n0Len-"<<n0Len<<std::endl;
    }

    while ( x1 < size ) {
        helparray[ index++ ] = S[ x1++ ];
//          std::cout<<" x0 "<<x0<<" x1 "<<x1<<" n0Len-"<<n0Len<<std::endl;
    }
}


uint findRecvCPU_S0( Quintuple* splitter, Quintuple element, uint size, int* sendcnt )
{
    uint left = 0, end = size, done = 0 ;
    uint middle = ( left + size / 2 );
    while ( !done && size > 0 ) {
/*      std::cout<<"element "<<element<<" splitter[middle] "<<splitter[middle]<<" size "<<size <<std::endl;*/
//  std::cout <<"element ";  element.print();
//  std::cout <<" + splitter "; splitter[middle].print();

        if ( element.index == splitter[ middle ].index ) {
            done = 1;
//          middle++;
        } else {
            if ( splitter[ middle ].index < MOD0 ) {                //mod0-splitter
                if ( element.cmpS0GreaterMod0( splitter[ middle ] ) )               left = middle + 1;
            } else  if ( splitter[ middle ].index < MOD1 ) {    //mod1-splitter
                if ( element.cmpS0GreaterMod1( splitter[ middle ] ) )               left = middle + 1;
            } else                                                          //mod2-splitter
                if ( element.cmpS0GreaterMod2( splitter[ middle ] ) )               left = middle + 1;

            size /= 2;
            middle = ( left + size / 2 ) >= end ? end - 1 : ( left + size / 2 );
        }
    }
//  std::cout<<"middle "<<middle<<" end "<<end<<std::endl;
    assert (middle<end);
    
    while ( !done && middle < end-1 ) {
        if ( splitter[ middle ].index < MOD0 ) {                //mod0-splitter
            if ( element.cmpS0GreaterMod0( splitter[ middle ] ) )           middle++;
            else done=1;
        } else  if ( splitter[ middle ].index < MOD1 ) {    //mod1-splitter
            if ( element.cmpS0GreaterMod1( splitter[ middle ] ) )           middle++;
            else done=1;
        } else{                                                             //mod2-splitter
            if ( element.cmpS0GreaterMod2( splitter[ middle ] ) )           middle++;
            else done=1;
        }
    }
    //  std::cout<<" array[middle] "<<array[middle]<<" middle nach while "<<middle<<std::endl;

    //  middle=(middle)%end;
//  std::cout <<"Eingabe ";  element.print();
//  std::cout <<"splitter "; splitter[middle].print();
//  std::cout<<"element "<<element.index<<" array[middle] "<<splitter[middle].index<<" middle "<<middle<<" "<<end<<std::endl;
    assert (middle<end);
    sendcnt[ middle ] ++;
    return middle;
}

uint findRecvCPU_S1( Quintuple* splitter, Quintuple element, uint size, int* sendcnt )
{
    uint left = 0, end = size, done = 0 ;
    uint middle = ( left + size / 2 );
    //  std::cout<<"splitter "<<splitter<<" array[middle] "<<array[middle]<<std::endl;
    while ( !done && size > 0 ) {
        if ( element.index == splitter[ middle ].index ) {
            done = 1;
//          middle++;
        } else {
            if ( splitter[ middle ].index < MOD0 ) {                //mod0-splitter
                if ( element.cmpS1GreaterMod0( splitter[ middle ] ) )               left = middle + 1;
            } else                                                          //mod12-splitter
                if ( element.cmpS1GreaterMod12( splitter[ middle ] ) )              left = middle + 1;

            //          std::cout<<" array[middle] "<<array[middle]<<" middle "<<middle<<std::endl;
            size /= 2;
            middle = ( left + size / 2 ) >= end ? end - 1 : ( left + size / 2 );
        }
    }
    while ( !done && middle < end -1) {
        if ( splitter[ middle ].index < MOD0 ) {                //mod0-splitter
            if ( element.cmpS1GreaterMod0( splitter[ middle ] ) )           middle++;
            else done=1;
        } else{                                                             //mod12-splitter
            if ( element.cmpS1GreaterMod12( splitter[ middle ] ) )          middle++;
            else done=1;
        }
    }
    assert (middle<end);
    sendcnt[ middle ] ++;
    return middle;
}

uint findRecvCPU_S2( Quintuple* splitter, Quintuple element, uint size, int* sendcnt )
{
    uint left = 0, end = size, done = 0 ;
    uint middle = ( left + size / 2 );
    while ( !done && size > 0 ) {
        if ( element.index == splitter[ middle ].index ) {
            done = 1;
//          middle++;
        } else {
            if ( splitter[ middle ].index < MOD0 ) {                //mod0-splitter
                if ( element.cmpS2GreaterMod0( splitter[ middle ] ) )
                    left = middle + 1;
            } else                                                          //mod12-splitter
                if ( element.cmpS2GreaterMod12( splitter[ middle ] ) )
                    left = middle + 1;

            size /= 2;
            middle = ( left + size / 2 ) >= end ? end - 1 : ( left + size / 2 );
        }
    }

//  std::cout <<"VOR 2. WHILE element "; element.print();
//  std::cout <<" + splitter "; splitter[middle].print();
//  std::cout<<"element "<<element.index<<" array[middle] "<<splitter[middle].index<<" middle "<<middle<<" "<<end<<std::endl;

    while ( !done && middle < end -1) {
        if ( splitter[ middle ].index < MOD0 ) {                //mod0-splitter
            if ( element.cmpS2GreaterMod0( splitter[ middle ] ) )       middle++;
            else done=1;
        } else{                                                             //mod12-splitter
            if ( element.cmpS2GreaterMod12( splitter[ middle ] ) )          middle++;
            else done=1;
        }
    }
/*  std::cout <<"NACH 2. WHILE element ";    element.print();
    std::cout <<" + splitter "; splitter[middle].print();
    std::cout<<"element "<<element.index<<" array[middle] "<<splitter[middle].index<<" middle "<<middle<<" "<<end<<std::endl;*/
    assert (middle<end);
    sendcnt[ middle ] ++;
    return middle;
}

/**
 * Array zum versenden erstellen (out)
 * sendcnt, pivpos out-Parameter
 */
void sort( Quintuple* S, Quintuple* out, Quintuple* pivbuf, uint* pivpos, int* sendcnt, uint* n, uint size )
{
    uchar * positions = new uchar[ n[ 0 ] + n[ 1 ] + n[ 2 ] ];
    uint i = 0;
    for ( ; i < n[ 0 ] ; i++ )                              //finde Empfänger von mod0
        positions[ i ] = ( uchar ) findRecvCPU_S0( pivbuf, S[ i ], size, sendcnt );
    for ( ; i < n[ 0 ] + n[ 1 ] ; i++ )                 //finde Empfänger von mod1
        positions[ i ] = ( uchar ) findRecvCPU_S1( pivbuf, S[ i ], size, sendcnt );
    for ( ; i < n[ 0 ] + n[ 1 ] + n[ 2 ] ; i++ )    //finde Empfänger von mod2
        positions[ i ] = ( uchar ) findRecvCPU_S2( pivbuf, S[ i ], size, sendcnt );
//  std::cout <<" sendcnt ";
//  for (int i=0 ; i<size ; i++)std::cout<<sendcnt[i]<<" ";
//  std::cout << std::endl;
    pivpos[ 0 ] = 0;
    uint counter[ size ];
    for ( i = 0 ; i < size ; i++ ) {
        pivpos[ i + 1 ] = pivpos[ i ] + sendcnt[ i ];
        counter[ i ] = 0;
//      std::cout<<pivpos[i+1];
    }

//  std::cout <<"alles ok  -------------------------------------------------- "<< std::endl;
    for ( i = 0 ; i < n[ 0 ] + n[ 1 ] + n[ 2 ] ; i++ ){
//      std::cout << pivpos[ positions[ i ] ] <<" "<< counter[ positions[ i ]]<<" "<<size<< std::endl;
        out[ pivpos[ positions[ i ] ] + counter[ positions[ i ] ] ++ ] = S[ i ];
    }
    delete[] positions;
}





/** used for testing SortS0S1S2 mit Samples aus S0,S1,S2 --- 23.11.05 */
/** Sx > Sy */
bool cmpS0S0( const Quintuple& a, const Quintuple& b ) {
    return ( a.name[ 0 ] > b.name[ 0 ] ) ||
            ( ( a.name[ 0 ] == b.name[ 0 ] ) && ( a.name[ 2 ] > b.name[ 2 ] ) );
}
bool cmpS0S1( const Quintuple& a, const Quadruple& b ) {
    return ( a.name[ 0 ] > b.name[ 1 ] ) ||
            ( ( a.name[ 0 ] == b.name[ 1 ] ) && ( a.name[ 2 ] > b.name[ 2 ] ) );
}
bool cmpS0S2( const Quintuple& a, const Quintuple& b ) {
    return  ( a.name[ 0 ] >  b.name[ 1 ] ) ||
            ( ( a.name[ 0 ] == b.name[ 1 ] ) &&
            ( ( a.name[ 1 ] >  b.name[ 2 ] ) || ( ( a.name[ 1 ] == b.name[ 2 ] ) && ( a.name[ 3 ] > b.name[ 3 ] ) ) ) );
}
bool cmpS1S0( const Quintuple& a, const Quintuple& b ) {
    return ( a.name[ 1 ] > b.name[ 0 ] ) ||
            ( ( a.name[ 1 ] == b.name[ 0 ] ) && ( a.name[ 2 ] > b.name[ 2 ] ) );
}
bool cmpS1S1( const Quintuple& a, const Quadruple& b ) {
    return ( a.name[ 0 ] > b.name[ 0 ] );
}
bool cmpS1S2( const Quintuple& a, const Quintuple& b ) {
    return ( a.name[ 0 ] > b.name[ 0 ] );
}
bool cmpS2S0( const Quintuple& a, const Quintuple& b ) {
    return  ( a.name[ 1 ] >  b.name[ 0 ] ) ||
            ( ( a.name[ 1 ] == b.name[ 0 ] ) &&
            ( ( a.name[ 2 ] >  b.name[ 1 ] ) || ( ( a.name[ 2 ] == b.name[ 1 ] ) && ( a.name[ 3 ] > b.name[ 3 ] ) ) ) );
}
bool cmpS2S1( const Quintuple& a, const Quadruple& b ) {
    return ( a.name[ 0 ] > b.name[ 0 ] );
}
bool cmpS2S2( const Quintuple& a, const Quintuple& b ) {
    return ( a.name[ 0 ] > b.name[ 0 ] );
}

/** Sx < Sy */
bool cmpS0S1l( const Quintuple& a, const Quintuple& b ) {
    return ( a.name[ 0 ] < b.name[ 1 ] ) ||
            ( ( a.name[ 0 ] == b.name[ 1 ] ) && ( a.name[ 2 ] < b.name[ 2 ] ) );
}

