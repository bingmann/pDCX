#include <cassert>

#define MOD0 ~0U/3
#define MOD1 2*(~0U/3)

typedef unsigned int uint;

class Pair {
public:
	uint name;
	uint index;

	bool less ( const Pair& a ) const {
		return *this < a;
	}
	bool operator< ( const Pair& a ) const {
		return index < a.index;
		//return (index%3 < a.index%3) || ((index%3 == a.index%3) &&     (index/3 < a.index/3));
	}
	void print() const {
		std::cout << name << " " << index << " ||";
	}
};

/**
 * splitter '% /' 3 < P '%/' 3  
 */
static inline
bool cmpIndexModDiv( const Pair& splitter, const Pair& P ) {
    return ( splitter.index % 3 < P.index % 3 ) ||
           ( ( splitter.index % 3 == P.index % 3 ) && ( splitter.index / 3 < P.index / 3 ) );
}
static inline
bool cmpGreaterIndex( const Pair& a, const Pair& b ) {
    return ( a.index > b.index );
}
static inline
bool cmpLessIndex( const Pair& a, const Pair& b ) {
    return ( a.index < b.index );
}
static inline
bool cmpLeqIndex( const Pair& a, const Pair& b ) {
    return ( a.index < b.index );
}
static inline
bool cmpSplitterPair( const uint& a, const Pair& b ) {
    return ( a < b.index );
}

class Triple {
public:
	uint name[2];
	uint index;

	bool less ( const Triple& a ) const {
		return * this < a;
	}
	bool operator< ( const Triple& a ) const {
		return name[ 0 ] < a.name[ 0 ];
		//return (index%3 < a.index%3) || ((index%3 == a.index%3) &&     (index/3 < a.index/3));
	}
	void print() const {
		std::cout << name[0] << " " << name[1] << " " << index << " ||";
	}
};

class Quadruple {
public:
	uint name[ 3 ];
	uint index;

	bool less ( const Quadruple& a ) const;
	bool equal ( const Quadruple& a ) const;
	bool operator< ( const Quadruple& a ) const;
	bool operator> ( const Quadruple& a ) const;
	bool operator== ( const Quadruple& a ) const;
	bool operator<= ( const Quadruple& a ) const;

	void print() const;
	void printchar() const;
	void printcharS1() const;
};

class Quintuple {
public:
	uint name[ 4 ];
	uint index;

	bool less ( const Quintuple& a ) const;
	bool equal ( const Quintuple& a ) const;
	bool operator< ( const Quintuple& a ) const;
	bool operator== ( const Quintuple& a ) const;
	bool operator<= ( const Quintuple& a ) const;

	//bool operator>(const Quintuple& a) const;
/** SO greater ModX-splitter */
	inline bool cmpS0GreaterMod0( const Quintuple& mod0 ) const{
		return ( name[ 0 ] > mod0.name[ 0 ] ) ||
			( ( name[ 0 ] == mod0.name[ 0 ] ) && ( name[ 2 ] > mod0.name[ 2 ] ) );
	}

	inline bool cmpS0GreaterMod1( const Quintuple& mod1 ) const{
		return ( name[ 0 ] > mod1.name[ 1 ] ) ||
			( ( name[ 0 ] == mod1.name[ 1 ] ) && ( name[ 2 ] > mod1.name[ 2 ] ) );
	}
	inline bool cmpS0GreaterMod2( const Quintuple& mod2 ) const{
		return ( name[ 0 ] > mod2.name[ 1 ] ) ||
			( ( name[ 0 ] == mod2.name[ 1 ] ) &&
			  ( ( name[ 1 ] > mod2.name[ 2 ] ) || ( ( name[ 1 ] == mod2.name[ 2 ] ) && ( name[ 3 ] > mod2.name[ 3 ] ) ) ) );
	}

	inline bool cmpS0Mod0(const Quintuple& mod0 ) const{
		return ( name[ 0 ] < mod0.name[ 0 ] ) ||
			( ( name[ 0 ] == mod0.name[ 0 ] ) && ( name[ 2 ] < mod0.name[ 2 ] ) );
	}
	inline bool cmpS0Mod1( const Quintuple& mod1 ) const{
		return ( name[ 0 ] < mod1.name[ 1 ] ) ||
			( ( name[ 0 ] == mod1.name[ 1 ] ) && ( name[ 2 ] < mod1.name[ 2 ] ) );
	}
	inline bool cmpS0Mod2( const Quintuple& mod2 ) const{
		return ( name[ 0 ] < mod2.name[ 1 ] ) ||
			( ( name[ 0 ] == mod2.name[ 1 ] ) &&
			  ( ( name[ 1 ] < mod2.name[ 2 ] ) || ( ( name[ 1 ] == mod2.name[ 2 ] ) && ( name[ 3 ] < mod2.name[ 3 ] ) ) ) );
	}

/** S1 greater ModX-splitter */
	inline bool cmpS1GreaterMod0( const Quintuple& mod0 ) const{
		return ( name[ 1 ] > mod0.name[ 0 ] ) ||
			( ( name[ 1 ] == mod0.name[ 0 ] ) && ( name[ 2 ] > mod0.name[ 2 ] ) );
	}
	inline bool cmpS1GreaterMod12( const Quintuple& mod1 ) const{
		return ( name[ 0 ] > mod1.name[ 0 ] );
	}

	inline bool cmpS1Mod0( const Quintuple& mod0 ) const{
		return ( name[ 1 ] < mod0.name[ 0 ] ) ||
			( ( name[ 1 ] == mod0.name[ 0 ] ) && ( name[ 2 ] < mod0.name[ 2 ] ) );
	}
	inline bool cmpS1Mod12( const Quintuple& mod1 ) const{
		return ( name[ 0 ] < mod1.name[ 0 ] );
	}
/** S2 greater ModX-splitter */
	inline bool cmpS2GreaterMod0( const Quintuple& mod0 ) const{
		return ( name[ 1 ] > mod0.name[ 0 ] ) ||
			( ( name[ 1 ] == mod0.name[ 0 ] ) &&
			  ( ( name[ 2 ] > mod0.name[ 1 ] ) || ( ( name[ 2 ] == mod0.name[ 1 ] ) && ( name[ 3 ] > mod0.name[ 3 ] ) ) ) );
	}
	inline bool cmpS2GreaterMod12( const Quintuple& mod1 )const {
		return ( name[ 0 ] > mod1.name[ 0 ] );
	}

	inline bool cmpS2Mod0( const Quintuple& mod0 ) const{
		return ( name[ 1 ] < mod0.name[ 0 ] ) ||
			( ( name[ 1 ] == mod0.name[ 0 ] ) &&
			  ( ( name[ 2 ] < mod0.name[ 1 ] ) || ( ( name[ 2 ] == mod0.name[ 1 ] ) && ( name[ 3 ] < mod0.name[ 3 ] ) ) ) );
	}
	inline bool cmpS2Mod12( const Quintuple& mod1 )const {
		return ( name[ 0 ] < mod1.name[ 0 ] );
	}


    inline bool operator> (const Quintuple& a) const {
        if(index < MOD0){
            if(a.index<MOD0){
                return cmpS0GreaterMod0(a);
            }else if(a.index<MOD1){
                return cmpS0GreaterMod1(a);
            }else 
                return cmpS0GreaterMod2(a);
        }else if(index < MOD1){
            if(a.index<MOD0){
                return cmpS1GreaterMod0(a);
            }else
                return cmpS1GreaterMod12(a);
        }else{
            if(a.index<MOD0){
                return cmpS2GreaterMod0(a);
            }else 
                return cmpS2GreaterMod12(a);
        }
    }

	void print() const;
	void printcharS0() const;
	void printcharS2() const;
	void printS() const;
	void printS0() const;
	void printS1() const;
	void printS2() const;

};


inline bool cmpS0S1S2(const Quintuple& a, const Quintuple& b )
{
	if (a.index < MOD0) {
		if (b.index<MOD0) {
			return a.cmpS0Mod0(b);
		} else if (b.index<MOD1) {
			return a.cmpS0Mod1(b);
		} else 
			return a.cmpS0Mod2(b);
	} else if (a.index < MOD1) {
		if (b.index<MOD0) {
			return a.cmpS1Mod0(b);
		} else
			return a.cmpS1Mod12(b);
	} else {
		if (b.index<MOD0) {
			return a.cmpS2Mod0(b);
		} else 
			return a.cmpS2Mod12(b);
	}
}
/**
 * sorts elements of type T1 by cmp(type T2,type T1)-function
 * @param sortarray     :array of sorted subarrays
 * @param left          :first element (first offset-element)
 * @param right         :last element (last offset-element)
 * @param offset        :elements from offset[i] to offset[i+1] are already sorted
 * @param cmp       :T2 and T1 may be different types
 */
template <class T1, class T2> 
void mergesort( T1* sortarray, int left, int right, int* offset, bool (*cmp)(const T2&, const T1&) )
{
    if ( ( right - left ) > 1 ) {
        uint i, j, k = offset[ left ], middle;
        middle = ( left + right ) / 2;

        mergesort( sortarray, left, middle, offset, cmp );
        mergesort( sortarray, middle, right, offset, cmp );
        T1 helparray[ offset [ middle ] ];
        for ( i = offset[ left ]; i < offset[ middle ]; i++ )
            helparray[ i ] = sortarray[ i ];

        j = offset[ middle ];
        i = offset[ left ];
        while ( k < j && j < offset[ right ] ) {
            if ( cmp(helparray[ i ],sortarray[ j ]) ) sortarray[ k ] = helparray[ i++ ];
            else sortarray[ k ] = sortarray[ j++ ];
            k++;
        }

        while ( k < j ) sortarray[ k++ ] = helparray[ i++ ];
    }
}

template <class T1, class T2, class T3> 
void mergesort( T1* sortarray, T1* helparray, T3 left, T3 right, T3* offset, bool (*cmp)(const T2&, const T1&) )
{
    if ( ( right - left ) > 1 ) {
        T3 i, j, k = offset[ left ], middle;
        middle = ( left + right ) / 2;

        mergesort( sortarray, helparray,left, middle, offset, cmp );
        mergesort( sortarray, helparray, middle, right, offset, cmp );
        for ( i = offset[ left ]; i < offset[ middle ]; i++ )
            helparray[ i ] = sortarray[ i ];

        j = offset[ middle ];
        i = offset[ left ];
        while ( k < j && j < offset[ right ]) {
            if ( cmp(helparray[ i ],sortarray[ j ]) ) sortarray[ k ] = helparray[ i++ ];
            else sortarray[ k ] = sortarray[ j++ ];
            k++;
        }
        while ( k < j ) sortarray[ k++ ] = helparray[ i++ ];
    }

}

/**
  * Function finds position of splitter in array. Elements of array and splitter are compared by *cmp-function
  * @param array        : sorted array  
  * @param splitter     :
  * @param start        : first element of array which has to be compared with splitter
  * @param end      : last element of array 
  * @param cmp      : compare function
  */
template <class T1, class T2>
inline void findPos( T1* array, T2 splitter, uint& start, uint end, bool (*cmp)(const T2&, const T1&)) {
    for ( ; start < end; start++ ) 
        if(cmp(splitter,array[start])) return;
}

template <class T1, class T2>
uint binsearch(T1* array, T2 element, int right, bool (*cmp)(const T2&, const T1&)) {
    int left=0;
    uint  middle;
    while (left < right){
        middle= (left+right) / 2;
        if (cmp(element,array[middle])) left= middle+1;
        else    right= middle; 
    }
    return left;
}

/**
 * cmp-Fct is true if T2 > T1
 */
template <class T1, class T2>
inline uint findPos( T1* array, T2 splitter, uint size, bool (*cmp)(const T2&, const T1&)) {
    uint left = 0, end = size, middle, done=0 ;
    middle = size/2;//(left + size/2)>=end? end-1 :(left + size/2);
    while (!done && size > 0) {
        if(splitter.index==array[middle].index)
            done=1;
        else{
            if ( cmp(splitter,array[middle]) )
                left = middle+1;            
            size /= 2;
            middle = (left + size/2)>=end? end-1 :(left + size/2);
        }
    }

    while (!done && middle < end && cmp(splitter,array[middle]) )
        ++middle;
    return (middle>end)?end:middle;
}

template <class T1, class T2>
inline uint findpostest( T1* array, T2 splitter, uint size, bool (*cmp)(const T2&, const T1&)) {
    /*  find position of key (or next larger) in sorted array */
    std::cout<<"splitter ";
    splitter.print();

    uint left = 0, end = size, middle, done=0 ;
    middle = size/2;//(left + size/2)>=end? end-1 :(left + size/2);

    while (!done && size > 0) {
        if(splitter.index==array[middle].index){
            std::cout<<"done "<<middle<<" "<<std::endl;
            done=1;
        }else{
            std::cout<<"middle "<<middle<<" ";
            array[middle].print();
            
            if ( cmp(splitter,array[middle]) )
                left = middle+1;            
            size /= 2;
            middle = (left + size/2)>=end? end-1 :(left + size/2);
        }
    }
    std::cout <<"cmp(splitter,array[middle]) "<<cmp(splitter,array[middle])<< std::endl;
    while (!done && middle < end && cmp(splitter,array[middle]) ){
            std::cout<<"hip-hip hurra "<<middle+1<<std::endl;
        ++middle;
    }
    return (middle>end)?end:middle;
}


uint findPosSA( uint* array, uint splitter, uint size, int* sendcnt);
uint findPosSATest( uint* array, uint splitter, uint size, int* sendcnt);
void findPosPair( Pair* P, uint splitter, uint& start, uint size, uint names);

bool cmpIndexModDiv( const Pair& splitter, const Pair& P);
bool cmpSplitterPair( const uint& a, const Pair& b);

/** findPos */
bool cmpGreaterUint( const uint& a, const uint& b );
bool cmpSplitterGreaterS12( const Quadruple& splitter, const Quadruple& S12 );
bool cmpSplitterGreaterS0( const Quintuple& splitter, const Quintuple& S0 );
bool cmpS0GreaterS1( const Quintuple& S0, const Quadruple& S1 ) ;
bool cmpS0GreaterS2( const Quintuple& S0, const Quintuple& S2 ) ;

bool cmpGreaterIndex( const Pair& splitter, const Pair& P);

/** mergesort */
bool cmpSplitterLeqS12( const Quadruple& splitter, const Quadruple& S12 );
bool cmpSplitterLeqS0( const Quintuple& splitter, const Quintuple& S0 );
bool cmpS0LessS1( const Quintuple& S0, const Quadruple& S1 ) ;
bool cmpS0LessS2( const Quintuple& S0, const Quintuple& S2 ) ;

bool cmpLessIndex( const Pair& splitter, const Pair& P);
bool cmpLeqIndex( const Pair& a, const Pair& b );

bool cmpFirstNameS1( const Quadruple& a, const Quadruple& b ) ;
bool cmpFirstNameS2( const Quintuple& a, const Quintuple& b ) ;
bool cmpNameRank( const Quintuple& a, const Quintuple& b ) ;

void merge( Quintuple* S0, Quadruple* S1, Quintuple* S2, uint* SA, uint* n );

void merge(Quintuple* S, Quintuple* helparray, uint n0Len, uint size);
void sort(Quintuple* S, Quintuple* out, Quintuple* pivbuf, uint* pivpos, int* sendcnt, uint* n, uint size);

/** used for testing SortS0S1S2 mit Samples aus S1 --- 12.11.05 */
bool cmpSplitterS0( const Quadruple& a, const Quintuple& b );
bool cmpSplitterS1( const Quadruple& a, const Quadruple& b );
bool cmpSplitterS2( const Quadruple& a, const Quintuple& b ) ;


/** used for testing SortS0S1S2 mit Samples aus S0,S1,S2 --- 23.11.05 */
/** Sx > Sy */
bool cmpS0S0( const Quintuple& a, const Quintuple& b ) ;
bool cmpS0S1( const Quintuple& a, const Quadruple& b ) ;
bool cmpS0S2( const Quintuple& a, const Quintuple& b ) ;
bool cmpS1S0( const Quintuple& a, const Quintuple& b ) ;
bool cmpS1S1( const Quintuple& a, const Quadruple& b ) ;
bool cmpS1S2( const Quintuple& a, const Quintuple& b ) ;
bool cmpS2S0( const Quintuple& a, const Quintuple& b ) ;
bool cmpS2S1( const Quintuple& a, const Quadruple& b ) ;
bool cmpS2S2( const Quintuple& a, const Quintuple& b ) ;
/** Sx < Sy */
bool cmpS0S1l( const Quintuple& a, const Quintuple& b ) ;

template <class T1>
void merge2( T1* S0, T1* S1, T1* S2, T1* Samples, uint size )
{
    uint x0 = 0, x1 = 0, x2 = 0;
    uint nsa = 0;
    //  std::cout<<"x0:"<<x0<<" n0:"<<n[0]<<std::endl;
    while ( x0 != size ) {
        /*  S0[x0].print();
            S1[x1].print();
            S2[x2].print();
        */
        if ( cmpS0S1l( S0[ x0 ], S1[ x1 ] ) )                                                                   // S0 < S1
            if ( cmpS0LessS2( S0[ x0 ], S2[ x2 ] ) )
                Samples[ nsa++ ] = S0[ x0++ ];                  // S0 < S2
            else
                Samples[ nsa++ ] = S2[ x2++ ];                                                          // S2 < S0 < S1
        else if ( S1[ x1 ].name[ 0 ] < S2[ x2 ].name[ 0 ] )
            Samples[ nsa++ ] = S1[ x1++ ]; // S1 < S0 und S1 < S2
        else
            Samples[ nsa++ ] = S2[ x2++ ];                                                              // S2 < S1 < S0
        //      std::cout<<"SA["<<nsa-1<<"]="<<SA[nsa-1 ]<<std::endl;
        if (x1 > size+1 || x2 > size+1) std::cout <<x0<<" "<< x1 <<" <= "<<size+1<<" "<<x2<<" <= "<<size+1<<std::endl;
        assert( x1 <= size+1 && x2 <= size+1 );
    }
    //  std::cout<<"x1:"<<x1<<" n1:"<<n[1]<<std::endl;
    while ( x1 !=size ) { //sind noch S1-Tuple übrig
        if ( S1[ x1 ].name[ 0 ] < S2[ x2 ].name[ 0 ] )
            Samples[ nsa++ ] = S1[ x1++ ];
        else
            Samples[ nsa++ ] = S2[ x2++ ];
        //      std::cout<<"SA["<<nsa-1<<"]="<<SA[nsa-1 ]<<std::endl;
        assert( x2 <= size+1 );
    }
    //  std::cout<<"x2:"<<x2<<" n2:"<<n[2]<<std::endl;
    while ( x2 != size ) { //sind noch S2-Tuple übrig
        Samples[ nsa++ ] = S2[ x2++ ];
        //      std::cout<<"SA["<<nsa-1<<"]="<<SA[nsa-1 ]<<std::endl;
    }
}

