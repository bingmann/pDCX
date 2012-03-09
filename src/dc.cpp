/************************************************************************************
 *      Program:    pDC3
 *      Author:     Fabian Kulla
 *      Mail:       @
 *      Description: parallel suffix array construction
**************************************************************************************/

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

#include <cassert>
#include <cmath>
#include <vector>
#include <limits>
#include <sstream>
#include <errno.h>

#include "yuta-sais-lite.h"
#include "sachecker.h"

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

// **********************************************************************
// * Loser tree implementation

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
	    m_tree[levelput + i] = m_less.done(i) ? -1 : i;

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
		else if ( m_tree[levelget + 2*i] < 0 )
		    m_tree[levelput + i] = m_tree[levelget + 2*i + 1];
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

    void replay()
    {
	int top = m_tree[0];

	int p = (m_tree.size() / 2) + top;

	if (m_less.done(top)) top = -1;	// mark sequence as done

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

template <typename Type, class Comparator>
struct MergeAreasLTCmp
{
    std::vector<Type>& 	v;
    std::vector<int>&	pos;
    const int*		endpos;
    const Comparator&	cmp;

    MergeAreasLTCmp(std::vector<Type>& _v, std::vector<int>& _pos, int* _endpos, const Comparator& _cmp)
	: v(_v), pos(_pos), endpos(_endpos), cmp(_cmp)
    {
    }

    bool done(int v) const
    {
	return pos[v] >= endpos[v+1];
    }
	
    bool operator()(int a, int b) const
    {
	return cmp( v[pos[a]], v[pos[b]] );
    }
};

/**
 * Merge an area of ordered sequences as received from other processors:
 *
 * @param v	the complete vector
 * @param area  array of indexes to the adjacent areas first position - size arealen+1 (!)
 * @param arealen number of areas.
 * @param cmp	comparator
 */

template < typename Type, class Comparator >
void merge_areas(std::vector<Type>& v, int* area, int areanum, const Comparator& cmp = std::less<Type>())
{
    std::vector<int> pos (&area[0], &area[areanum+1]);

    MergeAreasLTCmp<Type,Comparator> ltcmp(v, pos, area, cmp);
    LoserTree< MergeAreasLTCmp<Type,Comparator> > LT (areanum, ltcmp);

    std::vector<Type> out (v.size());

    int top, j = 0;

    while( (top = LT.top()) >= 0 )
    {
	out[j++] = v[ pos[top] ];

	++pos[top];

	LT.replay();
    }

    std::swap(v, out);
}

template < typename Type >
void merge_areas(std::vector<Type>& v, int* area, int areanum)
{
    return merge_areas(v, area, areanum, std::less<Type>());
}

// **********************************************************************
// * MPI Datatype Templates

template <typename Type>
struct MPI_DatatypeTemplate
{
    static MPI_Datatype getType();
};

#define MAKE_MPI_DATATYPE_TEMPLATE(Type,MpiConst)	\
    template <>						\
    MPI_Datatype MPI_DatatypeTemplate<Type>::getType()	\
    {							\
        return MpiConst;				\
    }

MAKE_MPI_DATATYPE_TEMPLATE(char, MPI_CHAR);
MAKE_MPI_DATATYPE_TEMPLATE(unsigned char, MPI_BYTE);

MAKE_MPI_DATATYPE_TEMPLATE(int, MPI_INT);
MAKE_MPI_DATATYPE_TEMPLATE(unsigned int, MPI_UNSIGNED);

template <typename Type>
MPI_Datatype getMpiDatatype(const Type& t)
{
    return MPI_DatatypeTemplate<Type>::getType();
}

#define GET_MPI_DATATYPE(Type)	MPI_DatatypeTemplate<Type>::getType()

// **********************************************************************
// * pDCX base class

struct DC3Param {
    static const unsigned int X = 3;
    static const unsigned int D = 2;

    static const unsigned int DC[D];

    static const int 	      cmpDepthRanks[X][X][3];
};

const unsigned int DC3Param::DC[] = { 1, 2 };

const int DC3Param::cmpDepthRanks[3][3][3] =
{
    { { 1,0,0 }, { 1,0,1 }, { 2,1,1 } },
    { { 1,1,0 }, { 0,0,0 }, { 0,0,0 } },
    { { 2,1,1 }, { 0,0,0 }, { 0,0,0 } },
};

struct DC7Param {
    static const unsigned int X = 7;
    static const unsigned int D = 3;

    static const unsigned int DC[D];

    static const int 	      cmpDepthRanks[X][X][3];
};

const unsigned int DC7Param::DC[] = { 0, 1, 3 };

const int DC7Param::cmpDepthRanks[7][7][3] =
{
    { { 0,0,0 }, { 0,0,0 }, { 1,1,0 }, { 0,0,0 }, { 3,2,0 }, { 3,2,1 }, { 1,1,0 } },
    { { 0,0,0 }, { 0,0,0 }, { 6,2,2 }, { 0,0,0 }, { 6,2,2 }, { 2,1,0 }, { 2,1,1 } },
    { { 1,0,1 }, { 6,2,2 }, { 1,0,0 }, { 5,1,2 }, { 6,2,2 }, { 5,1,2 }, { 1,0,0 } },
    { { 0,0,0 }, { 0,0,0 }, { 5,2,1 }, { 0,0,0 }, { 4,1,1 }, { 5,2,2 }, { 4,1,2 } },
    { { 3,0,2 }, { 6,2,2 }, { 6,2,2 }, { 4,1,1 }, { 3,0,0 }, { 3,0,1 }, { 4,1,2 } },
    { { 3,1,2 }, { 2,0,1 }, { 5,2,1 }, { 5,2,2 }, { 3,1,0 }, { 2,0,0 }, { 2,0,1 } },
    { { 1,0,1 }, { 2,1,1 }, { 1,0,0 }, { 4,2,1 }, { 4,2,1 }, { 2,1,0 }, { 1,0,0 } },
};

struct DC13Param {
    static const unsigned int X = 13;
    static const unsigned int D = 4;

    static const unsigned int DC[D];

    static const int 	      cmpDepthRanks[X][X][3];
};

const unsigned int DC13Param::DC[] = { 0, 1, 3, 9 };

const int DC13Param::cmpDepthRanks[13][13][3] =
{
    { { 0,0,0 }, { 0,0,0 }, { 1,1,0 }, { 0,0,0 }, { 9,3,1 }, { 9,3,2 }, { 3,2,0 }, { 9,3,3 }, { 1,1,0 }, { 0,0,0 }, { 3,2,0 }, { 3,2,1 }, { 1,1,0 } },
    { { 0,0,0 }, { 0,0,0 }, {12,3,3 }, { 0,0,0 }, {12,3,3 }, { 8,2,1 }, { 8,2,2 }, { 2,1,0 }, { 8,2,3 }, { 0,0,0 }, {12,3,3 }, { 2,1,0 }, { 2,1,1 } },
    { { 1,0,1 }, {12,3,3 }, { 1,0,0 }, {11,2,3 }, {12,3,3 }, {11,2,3 }, { 7,1,1 }, { 7,1,2 }, { 1,0,0 }, { 7,1,3 }, {12,3,3 }, {11,2,3 }, { 1,0,0 } },
    { { 0,0,0 }, { 0,0,0 }, {11,3,2 }, { 0,0,0 }, {10,2,2 }, {11,3,3 }, {10,2,3 }, { 6,1,1 }, { 6,1,2 }, { 0,0,0 }, { 6,1,2 }, {11,3,3 }, {10,2,3 } },
    { { 9,1,3 }, {12,3,3 }, {12,3,3 }, {10,2,2 }, { 5,0,0 }, { 9,1,2 }, {10,2,3 }, { 9,1,3 }, { 5,0,1 }, { 5,0,2 }, {12,3,3 }, { 5,0,2 }, {10,2,3 } },
    { { 9,2,3 }, { 8,1,2 }, {11,3,2 }, {11,3,3 }, { 9,2,1 }, { 4,0,0 }, { 8,1,2 }, { 9,2,3 }, { 8,1,3 }, { 4,0,1 }, { 4,0,1 }, {11,3,3 }, { 4,0,2 } },
    { { 3,0,2 }, { 8,2,2 }, { 7,1,1 }, {10,3,2 }, {10,3,2 }, { 8,2,1 }, { 3,0,0 }, { 7,1,2 }, { 8,2,3 }, { 7,1,3 }, { 3,0,0 }, { 3,0,1 }, {10,3,3 } },
    { { 9,3,3 }, { 2,0,1 }, { 7,2,1 }, { 6,1,1 }, { 9,3,1 }, { 9,3,2 }, { 7,2,1 }, { 2,0,0 }, { 6,1,2 }, { 7,2,3 }, { 6,1,2 }, { 2,0,0 }, { 2,0,1 } },
    { { 1,0,1 }, { 8,3,2 }, { 1,0,0 }, { 6,2,1 }, { 5,1,0 }, { 8,3,1 }, { 8,3,2 }, { 6,2,1 }, { 1,0,0 }, { 5,1,2 }, { 6,2,2 }, { 5,1,2 }, { 1,0,0 } },
    { { 0,0,0 }, { 0,0,0 }, { 7,3,1 }, { 0,0,0 }, { 5,2,0 }, { 4,1,0 }, { 7,3,1 }, { 7,3,2 }, { 5,2,1 }, { 0,0,0 }, { 4,1,1 }, { 5,2,2 }, { 4,1,2 } },
    { { 3,0,2 }, {12,3,3 }, {12,3,3 }, { 6,2,1 }, {12,3,3 }, { 4,1,0 }, { 3,0,0 }, { 6,2,1 }, { 6,2,2 }, { 4,1,1 }, { 3,0,0 }, { 3,0,1 }, { 4,1,2 } },
    { { 3,1,2 }, { 2,0,1 }, {11,3,2 }, {11,3,3 }, { 5,2,0 }, {11,3,3 }, { 3,1,0 }, { 2,0,0 }, { 5,2,1 }, { 5,2,2 }, { 3,1,0 }, { 2,0,0 }, { 2,0,1 } },
    { { 1,0,1 }, { 2,1,1 }, { 1,0,0 }, {10,3,2 }, {10,3,2 }, { 4,2,0 }, {10,3,3 }, { 2,1,0 }, { 1,0,0 }, { 4,2,1 }, { 4,2,1 }, { 2,1,0 }, { 1,0,0 } },
};

template <typename DCParam, typename alphabet_type>
class pDCX
{
public:

    typedef unsigned int	uint;

    // **********************************************************************
    // * global parameters

    static const unsigned int X		= DCParam::X;
    static const unsigned int D		= DCParam::D;

//static const unsigned int DCH[X] = { 3, 6, 6, 5, 6, 5, 4 };	// additional chars in tuple
//static const unsigned int DCD[X] = { 0, 0, 1, 0, 3, 2, 1 };	// depth to sort chars before using first rank

//static const unsigned int inDC[X] = { 1, 1, 0, 1, 0, 0, 0 };

    static const bool debug		= true;
    static const bool debug_sortsample	= true;
    static const bool debug_nameing	= true;
    static const bool debug_recursion	= true;
    static const bool debug_finalsort	= false;

    static const bool debug_compare	= false;
    

    // **********************************************************************
    // * tuple types

    class Pair {
    public:
	uint		name;
	uint		index;

	bool operator< ( const Pair& a ) const {
	    return index < a.index;
	}

	static inline
	bool cmpIndexModDiv( const Pair& a, const Pair& b ) {
	    return ( a.index % X < b.index % X ) ||
		( ( a.index % X == b.index % X ) && ( a.index / X < b.index / X ) );
	}

	friend std::ostream& operator<< (std::ostream& os, const Pair& p)
	{
	    return (os << "(" << p.name << "," << p.index << ")");
	}
    };

    class TupleS
    {
    public:
	alphabet_type	chars[X];
	uint		index;

	bool operator< (const TupleS& o) const
	{
	    for (unsigned int i = 0; i < X; ++i)
	    {
		if (chars[i] == o.chars[i]) continue;
		return chars[i] < o.chars[i];
	    }
	    return false;
	}

	bool operator== (const TupleS& o) const
	{
	    for (unsigned int i = 0; i < X; ++i)
	    {
		if (chars[i] != o.chars[i]) return false;
	    }
	    return true;
	}

	friend std::ostream& operator<< (std::ostream& os, const TupleS& t)
	{
	    os << "([";
	    for (unsigned int i = 0; i < X; ++i) {
		if (i != 0) os << " ";
		os << t.chars[i];
	    }
	    os << "]," << t.index << ")";
	    return os;
	}
    };

    struct TupleN
    {
	alphabet_type	chars[X-1];
	uint		ranks[D];
	uint		index;

	friend std::ostream& operator<< (std::ostream& os, const TupleN& t)
	{
	    os << "(c[";
	    for (unsigned int i = 0; i < X-1; ++i) {
		if (i != 0) os << " ";
		os << t.chars[i];
	    }
	    os << "],r[";
	    for (unsigned int i = 0; i < D; ++i) {
		if (i != 0) os << " ";
		os << t.ranks[i];
	    }
	    os << "]," << t.index << ")";
	    return os;
	}
    };

    template <int D>
    static inline bool cmpTupleNdepth(const TupleN& a, const TupleN& b)
    {
	for (unsigned int d = 0; d < D; ++d)
	{
	    if (a.chars[d] == b.chars[d]) continue;
	    return (a.chars[d] < b.chars[d]);
	}

	assert(a.ranks[0] != b.ranks[0]);

	return (a.ranks[0] < b.ranks[0]);
    }

    // **********************************************************************
    // *** MPI variables

    static const unsigned int MSGTAG	= 42;	// arbitrary number

    int				myproc;		// MPI my process number
    int			 	nprocs;		// MPI total processes

    static const int		ROOT = 0;	// arbitray master process number


    MPI_Datatype		MPI_PAIR;
    MPI_Datatype		MPI_TUPLE_SAMPLE;
    MPI_Datatype		MPI_TUPLE_NONSAMPLE;

    MPI_Status			status;

    int				samplefactor;

    void init_mpi_datatypes()
    {
#define DISP(obj,attr)	((char*)&obj.attr - (char*)&obj)

	// type: pair of indexes
	{
	    Pair 	p;

	    MPI_Datatype typelist[2] = { getMpiDatatype(p.name), getMpiDatatype(p.index) };
	    int 	blocklen[2] = { 1, 1 };
	    MPI_Aint	displace[2] = { DISP(p,name), DISP(p,index) };
	       
	    MPI_Type_struct( 2, blocklen, displace, typelist, &MPI_PAIR );
	    MPI_Type_commit( &MPI_PAIR );
	}

	// type: tuple with X chars and index
	{
	    TupleS 	t;

	    MPI_Datatype typelist[2] = { getMpiDatatype(t.chars[0]), getMpiDatatype(t.index) };
	    int 	blocklen[2] = { X, 1 };
	    MPI_Aint	displace[2] = { DISP(t,chars), DISP(t,index) };
	       
	    MPI_Type_struct( 2, blocklen, displace, typelist, &MPI_TUPLE_SAMPLE );
	    MPI_Type_commit( &MPI_TUPLE_SAMPLE );
	}

	// type: tuple with X-1 chars, D ranks and index
	{
	    TupleN 	t;

	    MPI_Datatype typelist[3] = { getMpiDatatype(t.chars[0]), getMpiDatatype(t.ranks[0]), getMpiDatatype(t.index) };
	    int 	blocklen[3] = { X-1, D, 1 };
	    MPI_Aint	displace[3] = { DISP(t,chars), DISP(t,ranks), DISP(t,index) };
	       
	    MPI_Type_struct( 3, blocklen, displace, typelist, &MPI_TUPLE_NONSAMPLE );
	    MPI_Type_commit( &MPI_TUPLE_NONSAMPLE );
	}

#undef DISP
    }

    void deinit_mpi_datatypes()
    {
	MPI_Type_free( &MPI_TUPLE_SAMPLE );
	MPI_Type_free( &MPI_PAIR );
	MPI_Type_free( &MPI_TUPLE_NONSAMPLE );
    }

    pDCX()
    {
	MPI_Comm_rank( MPI_COMM_WORLD, &myproc );
	MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

	samplefactor = nprocs;	// TODO

	init_mpi_datatypes();
    }

    ~pDCX()
    {
	deinit_mpi_datatypes();
    }

    // **********************************************************************
    // *** MPI variables

    static inline
    bool cmpTupleNCompare(const TupleN& t1, const TupleN& t2)
    {
	unsigned int v1 = t1.index % X, v2 = t2.index % X;

	const int* deprank = DCParam::cmpDepthRanks[v1][v2];

	if (debug_compare)
	    std::cout << "cmp " << v1 << "(" << t1.index << ") ? " << v2 << "(" << t2.index << ") - depth " << deprank[0] << "\n";

	for (int d = 0; d < deprank[0]; ++d)
	{
	    if (t1.chars[d] == t2.chars[d]) continue;
	    return (t1.chars[d] < t2.chars[d]);
	}

	if (debug_compare)
	    std::cout << "break tie using ranks " << deprank[1] << " - " << deprank[2] << " = " << t1.ranks[ deprank[1] ] << " - " << t2.ranks[ deprank[2] ] << "\n";

	//assert (t1.ranks[ deprank[1] ] != t2.ranks[ deprank[2] ]);

	return (t1.ranks[ deprank[1] ] < t2.ranks[ deprank[2] ]);
    }

    struct TupleNMerge
    {
	const std::vector<TupleN>*	 m_S;

	std::vector<unsigned int>	 m_ptr;

	TupleNMerge(const std::vector<TupleN>* S)
	    : m_S(S), m_ptr(X, 0)
	{
	}

	inline bool done(int v) const
	{
	    return (m_ptr[v] >= m_S[v].size());
	}

	inline bool operator()(int v1, int v2) const
	{
	    assert( v1 < v2 );

	    const int* deprank = DCParam::cmpDepthRanks[v1][v2];

	    const TupleN& t1 = m_S[v1][ m_ptr[v1] ];
	    const TupleN& t2 = m_S[v2][ m_ptr[v2] ];

	    assert( t1.index % X == (unsigned int)v1 );
	    assert( t2.index % X == (unsigned int)v2 );

	    if (debug_compare)
		std::cout << "cmp " << v1 << "(" << t1.index << ") ? " << v2 << "(" << t2.index << ") - depth " << deprank[0] << "\n";

	    for (int d = 0; d < deprank[0]; ++d)
	    {
		if (t1.chars[d] == t2.chars[d]) continue;
		return (t1.chars[d] < t2.chars[d]);
	    }

	    if (debug_compare)
		std::cout << "break tie using ranks " << deprank[1] << " - " << deprank[2] << " = " << t1.ranks[ deprank[1] ] << " - " << t2.ranks[ deprank[2] ] << "\n";

	    assert (t1.ranks[ deprank[1] ] != t2.ranks[ deprank[2] ]);

	    return (t1.ranks[ deprank[1] ] < t2.ranks[ deprank[2] ]);
	}
    };


    bool dcx( alphabet_type* input, uint globalSize, uint localStride )
    {
	uint samplesize = (uint)sqrt(localStride * D / X / nprocs) * samplefactor;

	if ( samplesize >= D * (localStride / X) ) samplesize = D * (localStride / X) - 1;

	const unsigned int* DC = DCParam::DC;

	const unsigned int globalMultipleOfX = (globalSize + X - 1) / X;	// rounded up number of global multiples of X

	const unsigned int localOffset = myproc * localStride;
	const unsigned int localSize = (myproc != nprocs-1) ? localStride : globalSize - localOffset;
	const unsigned int localSizeReal = (myproc != nprocs-1) ? localStride + (X-1) : globalSize - localOffset;

	const unsigned int M = (localSize + X - 1) / X;			// number of incomplete X chars in local area size

	if (debug)
	{
	    std::cout << "******************** DCX ********************" << std::endl;

	    std::cout << "Parameters:\n"
		      << "  globalSize = " << globalSize << "\n"
		      << "  localStride = " << localStride << "\n"
		      << "  localSize = " << localSize << "\n"
		      << "  localSizeReal = " << localSizeReal << "\n"
		      << "  globalMultipleOfX = " << globalMultipleOfX << "\n"
		      << "  localMultipleOfX (aka M) = " << M << "\n"
		      << "  samplesize = " << samplesize << "\n";
	}

	DBG_ARRAY2(debug, "Input", input, localSize);

	DBG_ARRAY2(debug, "Input (extra tuples)", (input + localSize), localSizeReal - localSize);

	// **********************************************************************
	// * calculate build DC-tuple array and sort locally

	std::vector<TupleS> R (D * M);			// create D * M tuples which might include up to D-1 dummies

	uint j = 0;
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

	assert( j == D * M );

	DBG_ARRAY(debug_sortsample, "Sample suffixes", R);

	// **********************************************************************
	// {{{ Sample sort of array R
	{
	    std::sort(R.begin(), R.end());		// sort locally

	    DBG_ARRAY(debug_sortsample, "Locally sorted sample suffixes", R);

	    // **********************************************************************
	    // * select equidistance samples and redistribute sorted DC-tuples

	    // select samples
	    TupleS* samplebuf = new TupleS[ samplesize ];

	    double dist = (double) R.size() / samplesize;
	    for ( uint i = 0; i < samplesize; i++ )
		samplebuf[i] = R[ int( i * dist ) ];

	    // root proc collects all samples
	    TupleS* samplebufall;
	    if (myproc == ROOT) samplebufall = new TupleS[ nprocs * samplesize ];

	    MPI_Gather( samplebuf, samplesize, MPI_TUPLE_SAMPLE,
			samplebufall, samplesize, MPI_TUPLE_SAMPLE, ROOT, MPI_COMM_WORLD );

	    delete[] samplebuf;

	    // root proc sorts samples as splitters

	    TupleS* splitterbuf = new TupleS[ nprocs ];

	    if ( myproc == ROOT )
	    {
		std::sort( samplebufall, samplebufall + samplesize * nprocs );

		DBG_ARRAY2(debug_sortsample, "Sample splitters", samplebufall, nprocs * samplesize);

		for ( int i = 0; i < nprocs; i++ )		// pick splitters
		    splitterbuf[i] = samplebufall[ i * samplesize ];

		DBG_ARRAY2(debug_sortsample, "Selected splitters", splitterbuf, nprocs);

		delete[] samplebufall;
	    }

	    // distribute splitters

	    MPI_Bcast( splitterbuf, nprocs, MPI_TUPLE_SAMPLE, ROOT, MPI_COMM_WORLD );

	    // find nearest splitters in locally sorted tuple list

	    uint* splitterpos = new uint[ nprocs + 1 ];

	    splitterpos[0] = 0;
	    for ( int i = 1; i < nprocs; i++ )
	    {
		typename std::vector<TupleS>::const_iterator it = std::lower_bound(R.begin(), R.end(), splitterbuf[i]);

		splitterpos[i] = it - R.begin();
	    }
	    splitterpos[ nprocs ] = R.size();

	    DBG_ARRAY2(debug_sortsample, "Splitters positions", splitterpos, nprocs+1);

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

	    std::vector<TupleS> recvbuf ( recvoff[ nprocs ] );

	    MPI_Alltoallv( R.data(), sendcnt, sendoff, MPI_TUPLE_SAMPLE,
			   recvbuf.data(), recvcnt, recvoff, MPI_TUPLE_SAMPLE, MPI_COMM_WORLD );

	    std::swap(R, recvbuf);

	    delete [] sendcnt;
	    delete [] recvcnt;
	    delete [] sendoff;

	    merge_areas(R, recvoff, nprocs);

	    delete [] recvoff;
	}
	// }}} end Sample sort of array R

	DBG_ARRAY(debug_sortsample, "Sorted sample suffixes", R);

	// R contains DC-sample tuples in sorted order

	// **********************************************************************
	// * Lexicographical naming

	std::vector<Pair> P ( R.size() );

	int recursion;

	{
	    // naming with local names

	    TupleS temp;	// get last tuple from previous process as basis for name comparison (cyclicly)
	    MPI_Sendrecv( &(R[R.size()-1]), 1, MPI_TUPLE_SAMPLE, (myproc + 1) % nprocs, MSGTAG,
			  &temp, 1, MPI_TUPLE_SAMPLE, (myproc - 1 + nprocs) % nprocs, MSGTAG,
			  MPI_COMM_WORLD, &status );

	    uint name = 0;
	    for ( uint i = 0; i < R.size(); i++ ) {
		if ( !( R[i] == temp ) ) {
		    name++;
		    if (debug_nameing)
			std::cout << "Giving name " << name << " to " << R[i] << "\n";
		    temp = R[i];
		}
		P[i].name = name;
		P[i].index = R[i].index;
	    }
	    R.clear();		// Why?: because it is easier to recreate the tuples later on

	    DBG_ARRAY(debug_nameing, "Local Names", P);

	    // renaming with global names: calculate using prefix sum

	    uint namesglob;
	    MPI_Scan( &name, &namesglob, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );

	    // update local names
	    for ( uint i = 0; i < P.size(); i++ )
		P[i].name += (namesglob - name);

	    DBG_ARRAY(debug_nameing, "Global Names", P);

	    // determine whether recursion is necessary: last proc checks its highest name

	    if (myproc == nprocs - 1)
	    {
		if (debug_nameing)
		    std::cout << "last name: " << P.back().name << " -? " << D * globalMultipleOfX << "\n";

		recursion = (P.back().name != D * globalMultipleOfX);

		if (debug_nameing)
		    std::cout << "recursion: " << recursion << "\n";
	    }

	    MPI_Bcast( &recursion, 1, MPI_INT, nprocs - 1, MPI_COMM_WORLD );
	}

	if ( recursion )
	{
	    uint namesGlobalSize = D * globalMultipleOfX;
	    uint namesLocalStride = ( namesGlobalSize + nprocs - 1 ) / nprocs;		// rounded up
	    namesLocalStride += X - namesLocalStride % X;				// round up to nearest multiple of X

	    if (namesGlobalSize > 2 * X * nprocs)
	    {
		if (debug_recursion)
		    std::cout << "---------------------   RECURSION pDCX ---------------- " << localSize << std::endl;

		// **********************************************************************
		// {{{ Sample sort of array P by (i mod X, i div X)

		std::sort(P.begin(), P.end(), Pair::cmpIndexModDiv);		// sort locally

		DBG_ARRAY(debug_recursion, "Global Names sorted cmpModDiv", P);

		uint* splitterpos = new uint[nprocs+1];
		int* sendcnt = new int[nprocs];
		int* sendoff = new int[nprocs+1];
		int* recvcnt = new int[nprocs];
		int* recvoff = new int[nprocs+1];

		// use equidistance splitters from 0..namesGlobalSize (because indexes are known in advance)
		splitterpos[0] = 0;
		Pair ptemp;
		ptemp.name=0;
		for ( int i = 1; i < nprocs; i++ ) {
		    ptemp.index = i * namesLocalStride;

		    unsigned int x = ptemp.index;

		    unsigned int divM = ptemp.index / globalMultipleOfX;
		    ptemp.index = DC[divM] + X * (ptemp.index - divM * globalMultipleOfX);

		    if (debug_recursion)
			std::cout << "splitter: " << ptemp.index << " = " << x << " - " << divM << "\n";

		    typename std::vector<Pair>::const_iterator it = std::lower_bound(P.begin(), P.end(), ptemp, Pair::cmpIndexModDiv);
		    splitterpos[i] = it - P.begin();
		}
		splitterpos[ nprocs ] = P.size();

		DBG_ARRAY2(1, "Splitters positions", splitterpos, nprocs+1);

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

		std::vector<Pair> recvBufPair ( recvoff[ nprocs ] + X-1 );
		unsigned int recvBufPairSize = recvoff[ nprocs ];

		MPI_Alltoallv( P.data(), sendcnt, sendoff, MPI_PAIR, recvBufPair.data(), recvcnt, recvoff, MPI_PAIR, MPI_COMM_WORLD );

		P.clear();

		// final X-1 tuples should be ignored due to recvoff areas
		merge_areas(recvBufPair, recvoff, nprocs, Pair::cmpIndexModDiv);

		// transfer additional X-1 names for last tuple

		MPI_Sendrecv( recvBufPair.data(), X-1, MPI_PAIR, ( myproc - 1 + nprocs ) % nprocs, MSGTAG,
			      recvBufPair.data() + recvBufPairSize, X-1, MPI_PAIR, ( myproc + 1 ) % nprocs, MSGTAG,
			      MPI_COMM_WORLD, &status );

		// TODO: merge and reduce at once

		uint* namearray = new uint[ recvBufPair.size() ];
		for (unsigned int i = 0; i < recvBufPair.size(); ++i)
		    namearray[i] = recvBufPair[i].name;

		DBG_ARRAY2(debug_recursion, "Pairs P (globally sorted by indexModDiv)", recvBufPair.data(), recvBufPairSize);

		DBG_ARRAY2(debug_recursion, "Extra pairs", (recvBufPair.data() + recvBufPairSize), recvBufPair.size() - recvBufPairSize);

		// }}} end Sample sort of array P

		recvBufPair.clear();

		pDCX<DCParam, uint> rdcx;

		rdcx.dcx( namearray, namesGlobalSize, namesLocalStride );

		//delete [] namearray;

		DBG_ARRAY(1, "Recursive localSA", rdcx.localSA);

		uint SAsize = rdcx.localSA.size();
		uint allSAsize[nprocs+1];

		MPI_Allgather( &SAsize, 1, MPI_UNSIGNED, allSAsize, 1, MPI_UNSIGNED, MPI_COMM_WORLD );

		uint sum = 0;
		for (unsigned int i = 0; i < nprocs; ++i)
		{
		    uint newsum = sum + allSAsize[i];
		    allSAsize[i] = sum;
		    sum = newsum;
		}
		allSAsize[nprocs] = sum;

		DBG_ARRAY2(1, "allSAsize", allSAsize, nprocs+1);

		// generate array of pairs (index,rank) from localSA

		P.resize( rdcx.localSA.size() );

		for (unsigned int i = 0; i < rdcx.localSA.size(); ++i)
		{
		    // generate index in ModDiv sorted input sequence

		    uint saidx = rdcx.localSA[i];

		    unsigned int divM = saidx / globalMultipleOfX;

		    uint index = DC[divM] + X * (saidx - divM * globalMultipleOfX);

		    P[i].index = index;
		    P[i].name = allSAsize[myproc] + i + 1;
		}

		goto SortP;
	    }
	    else // use sequential suffix sorter
	    {
		if (debug_recursion)
		    std::cout << "---------------------   RECURSION local sais ---------------- " << localSize << std::endl;

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
		    assert( recvoff[nprocs] == (int)namesGlobalSize );
		}

		std::vector<Pair> Pall;

		if (myproc == ROOT)
		    Pall.resize( namesGlobalSize );

		MPI_Gatherv(P.data(), P.size(), MPI_PAIR,
			    Pall.data(), recvcnt, recvoff, MPI_PAIR, ROOT, MPI_COMM_WORLD);

		if (myproc == ROOT)
		{
		    uint maxname = Pall.back().name;

		    std::sort(Pall.begin(), Pall.end(), Pair::cmpIndexModDiv);		// sort locally

		    DBG_ARRAY(debug_recursion, "Global Names sorted cmpModDiv", Pall);

		    uint* namearray = new uint[ Pall.size() ];
		    for (unsigned int i = 0; i < Pall.size(); ++i)
			namearray[i] = Pall[i].name;

		    DBG_ARRAY2(debug_recursion, "Recursion input", namearray, Pall.size());

		    int* rSA = new int [ Pall.size() ];

		    yuta_sais_lite::saisxx< uint*, int*, int >( namearray, rSA, Pall.size(), maxname+1 );

		    delete [] namearray;

		    DBG_ARRAY2(debug_recursion, "Recursive SA", rSA, Pall.size());

		    // generate rank array - same as updating pair array with correct names
		    for (unsigned int i = 0; i < Pall.size(); ++i)
		    {
			Pall[ rSA[i] ].name = i+1;
		    }

		    DBG_ARRAY(debug_recursion, "Fixed Global Names sorted cmpModDiv", Pall);

		    std::swap(P, Pall);
		}

		// **********************************************************************
		// *** redistribute pairs partitioned by localSize

		if (myproc == ROOT)
		{
		    std::sort(P.begin(), P.end());			// sort locally by index

		    DBG_ARRAY(debug_recursion, "Fixed Global Names sorted index", P);

		    uint* splitterpos = new uint[nprocs+1];

		    // use equidistance splitters from 0..globalSize (because indexes are fixed)
		    splitterpos[ 0 ] = 0;
		    Pair ptemp; ptemp.name=0;
		    for ( int i = 1; i < nprocs; i++ ) {
			ptemp.index = i * localStride;

			typename std::vector<Pair>::const_iterator it = std::lower_bound(P.begin(), P.end(), ptemp);
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

		    DBG_ARRAY2(debug_recursion, "recvcnt", recvcnt, nprocs);
		    DBG_ARRAY2(debug_recursion, "recvoff", recvoff, nprocs+1);
		}

		MPI_Bcast( recvcnt, nprocs, MPI_INT, ROOT, MPI_COMM_WORLD );

		std::vector<Pair> recvBufPair ( recvcnt[ myproc ] + D );
		unsigned int recvBufPairSize = recvcnt[ myproc ];

		MPI_Scatterv( P.data(), recvcnt, recvoff, MPI_PAIR, recvBufPair.data(), recvBufPairSize, MPI_PAIR, ROOT, MPI_COMM_WORLD );

		//delete[] P;
		P.clear();

		DBG_ARRAY2(debug_recursion, "Pairs P (globally sorted by index)", recvBufPair.data(), recvBufPairSize);

		// **********************************************************************
		// *** every PE needs D additional sample suffix ranks to complete the final tuple

		Pair temp[ D ];

		MPI_Sendrecv( recvBufPair.data(), D, MPI_PAIR, ( myproc - 1 + nprocs ) % nprocs, MSGTAG,
			      &temp, D, MPI_PAIR, ( myproc + 1 ) % nprocs, MSGTAG,
			      MPI_COMM_WORLD, &status );

		DBG_ARRAY2(debug_recursion, "Pairs P (extra tuples)", temp, D);

		if ( myproc == nprocs - 1 )	// last processor gets sentinel tuples with maximum ranks
		{
		    for (unsigned int i = 0; i < D; ++i)
		    {
			recvBufPair[ recvBufPairSize + i ].name = INT_MAX - D + i;
			recvBufPair[ recvBufPairSize + i ].index = recvBufPair[ recvBufPairSize - D ].index - DC[0] + X + DC[i];
		    }
		}
		else	// other processors get D following tuples with indexes from the DC
		{
		    for (unsigned int i = 0; i < D; ++i)
		    {
			recvBufPair[ recvBufPairSize + i ].name = temp[i].name;
			recvBufPair[ recvBufPairSize + i ].index = recvBufPair[ recvBufPairSize - D ].index - DC[0] + X + DC[i];
			assert( recvBufPair[ recvBufPairSize + i ].index == temp[i].index );
		    }
		}

		std::swap(recvBufPair, P);

		DBG_ARRAY(debug_recursion, "Pairs P (globally sorted by index + extra tuples)", P);

		delete[] recvcnt;
		delete[] recvoff;

	    } // end use sequential suffix sorter

	    if (debug_recursion)
		std::cout << "---------------------   END  RECURSION  --------------- " << std::endl;
	}
	else {
	    if (debug_recursion)
		std::cout << "---------------------   keine  Recursion---------------- " << localSize << std::endl;

	SortP:
	    // **********************************************************************
	    // *** sample sort pairs P by index

	    std::sort( P.begin(), P.end() );

	    DBG_ARRAY(debug_recursion, "Pairs P (sorted by index)", P);

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
		ptemp.index = i * localStride;	// TODO: check range (maybe index doesnt start at 0)?

		typename std::vector<Pair>::const_iterator it = std::lower_bound(P.begin(), P.end(), ptemp);
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

	    merge_areas(recvBufPair, recvoff, nprocs);

	    DBG_ARRAY2(debug_recursion, "Pairs P (globally sorted by index)", recvBufPair.data(), recvBufPairSize);

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
		    recvBufPair[ recvBufPairSize + i ].name = INT_MAX - D + i;
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

	    DBG_ARRAY(debug_recursion, "Pairs P (globally sorted by index + extra tuples)", P);

	    delete[] sendcnt;
	    delete[] sendoff;
	    delete[] recvcnt;
	    delete[] recvoff;
	}

	// P contains pairs of index and global rank: sorted by index and partitioned by localSize

	// **********************************************************************
	// *** Generate tuple arrays of samples and non-samples

	std::vector<TupleN> S [X];

	for (unsigned int k = 0; k < X; ++k)
	    S[k].resize(M);

	unsigned int dp = 0;	// running index into P incremented when crossing DC-indexes

	for (unsigned int i = 0; i < M; ++i)
	{
	    for (unsigned int k = 0; k < X; ++k)
	    {
		S[k][i].index = myproc * localStride + i * X + k;

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
		std::sort(S[k].begin(), S[k].end(), cmpTupleNdepth<6>);	// TODO: sort less
	    }

	    // select equidistant samples

	    TupleN* samplebuf = new TupleN[ X * samplesize ];

	    double dist = ( double ) M / samplesize;
	    for ( uint i = 0, j = 0; i < samplesize; i++ )
	    {
		for (unsigned int k = 0; k < X; ++k)
		{
		    samplebuf[j++] = S[k][ int( i * dist ) ];
		}
	    }

	    // root proc collects all samples
	    TupleN* samplebufall;
	    if (myproc == ROOT) samplebufall = new TupleN[ nprocs * X * samplesize ];

	    MPI_Gather( samplebuf, X * samplesize, MPI_TUPLE_NONSAMPLE, samplebufall, X * samplesize, MPI_TUPLE_NONSAMPLE, ROOT, MPI_COMM_WORLD );

	    delete[] samplebuf;

	    // root proc sorts samples as splitters

	    TupleN* splitterbuf = new TupleN[ nprocs ];

	    if ( myproc == ROOT )
	    {
		std::sort( samplebufall, samplebufall + nprocs * X * samplesize, cmpTupleNCompare );

		DBG_ARRAY2(debug_finalsort, "Sample splitters", samplebufall, nprocs * X * samplesize);

		for ( int i = 0; i < nprocs; i++ )		// pick splitters
		    splitterbuf[i] = samplebufall[ i * X * samplesize ];

		DBG_ARRAY2(debug_finalsort, "Selected splitters", splitterbuf, nprocs);

		delete[] samplebufall;
	    }

	    // distribute splitters

	    MPI_Bcast( splitterbuf, nprocs, MPI_TUPLE_NONSAMPLE, ROOT, MPI_COMM_WORLD );

	    // find nearest splitters in each of the locally sorted tuple list

	    uint** splitterpos = new uint*[X];

	    for (unsigned int k = 0; k < X; ++k)
	    {
		splitterpos[k] = new uint[ nprocs + 1 ];

		splitterpos[k][0] = 0;
		for ( int i = 1; i < nprocs; i++ )
		{
		    typename std::vector<TupleN>::const_iterator it = std::lower_bound(S[k].begin(), S[k].end(), splitterbuf[i], cmpTupleNCompare );

		    splitterpos[k][i] = it - S[k].begin();
		}
		splitterpos[k][ nprocs ] = M;
	    }

	    for (unsigned int k = 0; k < X; ++k)
	    {
		DBG_ARRAY2(debug_finalsort, "Splitters S." << k, splitterpos[k], nprocs+1);
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
		for ( int i = 1; i <= nprocs; i++ )
		{
		    sendoff[k][i] = sendoff[k][i-1] + sendcnt[k][i-1];
		    recvoff[k][i] = recvoff[k][i-1] + recvcnt[k][i-1];

		    assert(sendoff[k][i] >= 0);
		    assert(recvoff[k][i] >= 0);
		}

		std::vector<TupleN> recvbuf ( recvoff[k][ nprocs ] );

		MPI_Alltoallv( S[k].data(), sendcnt[k], sendoff[k], MPI_TUPLE_NONSAMPLE, recvbuf.data(), recvcnt[k], recvoff[k], MPI_TUPLE_NONSAMPLE, MPI_COMM_WORLD );

		std::swap(S[k], recvbuf);

		delete [] sendcnt[k];
		delete [] recvcnt[k];
		delete [] sendoff[k];

		totalsize += S[k].size();
	    }

	    delete [] sendcnt;
	    delete [] sendoff;
	    delete [] recvcnt;

	    // merge received array parts

	    for (unsigned int k = 0; k < X; ++k)
	    {
		if ( S[k].size() )
		{
		    merge_areas(S[k], recvoff[k], nprocs, cmpTupleNdepth<6>);

		    delete [] recvoff[k];
		}
	    }

	    delete [] recvoff;

	    for (unsigned int k = 0; k < X; ++k)
	    {
		DBG_ARRAY(debug_finalsort, "After samplesort S" << k, S[k]);
	    }

	    std::vector<TupleN> suffixarray (totalsize);
	    localSA.resize( totalsize );
	    int j = 0;

	    TupleNMerge tuplecmp(S);
	    LoserTree<TupleNMerge> LT(X, tuplecmp);

	    int top;

	    while( (top = LT.top()) >= 0 )
	    {
		suffixarray[j] = S[top][ tuplecmp.m_ptr[top] ];
		localSA[j] = suffixarray[j].index;

		if (debug_finalsort)
		    std::cout << "Winning tuple: " << suffixarray[j] << "\n";

		if (suffixarray[j].index < globalSize) ++j;

		tuplecmp.m_ptr[top]++;

		LT.replay();
	    }

	    DBG_ARRAY2(debug_finalsort, "Suffixarray merged", suffixarray, j);

	    localSA.resize(j);
	}

	return true;
    }

    std::vector<uint>   localSA;

    bool run(const char* filename)
    {
	// **********************************************************************
	// * Read input file size

	uint globalSize;

	if (myproc == ROOT)
	{
	    std::ifstream infile( filename );

	    if (!infile.good()) {
		perror("Cannot read input file");
		return false;
	    }

	    // determine file size
	    infile.seekg( 0, std::ios::end );
	    globalSize = infile.tellg();
	}

	MPI_Bcast( &globalSize, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD );

	// **********************************************************************
	// * Calculate local input size (= general stride)

	uint localStride = ( globalSize + nprocs - 1 ) / nprocs;	// divide by processors rounding up
	localStride += X - localStride % X;				// round up to nearest multiple of X

	assert( localStride * nprocs >= globalSize );

	if ( myproc == ROOT )
	{
	    std::cout << "Total input size = " << globalSize << " bytes. localStride = " << localStride << std::endl;
	}

	// **********************************************************************
	// * Read input file and send to other processors - with an overlap of (X-1) characters for the final tuples

	const int overlap = X-1;
	alphabet_type* input = new alphabet_type[ localStride + overlap ];

	assert( sizeof(alphabet_type) == 1 );

	if (myproc == ROOT)
	{
	    std::ifstream infile( filename, std::ios::binary );

	    // read input for processors 1 to n-1
	    for (int p = 1; p < nprocs; ++p)
	    {
		infile.seekg( p * localStride, std::ios::beg );

		uint readsize = (p != nprocs-1) ? localStride + overlap : globalSize - (p * localStride);

		std::cout << "Read for process " << p << " from pos " << p * localStride << " of length " << readsize << std::endl;

		infile.read( (char*)input, readsize );
		MPI_Send( input, readsize, MPI_CHAR, p, MSGTAG, MPI_COMM_WORLD );
	    }

	    if (!infile.good()) {
		perror("Error reading file.");
		return false;
	    }

	    // read input for processor 0 (ROOT)

	    std::cout << "Read for process 0 from pos 0 of length " << localStride + overlap << std::endl;

	    infile.seekg( 0, std::ios::beg );
	    infile.read( (char*)input, localStride + overlap );

	    if (!infile.good()) {
		perror("Error reading file.");
		return false;
	    }
	}
	else // not ROOT: receive data
	{
	    MPI_Recv( input, localStride + overlap, MPI_CHAR, ROOT, MSGTAG, MPI_COMM_WORLD, &status );
	}

	MPI_Barrier( MPI_COMM_WORLD );

	// **********************************************************************
	// * Construct suffix array recursively

	dcx( input, globalSize, localStride );

	delete [] input;

	return true;
    }

    bool writeSA(const char* filename)
    {
	std::vector<uint> allSAsize (nprocs);

	uint localSAsize = localSA.size();
	MPI_Gather( &localSAsize, 1, MPI_UNSIGNED, allSAsize.data(), 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD );

	if (myproc == ROOT)
	{
	    int fd = open(filename, O_WRONLY | O_CREAT | O_TRUNC , S_IRUSR | S_IWUSR);

	    if (fd < 0) {
		std::cout << "Error opening file: " << strerror(errno) << std::endl;
		return false;
	    }

	    // write data portion from the ROOT process
	    write(fd, localSA.data(), allSAsize[0] * sizeof(uint));

	    std::cout << "Wrote data from process 0." << std::endl;

	    // receive data from other processes

	    uint maxsize = *std::max_element(allSAsize.begin(), allSAsize.end());

	    uint* buffer = new uint[maxsize];

	    for (int p = 1; p < nprocs; p++)
	    {
		MPI_Recv( buffer, allSAsize[p], MPI_UNSIGNED, p, MSGTAG, MPI_COMM_WORLD, &status );

		ssize_t wb = write(fd, buffer, allSAsize[p] * sizeof(uint));

		if ((uint)wb != allSAsize[p] * sizeof(uint)) {
		    std::cout << "Error writing to file: " << strerror(errno) << std::endl;
		    return false;
		}

		std::cout << "Wrote data from process " << p << "." << std::endl;
	    }

	    delete[] buffer;

	    if (close(fd) != 0) {
		std::cout << "Error writing to file: " << strerror(errno) << std::endl;
		return false;
	    }
	}
	else
	{
	    MPI_Send( localSA.data(), localSA.size() , MPI_UNSIGNED, ROOT, MSGTAG , MPI_COMM_WORLD );
	}
	
	return true;
    }

    bool checkSAlocal(const char* filename)
    {
	// **********************************************************************
	// * Read input file size

	uint globalSize;

	if (myproc == ROOT)
	{
	    std::ifstream infile( filename );

	    if (!infile.good()) {
		perror("Cannot read input file");
		return false;
	    }

	    // determine file size
	    infile.seekg( 0, std::ios::end );
	    globalSize = infile.tellg();
	}

	MPI_Bcast( &globalSize, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD );

	// **********************************************************************
	// * Collect and check suffix array

	int* recvcnt = new int[nprocs];
	int* recvoff = new int[nprocs+1];

	uint localSAsize = localSA.size();
	MPI_Gather( &localSAsize, 1, MPI_INT, recvcnt, 1, MPI_INT, ROOT, MPI_COMM_WORLD );

	if (myproc == ROOT)
	{
	    recvoff[0] = 0;
	    for ( int i = 1; i <= nprocs; i++ ) {
		recvoff[i] = recvoff[i-1] + recvcnt[i-1];
	    }
	    assert( recvoff[nprocs] == (int)globalSize );
	}

	std::vector<uint> gSA;

	if (myproc == ROOT)
	    gSA.resize( globalSize );

	MPI_Gatherv(localSA.data(), localSA.size(), MPI_INT,
		    gSA.data(), recvcnt, recvoff, MPI_INT, ROOT, MPI_COMM_WORLD);

	delete [] recvcnt;
	delete [] recvoff;

	if (myproc == ROOT)
	{
	    DBG_ARRAY(1, "Suffixarray collected", gSA);

	    std::ifstream infile( filename );
	    if (!infile.good()) {
		perror("Cannot read input file");
		return -1;
	    }

	    std::vector<char> string(globalSize);

	    infile.read((char*)string.data(), globalSize);
	    if (!infile.good()) {
		perror("Cannot read input file");
		return -1;
	    }

	    std::cout << "result suffix array: \n";

	    for (unsigned int i = 0; i < gSA.size(); ++i)
	    {
		std::cout << i << " : " << gSA[i] << " : ";

		for (unsigned int j = 0; gSA[i]+j < globalSize; ++j)
		{
		    std::cout << string[gSA[i]+j] << " ";
		}

		std::cout << "\n";
	    }

	    assert( sachecker::sa_checker(string, gSA) );
	}
	return true;
    }
};

int main( int argc, char **argv )
{
    MPI_Init( &argc, &argv );

    int myproc, nprocs;

    MPI_Comm_rank( MPI_COMM_WORLD, &myproc );
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

    if ( nprocs <= 1 )
    {
	std::cerr << "Error: requires more than one MPI processor (use -np 2)." << std::endl;
	return -1;
    }

    if ( argc < 2 ) {
	std::cout << "No input file! Call using mpdrun -np 4 ./pDCX <input-file> [output-file]" << std::endl;
	return 0;
    }

    {
	pDCX<DC7Param, uint8_t> dcx;

	dcx.run(argv[1]);

	dcx.checkSAlocal(argv[1]);

	if ( argc >= 3 ) {
	    dcx.writeSA(argv[2]);
	}
    }

    MPI_Finalize();
    return 0;
}
