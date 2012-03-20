/**
 *
 * pDCX
 *
 * MPI-distributed and parallel suffix sorter using difference cover.
 *
 * Written by Timo Bingmann in 2012 loosely based on the previous work
 * by Fabian Kulla in 2006.
 *
 */

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

#include <proc/readproc.h>

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

unsigned int getmemusage()
{
    struct proc_t usage;
    look_up_our_self(&usage);
    return usage.vsize;
}

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

template<typename Type>
std::string strC(const Type& t)
{
    std::ostringstream os;
    os << t;
    return os.str();
}

template<>
std::string strC(const char& c)
{
    std::ostringstream os;
    if (isprint(c) && !isspace(c)) os << (char)c;
    else os << (int)c;
    return os.str();
}

template<>
std::string strC(const unsigned char& c)
{
    std::ostringstream os;
    if (isprint(c) && !isspace(c)) os << (char)c;
    else os << (int)c;
    return os.str();
}

template <typename Type>
void vector_free(std::vector<Type>& v)
{
    std::vector<Type> v2;
    std::swap(v,v2);
}

template <typename Type>
void exclusive_prefixsum(Type* array, unsigned int size)
{
    uint sum = 0;
    for (unsigned int i = 0; i < size; ++i)
    {
	uint newsum = sum + array[i];
	array[i] = sum;
	sum = newsum;
    }
    array[size] = sum;
}

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
    static const bool debug_input	= false;
    static const bool debug_rebalance	= false;
    static const bool debug_sortsample	= false;
    static const bool debug_nameing	= false;
    static const bool debug_recursion	= false;
    static const bool debug_rec_selfcheck = false;
    static const bool debug_finalsort	= false;

    static const bool debug_compare	= false;

    static const bool debug_checker1	= false;
    static const bool debug_checker2	= false;

    static const bool debug_output	= false;

    // **********************************************************************
    // * tuple types

    class Pair {
    public:
	uint		index;
	uint		name;
	unsigned char	unique;

	bool operator< ( const Pair& a ) const {
	    return (index < a.index);
	}

	static inline
	bool cmpName( const Pair& a, const Pair& b ) {
	    return (a.name < b.name);
	}

	static inline
	bool cmpIndexModDiv( const Pair& a, const Pair& b ) {
	    return ( a.index % X < b.index % X ) ||
		( ( a.index % X == b.index % X ) && ( a.index / X < b.index / X ) );
	}

	friend std::ostream& operator<< (std::ostream& os, const Pair& p)
	{
	    return (os << "(" << p.index << "," << p.name << "," << int(p.unique) << ")");
	}
    };

    class Triple {
    public:
	uint		rank1;
	uint		rank2;
	alphabet_type	char1;

	bool operator< ( const Triple& a ) const {
	    return (rank1 < a.rank1);
	}

	friend std::ostream& operator<< (std::ostream& os, const Triple& p)
	{
	    return (os << "(" << p.rank1 << "," << p.rank2 << "," << strC(p.char1) << ")");
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
		os << strC(t.chars[i]);
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
		os << strC(t.chars[i]);
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

    template <typename Tuple>
    void radixsort_CE0(Tuple* array, uint size, size_t depth, size_t K)
    {
	if (size < 32) {
	    //insertion_sort(strings, n, depth);
	    std::sort(array, array + size);
	    return;
	}

	size_t bucketsize[K];
	memset(bucketsize, 0, sizeof(bucketsize));
	for (size_t i=0; i < size; ++i)
	    ++bucketsize[ array[i].chars[depth] ];

	size_t bucketindex[K];
	{
	    std::vector<Tuple> sorted (size);

	    bucketindex[0] = 0;
	    for (size_t i=1; i < K; ++i)
		bucketindex[i] = bucketindex[i-1] + bucketsize[i-1];

	    for (size_t i=0; i < size; ++i)
		sorted[ bucketindex[ array[i].chars[depth] ]++ ] = array[i];

	    memcpy(array, sorted.data(), size * sizeof(Tuple));
	}

	//if (depth == 6) return;
	//if (depth == sizeof(array[0].chars) / sizeof(array[0].chars[0])) return;

	size_t bsum = 0;
	for (size_t i=0; i < K; ++i) {
	    if (bucketsize[i] == 0) continue;
	    radixsort_CE0(array + bsum, bucketsize[i], depth+1, K);
	    bsum += bucketsize[i];
	}
    }

    template <int Depth>
    static inline bool cmpTupleNdepth(const TupleN& a, const TupleN& b)
    {
	for (unsigned int d = 0; d < Depth; ++d)
	{
	    if (a.chars[d] == b.chars[d]) continue;
	    return (a.chars[d] < b.chars[d]);
	}

	// ranks must always differ, however for some reason a == b is possible.
	assert(a.ranks[0] != b.ranks[0] || a.index == b.index);

	return (a.ranks[0] < b.ranks[0]);
    }

    // **********************************************************************
    // *** MPI variables

    static const unsigned int MSGTAG	= 42;	// arbitrary number

    int				myproc;		// MPI my process number
    int			 	nprocs;		// MPI total processes

    static const int		ROOT = 0;	// arbitray master process number

    MPI_Datatype		MPI_PAIR;
    MPI_Datatype		MPI_TRIPLE;
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

	    MPI_Datatype typelist[3] = { getMpiDatatype(p.name), getMpiDatatype(p.index), getMpiDatatype(p.unique) };
	    int 	blocklen[3] = { 1, 1, 1 };
	    MPI_Aint	displace[3] = { DISP(p,name), DISP(p,index), DISP(p,unique) };

	    MPI_Type_struct( 3, blocklen, displace, typelist, &MPI_PAIR );
	    MPI_Type_commit( &MPI_PAIR );
	}

	// type: triple of 2 indexes + 1 char
	{
	    Triple	t;

	    MPI_Datatype typelist[3] = { getMpiDatatype(t.rank1), getMpiDatatype(t.rank2), getMpiDatatype(t.char1) };
	    int 	blocklen[3] = { 1, 1, 1 };
	    MPI_Aint	displace[3] = { DISP(t,rank1), DISP(t,rank2), DISP(t,char1) };

	    MPI_Type_struct( 3, blocklen, displace, typelist, &MPI_TRIPLE );
	    MPI_Type_commit( &MPI_TRIPLE );
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
	MPI_Type_free( &MPI_PAIR );
	MPI_Type_free( &MPI_TRIPLE );
	MPI_Type_free( &MPI_TUPLE_SAMPLE );
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

    template <typename Type>
    void gather_vector(const std::vector<Type>& v, std::vector<Type>& out, unsigned int removelap = 0, MPI_Datatype mdt = GET_MPI_DATATYPE(Type))
    {
	int size = v.size() - (myproc != nprocs-1 ? removelap : 0);

	int* recvcnt = new int[nprocs];
	int* recvoff = new int[nprocs+1];

	MPI_Gather( &size, 1, MPI_INT, recvcnt, 1, MPI_INT, ROOT, MPI_COMM_WORLD );

	if (myproc == ROOT)
	{
	    recvoff[0] = 0;
	    for ( int i = 1; i <= nprocs; i++ ) {
		recvoff[i] = recvoff[i-1] + recvcnt[i-1];
	    }

	    out.resize( recvoff[nprocs] );
	}

	MPI_Gatherv((Type*)v.data(), size, mdt,
		    out.data(), recvcnt, recvoff, mdt, ROOT, MPI_COMM_WORLD);
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
    
    // functions for rebalancing the input
    static inline uint RangeFix(uint a, uint b, uint limit)
    {
	if (b >= a) return 0;
	return std::min<uint>(limit, a - b);
    }

    // functions for rebalancing the input
    inline uint Extra(int i)
    {
	return (i != nprocs-1) ? (X-1) : 0;
    }

    bool dcx( std::vector<alphabet_type>& input, uint depth, uint K )
    {
	const unsigned int* DC = DCParam::DC;

	// **********************************************************************
        // * analyze input and rebalance input to localStride, which is a multiple of p and X.

	// collect all localSizes and calc prefix sum

	unsigned int localSize = input.size();

	unsigned int localSizes[nprocs+1];

	MPI_Allgather( &localSize, 1, MPI_UNSIGNED, localSizes, 1, MPI_UNSIGNED, MPI_COMM_WORLD );

	exclusive_prefixsum(localSizes, nprocs);

	DBG_ARRAY2(debug_rebalance, "localSizes", localSizes, nprocs+1);

	// calculate localStride

	const uint globalSize = localSizes[nprocs];			// global size of input

	uint localStride = ( globalSize + nprocs - 1 ) / nprocs;	// divide by processors rounding up
	localStride += X - localStride % X;				// round up to nearest multiple of X

	const unsigned int localOffset = myproc * localStride;
	localSize = (myproc != nprocs-1) ? localStride : globalSize - localOffset;	// target localSize (without extra tuples)
	const unsigned int localSizeExtra = (myproc != nprocs-1) ? localStride + (X-1) : globalSize - localOffset;	// target localSize with extra tuples

	const unsigned int globalMultipleOfX = (globalSize + X - 1) / X;	// rounded up number of global multiples of X
	const unsigned int M = (localSize + X - 1) / X;			// number of incomplete X chars in local area size

	uint samplesize = (uint)sqrt(localStride * D / X / nprocs) * samplefactor;
	if ( samplesize >= D * (localStride / X) ) samplesize = D * (localStride / X) - 1;

	if (debug)
	{
	    std::cout << "******************** DCX (process " << myproc << ") depth " << depth << " ********************" << std::endl;

	    std::cout << "Parameters:\n"
		      << "  globalSize = " << globalSize << "\n"
		      << "  localStride = " << localStride << "\n"
		      << "  localSize = " << localSize << "\n"
		      << "  localSizeExtra = " << localSizeExtra << "\n"
		      << "  globalMultipleOfX = " << globalMultipleOfX << "\n"
		      << "  localMultipleOfX (aka M) = " << M << "\n"
		      << "  samplesize = " << samplesize << "\n"
		      << "  current memusage = mem " << getmemusage() << "\n";
	}

	// rebalance input

	{
	    int* sendcnt = new int[ nprocs ];
	    int* sendoff = new int[ nprocs ];
	    int* recvcnt = new int[ nprocs ];
	    int* recvoff = new int[ nprocs ];

	    sendoff[0] = 0;
	    for (int i = 1; i < nprocs; ++i)
	    {
		if (debug_rebalance)
		{
		    std::cout << "range sent " << myproc << " -> " << i << " is "
			      << RangeFix(i * localStride, localSizes[myproc], input.size()) << " - "
			      << RangeFix((i+1) * localStride + Extra(i), localSizes[myproc], input.size()) << "\n";
		}

		sendoff[i] = RangeFix(i * localStride, localSizes[myproc], input.size());
		sendcnt[i-1] = RangeFix(i * localStride + Extra(i-1), localSizes[myproc], input.size()) - sendoff[i-1];
	    }
	    sendcnt[nprocs-1] = input.size() - sendoff[nprocs-1];
	    
	    DBG_ARRAY2(debug_rebalance, "sendcnt", sendcnt, nprocs);
	    DBG_ARRAY2(debug_rebalance, "sendoff", sendoff, nprocs);

	    recvoff[0] = 0;
	    for (int i = 1; i < nprocs; ++i)
	    {
		if (debug_rebalance)
		{
		    std::cout << "range recv " << i << " -> " << myproc << " is "
			      << RangeFix(localSizes[i], myproc * localStride, localSizeExtra) << "\n"
			      << RangeFix(localSizes[i+1], myproc * localStride, localSizeExtra) << "\n";
		}

		recvoff[i] = RangeFix(localSizes[i], myproc * localStride, localSizeExtra);
		recvcnt[i-1] = RangeFix(localSizes[i], myproc * localStride, localSizeExtra) - recvoff[i-1];
	    }
	    recvcnt[nprocs-1] = localSizeExtra - recvoff[nprocs-1];

	    DBG_ARRAY2(debug_rebalance, "recvcnt", recvcnt, nprocs);
	    DBG_ARRAY2(debug_rebalance, "recvoff", recvoff, nprocs);

	    std::vector<alphabet_type> recvbuf (localSizeExtra);

	    MPI_Alltoallv( input.data(), sendcnt, sendoff, GET_MPI_DATATYPE(alphabet_type),
			   recvbuf.data(), recvcnt, recvoff, GET_MPI_DATATYPE(alphabet_type), MPI_COMM_WORLD );
	    
	    std::swap(input, recvbuf);
	}

	DBG_ARRAY2(debug_input, "Input (without extra tuples)", input.data(), localSize);

	DBG_ARRAY2(debug_input, "Input (extra tuples)", (input.data() + localSize), localSizeExtra - localSize);

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
		    R[j].chars[y] = (x < localSizeExtra) ? input[x] : 0;

		++j;
	    }
	}

	assert( j == D * M );

	std::cout << "done local sort sample suffixes - mem = " << getmemusage() << "\n";
	std::cout << "sizeof R = " << R.size() * sizeof(R[0]) << " - " << R.capacity() * sizeof(R[0]) << "\n";

	DBG_ARRAY(debug_sortsample, "Sample suffixes", R);

	// **********************************************************************
	// {{{ Sample sort of array R
	{
	    // sort locally
	    if (K < 4096 && 0)
		radixsort_CE0(R.data(), R.size(), 0, K);
	    else
		std::sort(R.begin(), R.end());
	    
	    std::cout << "done local sort sample suffixes\n";

	    //DBG_ARRAY(1, "Locally sorted sample suffixes", R);

	    // **********************************************************************
	    // * select equidistance samples and redistribute sorted DC-tuples

	    // select samples
	    TupleS* samplebuf = new TupleS[ samplesize ];

	    double dist = (double) R.size() / samplesize;
	    for ( uint i = 0; i < samplesize; i++ )
		samplebuf[i] = R[ int( i * dist ) ];

	    // root proc collects all samples
	    TupleS* samplebufall = NULL;
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
	std::cout << "done global sort sample suffixes - mem = " << getmemusage() << "\n";

	// R contains DC-sample tuples in sorted order

	// **********************************************************************
	// * Lexicographical naming

	std::vector<Pair> P ( R.size() );

	uint lastname, recursion;

	{
	    // naming with local names

	    unsigned int dupnames = 0;

	    TupleS temp;	// get last tuple from previous process as basis for name comparison (cyclicly)
	    MPI_Sendrecv( &(R[R.size()-1]), 1, MPI_TUPLE_SAMPLE, (myproc + 1) % nprocs, MSGTAG,
			  &temp, 1, MPI_TUPLE_SAMPLE, (myproc - 1 + nprocs) % nprocs, MSGTAG,
			  MPI_COMM_WORLD, &status );

	    uint name = 0, unique = 0;
	    for ( uint i = 0; i < R.size(); i++ ) {
		if ( !( R[i] == temp ) ) {
		    name++;
		    if (debug_nameing)
			std::cout << "Giving name " << name << " to " << R[i] << "\n";
		    temp = R[i];
		    unique = 1;
		}
		else {
		    dupnames++;
		    if (i != 0) P[i-1].unique = 0;
		    unique = 0;
		}
		P[i].name = name;
		P[i].index = R[i].index;
		P[i].unique = unique;
	    }
	    vector_free(R);		// Why?: because it is easier to recreate the tuples later on

	    std::cout << "given: dupnames " << dupnames << " names given: " << name << " total: " << P.size() << "\n";

	    DBG_ARRAY(debug_nameing, "Local Names", P);

	    // renaming with global names: calculate using prefix sum

	    uint namesglob;
	    MPI_Scan( &name, &namesglob, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );

	    // update local names - and free up first D names for sentinel ranks
	    for ( uint i = 0; i < P.size(); i++ )
		P[i].name += (namesglob - name) + D;

	    DBG_ARRAY(debug_nameing, "Global Names", P);

	    // determine whether recursion is necessary: last proc broadcasts highest name

	    if (myproc == nprocs-1)
		lastname = P.back().name;

	    MPI_Bcast( &lastname, 1, MPI_INT, nprocs - 1, MPI_COMM_WORLD );

	    if (1 || debug_nameing)
		std::cout << "last name: " << lastname << " =? " << D * globalMultipleOfX + D << "\n";
	    
	    recursion = (lastname != D * globalMultipleOfX + D);

	    if (1 || debug_nameing)
		std::cout << "recursion: " << recursion << "\n";
	}

	std::cout << "done naming - mem = " << getmemusage() << "\n";

	if ( recursion )
	{
	    uint namesGlobalSize = D * globalMultipleOfX + D;				// add D dummies separating mod-X areas
	    uint namesLocalStride = ( namesGlobalSize + nprocs - 1 ) / nprocs;		// rounded up
	    namesLocalStride += X - namesLocalStride % X;				// round up to nearest multiple of X
	    uint namesGlobalMultipleOfX = globalMultipleOfX + 1;			// account one extra X-tuple for D separation dummies

	    std::cout << "namesGlobalSize = " << namesGlobalSize << "\n"
		      << "namesLocalStride = " << namesLocalStride << "\n";

	    if (myproc == nprocs-1)
	    {
		for (unsigned int i = 0; i < D; ++i)
		{
		    Pair x; x.index = globalMultipleOfX * X + DC[i]; x.name = D-1-i;
		    P.push_back(x);
		}
	    }

	    if (debug_recursion)
	    {
		std::vector<Pair> Pall;
		gather_vector(P, Pall, 0, MPI_PAIR);

		if (myproc == ROOT)
		{
		    std::sort(Pall.begin(), Pall.end(), Pair::cmpIndexModDiv);		// sort locally

		    DBG_ARRAY(debug_recursion, "Pall", Pall);
		}
	    }

	    if (namesGlobalSize > 2 * X * nprocs)
	    {
		if (debug_recursion)
		    std::cout << "---------------------   RECURSION pDCX ----------------  - mem = " << getmemusage() << std::endl;

		// **********************************************************************
		// {{{ Sample sort of array P by (i mod X, i div X)

		std::sort(P.begin(), P.end(), Pair::cmpIndexModDiv);		// sort locally

		DBG_ARRAY(debug_recursion, "Names locally sorted by cmpModDiv", P);

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

		    unsigned int divM = ptemp.index / namesGlobalMultipleOfX;
		    ptemp.index = DC[divM] + X * (ptemp.index - divM * namesGlobalMultipleOfX);

		    if (debug_recursion)
			std::cout << "splitter: " << ptemp.index << " = " << x << " - " << divM << "\n";

		    typename std::vector<Pair>::const_iterator it = std::lower_bound(P.begin(), P.end(), ptemp, Pair::cmpIndexModDiv);
		    splitterpos[i] = it - P.begin();
		}
		splitterpos[ nprocs ] = P.size();

		DBG_ARRAY2(debug_recursion, "Splitters positions", splitterpos, nprocs+1);

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

		MPI_Alltoallv( P.data(), sendcnt, sendoff, MPI_PAIR, recvBufPair.data(), recvcnt, recvoff, MPI_PAIR, MPI_COMM_WORLD );

		vector_free(P);

		// final X-1 tuples should be ignored due to recvoff areas
		merge_areas(recvBufPair, recvoff, nprocs, Pair::cmpIndexModDiv);

		// TODO: merge and reduce at once

		uint uniqueseq = 0;

		std::vector<uint> namearray ( recvBufPair.size() );
		for (unsigned int i = 0; i < recvBufPair.size(); ++i)
		{
		    if (i != 0) {
			if (recvBufPair[i-1].unique && recvBufPair[i].unique)
			    uniqueseq++;
		    }

		    namearray[i] = recvBufPair[i].name;
		}

		DBG_ARRAY(debug_recursion, "Pairs P (globally sorted by indexModDiv)", recvBufPair);

		std::cout << "uniques in sequence: " << uniqueseq << " - " << recvBufPair.size() / 2 << "\n";

		// }}} end Sample sort of array P

		if (uniqueseq > recvBufPair.size() / 2 && 0)
		{
		    // **********************************************************************
		    // * recurse on compressed sequence of duplicates and uniques

		    // reuse name array's second half for indexes
		    uint* indexarray = namearray.data() + recvBufPair.size() / 2;

		    uint j = 0;
		    for (unsigned int i = 0; i < recvBufPair.size(); ++i)
		    {
			if (i != 0) {
			    if (recvBufPair[i-1].unique && recvBufPair[i].unique)
				continue;
			}

			namearray[j] = recvBufPair[i].name;
			
			unsigned int divM = i / namesGlobalMultipleOfX;
			uint index = DC[divM] + X * (i - divM * namesGlobalMultipleOfX);

			indexarray[j] = index;
			std::cout << "dup/firstunique name: " << namearray[j] << " - index " << indexarray[j] << "\n";
			++j;
		    }

		    uint oldNamesGlobalSize = namesGlobalSize;
		    namesGlobalSize = j;
		    namesLocalStride = ( namesGlobalSize + nprocs - 1 ) / nprocs;		// rounded up
		    namesLocalStride += X - namesLocalStride % X;				// round up to nearest multiple of X

		    assert(j < recvBufPair.size() / 2);

		    vector_free(recvBufPair);

		    pDCX<DCParam, uint> rdcx;
		    //rdcx.dcx( namearray.data(), namesGlobalSize - (X-1), namesLocalStride, depth+1, oldNamesGlobalSize+1 );

		    std::cout << "SAlocal: " << rdcx.localSA.size() << " - indexes " << j << "\n";
		    std::cout << "SAlocal: " << rdcx.localSA.size() << " - indexes " << j << "\n";
		    std::cout << "SAlocal: " << rdcx.localSA.size() << " - indexes " << j << "\n";
		    std::cout << "SAlocal: " << rdcx.localSA.size() << " - indexes " << j << "\n";

		    assert(0);
		}
		else
		{
		    // recurse on full sequence of names

		    vector_free(recvBufPair);

		    DBG_ARRAY(debug_recursion, "namearray", namearray);

		    assert( namearray.size() == namesLocalStride || myproc == nprocs-1 );

		    pDCX<DCParam, uint> rdcx;
		    rdcx.dcx( namearray, depth+1, lastname+1 );

		    if (debug_rec_selfcheck)
		    {
			if (debug)
			    std::cout << "---------------------   RECURSION local checkSA ---------------- " << localSize << std::endl;

			std::vector<uint> nameAll;
			std::vector<uint> SAall;

			gather_vector(namearray, nameAll, X-1);
			gather_vector(rdcx.localSA, SAall, 0);

			DBG_ARRAY(debug_recursion, "nameAll", nameAll);
			DBG_ARRAY(debug_recursion, "SAall", SAall);

			if (myproc == ROOT)
			{
			    assert( sachecker::sa_checker(nameAll, SAall) );			    
			}

			MPI_Barrier( MPI_COMM_WORLD );
		    }

		    vector_free(namearray);

		    DBG_ARRAY(debug_recursion, "Recursive localSA", rdcx.localSA);

		    uint SAsize = rdcx.localSA.size();
		    uint allSAsize[nprocs+1];

		    MPI_Allgather( &SAsize, 1, MPI_UNSIGNED, allSAsize, 1, MPI_UNSIGNED, MPI_COMM_WORLD );

		    exclusive_prefixsum(allSAsize, nprocs);

		    DBG_ARRAY2(debug_recursion, "allSAsize", allSAsize, nprocs+1);

		    // generate array of pairs (index,rank) from localSA

		    P.resize( rdcx.localSA.size() );

		    for (unsigned int i = 0; i < rdcx.localSA.size(); ++i)
		    {
			// generate index in ModDiv sorted input sequence

			uint saidx = rdcx.localSA[i];

			unsigned int divM = saidx / namesGlobalMultipleOfX;

			uint index = DC[divM] + X * (saidx - divM * namesGlobalMultipleOfX);

			P[i].index = index;
			P[i].name = allSAsize[myproc] + i + 1;
		    }
		}
	    }
	    else // use sequential suffix sorter
	    {
		if (debug)
		    std::cout << "---------------------   RECURSION local sais ---------------- " << localSize << std::endl;

		std::vector<Pair> Pall;
		gather_vector(P, Pall, 0, MPI_PAIR);

		if (myproc == ROOT)
		{
		    assert( Pall.size() == (int)namesGlobalSize );

		    DBG_ARRAY(debug_recursion, "Global Names sorted index", Pall);

		    std::sort(Pall.begin(), Pall.end(), Pair::cmpIndexModDiv);		// sort locally

		    DBG_ARRAY(debug_recursion, "Global Names sorted cmpModDiv", Pall);

		    uint* namearray = new uint[ Pall.size() ];
		    for (unsigned int i = 0; i < Pall.size(); ++i)
			namearray[i] = Pall[i].name;

		    DBG_ARRAY2(debug_recursion, "Recursion input", namearray, Pall.size());

		    int* rSA = new int [ Pall.size() ];

		    yuta_sais_lite::saisxx< uint*, int*, int >( namearray, rSA, Pall.size(), lastname+1 );

		    delete [] namearray;

		    DBG_ARRAY2(debug_recursion, "Recursive SA", rSA, Pall.size());

		    // generate rank array - same as updating pair array with correct names
		    for (uint i = D; i < Pall.size(); ++i)
		    {
			Pall[ rSA[i] ].name = i + D;
		    }

		    DBG_ARRAY(debug_recursion, "Fixed Global Names sorted cmpModDiv", Pall);

		    std::swap(P, Pall);
		}
		else
		{
		    vector_free(P);
		}

	    } // end use sequential suffix sorter

	    if (debug_recursion)
		std::cout << "---------------------   END  RECURSION  ---------------  - mem = " << getmemusage() << std::endl;
	}
	else {
	    if (debug_recursion)
		std::cout << "---------------------   keine  Recursion---------------- - mem = " << getmemusage() << std::endl;
	}

	// in any outcome: here P contains pairs of index and unique rank.

	{
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
		ptemp.index = i * localStride;

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

	    vector_free(P);

	    merge_areas(recvBufPair, recvoff, nprocs);

	    DBG_ARRAY2(debug_recursion, "Pairs P (globally sorted by index)", recvBufPair.data(), recvBufPairSize);

	    // **********************************************************************
	    // *** every PE needs D additional sample suffix ranks to complete the final tuple

	    Pair temp[ D ];

	    MPI_Sendrecv( recvBufPair.data(), D, MPI_PAIR, ( myproc - 1 + nprocs ) % nprocs, MSGTAG,
			  &temp, D, MPI_PAIR, ( myproc + 1 ) % nprocs, MSGTAG,
			  MPI_COMM_WORLD, &status );

	    DBG_ARRAY2(debug_recursion, "Pairs P (extra tuples)", temp, D);

	    if ( myproc == nprocs - 1 )	// last processor gets sentinel tuples with lowest ranks
	    {
		for (unsigned int i = 0; i < D; ++i)
		{
		    // the first D ranks are freed up (above) for the following sentinel ranks 0..D-1:
		    recvBufPair[ recvBufPairSize + i ].name = D - i - 1;
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
		    S[k][i].chars[c] = (i*X + k + c < localSizeExtra) ? input[ i*X + k + c ] : 0;

		for (unsigned int d = 0; d < D; ++d)
		    S[k][i].ranks[d] = P[ dp + d ].name;

		if (DC[dp % D] == k) ++dp;
	    }
	}

	std::cout << "done creating S_i - mem = " << getmemusage() << "\n";

	// **********************************************************************
	// *** Sample sort tuple arrays
	{
	    for (unsigned int k = 0; k < X; ++k)
	    {
		// TODO: sort less - not all S[k] must be sorted to depth X-1 (needs additional lookup table)
		std::sort(S[k].begin(), S[k].end(), cmpTupleNdepth<X-1>);
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
	    TupleN* samplebufall = NULL;
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
		    merge_areas(S[k], recvoff[k], nprocs, cmpTupleNdepth<X-1>);

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

	    std::cout << "done merging suffixarray - mem = " << getmemusage() << "\n";

	    localSA.resize(j);
	}

	if (debug)
	{
	    std::cout << "******************** finished DCX (process " << myproc << ") depth " << depth << " ********************" << std::endl;
	}

	return true;
    }

    uint			globalSize;

    uint			localStride;

    std::vector<uint>   	localSA;

    std::vector<uint8_t>	localInput;

    bool run(const char* filename)
    {
	// **********************************************************************
	// * Read input file size

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
	    
	    char* endptr = NULL;
	    unsigned int reducesize = getenv("SIZE") ? strtoul(getenv("SIZE"),&endptr,10) : 0;
	    if (!endptr || *endptr != '\0') reducesize = 0;

	    if (reducesize && globalSize > reducesize)
		globalSize = reducesize;
	}

	MPI_Bcast( &globalSize, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD );

	// **********************************************************************
	// * Calculate local input size (= general stride)

	localStride = ( globalSize + nprocs - 1 ) / nprocs;	// divide by processors rounding up
	localStride += X - localStride % X;			// round up to nearest multiple of X

	assert( localStride * nprocs >= globalSize );

	if ( myproc == ROOT )
	{
	    std::cout << "Total input size = " << globalSize << " bytes. localStride = " << localStride << std::endl;
	}

	// **********************************************************************
	// * Read input file and send to other processors

	localInput.resize( localStride );

	assert( sizeof(alphabet_type) == 1 );

	if (myproc == ROOT)
	{
	    std::ifstream infile( filename, std::ios::binary );

	    // read input for processors 1 to n-1
	    for (int p = 1; p < nprocs; ++p)
	    {
		infile.seekg( p * localStride, std::ios::beg );

		uint readsize = (p != nprocs-1) ? localStride : globalSize - (p * localStride);

		std::cout << "Read for process " << p << " from pos " << p * localStride << " of length " << readsize << std::endl;

		infile.read( (char*)localInput.data(), readsize );
		MPI_Send( localInput.data(), readsize, MPI_CHAR, p, MSGTAG, MPI_COMM_WORLD );
	    }

	    if (!infile.good()) {
		perror("Error reading file.");
		return false;
	    }

	    // read input for processor 0 (ROOT)

	    std::cout << "Read for process 0 from pos 0 of length " << localStride << std::endl;

	    infile.seekg( 0, std::ios::beg );
	    infile.read( (char*)localInput.data(), localStride );

	    if (!infile.good()) {
		perror("Error reading file.");
		return false;
	    }
	}
	else // not ROOT: receive data
	{
	    MPI_Recv( localInput.data(), localStride, MPI_CHAR, ROOT, MSGTAG, MPI_COMM_WORLD, &status );

	    int recvcnt;
	    MPI_Get_count(&status, MPI_CHAR, &recvcnt);

	    uint readsize = (myproc != nprocs-1) ? localStride : globalSize - (myproc * localStride);

	    assert( (int)readsize == recvcnt );

	    localInput.resize(readsize);
	}

	MPI_Barrier( MPI_COMM_WORLD );

	// **********************************************************************
	// * Construct suffix array recursively

	dcx( localInput, 0, 256 );

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

    bool checkSA()
    {
	if (debug)
	{
	    std::cout << "******************** SAChecker (process " << myproc << ") ********************" << std::endl;
	    std::cout << "localStride = " << localStride << "\n";
	    std::cout << "localSA.size() = " << localSA.size() << "\n";
	    std::cout << "localInput.size() = " << localInput.size() << "\n";
	}

	assert( localStride + (X-1) == localInput.size() || myproc == nprocs-1 );

	// **********************************************************************
	// * Generate pairs (SA[i],i)

	uint localSAsize = localSA.size();

	MPI_Scan( MPI_IN_PLACE, &localSAsize, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );

	uint indexStart = localSAsize - localSA.size();
	
	std::vector<Pair> P ( localSA.size() );

	for (uint i = 0; i < localSA.size(); ++i)
	{
	    P[i].index = localSA[i];
	    P[i].name = indexStart + i;
	}

	DBG_ARRAY(debug_checker1, "(SA[i],i)", P);

	// **********************************************************************
	// * Sample sort of array P by (SA[i])
	{
	    std::sort(P.begin(), P.end());

	    uint* splitterpos = new uint[nprocs+1];
	    int* sendcnt = new int[nprocs];
	    int* sendoff = new int[nprocs+1];
	    int* recvcnt = new int[nprocs];
	    int* recvoff = new int[nprocs+1];

	    // use equidistance splitters from 0..globalSize (because indexes are known in advance)
	    splitterpos[0] = 0;
	    Pair ptemp;
	    ptemp.name = 0;
	    for ( int i = 1; i < nprocs; i++ ) {
		ptemp.index = i * localStride;

		typename std::vector<Pair>::const_iterator it = std::lower_bound(P.begin(), P.end(), ptemp);
		splitterpos[i] = it - P.begin();
	    }
	    splitterpos[ nprocs ] = P.size();

	    DBG_ARRAY2(debug_checker1, "Splitters positions", splitterpos, nprocs+1);

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

	    std::vector<Pair> recvBufPair ( recvoff[ nprocs ] + 1 );
	    unsigned int recvBufPairSize = recvoff[ nprocs ];

	    MPI_Alltoallv( P.data(), sendcnt, sendoff, MPI_PAIR, recvBufPair.data(), recvcnt, recvoff, MPI_PAIR, MPI_COMM_WORLD );

	    vector_free(P);

	    merge_areas(recvBufPair, recvoff, nprocs);

	    // **********************************************************************
	    // *** every P needs 1 additional pair

	    Pair temp;

	    MPI_Sendrecv( recvBufPair.data(), 1, MPI_PAIR, ( myproc - 1 + nprocs ) % nprocs, MSGTAG,
			  &temp, 1, MPI_PAIR, ( myproc + 1 ) % nprocs, MSGTAG,
			  MPI_COMM_WORLD, &status );

	    if ( myproc == nprocs - 1 )	// last processor gets sentinel pair: virtual pair of '$' position after string
	    {
		recvBufPair[ recvBufPairSize ].name = globalSize;
		recvBufPair[ recvBufPairSize ].index = INT_MAX;
	    }
	    else	// other processors get 1 following pair with indexes from the DC
	    {
		recvBufPair[ recvBufPairSize ] = temp;
	    }

	    std::swap(recvBufPair, P);
	}

	// now consider P as [ (i,ISA[i]) ]_{i=0..n-1} (by substituting i -> ISA[i])
	
	DBG_ARRAY(debug_checker1, "(SA[i],i) sorted by SA[i] equiv: (i,ISA[i]) including 1 extra pair", P);

	// **********************************************************************
	// * First check: is [P.name] the sequence [0..n)

	int error = false;

	for (uint i = 0; i < P.size()-1; ++i)		// -1 due to extra pair at end
	{
	    if (P[i].index != myproc * localStride + i)
	    {
		std::cout << "SA is not a permutation of [0,n) at position " << P[i].name << "\n";
		error = true;
		break;
	    }
	}

	MPI_Allreduce( MPI_IN_PLACE, &error, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD );

	if (error) return false;

	// **********************************************************************
	// * Generate triples (ISA[i], ISA[i+1], S[i])

	std::vector<Triple> S ( P.size()-1 );	// -1 due to extra pair at end

	for (uint i = 0; i < P.size()-1; ++i)
	{
	    S[i].rank1 = P[i].name;
	    S[i].rank2 = P[i+1].name;

	    S[i].char1 = localInput[i];
	}

	DBG_ARRAY(debug_checker2, "(ISA[i], ISA[i+1], S[i])", S);

	// **********************************************************************
	// * Sample sort of array S by (S[].rank1)
	{
	    std::sort(S.begin(), S.end());

	    uint* splitterpos = new uint[nprocs+1];
	    int* sendcnt = new int[nprocs];
	    int* sendoff = new int[nprocs+1];
	    int* recvcnt = new int[nprocs];
	    int* recvoff = new int[nprocs+1];

	    // use equidistance splitters from 0..globalSize (because indexes are known in advance)
	    splitterpos[0] = 0;
	    Triple ptemp;
	    for ( int i = 1; i < nprocs; i++ ) {
		ptemp.rank1 = i * localStride;

		typename std::vector<Triple>::const_iterator it = std::lower_bound(S.begin(), S.end(), ptemp);
		splitterpos[i] = it - S.begin();
	    }
	    splitterpos[ nprocs ] = S.size();

	    DBG_ARRAY2(debug_checker2, "Splitters positions", splitterpos, nprocs+1);

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

	    std::vector<Triple> recvBuf ( recvoff[ nprocs ] + 1 );
	    unsigned int recvBufSize = recvoff[ nprocs ];

	    MPI_Alltoallv( S.data(), sendcnt, sendoff, MPI_TRIPLE, recvBuf.data(), recvcnt, recvoff, MPI_TRIPLE, MPI_COMM_WORLD );

	    vector_free(S);

	    merge_areas(recvBuf, recvoff, nprocs);

	    // **********************************************************************
	    // *** every P needs 1 additional triple

	    Triple temp;

	    MPI_Sendrecv( recvBuf.data(), 1, MPI_TRIPLE, ( myproc - 1 + nprocs ) % nprocs, MSGTAG,
			  &temp, 1, MPI_TRIPLE, ( myproc + 1 ) % nprocs, MSGTAG,
			  MPI_COMM_WORLD, &status );

	    if ( myproc == nprocs - 1 )	// last processor gets sentinel triple - which shouldnt be compared later on.
	    {
		recvBuf[ recvBufSize ].rank1 = INT_MAX;
		recvBuf[ recvBufSize ].rank2 = INT_MAX;
		recvBuf[ recvBufSize ].char1 = 0;
	    }
	    else	// other processors get 1 following pair with indexes from the DC
	    {
		recvBuf[ recvBufSize ] = temp;
	    }

	    std::swap(recvBuf, S);
	}
	
	DBG_ARRAY(debug_checker2, "(ISA[i], ISA[i+1], S[i]) sorted by ISA[i]\nequiv: (ISA[SA[i]], ISA[SA[i]+1], S[SA[i]]) sorted by i", S);

	// now consider S as [ (i, ISA[SA[i]+1], S[SA[i]]) ]_{i=0..n-1} (by substituting i -> SA[i])

	// **********************************************************************
	// * Second check: use ISA to check suffix of suffixes for correct order

	unsigned int iend = S.size()-1 - (myproc == nprocs-1 ? 1 : 0);

	for (uint i = 0; !error && i < iend; ++i)		// -1 due to extra pair at end
	{
	    if (S[i].char1 > S[i+1].char1) {
		// simple check of first character of suffix
		std::cout << "Error: suffix array position " << i + myproc * localStride  << " ordered incorrectly.\n";
		error = true;
	    }
	    else if (S[i].char1 == S[i+1].char1)
	    {
		if ( S[i+1].rank2 == globalSize ) {
		    // last suffix of string must be first among those
		    // with same first character
		    std::cout << "Error: suffix array position " << i + myproc * localStride << " ordered incorrectly.\n";
		    error = true;
		}
		if ( S[i].rank2 != globalSize &&
		     S[i].rank2 > S[i+1].rank2 )
		{
		    // positions SA[i] and SA[i-1] has same first
		    // character but their suffixes are ordered
		    // incorrectly: the suffix position of SA[i] is given
		    // by ISA[SA[i]]
		    std::cout << "Error: suffix array position " << i + myproc * localStride << " ordered incorrectly.\n";
		    error = true;
		}
	    }
	}

	MPI_Allreduce( MPI_IN_PLACE, &error, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD );

	return (error == false);
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

	    char* endptr = NULL;
	    unsigned int reducesize = getenv("SIZE") ? strtoul(getenv("SIZE"),&endptr,10) : 0;
	    if (!endptr || *endptr != '\0') reducesize = 0;

	    if (reducesize && globalSize > reducesize)
		globalSize = reducesize;
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
	    //DBG_ARRAY(debug_output, "Suffixarray collected", gSA);

	    std::ifstream infile( filename );
	    if (!infile.good()) {
		perror("Cannot read input file");
		return -1;
	    }

	    std::vector<uint8_t> string(globalSize);

	    infile.read((char*)string.data(), globalSize);
	    if (!infile.good()) {
		perror("Cannot read input file");
		return -1;
	    }

	    if (debug_output)
	    {
		std::cout << "result suffix array: \n";

		for (unsigned int i = 0; i < gSA.size(); ++i)
		{
		    std::cout << i << " : " << gSA[i] << " : ";

		    for (unsigned int j = 0; gSA[i]+j < globalSize && j < 32; ++j)
		    {
			std::cout << strC(string[gSA[i]+j]) << " ";
		    }

		    std::cout << "\n";
		}
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

	if ( argc >= 3 ) {
	    dcx.writeSA(argv[2]);
	}

	MPI_Barrier( MPI_COMM_WORLD );

	std::cout << "Suffix array checker: " << dcx.checkSA() << "\n";

	dcx.checkSAlocal(argv[1]);
    }

    MPI_Finalize();
    return 0;
}
