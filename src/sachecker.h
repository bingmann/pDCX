/*
 * Algorithm to check a suffix array. Loosely based on the ideas of
 * Kaerkkaeinen.
 */

namespace sachecker {

template <typename alphabet_type>
bool string_leq(const std::vector<alphabet_type>& string, unsigned int p1, unsigned int p2)
{
    assert(p1 != p2 && "Bad string comparison");

    for (unsigned int i = 0; ; ++i)
    {
        if (p1+i >= string.size() && p2+i < string.size()) return true;
	if (p2+i < string.size() && p2+i >= string.size()) return false;

	if (p1+i >= string.size() && p2+i >= string.size()) {
	    assert(0 && "Bad string comparison");
	    return true;
	}

	if (string[p1+i] < string[p2+i]) return true;
	if (string[p1+i] > string[p2+i]) return false;
    }

    assert(0 && "Bad string comparison");
} 

template <typename alphabet_type, typename offset_type>
bool sa_checker_full(const std::vector<alphabet_type>& string, const std::vector<offset_type>& SA)
{
    if (string.size() != SA.size())
	return false;

    // *** check that SA contains a permutation of [0,n)
    // done via polynomial evaluation of (z-0)*(z-1)*(z-2)...(z-(n-1))
    uint64_t p = 4294967291;	// last prime < 2^32
    uint64_t z = SA.size() + 2;
    assert(p > z);

    uint64_t val1 = 1, val2 = 1;

    for (unsigned int i = 0; i < SA.size(); ++i)
    {
	val1 *= (z - SA[i]);
	val2 *= (z - i);

	val1 %= p; val2 %= p;
    }

    if (val1 != val2)
	return false;
    
    // *** check naively that SA is a sorted suffix array for string

    for (unsigned int i = 0;  i < SA.size()-1;  i++)
    {
	if (!string_leq(string, SA[i], SA[i+1])) {
	    std::cout << "Error: suffix " << i << " ordered incorrectly.\n";
	    return false;
	}
    }

    return true;
}

template <typename alphabet_type, typename offset_type>
bool sa_checker(const std::vector<alphabet_type>& string, const std::vector<offset_type>& SA)
{
    if (string.size() != SA.size()) {
	std::cout << "String length " << string.size() << " != SA length " << SA.size() << "\n";
	return false;
    }

    std::vector<offset_type> ISA(SA.size());

    // *** check that SA contains a permutation of [0,n)
    // done via polynomial evaluation of (z-0)*(z-1)*(z-2)...(z-(n-1))
    // and build ISA simultaneously
    uint64_t p = 4294967291;	// greatest prime < 2^32
    uint64_t z = SA.size() + 2;
    assert(p > z);

    uint64_t val1 = 1, val2 = 1;

    for (unsigned int i = 0; i < SA.size(); ++i)
    {
	val1 *= (z - SA[i]);
	val2 *= (z - i);

	val1 %= p; val2 %= p;
	
	if ((int)SA[i] < 0 || SA[i] >= (offset_type)SA.size()) {
	    std::cout << "Error: suffix array contents exceeds n-1.\n";
	    return false;
	}

	ISA[ SA[i] ] = i;
    }

    if (val1 != val2) {
	std::cout << "Error: suffix array is not a permutation of 0..n-1.\n";
	return false;
    }
    
    // *** check naively that SA is a sorted suffix array for string

    for (unsigned int i = 1;  i < SA.size();  i++)
    {
	if (string[ SA[i-1] ] > string[ SA[i] ]) {
	    // simple check of first character
	    std::cout << "Error: suffix " << i << " ordered incorrectly.\n";
	    return false;
	}
	else if (string[ SA[i-1] ] == string[ SA[i] ])
	{
	    if ( SA[i] + 1 == (offset_type)string.size() ) {
		 // last suffix of string must be first among those
		 // with same first character
		std::cout << "Error: suffix " << i << " ordered incorrectly.\n";
		return false;
	    }
	    if ( SA[i-1] + 1 != (offset_type)string.size() &&
		 ISA[ SA[i]+1 ] < ISA[ SA[i-1]+1 ] )
	    {
		// positions SA[i] and SA[i-1] has same first
		// character but their suffixes are ordered
		// incorrectly: the suffix position of SA[i] is given
		// by ISA[SA[i]]
		std::cout << "Error: suffix " << i << " ordered incorrectly.\n";
		return false;
	    }
	}
    }

    return true;
}

} // namespace sachecker
