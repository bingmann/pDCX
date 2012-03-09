#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <assert.h>
#include <inttypes.h>

bool isdiffcover(int n, const std::vector<int>& DC)
{
    unsigned int x = 0;
    
    for (unsigned int i = 0; i < DC.size(); ++i)
    {
	for (unsigned int j = 0; j < DC.size(); ++j)
	{
	    x |= 1 << ((DC[i] - DC[j] + n) % n);
	}
    }

    return (x + 1 == (1 << n));
}

void find_diffcover(int n)
{
    // enumerate all 2^n vectors
    // but skip those equivalent via cyclic shifts by fixing the first bit = 1

    // second version: fix the second bit = 1, as 0 is always included in a diffcover.

    unsigned int minsize = n;

    unsigned int mindelta = n*n;

    std::vector< std::vector<int> > minlist;

    for (uint64_t x = 0; x < (1 << n); x += 1)
    {
	// generate DC from x regarded as a bitfield

	std::vector<int> DC;

	for (unsigned int i = 0; i < n; ++i)
	{
	    if (x & (1LLU << i)) {
		DC.push_back(i);
	    }
	}

	if (DC.size() > minsize) continue;

#if 0
	std::cout << "DC" << n << " : ";
	for (unsigned int i = 0; i < DC.size(); ++i)
	{
	    std::cout << DC[i] << " ";
	}
	std::cout << "\n";
#endif

	if (!isdiffcover(n,DC)) continue;

	if (DC.size() < minsize)
	{
	    minsize = DC.size();
	    
	    mindelta = n*n;
	    minlist.clear();
	}

#if 0
	std::cout << "DC" << n << " : ";
	for (unsigned int i = 0; i < DC.size(); ++i)
	{
	    std::cout << DC[i] << " ";
	}
	std::cout << "\n";
#endif
	// calculate DX complement

	std::vector<int> DX (n - DC.size(), -1);

	{
	    unsigned int e = 0, j = 0;
	    for (unsigned int d = 0; d < DC.size(); ++d)
	    {
		while (e++ < DC[d])
		{
		    DX.at(j++) = e-1;
		}
	    }
	    while (e++ < n) {
		DX.at(j++) = e-1;
	    }
	    assert(j == n - DC.size());
#if 0
	    std::cout << "DX (complement): ";
	    for (unsigned int dc = 0; dc < n - DC.size(); ++dc) {
		std::cout << DX[dc] << " ";
	    }
	    std::cout << "\n";
#endif
	}
	
	const unsigned int X = n;
	const unsigned int D = DC.size();

	// calculate character tuple length in comparison of B_k in complement.
	std::vector<int> DXL (X - D, -1);
	
	for (unsigned int i = 0; i < X - D; ++i)
	{
	    unsigned int p = DX[i];

	    unsigned int delta = X;

	    // find smallest increment so that p+delta hits element in D
	    for (unsigned int j = 0; j < D; ++j)
	    {
		delta = std::min(delta, (DC[j] - p + X) % X);
	    }

	    DXL[i] = delta;
	}

	unsigned int dxlsum = 0;
#if 1
	for (unsigned int dc = 0; dc < X-D; ++dc) {
	    dxlsum += DXL[dc];
	}
#else
	std::cout << "DXL (delta): ";
	for (unsigned int dc = 0; dc < X-D; ++dc) {
	    std::cout << DXL[dc] << " ";
	    dxlsum += DXL[dc];
	}
	std::cout << " sum = " << dxlsum << "\n";
#endif
	if (dxlsum < mindelta)
	{
	    mindelta = dxlsum;
	    minlist.clear();
	}
	if (dxlsum == mindelta)
	{
	    minlist.push_back( DC );
	}
    }

    std::sort(minlist.begin(),minlist.end());

    std::ostringstream oss;
    oss << "diffcovers-" << n << ".txt";
    //std::ofstream of(oss.str().c_str());
    std::ostream& of = std::cout;

    for (unsigned int l = 0; l < minlist.size(); ++l)
    {
	const std::vector<int> &DC = minlist[l];

	of << "DC" << n << " : ";
	for (unsigned int i = 0; i < DC.size(); ++i)
	{
	    of << DC[i] << " ";
	}
	of << "\n";

	// calculate DX complement

	std::vector<int> DX (n - DC.size(), -1);

	{
	    unsigned int e = 0, j = 0;
	    for (unsigned int d = 0; d < DC.size(); ++d)
	    {
		while (e++ < DC[d])
		{
		    DX.at(j++) = e-1;
		}
	    }
	    while (e++ < n) {
		DX.at(j++) = e-1;
	    }
	    assert(j == n - DC.size());

	    of << "DX (complement): ";
	    for (unsigned int dc = 0; dc < n - DC.size(); ++dc) {
		of << DX[dc] << " ";
	    }
	    of << "\n";
	}
	
	const unsigned int X = n;
	const unsigned int D = DC.size();

	// calculate character tuple length in comparison of B_k in complement.
	std::vector<int> DXL (X - D, -1);
	
	for (unsigned int i = 0; i < X - D; ++i)
	{
	    unsigned int p = DX[i];

	    unsigned int delta = X;

	    // find smallest increment so that p+delta hits element in D
	    for (unsigned int j = 0; j < D; ++j)
	    {
		delta = std::min(delta, (DC[j] - p + X) % X);
	    }

	    DXL[i] = delta;
	}

	unsigned int dxlsum = 0;

	of << "DXL (delta): ";
	for (unsigned int dc = 0; dc < X-D; ++dc) {
	    of << DXL[dc] << " ";
	    dxlsum += DXL[dc];
	}
	of << " sum = " << dxlsum << "\n";

	std::vector<bool> inDC (X, false);
	
	for (unsigned int i = 0; i < DC.size(); ++i)
	    inDC[DC[i]] = true;

	of << "const int DC" << X << "Param::cmpDepthRanks[" << X << "][" << X << "][3] =\n"
	   << "{\n";
	for (unsigned int a = 0; a < X; ++a)
	{
	    of << "    {";
	    for (unsigned int b = 0; b < X; ++b)
	    {
		if (b != 0) of << ",";

		int r1 = 0, r2 = 0;

		unsigned int j = 0;
		while (j < X && (!inDC[(a+j) % X] || !inDC[(b+j) % X]))
		{
		    if (inDC[(a+j) % X]) ++r1;
		    if (inDC[(b+j) % X]) ++r2;

		    ++j;
		}
		       
		std::cout << " {" << std::setw(2) << j << "," << r1 << "," << r2 << " }" ;
	    }
	    of << " },\n";
	}
	of << "};\n";
    }

    std::cout.flush();
}

int main()
{
//#pragma omp parallel for schedule(dynamic,1)
    for (int i = 3; i <= 63; ++i)
    {
	find_diffcover(i);
    }
}
