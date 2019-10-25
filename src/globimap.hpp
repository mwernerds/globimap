/*
Global Binary Map Core Implementation

Some remarks on the hashing trick performance:

Interesting: the assembly for the hashing trick:
32 bit, gcc, fast, native leads to
	mov     eax, DWORD PTR i[rip]
        inc     eax
        imul    eax, DWORD PTR h2[rip]
        add     eax, DWORD PTR h1[rip]
        and     eax, DWORD PTR mask[rip]


        mov     rax, QWORD PTR i[rip]
        inc     rax
        imul    rax, QWORD PTR h2[rip]
        add     rax, QWORD PTR h1[rip]
        and     rax, QWORD PTR mask[rip]
        ret

imul has 3 cycle latency


Functions:
- hash (either murmur or DJB64, murmur to be preferred as it well-documented (see paper) to work with hashing trick and it would  allow for prefixing)


class Globimap:
    void clear()
	delete the image and the correction information. release memory

    void setNamespace()
	set the namespace for the murmur hash multi-namespace support

    void add_error(std::vector<uint32_t> a)   
	register an error at (a[0], a[1]) in the error correction engine
    void put(std::vector<uint32_t> a)
	set the pixel (a[0],a[1])
    bool get(std::vector<uint32_t> a)
	get the pixel (a[0],a[1])
    void configure (size_t _d, size_t logm)
	configure the filter with _d hash functions and 2^logm bit)
    void summary()
	give a summary (compute-intensive) of the data structure
    std::vector<double> &rasterize(uint32_t x, uint32_t y, uint32_t s0, uint32_t s1)
	rasterize the rectangle (x,y) -> (x+s0, y+s1) with stride of s1. OMP loop parallel

    std::vector<double> & apply_correction(uint32_t x, uint32_t y, uint32_t s0, uint32_t s1)
	apply correction information to suppress false-positives (OMP task parallel)

    void tobuffer(std::string &buf)
	serialize the buffer into a string for writing/storing/communicating

    void frombuffer(std::string &buf, size_t n)
	deserialize the buffer from a string, n is filter size
    
*/

#ifndef GLOBIMAP_HPP_INC
#define GLOBIMAP_HPP_INC
#include<limits>
#include<vector>
#include<list>
#include<unordered_set>
#include<set>
#include<sstream>

#define GLOBIMAP_USE_MURMUR

//#define GLOBIMAP_USE_MURMUR_PREFIX
//#define GLOBIMAP_USE_DJB64


#include <string.h>

#ifdef GLOBIMAP_USE_MURMUR
#include "murmur.hpp"

#ifdef GLOBIMAP_USE_MURMUR_PREFIX
inline void hash(uint64_t *data, size_t len, uint64_t *v1, uint64_t *v2, const char *prefix="")
{
    uint64_t hash[2];
    char buffer[1024];
    snprintf(buffer, 1024,"%s-%lu-%lu", prefix,data[0],data[1]);
    trajcomp::murmur::MurmurHash3_x86_128 ( data, strlen(buffer), *v1, (void*) hash ); // len * sizeof(uint64_t)
    *v1 = hash[0];
    *v2 = hash[1];
}
#else
inline void hash(const uint64_t *data, const size_t len, uint64_t *v1, uint64_t *v2)
{
    uint64_t hash[2];
/*    char buffer[1024];
    snprintf(buffer, 1024,"%s-%lu-%lu", prefix,data[0],data[1]);*/
    trajcomp::murmur::MurmurHash3_x86_128 ( data, 16, *v1, (void*) hash ); // len*sizeof(uint64_t) 
    *v1 = hash[0];
    *v2 = hash[1];
}


#endif
#endif


#ifdef GLOBIMAP_USE_DJB64

// This should not be used, it would need careful evaluation on many layers. Stays here for completeness!!!

inline void hash(uint32_t *data, size_t len, uint32_t *v1, uint32_t *v2)
{
    const char *prefix1 = "p1";
    const char *prefix2 = "x8";
    char *str = reinterpret_cast<char *> (data);
    uint32_t hash1 = 5381;
    uint32_t hash2 = 5381;
    for (size_t i=0; i < 2; i++)
    {
        hash1 = ((hash1 << 5) + hash1) + prefix1[i]; /* hash * 33 + c */
        hash2 = ((hash2 << 5) + hash2) + prefix2[i]; /* hash * 33 + c */
    }
    
    for (size_t i=0; i < len * sizeof(uint32_t); i++)
    {
        hash1 = ((hash1 << 5) + hash1) + *str; /* hash * 33 + c */
        hash2 = ((hash2 << 5) + hash2) + *str; /* hash * 33 + c */
	str++;
    }
//    *v1 = static_cast<uint32_t>((hash & 0xFFFFFFFF00000000LL) >> 32);
//    *v2 = static_cast<uint32_t>(hash & 0xFFFFFFFFLL) ;
   *v1 = hash1;
   *v2 = hash2;
}
#endif


#ifdef GLOBIMAP_USE_LOOKUP3
// This should not be used, it would need careful evaluation on many layers. Not published, google for lookup3.h ;-)
#include "lookup3.hpp"

inline void hash(uint32_t *data, size_t len, uint32_t *v1, uint32_t *v2)
{
    hashword2 (data, len ,v1, v2);
 }
#endif





template<typename element_type=bool>
class GloBiMap {
    public:
	uint64_t maxhash=0; ///< was used for debugging that the hash numbers actually are large enough
	typedef std::set<std::pair<uint32_t, uint32_t>> error_container_t; // could be unordered_set dep. on your situation.
    private:
	std::vector<element_type> filter;
	int d;
	uint64_t mask;

protected:
       std::vector<double> storage;
       error_container_t errors; 

    public:
void clear()
{
    filter.clear();
    errors.clear();
}

void add_error(std::vector<uint32_t> a)
{
//    std::cout <<"Adding error information for " << a[0]<< "/" << a[1] << std::endl; 
    errors.emplace(std::make_pair(a[0],a[1]));
}



void put(std::vector<uint64_t> a)
{
//    std::cout << "put" << a[0] << ";" << a[1] <<";";

// get the two hashs:
    uint64_t h1=8589845122,h2=8465418721;
    double maxp = 0;
    hash (&a[0], 2 ,&h1, &h2);
    for (size_t i=0; i < static_cast<size_t>(d); i++)
    {
       uint64_t k = (h1 + (i+1)*h2)  & mask;

#ifdef GLOBIMAP_COMPUTE_MAXHASH
       if (k > maxhash)
         maxhash = k;
#endif
//       std::cout << "put: k=" << k << std::endl;
#ifdef DEBUG_HASH_PUT
	if (static_cast<long double>(k) / mask > maxp)
	   maxp = static_cast<long double>(k) / mask;
	#pragma omp critical
	std::cout << k << "(;" << static_cast<long double>(k) / mask <<") => " << maxp << std::endl;
#endif
	#pragma omp critical
	filter[k] = 1;
    }
#ifdef DEBUG_HASH_PUT
    std::cout << std::endl;
#endif

}
bool get(std::vector<uint64_t> a)
{
//    std::cout << "GET for " << a[0] << "/" << a[1] << std::endl;
    uint64_t h1=8589845122,h2=8465418721;
    hash (&a[0], 2 ,&h1, &h2);
    for (size_t i=0; i < static_cast<size_t>(d); i++)
    {
       uint64_t k = (h1 + (i+1)*h2)  & mask;
       if (filter[k] == 0)
          return false;
    }
    return true;
}

    
    void configure (size_t _d, size_t logm)
    {
        d = _d;
	mask = (static_cast<uint64_t>(1) << logm)-1;
	std::cout << "logm:" << logm << "mask=" << mask << std::hex << "0x" << mask << std::dec  << std::endl;
	filter.resize(mask +1);
	std::cout << "filter.size=" << filter.size() << std::endl;
    }
    


    std::string summary()
    {
       size_t ones = 0;
       #pragma omp parallel for
       for (size_t i=0; i < filter.size(); i++)
       {
	  if (filter[i] == 1)
          #pragma omp atomic
	      ones++;
       }
       std::stringstream ss;
//       ss << std::hex << t1 << std::endl << t2 << std::endl << t3 << std::endl << std::dec;
       ss << "{" << std::endl;
#ifdef GLOBIMAP_COMPUTE_MAXHASH
       ss << "\"maxhash\":" << maxhash << "," << std::endl;
#endif
       ss << "\"storage:\": " << static_cast<double>(filter.size()) / 8 / 1024 / 1024 << ","<< std::endl;
       ss << "\"ones:\": " << ones << ","<< std::endl;
       ss << "\"foz:\": " << static_cast<double>((filter.size() - ones)) / (double) filter.size() << ","<< std::endl;
       ss << "\"eci\": " <<        errors.size() << std::endl;
       ss << "}" << std::endl;
       return ss.str();
  }
    

    std::vector<double> &rasterize(uint64_t x, uint64_t y, uint32_t s0, uint32_t s1)
    {
	storage.resize(s0*s1);
	#pragma omp parallel for
	for (uint32_t i = 0; i < s0; i++)
	for (uint32_t j = 0; j < s1; j++)
	{
	   storage[i*s1+j] = get({x+i,y+j});
	}
	return storage;
    }

    std::vector<double> & apply_correction(uint32_t x, uint32_t y, uint32_t s0, uint32_t s1)
    {
	// apply corrections over storage
	if (storage.size() != s0 * s1)
	  throw(std::runtime_error("corrections can only be applied after rasterize with same extends (parameters!)"));
	/*
	task-based parallelism: the loop through the list is on a single scope, yet the processing is parallel
	*/
	#pragma omp parallel
	#pragma omp single
	{
	  for(auto it = errors.begin(); it != errors.end(); ++it)
	     #pragma omp task firstprivate(it)
	     {
	       auto first = (it->first-x);
	       auto second = (it->second-y);
	       if (first >= 0 && second >= 0 && first < s0 && second < s1)
		  storage[first*s1+second] = 0;
	     }
	  #pragma omp taskwait
	}
	return storage;
    }

    void tobuffer(std::string &buf)
    {
	buf.clear();
	char ch=0;
	for (size_t i=0; i<filter.size(); i++)
	{
	    auto bit = i %8;
	    ch |= (filter[i] << bit);
	    if (bit == 7){
		// full byte has been built.
		buf += ch;
		ch = 0;
	    }
	}
	if (filter.size() % 8 != 0)  // we have collected beyond the end.
	  buf += ch;
    }
    
    void frombuffer(std::string &buf, size_t n)
    {
	filter.resize(n);
	size_t k=0;
	for (size_t i=0; i < buf.size(); i++)
	{
	    char ch = buf[i];
	    for(size_t j=0; j < 8; j++)
	      if (k < n)
	        filter[k++] = ch & (1<<j);
	
	}
    }
    void _frombuffer(std::string &buf )
    {
	size_t k=0;
	auto n = filter.size();
	for (size_t i=0; i < buf.size(); i++)
	{
	    char ch = buf[i];
	    for(size_t j=0; j < 8; j++)
	      if (k < n)
	        filter[k++] = ch & (1<<j);
	
	}
    }

};


#endif
