/*
This is the Python Interface based on Boost::Python and Boost::Numpy


*/
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include<iostream>

#include<boost/python.hpp>
#include<boost/python/numpy.hpp>
#include<boost/python/numpy/ndarray.hpp>
#include <numpy/ndarrayobject.h> 
#include <numpy/ndarraytypes.h> 

// This holds the actual implementation. Copy this header to your projects (and a hasher, for example murmur.hpp)
#include "globimap.hpp"


namespace bp = boost::python;
namespace np = boost::python::numpy;

// Wrap a 1D C++ array (given as a pointer) to a numpy object.
static boost::python::object wrap(double* data, npy_intp size) {
  using namespace boost::python;

  npy_intp shape[1] = { size }; // array size
  PyObject* obj = PyArray_New(&PyArray_Type, 1, shape, NPY_DOUBLE, // data type
                              NULL, data, 
                              0, NPY_ARRAY_CARRAY, // NPY_ARRAY_CARRAY_RO for readonly
                              NULL);
  handle<> array( obj );
  return object(array);
}

// Wrap 2D C++ array (given as pointer) to a numpy object.
static boost::python::object wrap2D(double* data, npy_intp h, npy_intp w) {
  using namespace boost::python;

  npy_intp shape[2] = { h,w }; // array size
  PyObject* obj = PyArray_New(&PyArray_Type, 2, shape, NPY_DOUBLE, // data type
                              NULL, data, // data pointer
                              0, NPY_ARRAY_CARRAY, // NPY_ARRAY_CARRAY_RO for readonly
                              NULL);
  handle<> array( obj );
  return object(array);
}


// This template is a handy tool to call a function f(i,j,value) for each entry of a 2D matrix self.
template<typename func>
static void map_matrix(np::ndarray &self, func f)
{
   auto nd = self.get_nd();
   if (nd != 2) throw(std::runtime_error("2D array expected"));
   auto s1 = self.strides(0);
   auto s2 = self.strides(1);
   auto data = self.get_data();

   for (int i1=0; i1 < self.shape(0); i1++)
   {
     for(int i2=0; i2 < self.shape(1); i2++)
     {
         auto offset = i1 * s1 + i2 * s2;
	 double *d = reinterpret_cast<double *> (data + offset);
	 f(i1,i2,*d);
     }
   }
   

}

// This wil be our implementation in C++ of a Python class globimap.
typedef GloBiMap<bool> globimap_t ;

// The module begins
BOOST_PYTHON_MODULE(globimap) {
    import_array();
    np::initialize();
    using namespace boost::python;


    // It exports a class (named globimap) with chained functions, see README.md
   class_<globimap_t>("globimap")      
       .def("rasterize", +[](globimap_t& self, size_t x, size_t y, size_t s0, size_t s1) -> object {
	  auto &data = self.rasterize(x,y,s0,s1);
          return wrap2D(&data[0],s0,s1);
	})
       .def("correct", +[](globimap_t& self, size_t x, size_t y, size_t s0, size_t s1) -> object {
	  //auto &data = self.rasterize(x,y,s0,s1);
	  auto &data = self.apply_correction(x,y,s0,s1);	  
          return wrap2D(&data[0],s0,s1);
	})
       .def("put", +[](globimap_t &self, uint32_t x, uint32_t y) { self.put({x,y});})
       .def("get", +[](globimap_t &self, uint32_t x, uint32_t y) -> bool { return self.get({x,y});})
       .def("configure", +[](globimap_t &self, size_t k, size_t m)  {  self.configure(k,m);})
       .def("clear", +[](globimap_t &self){self.clear();})
       .def("summary", +[](globimap_t &self)->std::string { return self.summary();})
       .def("map", +[](globimap_t &self,object mat, int o0, int o1) {
		np::ndarray a = np::from_object(mat, np::dtype::get_builtin<double>());
		map_matrix(a, [&](int i0, int i1, double v){
		if (v != 0 && v!= 1)
		 throw(std::runtime_error("data is not binary."));
		if (v == 1)
   		    self.put({static_cast<uint32_t>(o0+i0),static_cast<uint32_t>(o1+i1)});
		});
	})
	.def("enforce",+[](globimap_t &self, object mat, int o0, int o1){
		np::ndarray a = np::from_object(mat, np::dtype::get_builtin<double>());
		map_matrix(a, [&](int i0, int i1, double v){
		if (v==0 && self.get({static_cast<uint32_t>(o0+i0),static_cast<uint32_t>(o1+i1)}))
		{// this is a false positive
		    self.add_error({static_cast<uint32_t>(o0+i0),static_cast<uint32_t>(o1+i1)});
		}
		});
	});
}
