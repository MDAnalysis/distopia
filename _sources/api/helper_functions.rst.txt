Helper Functions
****************

.. doxygenfunction:: genericload(const VectorToScalarT<VectorT> *source)
.. doxygenfunction:: genericload(const T *source)
.. doxygenfunction:: genericidxload(const T *source, const std::size_t *idxs, T &x, T &y, T &z)
.. doxygenfunction:: genericidxload(const VectorToScalarT<VectorT> *source, const std::size_t *idxs, VectorT &x, VectorT &y, VectorT &z)
.. doxygenfunction:: genericstore(VectorToScalarT<VectorT> *target, const VectorT val)
.. doxygenfunction:: genericstore(T *target, const T val)
.. doxygenfunction:: genericstream(VectorToScalarT<VectorT> *target, const VectorT val)
.. doxygenfunction:: genericstream(T *target, const T val)