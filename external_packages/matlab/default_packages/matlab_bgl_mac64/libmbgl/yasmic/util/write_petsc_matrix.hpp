#ifndef YASMIC_UTIL_WRITE_PETSC_MATRIX
#define YASMIC_UTIL_WRITE_PETSC_MATRIX

namespace impl
{
    namespace endian
    {
        void swap_int_4(int *tni4)                  /* 4 byte signed integers   */
        {
          *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
                 ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));  
        }

        void swap_double_8(double *tndd8)          /* 8 byte double numbers          */
        {
          char *tnd8=(char *)tndd8;
          char tnc;

          tnc=*tnd8;
          *tnd8=*(tnd8+7);
          *(tnd8+7)=tnc;

          tnc=*(tnd8+1);
          *(tnd8+1)=*(tnd8+6);
          *(tnd8+6)=tnc;

          tnc=*(tnd8+2);
          *(tnd8+2)=*(tnd8+5);
          *(tnd8+5)=tnc;

          tnc=*(tnd8+3);
          *(tnd8+3)=*(tnd8+4);
          *(tnd8+4)=tnc;
        }
    }
}

/**
 * Warning: This method modifies the data in rows, cols,
 * and vals.  It tries to restore it afterwards, but no 
 * guarantees.
 */
template <class VecRows, class VecCols, class VecVals>
bool write_petsc_matrix(std::ostream& f,
                 VecRows& rows, VecCols& cols, VecVals& vals,
                 int nr, int nc, std::size_t nz)
{	
    // convert to differences
	adjacent_difference(++rows.begin(), rows.end(), rows.begin());

    // convert to big endian...
	{
		int* intptr = &rows[0];
		unsigned int maxi = rows.size();
		for (unsigned int i = 0; i < maxi; ++i)
		{
            impl::endian::swap_int_4(intptr);
			++intptr;
		}
	}
	
	if (cols.size() > 0)
	{
		int* intptr = &cols[0];
		unsigned int maxi = cols.size();
		for (unsigned int i = 0; i < maxi; ++i)
		{
			impl::endian::swap_int_4(intptr);
			++intptr;
		}
	}

	if (vals.size() > 0)
	{
		double* doubleptr = &vals[0];
		unsigned int maxi = vals.size();
		for (unsigned int i = 0; i < maxi; ++i)
		{
			impl::endian::swap_double_8(doubleptr);
			++doubleptr;
		}
	}

	const int PETSC_COOKIE = 1211216;

	int header[4];

	header[0] = PETSC_COOKIE;
	header[1] = nr;
	header[2] = nc;
	header[3] = 0;

	impl::endian::swap_int_4(&header[0]);
	impl::endian::swap_int_4(&header[1]);
	impl::endian::swap_int_4(&header[2]);

	f.write((char *)&header, sizeof(int)*4);
	f.write((char *)&rows[0], sizeof(int)*nr); 

	if (cols.size() > 0)
	{
		f.write((char *)&cols[0], sizeof(int)*nz);
	}

	if (vals.size() > 0)
	{
		f.write((char *)&vals[0], sizeof(double)*nz);
	}

	return (true);
}

/**
 * Warning: This method makes a copy of the data.  
 *
 * Todo: Dispatch this to a method which will use
 * a non-copying routine if the matrix has a row_iterator.
 */
template <class Matrix>
bool write_petsc_matrix(std::ostream& f, const Matrix& m)
{
	/* 
	 * First, we'll pack the data; this isn't the world's
	 * most efficient code, but it should be sufficient...
	 */
	using namespace yasmic;
	using namespace std;

	typedef typename smatrix_traits<Matrix>::size_type size_type;

	size_type nr = nrows(m);
	size_type nc = ncols(m);
	size_type nz = nnz(m);

	vector<int> rows(nr+1);
	vector<int> cols(nz);
	vector<double> vals(nz);

	load_matrix_to_crm(m, rows.begin(), cols.begin(), vals.begin(), false);

    return (write_petsc(f, rows, cols, vals, nr, nc, nz));
}


#endif //YASMIC_UTIL_WRITE_PETSC_MATRIX

