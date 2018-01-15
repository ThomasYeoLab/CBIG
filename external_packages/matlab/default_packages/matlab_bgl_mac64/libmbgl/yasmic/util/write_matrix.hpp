#ifndef YASMIC_UTIL_WRITE_MATRIX
#define YASMIC_UTIL_WRITE_MATRIX

#include <fstream>

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

    namespace write
    {
        template <class index_type, class nz_index_type, class value_type, bool header>
        struct custom_smat_writer
        {
            template <class Matrix>
            bool write_matrix(std::ostream& f, Matrix& m)
            {
                using namespace yasmic;

                if (header)
                {
                    f << (index_type)nrows(m) << " " 
                            << (index_type)ncols(m) << " " 
                            << (nz_index_type)nnz(m) << std::endl;
                }

                typename smatrix_traits<Matrix>::nonzero_iterator nzi, nziend;
                boost::tie(nzi,nziend) = nonzeros(m);
                for (; nzi != nziend; ++nzi)
                {
                    f << (index_type)row(*nzi, m) << " "
                            << (index_type)column(*nzi, m) << " "
                            << (value_type)value(*nzi, m) << std::endl;
                }

                return (true);
            }
        };

        template <class index_type, class nz_index_type, class value_type, bool header>
        struct custom_bsmat_writer
        {
            template <class Matrix>
            bool write_matrix(std::ostream& f, Matrix& m)
            {
                using namespace yasmic;

                if (header)
                {
                	index_type nr = (index_type)nrows(m);
                	index_type nc = (index_type)ncols(m);
                	nz_index_type nz = (index_type)nnz(m);
                	
                	f.write((char*)&nr, sizeof(index_type));
                	f.write((char*)&nc, sizeof(index_type));
                	f.write((char*)&nz, sizeof(nz_index_type));
                }

                typename smatrix_traits<Matrix>::nonzero_iterator nzi, nziend;
                boost::tie(nzi,nziend) = nonzeros(m);
                for (; nzi != nziend; ++nzi)
                {
                	index_type r = row(*nzi, m);
                	index_type c = column(*nzi, m);
                	value_type v = value(*nzi, m);
                	
                	f.write((char*)&r, sizeof(index_type));
                	f.write((char*)&c, sizeof(index_type));
                	f.write((char*)&v, sizeof(value_type));
                }

                return (true);
            }
        };

        struct smat_writer
            : public custom_smat_writer<int, unsigned int, double, true>
        {
        };

        struct bsmat_writer
            : public custom_bsmat_writer<unsigned int, unsigned int, double, true>
        {
        };

        struct cluto_writer
        {
            template <class RowAccessMatrix>
            bool write_matrix(std::ostream& f, RowAccessMatrix& m)
            {
                using namespace yasmic;
                f << nrows(m) << " " 
                  << ncols(m) << " " 
                  << nnz(m) << std::endl;
                  
                typename smatrix_traits<RowAccessMatrix>::row_iterator ri, riend;
                typename smatrix_traits<RowAccessMatrix>::row_nonzero_iterator rnzi, rnziend;
                for (boost::tie(ri,riend) = rows(m); ri!=riend; ++ri) {
                    for (boost::tie(rnzi,rnziend)=row_nonzeros(*ri,m); rnzi!=rnziend; ++rnzi) {
                        f << (column(*rnzi,m)+1) << " " 
                          << (value(*rnzi,m)) << " ";
                    }
                    f << std::endl;
                }
                
                return (true);
            }
        };
        
        

        struct petsc_writer
        {
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
					unsigned int maxi = (unsigned int)rows.size();
					for (unsigned int i = 0; i < maxi; ++i)
					{
			            impl::endian::swap_int_4(intptr);
						++intptr;
					}
				}
				
				if (cols.size() > 0)
				{
					int* intptr = &cols[0];
					unsigned int maxi = (unsigned int)cols.size();
					for (unsigned int i = 0; i < maxi; ++i)
					{
						impl::endian::swap_int_4(intptr);
						++intptr;
					}
				}
			
				if (vals.size() > 0)
				{
					double* doubleptr = &vals[0];
					unsigned int maxi = (unsigned int)vals.size();
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
                    f.write((char *)&cols[0], sizeof(int)*(std::streamsize)nz);
				}
			
				if (vals.size() > 0)
				{
					f.write((char *)&vals[0], sizeof(double)*(std::streamsize)nz);
				}
			
				return (true);
			}
			
            template <class NonzeroAccessMatrix>
            bool write_matrix(std::ostream& f, NonzeroAccessMatrix& m)
            {
                /* 
				 * First, we'll pack the data; this isn't the world's
				 * most efficient code, but it should be sufficient...
				 */
				using namespace yasmic;
				using namespace std;
			
				typedef typename smatrix_traits<NonzeroAccessMatrix>::size_type size_type;
			
				size_type nr = nrows(m);
				size_type nc = ncols(m);
				size_type nz = nnz(m);
			
				vector<int> rows(nr+1);
				vector<int> cols(nz);
				vector<double> vals(nz);
			
				load_matrix_to_crm(m, rows.begin(), cols.begin(), vals.begin(), false);
			
			    return (write_petsc_matrix(f, rows, cols, vals, nr, nc, nz));
            }
        };

    }
}

typedef struct impl::write::smat_writer smat_writer;
typedef struct impl::write::bsmat_writer bsmat_writer;
typedef struct impl::write::petsc_writer petsc_writer;
typedef struct impl::write::cluto_writer cluto_writer;

template <class index_type, class nz_index_type, class value_type>
struct parametrized_writers
{
    typedef struct impl::write::custom_bsmat_writer<index_type,nz_index_type,value_type,true>
        bsmat_writer;
    typedef struct impl::write::custom_smat_writer<index_type,nz_index_type,value_type,true>
        smat_writer;
};

template <class Matrix, class Tag>
void write_matrix(std::ostream& f, Matrix& m, Tag type)
{
    type.write_matrix(f,m);
}


#endif //YASMIC_UTIL_WRITE_MATRIX

