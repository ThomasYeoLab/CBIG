/**
 * @file bvgraph_matrix.hpp
 * Interprete a Boldi-Vigna Graph as a matrix.
 */
 
/*
 * David Gleich
 * 12 January 2007
 * Copyright, Stanford University, 2007
 */

#ifndef YASMIC_BVGRAPH_MATRIX
#define YASMIC_BVGRAPH_MATRIX

#ifdef BOOST_MSVC
#if _MSC_VER >= 1400
    // disable the warning for ifstream::read
    #pragma warning( push )
    #pragma warning( disable : 4996 )
#endif // _MSC_VER >= 1400
#endif // BOOST_MSVC

/**
 * @file bvgraph_matrix.hpp
 * Interprete a Boldi-Vigna Graph as a matrix.
 */
 
#include <iostream>
#include <istream>
#include <string>
#include <algorithm>
#include <map>
#include <vector>
#include <util/file.hpp>
#include <boost/lexical_cast.hpp>
#include <iterator>

namespace yasmic 
{
    // prototype the internal classes
    namespace impl {
        class bit_istream;
        class bvgraph_sequential_iterator;
    }

    class bvgraph_matrix 
    {
    private:
        // number of nodes
        int n;
        
        // number of arcs
        int m;
        
        // the base filename of the graph
        std::string basename;
        
        // reference counts
        static const int DEFAULT_MAX_REF_COUNT = 3;
        int _max_ref_count;
               
        static const int DEFAULT_WINDOW_SIZE = 7;
        int _window_size;
        
        static const int DEFAULT_MIN_INTERVAL_LENGTH = 3;
        int _min_interval_length;
        
        static const int DEFAULT_ZETA_K = 3;
        int _zeta_k;
      
        void load_internal()
        {
            // read properties
            std::string propfilename = basename + ".properties";
            std::string graphfilename = basename + ".graph";
            std::ifstream propif(propfilename.c_str());
            std::cout << propfilename << " " << util::file_exists(propfilename.c_str()) << std::endl;
            std::cout << propfilename << " " << !propif << std::endl;
            std::cout << graphfilename << " " << util::file_exists(graphfilename.c_str()) << std::endl; 
            if (!propif || !util::file_exists(graphfilename.c_str())) {
                assert( ("file does not exist", 0) );
                return;
            }
            
            std::map<std::string, int> options;
            
            // initialize the options we want
            options["nodes"];
            options["arcs"];
            options["windowsize"];
            options["minintervallength"];
            options["maxrefcount"];
            options["zetak"];
            
            typedef std::string::size_type position;
            std::string property_line;
            while (std::getline(propif, property_line))
            {
                // check if there is a property listed
                position eq = property_line.find('=');
                if (eq == std::string::npos) { continue; }
                std::string key = property_line.substr(0,eq);
                // TODO throw error on non-trivial compressionflags
                if (options.find(key) != options.end())
                {
                    // we want to save this key
                    std::string value = property_line.substr(eq+1);
                    options[key] = boost::lexical_cast<int>(value);
                }
            }
            
            n = options["nodes"];
            m = options["arcs"];
            _window_size = options["windowsize"];
            _min_interval_length = options["minintervallength"];
            _max_ref_count = options["maxrefcount"];
            _zeta_k = options["zetak"];
            
            // note that if compressionflags is specified, boost::lexical_cast 
            // will through an exception because it is an invalid type.
            if (n == 0 || m == 0 || _min_interval_length <= 1) {
                assert( ("property errors", 0) );
                // TODO throw an error
            }
        }
        
        // disable copy construction
        bvgraph_matrix(const bvgraph_matrix&);
        bvgraph_matrix& operator= (const bvgraph_matrix&);
        
    public:
        // constructor
        bvgraph_matrix(const char* filename)
        : n(0), m(0), basename(filename), 
          _max_ref_count(DEFAULT_MAX_REF_COUNT),
          _window_size(DEFAULT_WINDOW_SIZE),
          _min_interval_length(DEFAULT_MIN_INTERVAL_LENGTH),
          _zeta_k(DEFAULT_ZETA_K)
        {
            load_internal(); 
        }
        int num_nodes() const  { return (n); }
        int num_arcs() const { return (m); }
        
        std::string graph_filename() const { return (basename + ".graph"); }
        int max_ref_count() const  { return (_max_ref_count); }
        int window_size() const { return (_window_size); }
        int min_interval_length() const { return (_min_interval_length); }
        int zeta_k() const { return (_zeta_k); }
        
        typedef impl::bvgraph_sequential_iterator sequential_iterator;
        
        
    }; // class bvgraph_matrix
    
    namespace impl 
    {
        const int BYTELSB[] = {
            -1, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0
        };
            
    
        // Precomputed most significant bits for bytes (-1 for 0 ).
        const int BYTEMSB[] = {
            -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
        };
        
        
        class bit_istream
        {
        public:
            bit_istream(std::istream& is, const int buffer_size)
            : f(is), bufsize(buffer_size), fill(0), pos(0), avail(0), buffer(NULL), position(0),
              read_bits(0), current(0)
            {
                assert ( bufsize > 0 );
                buffer = new unsigned char[bufsize];
            }
            
            ~bit_istream() { if (buffer) { delete[] buffer; } }
            
            /**
             * Completely reset the internal buffers so that the 
             * underlying istream is positionable.
             */
            void flush() {
                position += pos;
                avail = 0;
                pos = 0;
                fill = 0;
            }
            
            /**
             * Close the bitstream
             */
            void close() {
                if (buffer) { delete[] buffer; }
                buffer = NULL;
            }
            
            int read_bit()
            {
                return read_from_current(1); 
            }
            
            /**
             * Read a fixed number of bits into an integer.
             */
            int read_int(int len) 
            {
                int i, x = 0;
                assert (len >= 0);
                assert (len <= 32);
                

                if (len <= fill) return read_from_current( len );

                len -= fill;
                x = read_from_current( fill );
                i = len >> 3;
                while (i-- != 0) { x = x << 8 | read(); }
                read_bits += len & ~7;
                
                len &= 7;

                return ( x << len ) | read_from_current( len );
            }
            
            int read_unary() 
            {
                assert ( fill < 24 );
                int x;
                
                const unsigned int current_left_aligned = current << (24 - fill) & 0xFFFFFF;
                if (current_left_aligned != 0)
                {
                    if ((current_left_aligned & 0xFF0000) != 0) { x = 7 - BYTEMSB[current_left_aligned >> 16]; }
                    else if ((current_left_aligned & 0xFF00) != 0) { x = 15 - BYTEMSB[current_left_aligned >> 8]; }
                    else { x = 23 - BYTEMSB[current_left_aligned & 0xFF]; }
                    read_bits += x + 1;
                    fill -= x + 1;
                    return (x);
                }
                
                x = fill;
                while ( (current = read()) == 0) { x += 8; }
                x += 7 - ( fill = BYTEMSB[current] );
                read_bits += x + 1;
                return (x);
            }
            
            int read_gamma()
            {
                const int msb = read_unary();
                return ( ( 1 << msb ) | read_int(msb) ) - 1;
            }
            
            int read_zeta(const int k)
            {
                const int h = read_unary();
                const int left = 1 << h * k;
                const int m = read_int(h*k + k - 1);
                if (m < left) { return (m + left - 1); }
                else { return (m << 1) + read_bit() - 1; }
            }
            
        private:
            // the underlying file stream
            std::istream &f;
            // the number of bits actually read
            long read_bits;
            // current bit buffer, the lowest fill bits represent the current contents
            int current;
            // the stream buffer
            unsigned char *buffer;
            // the buffer size
            unsigned int bufsize;
            // current number of bits in the bit buffer (stored low)
            int fill;
            // current position in the byte buffer
            int pos;
            // current number of of bytes available in the byte buffer
            int avail;
            // current position of the first byte in the byte buffer
            long position;
            
            /**
             * Read the next byte from the underlying stream.
             * 
             * This method does not update read_bits.
             * 
             * @return the byte
             */
            int read() 
            {
                if (avail == 0) 
                {
                    f.read((char*)buffer, bufsize);
                    avail = f.gcount(); 
                    if (avail <= 0) {
                        // throw an exception
                    }
                    else {
                        position += pos;
                        pos = 0;
                    }
                }
                
                avail--;
                return buffer[pos++] & 0xFF;
            }
            
            /** 
             * Fills {@link #current} to 16 bits.
             */
            int refill16() 
            {
                assert ( fill >= 8 );
                assert ( fill < 16 );
                
                if (avail > 0) {
                    // if there is a current byte in the buffer, use it directly.
                    avail--;
                    current = (current << 8) | (buffer[pos++] & 0xFF);
                    return (fill += 8);
                }
                
                current = (current << 8) | read();
                return (fill += 8);
            }
            
            int refill() 
            {
                if (fill == 0) {
                    current = read();
                    return (fill = 8);
                }
                
                if (avail > 0) {
                    avail--;
                    current = (current << 8) | (buffer[pos++] & 0xFF);
                    return (fill += 8);
                }
                
                current = (current << 8) | read();
                return (fill += 8);
            }
            
            /**
             * Read bits from the buffer, possibly refilling it.
             */
            int read_from_current(const int len) 
            {
                if (len == 0) { return 0; }
                if (fill == 0) { current = read(); fill = 8; }
                read_bits += len;
                unsigned int rval = (unsigned)current;
                return (rval >> ( fill -= len) & (1 << len) - 1);
            }
        }; // class bit_istream
        
        class bvgraph_sequential_iterator
        {
        private:
            // the graph size
            int n;
            
            // the underlying stream
            std::ifstream graph_stream;
            bit_istream bis;
            
            int max_ref_count;
            int window_size;
            int min_interval_length;
            int zeta_k;
            
            // variables for the internal iterators
            bool _row_arcs_end;
            bool _rows_end; 
            int cyclic_buffer_size;
            std::vector<int> outd;
            int curr;
            int curr_outd;
            int curr_arc;
            std::vector<std::vector<int> > window;
            
            std::vector<int> buffer1;
            std::vector<int> buffer2;
            std::vector<int> arcs;
            
            int read_offset() { return (bis.read_gamma()); }    
            int read_outdegree() { return (bis.read_gamma()); }
            int read_reference() { return (bis.read_unary()); }
            int read_block() { return (bis.read_gamma()); }
            int read_block_count() { return (bis.read_gamma()); }
            int read_residual() { return (bis.read_zeta(zeta_k)); }
            
            int nat2int(const int x) { return x % 2 == 0 ? x >> 1 : -( ( x + 1 ) >> 1 ); }
            
            void load_successors()
            {
                const int x = curr;
                int ref, ref_index;
                int i, extra_count, block_count = 0;
                std::vector<int> block, left, len;
                
                int d = outd[x%cyclic_buffer_size]=read_outdegree();
                curr_outd = d;
                if (d == 0) { return; }
                
                // std::cout << "** Start successors " << x << " outdegree " << curr_outd << std::endl;
                
                // we read the reference only if the actual window size is larger than one 
                // (i.e., the one specified by the user is larger than 0).
                if ( window_size > 0 ) {
                    ref = read_reference();
                }
                
                ref_index = (x - ref + cyclic_buffer_size) % cyclic_buffer_size;
                
                if (ref > 0)
                {
                    if ( (block_count = read_block_count()) != 0 ) {
                        block.resize(block_count);
                    }
                    
                    // std::cout << "block_count = " << block_count << std::endl;
                    
                    // the number of successors copied, and the total number of successors specified
                    // in some copy
                    int copied = 0, total = 0;
                    
                    for (i = 0; i < block_count; i++) {
                        block[i] = read_block() + (i == 0 ? 0 : 1);
                        total += block[i];
                        if (i % 2 == 0) {
                            copied += block[i];
                        }
                    }
                    if (block_count%2 == 0) {
                        copied += (outd[ref_index] - total);
                    }
                    extra_count = d - copied;
                }
                else {
                    extra_count = d;
                }
                
                int interval_count = 0;
                if (extra_count > 0)
                {
                    if (min_interval_length != 0 && (interval_count = bis.read_gamma()) != 0) 
                    {
                        int prev = 0;
                        left.resize(interval_count);
                        len.resize(interval_count);
                        
                        // now read the intervals
                        left[0] = prev = nat2int(bis.read_gamma()) + x;
                        len[0] = bis.read_gamma() + min_interval_length;
                        
                        prev += len[0];
                        extra_count -= len[0];
                        
                        for (i=1; i < interval_count; i++) {
                            left[i] = prev = bis.read_gamma() + prev + 1;
                            len[i] = bis.read_gamma() + min_interval_length;
                            prev += len[i];
                            extra_count -= len[i];
                        }
                    }
                }
                
                // allocate a sufficient buffer for the output
                if (arcs.size() < d) {
                    arcs.resize(d);
                    buffer1.resize(d);
                    buffer2.resize(d);
                }
                
                int buf1_index = 0;
                int buf2_index = 0;
                
                // std::cout << "extra_count = " << extra_count << std::endl;
                // std::cout << "interval_count = " << interval_count << std::endl;
                // std::cout << "ref = " << ref << std::endl;
                
                // read the residuals into a buffer
                {
                    int prev = -1;
                    int residual_count = extra_count;
                    while (residual_count > 0) {
                        residual_count--;
                        if (prev == -1) { buffer1[buf1_index++] = prev = x + nat2int(read_residual()); }
                        else { buffer1[buf1_index++] = prev = read_residual() + prev + 1; }
                    }
                    // std::cout << "residuals: buf1" << std::endl;
                    // std::copy(buffer1.begin(),buffer1.begin()+buf1_index,std::ostream_iterator<int>(std::cout, "\n"));
                }
                // std::cout << "buf1_index = " << buf1_index << std::endl;
                    
                if (interval_count == 0)
                {
                    // don't do anything
                }
                else
                {
                    // copy the extra interval data
                    for (i = 0; i < interval_count; i++)
                    {
                        int cur_left = left[i];
                        for (int j = 0; j < len[i]; j++) {
                            buffer2[buf2_index++] = cur_left + j;
                        }
                    }
                    
                    // std::cout << "sequences: buf2" << std::endl;
                    // std::copy(buffer2.begin(),buffer2.begin()+buf2_index,std::ostream_iterator<int>(std::cout, "\n"));
                    
                    if (extra_count > 0)
                    {
                        std::merge(
                            buffer1.begin(), buffer1.begin()+buf1_index,
                            buffer2.begin(), buffer2.begin()+buf2_index,
                            arcs.begin()
                            );
                        buf1_index = buf1_index + buf2_index;
                        buf2_index = 0;           
                        std::copy(arcs.begin(), arcs.end(),
                            buffer1.begin());
                    }
                    else
                    {
                        std::copy(buffer2.begin(), buffer2.begin()+buf2_index,
                            buffer1.begin());
                        buf1_index = buf2_index;
                        buf2_index = 0;
                    }
                }
                
                if (ref <= 0)
                {
                    // don't do anything except copy
                    // the data to arcs
                    if (interval_count == 0 || extra_count == 0) {
                        std::copy(buffer1.begin(), buffer1.end(),
                            arcs.begin()
                            );
                    }
                }
                else
                {          
                    // TODO clean this code up          
                    // copy the information from the masked iterator
                    
                    int mask_index = 0;
                    int len = 0;
                    for (i=0; i < outd[ref_index]; )
                    {
                        if (len <= 0)
                        {
                            if (block_count == mask_index) 
                            {
                                if (block_count % 2 == 0) {
                                    len = outd[ref_index] - i;
                                }
                                else {
                                    break;
                                }
                            }
                            else {
                                if (mask_index % 2 == 0) { len = block[mask_index++]; }
                                else { i += block[mask_index++]; continue; }
                            }
                            
                            // in the case that length is 0, we continue.
                            if (len == 0) { continue; }
                        }
                        buffer2[buf2_index++] = window[ref_index][i];
                        len--;
                        i++;
                    }
                    
                    // std::cout << "masked" << std::endl;
                    // std::copy(buffer2.begin(),buffer2.begin()+buf2_index,std::ostream_iterator<int>(std::cout, "\n"));
                    
                    std::merge(
                        buffer1.begin(), buffer1.begin()+buf1_index,
                        buffer2.begin(), buffer2.begin()+buf2_index,
                        arcs.begin() 
                        );
                    buf1_index = buf1_index + buf2_index;
                    buf2_index = 0;
                    
                }

                // std::cout << "arcs" << std::endl;
                // std::copy(arcs.begin(),arcs.begin()+d,std::ostream_iterator<int>(std::cout, "\n"));
                assert (buf1_index == d);
                // std::cout << "end arcs" << std::endl;
            }
            
        public:
            // constructor
            bvgraph_sequential_iterator(
                const bvgraph_matrix& m)
            : n(m.num_nodes()),
              graph_stream(m.graph_filename().c_str(),std::ios::binary),
              bis(graph_stream, 16*1024),
              max_ref_count(m.max_ref_count()),
              window_size(m.window_size()),
              min_interval_length(m.min_interval_length()),
              zeta_k(m.zeta_k()),
              _row_arcs_end(false),
              _rows_end(false),
              cyclic_buffer_size(window_size+1),
              outd(cyclic_buffer_size),
              curr(-1),
              window(cyclic_buffer_size)                
            {}
                    
            
            // reset the iterators associated with this graph,
            // the row pointer is set to the first row
            void reset()
            {
                // reset all the file pointers
                graph_stream.seekg(0, std::ios_base::beg);
                bis.flush();
                
                _row_arcs_end = false;
                _rows_end = false;
                curr = -1;
            }
            
            // step to the next row of the matrix
            void next_row()
            {
                // check if we are done
                curr++;
                if (curr > n-1) { _rows_end = true; return; }
                
                int curr_index = curr % cyclic_buffer_size;
                load_successors();
                curr_arc = -1;
                _row_arcs_end = false;
                if (window[curr_index].size() < outd[curr_index]) {
                    window[curr_index].resize(outd[curr_index]);
                }
                
                // unwrap the buffered output
                for (int i = 0; i < outd[curr_index]; i++) {
                    window[curr_index][i] = arcs[i];
                } 
                
                // check to make sure there is something to do
                if (curr_outd == 0) { _row_arcs_end = true; }
            }
            // get the index of the current row
            int cur_row() { return (curr); }
            // get the current outdegree
            int cur_row_outdegree() { return (curr_outd); }
            // returns true when there are no more rows
            bool rows_end() { return (_rows_end); }
            
            // step to the next arc for the current row
            void next_row_arc()
            {
                curr_arc++;
                if (curr_arc >= curr_outd - 1) { _row_arcs_end = true; return; }
            }
            // get the index of the target of the current arc
            int cur_row_arc_target() { return (arcs[curr_arc]); }
            // returns true when there are no more arcs for the current row
            bool row_arcs_end() { return (_row_arcs_end); }
        }; // class bvgraph_iterator
    } // namespace impl  
    
              
}// namespace yasmic

