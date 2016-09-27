//============================================================================
//
//  This file is part of GPSTk, the GPS Toolkit.
//
//  The GPSTk is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation; either version 3.0 of the License, or
//  any later version.
//
//  The GPSTk is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with GPSTk; if not, write to the Free Software Foundation,
//  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110, USA
//  
//  Copyright 2004, The University of Texas at Austin
//
//============================================================================

//============================================================================
//
//This software developed by Applied Research Laboratories at the University of
//Texas at Austin, under contract to an agency or agencies within the U.S. 
//Department of Defense. The U.S. Government retains all rights to use,
//duplicate, distribute, disclose, or release this software. 
//
//Pursuant to DoD Directive 523024 
//
// DISTRIBUTION STATEMENT A: This software has been approved for public 
//                           release, distribution is unlimited.
//
//=============================================================================

/**
 * @file ZStreamBuf.cpp
 * Stream buffer wrapper for reading gzip-compressed data
 */

#include <cstring>
#include <stdexcept>
#include "ZStreamBuf.hpp"

#define BUFSIZE      4
#define IN_BUFSIZE   (BUFSIZE + 64)
#define OUT_BUFSIZE  (BUFSIZE + 2048)
#define BITS         16
#define INIT_BITS    9			/* initial number of bits/code */
#define MAXCODE(n)   (1L << (n))
#define FIRST        257					/* first free entry */
#define CLEAR        256

namespace gpstk
{
    ZStreamBuf ::
    ZStreamBuf(std::streambuf* _source,
                  size_t stream_buffer_size)
        : source(_source), outbuffer(stream_buffer_size)
    {
        setg(&outbuffer[0], &outbuffer[0], &outbuffer[0]);
        const char* LZW_MAGIC = "\037\235";
        char buf[3];

        if (source->sgetn(buf, 3) != 3)
            return;

        if (strncmp(buf, LZW_MAGIC, 2) != 0 || buf[2] & 0x60)
            throw new std::runtime_error("Not a LZW-compressed stream");

        memset(&lzwstate, 0x00, sizeof(lzwstate));
        lzwstate.eof = 0;
        lzwstate.inbuf.resize(IN_BUFSIZE);
        lzwstate.outbuf.resize(OUT_BUFSIZE);
        lzwstate.stackp = NULL;
        lzwstate.insize = 3; /* we read three bytes above */
        lzwstate.outpos = 0;
        lzwstate.rsize = 0;

        lzwstate.flags = (unsigned char)buf[2];
        lzwstate.maxbits = lzwstate.flags & 0x1f;    /* Mask for 'number of compresssion bits' */
        lzwstate.block_mode = lzwstate.flags & 0x80;

        lzwstate.n_bits = INIT_BITS;
        lzwstate.maxcode = MAXCODE(INIT_BITS) - 1;
        lzwstate.bitmask = (1<<INIT_BITS)-1;
        lzwstate.oldcode = -1;
        lzwstate.finchar = 0;
        lzwstate.posbits = 3<<3;
        lzwstate.free_ent = ((lzwstate.block_mode) ? FIRST : 256);

        /* initialize the first 256 entries in the table */
        memset(lzwstate.codetab, 0x00, sizeof(lzwstate.codetab));
        for (lzwstate.code = 255; lzwstate.code >= 0; --lzwstate.code)
            lzwstate.htab[lzwstate.code] = (unsigned long int)lzwstate.code;

        if (lzwstate.maxbits > BITS) {
            throw new std::runtime_error("LZW bits flag exceeds maximum");
        }
    }
    
    ZStreamBuf ::
    ~ZStreamBuf() {}

    std::streambuf::int_type ZStreamBuf ::
    underflow()
    {
        // Refill the buffer
        ssize_t bytes_read = lzw_read(&outbuffer[0], outbuffer.size());
        if (bytes_read < 0)
        {
            throw new std::runtime_error("Error while decompressing");
        }
        setg(&outbuffer[0], &outbuffer[0], &outbuffer[0]+bytes_read);
        if (bytes_read == 0)
            return traits_type::eof();
        else
            return traits_type::to_int_type(*gptr());
    }



#ifndef	NOALIGN
#define NOALIGN 0
#endif

#if defined(WORDS_BIGENDIAN) && NOALIGN == 1
    # define input(b,o,c,n,m) \
	do { \
		(c) = (*(long *)(&(b)[(o)>>3])>>((o)&0x7))&(m); \
		(o) += (n); \
	} while (0)
#else
# define input(b,o,c,n,m) \
	do { \
		unsigned char *p = (unsigned char*)&(b)[(o)>>3]; \
		(c) = ((((long)(p[0]))|((long)(p[1])<<8)| \
		       ((long)(p[2])<<16))>>((o)&0x7))&(m); \
		(o) += (n); \
	} while (0)
#endif

    ssize_t ZStreamBuf::lzw_read(char *readbuf, size_t count)
    {
        char* de_stack = (char *)&(lzwstate.htab[HSIZE-1]);

        size_t count_left = count;
        char *inbuf = &lzwstate.inbuf[0];
        char *outbuf = &lzwstate.outbuf[0];

        long int maxmaxcode = MAXCODE(lzwstate.maxbits);

        if (!count || lzwstate.eof)
            return 0;

        if (lzwstate.stackp != NULL) {
            if (lzwstate.outpos) {
                if (lzwstate.outpos >= count) {
                    outbuf = lzwstate.unreadbuf;
                    goto empty_existing_buffer;
                } else /*if (lzwstate.outpos < count)*/ {
                    memcpy(readbuf, lzwstate.unreadbuf, lzwstate.outpos);
                    goto resume_partial_reading;
                }
            }
            goto resume_reading;
        }

        do {
            resetbuf:
            {
                size_t i, e, o;
                o = lzwstate.posbits >> 3;
                e = o <= lzwstate.insize ? lzwstate.insize - o : 0;

                for (i = 0; i < e; ++i)
                    inbuf[i] = inbuf[i+o];

                lzwstate.insize = e;
                lzwstate.posbits = 0;
            }

            if (lzwstate.insize < IN_BUFSIZE-BUFSIZE) {
                if ((lzwstate.rsize = source->sgetn(inbuf+lzwstate.insize, BUFSIZE)) < 0)
                    return -1;
                lzwstate.insize += lzwstate.rsize;
            }

            lzwstate.inbits = ((lzwstate.rsize > 0) ? (lzwstate.insize - lzwstate.insize%lzwstate.n_bits)<<3 :
                           (lzwstate.insize<<3) - (lzwstate.n_bits-1));

            while (lzwstate.inbits > lzwstate.posbits) {
                if (lzwstate.free_ent > lzwstate.maxcode) {
                    lzwstate.posbits = ((lzwstate.posbits-1) + ((lzwstate.n_bits<<3) -
                                                        (lzwstate.posbits-1 + (lzwstate.n_bits<<3)) % (lzwstate.n_bits<<3)));

                    ++lzwstate.n_bits;
                    if (lzwstate.n_bits == lzwstate.maxbits)
                        lzwstate.maxcode = maxmaxcode;
                    else
                        lzwstate.maxcode = MAXCODE(lzwstate.n_bits)-1;

                    lzwstate.bitmask = (1 << lzwstate.n_bits) - 1;
                    goto resetbuf;
                }

                input(inbuf,lzwstate.posbits,lzwstate.code,lzwstate.n_bits,lzwstate.bitmask);

                if (lzwstate.oldcode == -1) {
                    if (lzwstate.code >= 256) return -1; /* error("corrupt input."); */
                    outbuf[lzwstate.outpos++] = lzwstate.finchar = lzwstate.oldcode = lzwstate.code;
                    continue;
                }

                if (lzwstate.code == CLEAR && lzwstate.block_mode) {
                    memset(lzwstate.codetab, 0x00, sizeof(lzwstate.codetab));
                    lzwstate.free_ent = FIRST - 1;
                    lzwstate.posbits = ((lzwstate.posbits-1) + ((lzwstate.n_bits<<3) -
                                                        (lzwstate.posbits-1 + (lzwstate.n_bits<<3)) % (lzwstate.n_bits<<3)));
                    lzwstate.maxcode = MAXCODE(lzwstate.n_bits = INIT_BITS)-1;
                    lzwstate.bitmask = (1 << lzwstate.n_bits) - 1;
                    goto resetbuf;
                }

                lzwstate.incode = lzwstate.code;
                lzwstate.stackp = de_stack;

                /* Special case for KwKwK string.*/
                if (lzwstate.code >= lzwstate.free_ent) {
                    if (lzwstate.code > lzwstate.free_ent) {
                        return -1;
                    }

                    *--lzwstate.stackp = lzwstate.finchar;
                    lzwstate.code = lzwstate.oldcode;
                }

                /* Generate output characters in reverse order */
                while (lzwstate.code >= 256) {
                    *--lzwstate.stackp = lzwstate.htab[lzwstate.code];
                    lzwstate.code = lzwstate.codetab[lzwstate.code];
                }

                *--lzwstate.stackp = (lzwstate.finchar = lzwstate.htab[lzwstate.code]);

                /* And put them out in forward order */
                {
                    lzwstate.stackp_diff = de_stack - lzwstate.stackp;

                    if (lzwstate.outpos+lzwstate.stackp_diff >= BUFSIZE) {
                        do {
                            if (lzwstate.stackp_diff > BUFSIZE-lzwstate.outpos)
                                lzwstate.stackp_diff = BUFSIZE-lzwstate.outpos;

                            if (lzwstate.stackp_diff > 0) {
                                memcpy(outbuf+lzwstate.outpos, lzwstate.stackp, lzwstate.stackp_diff);
                                lzwstate.outpos += lzwstate.stackp_diff;
                            }

                            if (lzwstate.outpos >= BUFSIZE) {
                                if (lzwstate.outpos < count_left) {
                                    memcpy(readbuf, outbuf, lzwstate.outpos);
                                    resume_partial_reading:
                                    readbuf += lzwstate.outpos;
                                    count_left -= lzwstate.outpos;
                                } else {
                                    empty_existing_buffer:
                                    lzwstate.outpos -= count_left;
                                    memcpy(readbuf, outbuf, count_left);
                                    lzwstate.unreadbuf = outbuf + count_left;
                                    return count;
                                }
                                resume_reading:
                                lzwstate.outpos = 0;
                            }
                            lzwstate.stackp += lzwstate.stackp_diff;
                        } while ((lzwstate.stackp_diff = (de_stack-lzwstate.stackp)) > 0);
                    } else {
                        memcpy(outbuf+lzwstate.outpos, lzwstate.stackp, lzwstate.stackp_diff);
                        lzwstate.outpos += lzwstate.stackp_diff;
                    }
                }

                /* Generate the new entry. */
                if ((lzwstate.code = lzwstate.free_ent) < maxmaxcode) {
                    lzwstate.codetab[lzwstate.code] = (unsigned short)lzwstate.oldcode;
                    lzwstate.htab[lzwstate.code] = (unsigned long int)lzwstate.finchar;
                    lzwstate.free_ent = lzwstate.code+1;
                }

                lzwstate.oldcode = lzwstate.incode;	/* Remember previous code. */
            }
        } while (lzwstate.rsize != 0);

        if (lzwstate.outpos < count_left) {
            lzwstate.eof = 1;
            memcpy(readbuf, outbuf, lzwstate.outpos);
            count_left -= lzwstate.outpos;
            return (count - count_left);
        } else {
            goto empty_existing_buffer;
        }
    }

    std::streampos ZStreamBuf::seekoff(std::streamoff off, std::ios::seekdir dir, std::ios_base::openmode which) {
        if (egptr() - eback() == 0)
            // No data present, fill buffer
            underflow();
        if (dir == std::ios_base::cur)
            gbump(off);
        return gptr() - eback();
    }
}

