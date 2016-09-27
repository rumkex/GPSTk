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
 * @file ZStreamBuf.hpp
 * Stream buffer wrapper for reading LZW-compressed data
 * Public domain LZW decompression code is used here.
 */

#ifndef GPSTK_ZSTREAMBUF_HPP
#define GPSTK_ZSTREAMBUF_HPP

#include <ios>
#include <streambuf>
#include <vector>

#ifdef WIN32
#define ssize_t unsigned long
#endif

namespace gpstk
{

    class ZStreamBuf : public std::streambuf
    {
        static const size_t HSIZE = (1<<17);

    private:
        ssize_t lzw_read(char* readbuf, size_t count);

    protected:
        std::streambuf* source;
        std::vector<char> outbuffer;

        struct {
            int eof;

            std::vector<char> inbuf;
            std::vector<char> outbuf;
            char *stackp, *unreadbuf;
            size_t stackp_diff;
            size_t insize, outpos;
            ssize_t rsize;

            unsigned char flags;
            int maxbits, block_mode;

            unsigned long int htab[HSIZE];
            unsigned short codetab[HSIZE];

            int n_bits, posbits, inbits, bitmask, finchar;
            long int maxcode, oldcode, incode, code, free_ent;
        } lzwstate;

        virtual std::streambuf::int_type underflow();

        virtual std::streampos seekoff(std::streamoff off,
                                       std::ios::seekdir way,
                                       std::ios::openmode which = std::ios::in | std::ios::out);

    public:
        ZStreamBuf(std::streambuf* _source,
                      size_t stream_buffer_size = 1024);
        ~ZStreamBuf();
    };
}

#endif

