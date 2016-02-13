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
 * @file CRinexStreamBuf.hpp
 * Stream buffer wrapper for reading Hatanaka-compressed data
 */

#ifndef GPSTK_CRINEXSTREAMBUF_HPP
#define GPSTK_CRINEXSTREAMBUF_HPP

#include <streambuf>
#include <istream>
#include <vector>
#include <string>
#include <climits>

namespace gpstk {

    class CRinexStreamBuf : public std::streambuf {
    protected:
        std::vector<char> out_buffer;
        char* out_ptr;
        std::string in_buffer, epoch_buffer;
        std::istream input;

    private:
        static const int MAXSAT = 90;
        static const int MAXTYPE = 20;
        static const int MAX_DIFF_ORDER = 5;

        typedef struct data_format {
            long u[MAX_DIFF_ORDER + 1];
            /* upper X digits for each difference order */
            long l[MAX_DIFF_ORDER + 1];
            /* lower 5 digits */
            int order;
            int arc_order;
        } data_format;

        bool header_done;

        struct {
            int version;
            int rinex_version;
            char ep_top_from, ep_top_to;
            char* p_event, * p_nsat, * p_satlst;
            int shift_clk;
            size_t offset;
            int nsat1;

            data_format clk1;
            data_format dy1[MAXSAT][MAXTYPE], dy0[MAXSAT][MAXTYPE];
            char flag1[MAXSAT][MAXTYPE * 2 + 1], flag[MAXSAT][MAXTYPE * 2 + 1];

            int sattbl[MAXSAT];
            char sat_lst_old[MAXSAT * 3];
            char dflag[MAXSAT][MAXTYPE * 2];

            int ntype, ntype_gnss[UCHAR_MAX],
                ntype_record[MAXSAT];
        } crx;

    public:
        CRinexStreamBuf(std::streambuf* _source,
                        size_t stream_buffer_size = 131072,
                        size_t in_buffer_size = 1024);

        ~CRinexStreamBuf();

    protected:
        virtual std::streambuf::int_type underflow();

    private:
        void trim();

        void read_preamble();

        bool read_header();

        bool read_data();

        void repair(std::string& outstr, const std::string& instr);

        void skip_to_next();

        void set_sat_table(char* p_new, char p_old[], int nsat1, int sattbl[]);

        void write_data(int nsat, char* p_sat_lst, int* sattbl, char dflag[][MAXTYPE * 2]);

        void putfield(data_format* y, char* flag);

        bool put_event_data(std::string& instr, char* p_event);

        bool getdiff(int i);

        void read_clock(std::string& instr, data_format* clk);

        void print_clock(long yu, long yl, int shift_clk);

        void process_clock(const data_format& oldclk);
    };
}

#endif
