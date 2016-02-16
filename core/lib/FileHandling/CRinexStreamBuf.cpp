#include <stdexcept>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <algorithm>

#include "CRinexStreamBuf.hpp"

namespace gpstk {

    CRinexStreamBuf::
    CRinexStreamBuf(std::streambuf* _source,
                    size_t stream_buffer_size,
                    size_t in_buffer_size)
        : input(_source), header_done(false) {
        in_buffer.resize(in_buffer_size);
        epoch_buffer.resize(in_buffer_size);
        out_buffer.resize(stream_buffer_size);
        read_preamble();
        // force underflow
        char* start = &out_buffer[0];
        setg(start, start, start);
    }

    CRinexStreamBuf::
    ~CRinexStreamBuf() {
    }

    std::streambuf::int_type CRinexStreamBuf::
    underflow() {
        out_ptr = &out_buffer[0];
        if (!header_done) {
            if (read_header()) {
                setg(eback(), eback(), out_ptr);
                return traits_type::to_int_type(*eback());
            }
        }
        else {
            if (read_data()) {
                setg(eback(), eback(), out_ptr);
                return traits_type::to_int_type(*eback());
            }
        }
        return traits_type::eof();
    }

    void CRinexStreamBuf::
    read_preamble() {
        std::getline(input, in_buffer);
        bool isVersionCorrect = in_buffer.compare(0, 3, "1.0") == 0 ||
                                in_buffer.compare(0, 3, "3.0") == 0;
        bool isHeaderCorrect = in_buffer.compare(60, 20, "CRINEX VERS   / TYPE") == 0;
        if (!input ||
            !isVersionCorrect ||
            !isHeaderCorrect) {
            throw new std::runtime_error("Not a compressed RINEX file");
        }
        crx.version = std::atoi(&in_buffer[0]);
        crx.rinex_version = -1;
        std::getline(input, in_buffer);
    }

    bool CRinexStreamBuf::read_header() {
        if (!std::getline(input, in_buffer)) {
            return false;
        }
        if (crx.rinex_version == -1) {
            // First line has to be RINEX header
            if (std::string("RINEX VERSION / TYPE").compare(0, 20, &in_buffer[60]) != 0 ||
                (in_buffer[5] != '2' && in_buffer[5] != '3')) {
                    throw new std::runtime_error("Incorrect RINEX header");
            }
            crx.rinex_version = atoi(&in_buffer[0]);
        }
        else if (in_buffer.compare(60, 19, "# / TYPES OF OBSERV") == 0 && in_buffer[5] != ' ') { /** for RINEX2 **/
            crx.ntype = atoi(&in_buffer[0]);
        }
        else if (in_buffer.compare(60, 19, "SYS / # / OBS TYPES") == 0) { /** for RINEX3  **/
            if (in_buffer[0] != ' ') {
                crx.ntype_gnss[(size_t) in_buffer[0]] = atoi(&in_buffer[3]);
            }
            if (crx.ntype_gnss[(size_t) in_buffer[0]] > MAXTYPE) {
                throw new std::runtime_error("Too many observation types");
            }
        }
        else if (in_buffer.compare(60, 13, "END OF HEADER") == 0) {
            header_done = true;
            crx.clk1.order = 0;
            crx.clk1.arc_order = 0;
            if (crx.rinex_version == 2) {
                crx.ep_top_from = '&';
                crx.ep_top_to = ' ';
                crx.p_event = &in_buffer[28];  /** pointer to event flag **/
                crx.p_nsat = &epoch_buffer[29];  /** pointer to n_sat **/
                crx.p_satlst = &epoch_buffer[32];  /** pointer to address to add satellite list **/
                crx.shift_clk = 1;
                crx.offset = 3;
            } else {
                crx.ep_top_from = '>';
                crx.ep_top_to = '>';
                crx.p_event = &in_buffer[31];
                crx.p_nsat = &epoch_buffer[32];
                crx.p_satlst = &epoch_buffer[41];
                crx.shift_clk = 4;
                crx.offset = 6;
            }
        }
        out_ptr += sprintf(out_ptr, "%s\n", in_buffer.c_str());
        return true;
    }

    bool CRinexStreamBuf::read_data() {
        while (std::getline(input, in_buffer)) {
            SKIP:
            if (crx.version == 3 && in_buffer[0] == '&')
            {
                // Skip escape lines of CRINEX version 3
                continue;
            }
            if (in_buffer[0] == crx.ep_top_from) {
                in_buffer[0] = crx.ep_top_to;
                if (*crx.p_event != '0' && *crx.p_event != '1') {
                    if (!put_event_data(in_buffer, crx.p_event)) {
                        skip_to_next();
                    }
                    goto SKIP;
                }
                // Initialize new epoch
                epoch_buffer.clear();          /**** initialize arc for epoch data ***/
                crx.nsat1 = 0;               /**** initialize the all satellite arcs ****/
            }

            // Decode delta encoding
            repair(epoch_buffer, in_buffer);

            if (epoch_buffer[0] != crx.ep_top_to || epoch_buffer.length() < (crx.offset + 26) || epoch_buffer[crx.offset + 23] != ' '
                || epoch_buffer[crx.offset + 24] != ' ' || !isdigit(epoch_buffer[crx.offset + 25])) {
                // Incorrect epoch record
                skip_to_next();
                goto SKIP;
            }

            int nsat = atoi(crx.p_nsat);
            if (nsat > MAXSAT) {
                throw new std::runtime_error("Satellite number exceeds the maximum");
            }

            set_sat_table(crx.p_satlst, crx.sat_lst_old, nsat, crx.sattbl); /****  set satellite table  ****/
            if (!std::getline(input, in_buffer)) {
                skip_to_next();
                goto SKIP;
            }
            data_format oldclk = crx.clk1;
            read_clock(in_buffer, &crx.clk1);
            for (int i = 0; i < nsat; i++) {
                crx.ntype = crx.ntype_record[i];
                if (!getdiff(i)) {
                    skip_to_next();
                    goto SKIP;
                }
            }

            /*************************************/
            /**** print the recovered line(s) ****/
            /*************************************/
            if (in_buffer[0] != '\0') {
                process_clock(oldclk);
            }

            if (crx.rinex_version == 2) {
                if (crx.clk1.order >= 0) {
                    out_ptr += sprintf(out_ptr, "%-68.68s", epoch_buffer.c_str());
                    print_clock(crx.clk1.u[crx.clk1.order], crx.clk1.l[crx.clk1.order], crx.shift_clk);
                } else {
                    out_ptr += sprintf(out_ptr, "%.68s\n", epoch_buffer.c_str());
                }
                for (int i = 1; nsat - 12*i > 0; i++) {
                    out_ptr += sprintf(out_ptr, "%32.s%.36s\n", " ", &epoch_buffer[32 + 36*i]);
                }
            } else {
                if (crx.clk1.order >= 0) {
                    out_ptr += sprintf(out_ptr, "%.41s", epoch_buffer.c_str());
                    print_clock(crx.clk1.u[crx.clk1.order], crx.clk1.l[crx.clk1.order], crx.shift_clk);
                } else {
                    out_ptr += sprintf(out_ptr, "%.41s\n", epoch_buffer.c_str());
                }
            }

            write_data(nsat, crx.p_satlst, crx.sattbl, crx.dflag);

            /****************************/
            /**** save current epoch ****/
            /****************************/
            crx.nsat1 = nsat;
            std::strncpy(crx.sat_lst_old, crx.p_satlst, (size_t)nsat*3);
            for (int i = 0; i < nsat; i++) {
                strncpy(crx.flag1[i], crx.flag[i], (size_t)crx.ntype_record[i]*2);
                for (int j = 0; j < crx.ntype_record[i]; j++) {
                    crx.dy0[i][j] = crx.dy1[i][j];
                }
            }
            return true;
        }
        return false;
    }

    void CRinexStreamBuf::repair(std::string& outstr, const std::string& instr) {
        bool finalize = false;
        for (size_t n = 0; n < instr.length(); n++) {
            // Fill unitialized space
            if (outstr.length() == n) {
                outstr.push_back(' ');
                finalize = true;
            }
            if (instr[n] == ' ') {
                continue;
            }
            if (instr[n] == '&') {
                outstr[n] = ' ';
            } else {
                outstr[n] = instr[n];
            }
        }
        if (finalize) {
            outstr.push_back('\0');
        }
    }

    void CRinexStreamBuf::skip_to_next() {
        // pointer to the space between year and month
        const char* message_fmt = crx.rinex_version == 2 ? "%29d%3d\n%-60sCOMMENT\n": ">%31d%3d\n%-60sCOMMENT\n";

        while (std::getline(input, in_buffer))
        {
            if (in_buffer[0] == crx.ep_top_from &&
                in_buffer.length() >= 29 &&
                in_buffer[crx.offset] == ' ' &&
                in_buffer[crx.offset+3] == ' ' &&
                in_buffer[crx.offset+6] == ' ' &&
                in_buffer[crx.offset+9] == ' ' &&
                in_buffer[crx.offset+12] == ' ' &&
                in_buffer[crx.offset+23] == ' ' &&
                in_buffer[crx.offset+24] == ' ' &&
                std::isdigit(in_buffer[crx.offset+25]))
            {
                break;
            }
        }

        out_ptr += sprintf(out_ptr, message_fmt, 4, 1, "  *** Some epochs are skipped by CRX2RNX ***");
    }

    void CRinexStreamBuf::set_sat_table(char* p_new, char* p_old, int nsat, int* sattbl) {
        /***********************************************************************/
        /*  - Read number of satellites (nsat)                                 */
        /*  - Compare the satellite list at the epoch (*p_new) and that at the */
        /*    previous epoch(*p_old), and make index (*sattbl) for the         */
        /*    corresponding order of the satellites.                           */
        /*    *sattbl is set to -1 for new satellites.                         */
        /***********************************************************************/
        /*** set # of data types for each satellite ***/
        if (crx.rinex_version == 2) {             /** for RINEX2 **/
            for (int i = 0; i < nsat; i++) {
                crx.ntype_record[i] = crx.ntype;
            }
        } else {                                /** for RINEX3 **/
            for (int i = 0; i < nsat; i++) {
                crx.ntype_record[i] = crx.ntype_gnss[(unsigned int)p_new[3*i]];  /*** # of data type for the GNSS system ***/
                if (crx.ntype_record[i] < 0) {
                    throw new std::runtime_error("GNSS type not defined in header");
                }
            }
        }
        for (int i = 0; i < nsat; i++) {
            sattbl[i] = -1;
            for (int j = 0; j < crx.nsat1; j++) {
                if (strncmp(&p_new[3*i], &p_old[3*j], 3) == 0) {
                    sattbl[i] = j;
                    break;
                }
            }
        }
    }

    void CRinexStreamBuf::write_data(int nsat, char* p_sat_lst, int* sattbl, char dflag[][MAXTYPE * 2]) {
/********************************************************************/
/*  Functions                                                       */
/*      (1) compose the original data from 3rd order difference     */
/*      (2) repair the flags                                        */
/*  sattbl : previous column on which the satellites are set        */
/*           new satellites are set to -1                           */
/*       u : upper X digits of the data                             */
/*       l : lower 5 digits of the data                             */
/*            ( y = u*100 + l/1000)                                 */
/*   date of previous epoch are set to dy0                           */
/********************************************************************/
        for (int i = 0; i < nsat; i++) {
            /**** set # of data types for the GNSS type    ****/
            /**** and write satellite ID in case of RINEX3 ****/
            /**** ---------------------------------------- ****/
            if (crx.rinex_version == 3) {
                crx.ntype = crx.ntype_record[i];
                strncpy(out_ptr, &p_sat_lst[i * 3], 3);
                out_ptr += 3;
            }
            /**** repair the data flags ****/
            /**** ----------------------****/
            if (sattbl[i] < 0) {       /* new satellite */
                if (crx.rinex_version == 3) {
                    *crx.flag[i] = '\0';
                } else {
                    sprintf(crx.flag[i], "%-*s", crx.ntype * 2, dflag[i]);
                }
            } else {
                strncpy(crx.flag[i], crx.flag1[sattbl[i]], (size_t)crx.ntype * 2);
            }
            for (size_t n = 0; n < strlen(dflag[i]); n++) {
                if (dflag[i][n] == ' ') {
                    continue;
                }
                if (dflag[i][n] == '&') { // TODO: do we even need this?
                    crx.flag[i][n] = ' ';
                } else {
                    crx.flag[i][n] = dflag[i][n];
                }
            }

            /**** recover the date, and output ****/
            /**** ---------------------------- ****/
            for (int j = 0; j < crx.ntype; j++) {
                if (crx.dy1[i][j].arc_order >= 0) {
                    if (crx.dy1[i][j].order < crx.dy1[i][j].arc_order) {
                        (crx.dy1[i][j].order)++;
                        for (int k = 0; k < crx.dy1[i][j].order; k++) {
                            crx.dy1[i][j].u[k + 1] = crx.dy1[i][j].u[k] + crx.dy0[sattbl[i]][j].u[k];
                            crx.dy1[i][j].l[k + 1] = crx.dy1[i][j].l[k] + crx.dy0[sattbl[i]][j].l[k];
                            crx.dy1[i][j].u[k + 1] += crx.dy1[i][j].l[k + 1] / 100000;  /*** to avoid overflow of dy1.l ***/
                            crx.dy1[i][j].l[k + 1] %= 100000;
                        }
                    } else {
                        for (int k = 0; k < crx.dy1[i][j].order; k++) {
                            crx.dy1[i][j].u[k + 1] = crx.dy1[i][j].u[k] + crx.dy0[sattbl[i]][j].u[k + 1];
                            crx.dy1[i][j].l[k + 1] = crx.dy1[i][j].l[k] + crx.dy0[sattbl[i]][j].l[k + 1];
                            crx.dy1[i][j].u[k + 1] += crx.dy1[i][j].l[k + 1] / 100000;
                            crx.dy1[i][j].l[k + 1] %= 100000;
                        }
                    }
                    /* Signs of dy1[i][j].u and dy1[i][j].l can be different at this stage */
                    /*   and will be adjusted before outputting                 */
                    putfield(&crx.dy1[i][j], &crx.flag[i][j * 2]);
                } else {
                    if (crx.version == 1) {                       /*** CRINEX 1 assumes that flags are always ***/
                        out_ptr += sprintf(out_ptr,
                                           "                "); /*** blank if data field is blank           ***/
                        crx.flag[i][j * 2] = crx.flag[i][j * 2 + 1] = ' ';
                    } else {                                            /*** CRINEX 3 evaluate flags independently **/
                        out_ptr += sprintf(out_ptr, "              %c%c", crx.flag[i][j * 2], crx.flag[i][j * 2 + 1]);
                    }
                }
                if ((j + 1) == crx.ntype || (crx.rinex_version == 2 && (j + 1) % 5 == 0)) {
                    while (*--out_ptr == ' ') {
                    }
                    out_ptr++;  /*** cut spaces ***/
                    *out_ptr++ = '\n';
                }
            }
        }
    }

    void CRinexStreamBuf::putfield(CRinexStreamBuf::data_format* y, char* flag) {
        int i = y->order;

        if (y->u[i] < 0 && y->l[i] > 0) {
            y->u[i]++;
            y->l[i] -= 100000;
        } else if (y->u[i] > 0 && y->l[i] < 0) {
            y->u[i]--;
            y->l[i] += 100000;
        }
        /* The signs of y->u and y->l are the same (or zero) at this stage */

        if (y->u[i] != 0) {                                    /* ex) 123.456  -123.456 */
            out_ptr += sprintf(out_ptr, "%8ld %5.5ld%c%c", y->u[i], labs(y->l[i]), flag[0], flag[1]);
            out_ptr[-8] = out_ptr[-7];
            out_ptr[-7] = out_ptr[-6];
            if (y->u[i] > 99999999 || y->u[i] < -9999999) {
                throw std::runtime_error("Data record out of range");
            }
        } else {
            out_ptr += sprintf(out_ptr, "         %5.5ld%c%c", labs(y->l[i]), flag[0], flag[1]);
            if (out_ptr[-7] != '0') {                        /* ex)  12.345    -2.345 */
                out_ptr[-8] = out_ptr[-7];
                out_ptr[-7] = out_ptr[-6];
                if (y->l[i] < 0) {
                    out_ptr[-9] = '-';
                }
            } else if (out_ptr[-6] != '0') {                  /* ex)   1.234    -1.234 */
                out_ptr[-7] = out_ptr[-6];
                out_ptr[-8] = (y->l[i] < 0) ? '-' : ' ';
            } else {                                          /* ex)    .123     -.123 */
                out_ptr[-7] = (y->l[i] < 0) ? '-' : ' ';
            }
        }
        out_ptr[-6] = '.';
    }

    bool CRinexStreamBuf::put_event_data(std::string& instr, char* p_event) {
/***********************************************************************/
/*  - Put event data for one event.                                    */
/*  - This function is called when the event flag > 1.                 */
/***********************************************************************/
        do {
            instr[0] = crx.ep_top_to;
            out_ptr += sprintf(out_ptr, "%s\n", instr.c_str());
            if (instr.length() > 29) {
                int n = atoi((p_event + 1));
                for (int i = 0; i < n; i++) {
                    if (!std::getline(input, instr)) {
                        return false;
                    }
                    out_ptr += sprintf(out_ptr, "%s\n", instr.c_str());
                    if (instr.compare(60, 19, "# / TYPES OF OBSERV") == 0 && instr[5] != ' ') {
                        crx.ntype = atoi(&instr[0]);                                        /** for RINEX2 **/
                    } else if (instr.compare(60, 19, "SYS / # / OBS TYPES") == 0) { /** for RINEX3 **/
                        if (instr[0] != ' ') {
                            crx.ntype_gnss[(unsigned int) instr[0]] = atoi(&instr[3]);
                        }
                        if (crx.ntype_gnss[(unsigned int) instr[0]] > MAXTYPE) {
                            throw new std::runtime_error("GNSS type ID exceeds maximum");
                        }
                    }
                }
            }

            do {
                if (!std::getline(input, instr)) {
                    return false;
                }
            } while (crx.version == 3 && instr[0] == '&');

            if (instr[0] != crx.ep_top_from || instr.length() < 29 || !isdigit(*p_event)) {
                throw new std::runtime_error("Uninitialized epoch");
            }
        } while (*p_event != '0' && *p_event != '1');
        return true;
    }

    bool CRinexStreamBuf::getdiff(int sat) {
        int oldsat = crx.sattbl[sat];
        std::string line;
        if (!std::getline(input, line)) {
            return false;
        }
        line.resize(in_buffer.capacity(), '\0');

        /******************************************/
        /****  separate the fields with '\0'   ****/
        /******************************************/
        size_t p = 0;
        for (size_t nfields = 0; nfields < crx.ntype; p++) {
            if (line[p] == '\0') {
                nfields++;
                // make empty flag area
                line[p + 1] = '\0';
            }
            else if (line[p] == ' ') {
                // split fields
                nfields++;
                line[p] = '\0';
            }
        }
        strcpy(crx.dflag[sat], &line[p]);

        /************************************/
        /*     read the differenced data    */
        /************************************/
        size_t fbegin = 0, fend = 0;
        for (int j = 0; j < crx.ntype; j++) {
            if (line[fbegin] == '\0') {
                crx.dy1[sat][j].arc_order = -1;      /**** arc_order < 0 means that the field is blank ****/
                crx.dy1[sat][j].order = -1;
                fbegin++;
            } else {
                if (line[fbegin + 1] == '&') {     /**** arc initialization ****/
                    crx.dy1[sat][j].order = -1;
                    crx.dy1[sat][j].arc_order = atoi(&line[fbegin]);
                    fbegin += 2;
                    if (crx.dy1[sat][j].arc_order > MAX_DIFF_ORDER) {
                        throw new std::runtime_error("Difference order exceeds the maximum");
                    }
                } else if (crx.sattbl[sat] < 0) {
                    throw new std::runtime_error("Uninitialized new satellite arc");
                }
                else if (crx.dy0[oldsat][j].arc_order < 0) {
                    throw new std::runtime_error("Uninitialized data sequence");
                } else {
                    crx.dy1[sat][j].order = crx.dy0[oldsat][j].order;
                    crx.dy1[sat][j].arc_order = crx.dy0[oldsat][j].arc_order;
                }
                fend = line.find('\0', fbegin);
                size_t length = fend - fbegin;
                if (line[fbegin] == '-') {
                    length--;
                }
                if (length < 6) {
                    crx.dy1[sat][j].u[0] = 0;
                    crx.dy1[sat][j].l[0] = atol(&line[fbegin]);
                } else {
                    size_t lowstart = fend - 5;
                    crx.dy1[sat][j].l[0] = atol(&line[lowstart]);
                    line[lowstart] = '\0';
                    crx.dy1[sat][j].u[0] = atol(&line[fbegin]);
                    if (crx.dy1[sat][j].u[0] < 0) {
                        crx.dy1[sat][j].l[0] = -crx.dy1[sat][j].l[0];
                    }
                }
                fbegin = fend + 1;
            }
        }
        return true;
    }

    void CRinexStreamBuf::read_clock(std::string& instr, data_format* clk) {
        if (instr.length() == 0) {
            clk->order = -1;
        } else {
            size_t clock_start = 0;
            if (instr[1] == '&') {        /**** for the case of arc initialization ****/
                sscanf(&instr[0], "%d&", &clk->arc_order);
                if (clk->arc_order > MAX_DIFF_ORDER) {
                    throw new std::runtime_error("Difference order exceeds the maximum");
                }
                clk->order = -1;
                clock_start += 2;
            }
            size_t number_start = clock_start;
            if (instr[clock_start] == '-') {
                number_start++;
            }
            size_t len = instr.length();
            if ((len - number_start) < 9) {                /** s-p1 == strlen(p1) ***/
                clk->u[0] = 0;
                clk->l[0] = atol(&instr[clock_start]);
            } else {
                len -= 8;
                clk->l[0] = atol(&instr[len - 8]);
                instr[len - 8] = '\0';
                clk->u[0] = atol(&instr[clock_start]);
                if (clk->u[0] < 0) {
                    clk->l[0] = -clk->l[0];
                }
            }
        }
    }

    void CRinexStreamBuf::print_clock(long yu, long yl, int shift_clk) {
        char tmp[8], * p_tmp, * p;
        int n, sgn;

        if (yu < 0 && yl > 0) {
            yu++;
            yl -= 100000000;
        } else if (yu > 0 && yl < 0) {
            yu--;
            yl += 100000000;
        }
        /* The signs of yu and yl are the same (or zero) at this stage */

        /** add ond more digit to handle '-0'(RINEX2) or '-0000'(RINEX3) **/
        sgn = (yl < 0) ? -1 : 1;
        n = sprintf(tmp, "%.*ld", shift_clk + 1, yu * 10 + sgn); /** AT LEAST fractional parts are filled with 0 **/
        n--;                           /** n: number of digits excluding the additional digit **/
        p_tmp = &tmp[n];
        *p_tmp = '\0';
        p_tmp -= shift_clk;       /** pointer to the top of last "shift_clk" digits **/
        out_ptr += sprintf(out_ptr, "  .%s", p_tmp);  /** print last "shift_clk" digits.  **/
        if (n > shift_clk) {
            p_tmp--;
            p = out_ptr - shift_clk - 2;
            *p = *p_tmp;

            if (n > shift_clk + 1) {
                *(p - 1) = *(p_tmp - 1);
                if (n > shift_clk + 2) {
                    throw new std::runtime_error("Clock offset out of range");
                }
            }
        }

        out_ptr += sprintf(out_ptr, "%8.8ld\n", labs(yl));
    }

    void CRinexStreamBuf::process_clock(const data_format& oldclk) {
        /****************************************/
        /**** recover the clock offset value ****/
        /****************************************/
        if (crx.clk1.order < crx.clk1.arc_order)
        {
            crx.clk1.order++;
        }
        for (int i = 0, j = 1; i < crx.clk1.order; i++, j++) {
            int k = crx.clk1.order < crx.clk1.arc_order ? i: j;
            crx.clk1.u[j] = crx.clk1.u[i] + oldclk.u[k];
            crx.clk1.l[j] = crx.clk1.l[i] + oldclk.l[k];
            crx.clk1.u[j] += crx.clk1.l[j] / 100000000;
            crx.clk1.l[j] %= 100000000;
        }
        /* Signs of py1->u and py1->l can be different at this stage */
        /*   and will be adjustied before outputting */
    }

    void CRinexStreamBuf::trim() {
        // Remove trailing whitespace in current output line
        while (*out_ptr == ' ' && out_ptr >= &out_buffer[0]) {
            *out_ptr = '\0';
            out_ptr--;
        }
    }
}

