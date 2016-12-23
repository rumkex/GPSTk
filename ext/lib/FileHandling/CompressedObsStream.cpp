#include "CompressedObsStream.hpp"

#include <string>

#include "StringUtils.hpp"

using namespace gpstk::StringUtils;

namespace gpstk {
namespace experimental {
    CompressedObsStream::CompressedObsStream(const char* fn, std::ios_base::openmode mode) :
        Rinex3ObsStream(fn, mode), epochID(0)
    {
        readPreamble();
    }

    CompressedObsStream::CompressedObsStream(const std::string& fn, std::ios_base::openmode mode):
        Rinex3ObsStream(fn, mode), epochID(0)
    {
        readPreamble();
    }

    void CompressedObsStream::readPreamble()
    {
        std::string line;
        // Read main compressed RINEX header
        formattedGetLine(line, true);
        // It should end in "CRINEX VERS   / TYPE"
        if (line.length() < 80 || line.compare(60, 20, "CRINEX VERS   / TYPE") != 0) {
            FFStreamError e("Not a Hatanaka-compressed RINEX file");
            GPSTK_THROW(e);
        }
        version = asInt(line.substr(0, 20));

        // Setup version-specific parameters
        if (version == 1) {
            cfg.ep_from = '&';
            cfg.ep_to = ' ';
            cfg.event_offset = 28;
            cfg.nsat_offset = 29;
            cfg.satlist_offset = 32;
            cfg.offset = 3;
        }
        else {
            cfg.ep_from = '>';
            cfg.ep_to = '>';
            cfg.event_offset = 31;
            cfg.nsat_offset = 32;
            cfg.satlist_offset = 41;
            cfg.offset = 6;
        }

        // Read CRINEX PROG / DATE header line
        formattedGetLine(line, true);
        // TODO: actually parse it, if necessary
    }

    void CompressedObsStream::readData()
    {
        std::string epoch_update;

        if (epoch.length() > cfg.event_offset &&
            epoch[cfg.event_offset] != '0' && epoch[cfg.event_offset] != '1')
        {
            // We have backlog from parsing event data
            epoch_update = epoch_backlog;
            epoch_backlog.clear();
        }
        else
        {
            // Read epoch line
            formattedGetLine(epoch_update);
        }

        while (version == 3 && epoch_update[0] == '&')
        {
            // Skip escape lines of CRINEX version 3
            formattedGetLine(epoch_update);
        }

        // Check if this is a newly-initialized epoch
        if (epoch_update[0] == cfg.ep_from)
        {
            // Start from scratch
            epoch_update[0] = cfg.ep_to;
            epoch.clear();
            state.clear();
        }

        // Repair delta encoding
        if (epoch.length() < epoch_update.length())
        {
            epoch.resize(epoch_update.size(), ' ');
        }
        for (size_t i = 0; i < epoch_update.size(); i++)
        {
            if (epoch_update[i] == '&')
                epoch[i] = ' ';
            else if (epoch_update[i] != ' ')
                epoch[i] = epoch_update[i];
        }

        // Check the epoch flag
        if (epoch[cfg.event_offset] != '0' && epoch[cfg.event_offset] != '1')
        {
            // Event epoch, reset the state and read aux records
            while (true) {
                formattedGetLine(epoch_backlog);
                // Find next epoch line
                if (epoch_backlog[0] == cfg.ep_from &&
                    epoch_backlog.length() >= 29 &&
                    epoch_backlog[cfg.offset] == ' ' &&
                    epoch_backlog[cfg.offset + 3] == ' ' &&
                    epoch_backlog[cfg.offset + 6] == ' ' &&
                    epoch_backlog[cfg.offset + 9] == ' ' &&
                    epoch_backlog[cfg.offset + 12] == ' ' &&
                    epoch_backlog[cfg.offset + 23] == ' ' &&
                    epoch_backlog[cfg.offset + 24] == ' ' &&
                    std::isdigit(epoch_backlog[cfg.offset + 25]))
                {
                    break;
                }
                // Store aux records that will be later parsed by ObsData
                if (epoch[cfg.event_offset] >= '2' || epoch[cfg.event_offset] <= '5')
                {
                    auxRecords.push(epoch_backlog);
                }
            }
            state.clear();
            return;
        }

        // TODO: Read the epoch date
        epochID++;

        // Read clock offset, if present
        std::string clock_diff;
        formattedGetLine(clock_diff);
        if (clock_diff.length() > 0)
            rcv_clock.update(clock_diff);

        // Read the differences for each listed satellite
        size_t nsat = asInt(epoch.substr(cfg.nsat_offset, 3));
        for (size_t n = 0; n < nsat; n++)
        {
            RinexSatID sat(epoch.substr(cfg.satlist_offset + 3 * n, 3));
            size_t nobs = version == 3 ?
                header.mapObsTypes[asString(sat.systemChar())].size() :
                header.mapObsTypes["G"].size();
            // Prepare the storage if necessary
            if (state.find(sat) == state.end())
            {
                state[sat] = SatState();
                state[sat].obs = std::vector<ObsState>(nobs);
            }

            // Bump the epoch counter, so we know the satellite data is up-to-date
            state[sat].lastEpoch = epochID;
            state[sat].present = true;

            std::string satdiffs;
            formattedGetLine(satdiffs);

            // Tokenize the string, and parse the fields
            size_t f_begin = 0, f_end = 0;
            for (size_t i = 0; i < nobs && f_begin < satdiffs.length(); i++)
            {
                if (satdiffs[f_begin] == ' ')
                {
                    f_begin++;
                    state[sat].obs[i].reset();
                    continue;
                }
                size_t f_end = satdiffs.find(' ', f_begin);
                if (f_end == std::string::npos)
                {
                    f_end = satdiffs.length();
                }
                size_t f_length = f_end - f_begin;
                state[sat].obs[i].update(satdiffs.substr(f_begin, f_length));
                f_begin = f_end + 1;
            }
            if (f_begin < satdiffs.length())
                // Have to set the flags
                setFlags(state[sat], satdiffs.substr(f_begin));
        }

        // Invalidate all satellites not up-to-date
        for (State::iterator kv = state.begin(); kv != state.end(); kv++)
        {
            if (kv->second.lastEpoch != epochID)
            {
                state[kv->first].present = false;
                for (std::vector<ObsState>::iterator it = state[kv->first].obs.begin();
                     it != state[kv->first].obs.end(); it++)
                    it->reset();
            }
        }
    }

    void CompressedObsStream::setFlags(SatState obslist, const std::string & diff)
    {
        // TODO: Actually set the flags
    }

    void CompressedObsStream::ObsState::update(const std::string & diff)
    {
        if (diff.length() > 2 && diff[1] == '&')
        {
            // Initialize arc
            order = 0;
            arcOrder = diff[0] - '0';
            dy[0] = std::strtoll(&diff.c_str()[2], NULL, 10);
            for (size_t i = 1; i < arcOrder; i++)
            {
                dy[i] = 0;
            }
        } else if (diff.length() == 0) {
            order = -1;
        } else {
            if (order == -1)
            {
                FFStreamError e("Uninitialized satellite arc");
                GPSTK_THROW(e);
            }
            if (order < arcOrder)
                order++;
            dy[order] = strtoll(diff.c_str(), NULL, 10);
            for (size_t i = order-1; i >= 0; i--)
            {
                dy[i] += dy[i + 1];
                if (i == 0) break;
            }
        }
    }

    void CompressedObsStream::ObsState::reset()
    {
        std::fill(dy, &dy[MAX_ORDER], 0LL);
        lli = -1;
        ssi = -1;
        order = -1;
        arcOrder = 0;
    }

    void CompressedObsData::reallyGetRecord(FFStream & ffs)
        throw (std::exception, gpstk::FFStreamError, gpstk::StringUtils::StringException)
    {
        CompressedObsStream& strm = dynamic_cast<CompressedObsStream&>(ffs);

        // If the header hasn't been read, read it.
        if (!strm.headerRead) strm >> strm.header;

        // TODO: clear this ObsData

        try
        {
            strm.readData();
            epochFlag = strm.epoch[strm.cfg.event_offset] - '0';
            if (epochFlag < 0 || epochFlag > 6)
            {
                FFStreamError e("Invalid epoch flag: " + asString(epochFlag));
                GPSTK_THROW(e);
            }

            time = strm.parseTime();

            if (strm.rcv_clock.present())
                clockOffset = strm.rcv_clock.value() / 1000.0; // TODO: find out the coefficient
            else
                clockOffset = 0.0;

            if (epochFlag == 0 || epochFlag == 1 || epochFlag == 6)
            {
                for (CompressedObsStream::State::iterator kv = strm.state.begin(); kv != strm.state.end(); kv++)
                {
                    for (std::vector<CompressedObsStream::ObsState>::iterator it = kv->second.obs.begin();
                         it != kv->second.obs.end(); it++)
                    {
                        RinexDatum d;
                        d.lli = it->lli;
                        d.lliBlank = it->lli == -1;
                        d.ssi = it->ssi;
                        d.ssiBlank = it->ssi == -1;
                        d.data = it->present() ? it->value() / 1000.0 : 0.0;
                        d.dataBlank = !it->present();
                        obs[kv->first].push_back(d);
                    }
                }
            }
            else
            {
                // Write aux header
                while (strm.auxRecords.size() > 0)
                {
                    auxHeader.parseHeaderRecord(strm.auxRecords.front());
                    strm.auxRecords.pop();
                }
            }
        }
        catch (Exception& e)
        {
            GPSTK_RETHROW(e);
        }
    }


    CommonTime CompressedObsStream::parseTime() const
    {
        try
        {
            size_t offset;
            int yy = 0;
            if (epoch[0] == '>')
            {
                // RINEX 3 - "> YYYY MM DD..."
                offset = 4;
            }
            else
            {
                // RINEX 2 - " YY MM DD..."
                yy = ((CivilTime)(header.firstObs)).year / 100;
                yy *= 100;
                offset = 1;
            }
            // check if the spaces are in the right place - an easy
            // way to check if there's corruption in the file
            if ((epoch[offset+2] != ' ') || (epoch[offset+5] != ' ') || (epoch[offset+8] != ' ') ||
                (epoch[offset+11] != ' ') || (epoch[offset+14] != ' ') ||
                (epoch[offset+25] != ' ') || (epoch[offset+26] != ' '))
            {
                FFStreamError e("Invalid time format");
                GPSTK_THROW(e);
            }

            // if there's no time, just return a bad time
            if (epoch.substr(2, 27) == std::string(27, ' '))
                return CommonTime::BEGINNING_OF_TIME;

            int year, month, day, hour, min;
            double sec;

            year = offset == 1?
                yy + asInt(epoch.substr(1, 2)):
                asInt(epoch.substr(2, 4));
            month = asInt(epoch.substr(offset+3, 2));
            day = asInt(epoch.substr(offset+6, 2));
            hour = asInt(epoch.substr(offset+9, 2));
            min = asInt(epoch.substr(offset+12, 2));
            sec = asDouble(epoch.substr(offset+15, 11));

            // Real Rinex has epochs 'yy mm dd hr 59 60.0' surprisingly often.
            double ds = 0;
            if (sec >= 60.)
            {
                ds = sec;
                sec = 0.0;
            }

            CommonTime rv = CivilTime(year, month, day, hour, min, sec).convertToCommonTime();
            if (ds != 0) rv += ds;

            rv.setTimeSystem(timesystem);

            return rv;
        }
        // string exceptions for substr are caught here
        catch (std::exception &e)
        {
            FFStreamError err("std::exception: " + std::string(e.what()));
            GPSTK_THROW(err);
        }
        catch (gpstk::Exception& e)
        {
            FFStreamError err(e);
            GPSTK_THROW(err);
        }
    }

    CompressedObsStream& operator >> (CompressedObsStream& strm, gnssRinex& f)
    {
        // If the header hasn't been read, read it...
        if (!strm.headerRead) strm >> strm.header;

        // Clear out this object
        Rinex3ObsHeader& roh = strm.header;

        CompressedObsData rod;
        strm >> rod;

        // Fill data
        f.header.source.type = SatIDsystem2SourceIDtype(roh.fileSysSat);
        f.header.source.sourceName = roh.markerName;
        f.header.antennaType = roh.antType;
        f.header.antennaPosition = roh.antennaPosition;
        f.header.epochFlag = rod.epochFlag;
        f.header.epoch = rod.time;

        f.body = satTypeValueMapFromRinex3ObsData(roh, rod);

        return strm;
    }
}
}