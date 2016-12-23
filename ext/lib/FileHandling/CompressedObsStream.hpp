#pragma once
#include "Rinex3ObsStream.hpp"
#include "Rinex3ObsData.hpp"
#include "DataStructures.hpp"

#include <queue>
#include <stdint.h>

namespace gpstk {
namespace experimental {

    class CompressedObsStream :
        public Rinex3ObsStream
    {
    private:
        class ObsState {
        public:
            static const size_t MAX_ORDER = 10;
            int8_t order = -1, arcOrder = 0;
            int8_t lli = -1, ssi = -1;
            int64_t dy[MAX_ORDER];

            void update(const std::string& diff);
            int64_t value() { return dy[0]; }
            bool present() { return order >= 0; }
            void reset();
        };

        struct Config {
            char ep_from, ep_to;
            size_t offset, event_offset, nsat_offset,
                satlist_offset;
        };

        struct SatState {
            std::vector<ObsState> obs;
            uint32_t lastEpoch = 0;
            bool present = false;
        };

        typedef std::map<RinexSatID, SatState> State;

        // Compact RINEX version (so far only 1 or 3)
        int version;

        // Epoch counter (helps to find missing observations)
        int epochID = 0;

        // Current epoch
        std::string epoch;

        // Receiver clock offset
        ObsState rcv_clock;

        // Auxiliary header records waiting for processing
        std::queue<std::string> auxRecords;

        // Storage for prematurely-read epoch lines
        std::string epoch_backlog;

        // Reader state
        State state;

        // Compact RINEX 1.0 / 3.0 differences in form of flags
        Config cfg;

        void readPreamble();

        void readData();

        void setFlags(SatState obslist, const std::string& diff);

        CommonTime parseTime() const;

    public:
        CompressedObsStream(const char* fn, std::ios_base::openmode mode = std::ios::in);
        CompressedObsStream(const std::string& fn, std::ios_base::openmode mode = std::ios::in);

        friend class CompressedObsData;
    };

    // Wraps CompressedObsStream's state in Rinex3ObsData interface
    class CompressedObsData :
        public Rinex3ObsData {
    protected:
        virtual void reallyGetRecord(FFStream& ffs);
    };

    // Helper method for the Processing Framework
    CompressedObsStream& operator >> (CompressedObsStream& strm, gnssRinex& f);
}
}