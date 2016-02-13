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
 * @file Rinex3ObsStream.cpp
 * File stream for RINEX 3 observation file data.
 */

#include "Rinex3ObsStream.hpp"
#include "ZStreamBuf.hpp"
#include "CRinexStreamBuf.hpp"

namespace gpstk
{
   Rinex3ObsStream ::
   Rinex3ObsStream()
   {
      init(false);
   }


   Rinex3ObsStream ::
   Rinex3ObsStream( const char* fn,
                    std::ios::openmode mode )
         : FFTextStream(fn, mode)
   {
      init((mode & std::ios::in) == std::ios::in);
   }


   Rinex3ObsStream ::
   Rinex3ObsStream( const std::string fn,
                    std::ios::openmode mode )
         : FFTextStream(fn.c_str(), mode)
   {
      init((mode & std::ios::in) == std::ios::in);
   }


   Rinex3ObsStream ::
   ~Rinex3ObsStream()
   {
      // Cleanup allocated streambuffer filters
      for (std::vector<std::streambuf*>::iterator it = filters.begin(); it != filters.end(); it++)
      {
         delete *it;
      }
   }


   void Rinex3ObsStream ::
   open( const char* fn,
         std::ios::openmode mode )
   {
      FFTextStream::open(fn, mode);
      init((mode & std::ios::in) == std::ios::in);
   }


   void Rinex3ObsStream ::
   init(bool tryDecode)
   {
      headerRead = false;
      header = Rinex3ObsHeader();
      timesystem = TimeSystem::GPS;

      // If reading a file, detour from parent baseStream to provide unpacking/decompressing capabilities
      if (tryDecode)
      {
         // Read a few bytes to check the header for LZW signature
         char magic[2];
         if (!read(magic, 2))
            return;
         // Roll back
         seekg(-2, std::ios::cur);
         if (magic[0] == '\037' && magic[1] == '\235')
         {
            // Create a streambuffer to decompress LZW on-the-fly
            std::streambuf* filter = new ZStreamBuf(rdbuf());
            rdbuf(filter);
            filters.push_back(filter);
         }

         char crxHeader[6];
         // Skip to "CRINEX VERS   / TYPE" comment
         seekg(60, std::ios::cur);
         if (!read(crxHeader, 6))
            return;
         // Roll back
         seekg(-66, std::ios::cur);
         // Check if we have a correct header
         if (std::string("CRINEX").compare(0, 6, crxHeader) == 0)
         {
            // TODO: check if we're reading a Hatanaka-compressed file
            std::streambuf* filter = new CRinexStreamBuf(rdbuf());
            rdbuf(filter);
            filters.push_back(filter);
         }
      }
   }


   void Rinex3ObsStream ::
   open( const std::string& fn,
         std::ios::openmode mode )
   {
      open(fn.c_str(), mode);
   }


   bool Rinex3ObsStream ::
   isRinex3ObsStream(std::istream& i)
   {
      try
      {
         (void)dynamic_cast<Rinex3ObsStream&>(i);
      }
      catch(...)
      {
         return false;
      }

      return true;
   }

} // namespace gpstk
