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

// @file PositionSatStore.cpp
// Store a tabular list of ephemeris data (position, option velocity & acceleration)
// for several satellites, and compute values at any timetag from this table.
// Inherits TabularSatStore.

#include "PositionSatStore.hpp"
#include "MiscMath.hpp"
#include <vector>

using namespace std;

namespace gpstk
{
   /** @addtogroup ephemstore */
   //@{

   // Output stream operator is used by dump() in TabularSatStore
   ostream& operator<<(ostream& os, const PositionRecord& rec) throw()
   {
      os << "Pos" << fixed << setprecision(6)
         << " " << setw(13) << rec.Pos[0]
         << " " << setw(13) << rec.Pos[1]
         << " " << setw(13) << rec.Pos[2]
         << " sigP" << scientific << setprecision(2)
         << " " << setw(9) << rec.sigPos[0]
         << " " << setw(9) << rec.sigPos[1]
         << " " << setw(9) << rec.sigPos[2]
         << " Vel" << fixed << setprecision(6)
         << " " << setw(13) << rec.Vel[0]
         << " " << setw(13) << rec.Vel[1]
         << " " << setw(13) << rec.Vel[2]
         << " sigV" << scientific << setprecision(2)
         << " " << setw(9) << rec.sigVel[0]
         << " " << setw(9) << rec.sigVel[1]
         << " " << setw(9) << rec.sigVel[2]
         //<< " Acc" << fixed << setprecision(6)
         //<< " " << setw(13) << rec.Acc[0]
         //<< " " << setw(13) << rec.Acc[1]
         //<< " " << setw(13) << rec.Acc[2]
         //<< " sigA" << scientific << setprecision(2)
         //<< " " << setw(9) << rec.sigAcc[0]
         //<< " " << setw(9) << rec.sigAcc[1]
         //<< " " << setw(9) << rec.sigAcc[2]
         ;
      return os;
   }

   void PositionSatStore::updateInterval(const SatID& sat, const CommonTime& newTtag)
   {
       if (tables[sat].size() == 0)
           ;
       else
       {
           // Estimate interval from closest point

           auto it = tables[sat].lower_bound(newTtag);
           CommonTime ref;
           if (it != tables[sat].end())
               ref = it->first;
           else
               ref = tables[sat].rbegin()->first;            

           long day, msod, day1, msod1;
           double fsod;
           newTtag.getInternal(day1, msod1, fsod);
           ref.getInternal(day, msod, fsod);
           unsigned long long ms = std::abs((day1 - day) * MS_PER_DAY + (msod1 - msod));
           if (ms > dataInterval)
           {
              std::swap(ms, dataInterval);
           }

           // dataInterval = GCD(ms, dataInterval)
           while (ms != 0) {
               unsigned long long r = dataInterval % ms;
               dataInterval = ms;
               ms = r;
           }
       }
   }

   void PositionSatStore::updateWeights()
   {
       bWeights.clear();
       bWeights.resize(interpOrder, 1.0);
       for (std::size_t i = 0; i < interpOrder; i++)
       {
           for (std::size_t j = 0; j < interpOrder; j++)
           {
               if (i == j) continue;
               // Since a fixed interval is assumed, 
               // it's irrelevant (weights can be scaled by any number)
               bWeights[i] = bWeights[i] * (double(i) - double(j));
           }
           bWeights[i] = 1. / bWeights[i];
       }
   }


   void PositionSatStore::barycentricInterp(DataTableIterator it1, DataTableIterator it2,
       double t0, PositionRecord& rec) const
   {
       for (std::size_t k = 0; k < 3; k++) {
           double A(0), B(0), C(0), D(0);
           std::size_t i;
           auto it = it1;
           for (it = it1, i = 0; it != it2; it++, i++)
           {
               double y = it->second.Pos[k];
               double t(it->first - it1->first);
               double w = bWeights[i] / (t0 - t);
               A += y * w;
               B += w;
               C += y * w / (t0 - t);
               D += w / (t0 - t);
           }
           rec.Pos[k] = A / B;
           rec.Vel[k] = (D * rec.Pos[k] - C) / B;
       }
   }
   
   // Return value for the given satellite at the given time (usually via
   // interpolation of the data table). This interface from TabularSatStore.
   // @param[in] sat the SatID of the satellite of interest
   // @param[in] ttag the time (CommonTime) of interest
   // @return object of type PositionRecord containing the data value(s).
   // @throw InvalidRequest if data value cannot be computed, for example because
   //  a) the time t does not lie within the time limits of the data table
   //  b) checkDataGap is true and there is a data gap
   //  c) checkInterval is true and the interval is larger than maxInterval
   PositionRecord PositionSatStore::getValue(const SatID& sat, const CommonTime& ttag)
      const throw(InvalidRequest)
   {
      try {
         bool isExact;
         int i;
         PositionRecord rec;
         rec.sigAcc = rec.Acc = Triple(0, 0, 0);        // default
         DataTableIterator it1, it2, kt;        // cf. TabularSatStore.hpp

         isExact = getTableInterval(sat, ttag, Nhalf, it1, it2, haveVelocity);
         if(isExact && haveVelocity) {
            rec = it1->second;
            return rec;
         }

         size_t n, Nlow(Nhalf - 1), Nhi(Nhalf), Nmatch(Nhalf);
         DataTableIterator itmatch;
         CommonTime ttag0(it1->first);         
         int nIntervals = int(it2->first - it1->first) / (dataInterval / 1000);
         double dt(ttag - ttag0), err;           // dt in seconds
         
         for (kt = it1, n = 0; kt != it2; kt++, n++) {
             // find index matching ttag
             if (isExact && ABS(kt->first - ttag) < 1.e-8)
             {
                 itmatch = kt;
                 Nmatch = n;
             }
         };

         if (isExact && Nmatch == (int)(Nhalf - 1)) { Nlow++; Nhi++; }

         // Special case: use fast barycentric interpolation without data copying if applicable
         if (!haveVelocity && dataInterval > 1 && nIntervals == interpOrder - 1)
         {
             barycentricInterp(it1, it2, dt, rec);
             for (i = 0; i < 3; i++) {
                 rec.Vel[i] *= 10000.;         // km/sec -> dm/sec

                 if (isExact) {
                     rec.sigPos[i] = itmatch->second.sigPos[i];
                 }
                 else {
                     // TODO: How to fill this without scavenging through the iterator?
                     // rec.sigPos[i] = RSS(sigP[i][Nhi], sigP[i][Nlow]);
                     rec.sigPos[i] = 2.0;
                 }
                 // TD
                 rec.sigVel[i] = 0.0;
             }
             return rec;
         }

         // pull data out of the data table
         vector<double> times(interpOrder);

         vector<double> P[3] = { vector<double>(interpOrder), vector<double>(interpOrder), vector<double>(interpOrder) };
         vector<double> V[3] = { vector<double>(interpOrder), vector<double>(interpOrder), vector<double>(interpOrder) };
         vector<double> A[3] = { vector<double>(interpOrder), vector<double>(interpOrder), vector<double>(interpOrder) };
         vector<double> sigP[3] = { vector<double>(interpOrder), vector<double>(interpOrder), vector<double>(interpOrder) };
         vector<double> sigV[3] = { vector<double>(interpOrder), vector<double>(interpOrder), vector<double>(interpOrder) };
         vector<double> sigA[3] = { vector<double>(interpOrder), vector<double>(interpOrder), vector<double>(interpOrder) };

         for (kt = it1, n = 0; kt != it2; kt++, n++) {
            times[n] = kt->first - ttag0;          // sec
            for(i=0; i<3; i++) {
               P[i][n] = kt->second.Pos[i];
               V[i][n] = kt->second.Vel[i];
               A[i][n] = kt->second.Acc[i];
               sigP[i][n] = kt->second.sigPos[i];
               sigV[i][n] = kt->second.sigVel[i];
               sigA[i][n] = kt->second.sigAcc[i];
            }
         };
         
         // Lagrange interpolation
         if(haveVelocity) {
            for(i=0; i<3; i++) {
               // interpolate the positions
               rec.Pos[i] = LagrangeInterpolation(times,P[i],dt,err);
               if(haveAcceleration) {
                  // interpolate velocities and acclerations
                  rec.Vel[i] = LagrangeInterpolation(times,V[i],dt,err);
                  rec.Acc[i] = LagrangeInterpolation(times,A[i],dt,err);
               }
               else {
                  // interpolate velocities(dm/s) to get V and A
                  LagrangeInterpolation(times,V[i],dt,rec.Vel[i],rec.Acc[i]);
                  rec.Acc[i] *= 0.1;      // dm/s/s -> m/s/s
               }

               if(isExact) {
                  rec.sigPos[i] = sigP[i][Nmatch];
                  rec.sigVel[i] = sigV[i][Nmatch];
                  if(haveAcceleration) rec.sigAcc[i] = sigA[i][Nmatch];
               }
               else {
                  // TD is this sigma related to 'err' in the Lagrange call?
                  rec.sigPos[i] = RSS(sigP[i][Nhi],sigP[i][Nlow]);
                  rec.sigVel[i] = RSS(sigV[i][Nhi],sigV[i][Nlow]);
                  if(haveAcceleration)
                     rec.sigAcc[i] = RSS(sigA[i][Nhi],sigA[i][Nlow]);
               }
               // else Acc=sig_Acc=0   // TD can we do better?
            }
         }
         else {               // no V data - must interpolate position to get velocity
            for(i=0; i<3; i++) {
               // interpolate positions(km) to get P and V
               LagrangeInterpolation(times,P[i],dt,rec.Pos[i],rec.Vel[i]);
               rec.Vel[i] *= 10000.;         // km/sec -> dm/sec

               if(isExact) {
                  rec.sigPos[i] = sigP[i][Nmatch];
               }
               else {
                  rec.sigPos[i] = RSS(sigP[i][Nhi],sigP[i][Nlow]);
               }
               // TD
               rec.sigVel[i] = 0.0;
            }
         }
         return rec;
      }
      catch(InvalidRequest& e) { GPSTK_RETHROW(e); }
   }

   // Return the position for the given satellite at the given time
   // @param[in] sat the SatID of the satellite of interest
   // @param[in] ttag the time (CommonTime) of interest
   // @return Triple containing the position ECEF XYZ meters
   // @throw InvalidRequest if result cannot be computed, for example because
   //  a) the time t does not lie within the time limits of the data table
   //  b) checkDataGap is true and there is a data gap
   //  c) checkInterval is true and the interval is larger than maxInterval
   Triple PositionSatStore::getPosition(const SatID& sat, const CommonTime& ttag)
      const throw(InvalidRequest)
   {
      try {
         int i;
         DataTableIterator it1, it2, kt;

         if(getTableInterval(sat, ttag, Nhalf, it1, it2, true)) {
            // exact match
            for(unsigned int i=0; i<Nhalf; i++) ++it1;
            PositionRecord rec(it1->second);
            return rec.Pos;
         }

         // pull data out of the data table
         vector<double> times,P[3];
         CommonTime ttag0(it1->first);
         kt = it1;
         while(1) {
            times.push_back(kt->first - ttag0);    // sec
            for(i=0; i<3; i++)
               P[i].push_back(kt->second.Pos[i]);
            if(kt == it2) break;
            ++kt;
         };

         // interpolate
         Triple pos;
         double dt(ttag-ttag0), err;
         for(i=0; i<3; i++)
            pos[i] = LagrangeInterpolation(times,P[i],dt,err);

         return pos;
      }
      catch(InvalidRequest& e) { GPSTK_RETHROW(e); }
   }

   // Return the velocity for the given satellite at the given time
   // @param[in] sat the SatID of the satellite of interest
   // @param[in] ttag the time (CommonTime) of interest
   // @return Triple containing the velocity ECEF XYZ meters/second
   // @throw InvalidRequest if result cannot be computed, for example because
   //  a) the time t does not lie within the time limits of the data table
   //  b) checkDataGap is true and there is a data gap
   //  c) checkInterval is true and the interval is larger than maxInterval
   Triple PositionSatStore::getVelocity(const SatID& sat, const CommonTime& ttag)
      const throw(InvalidRequest)
   {
      try {
         int i;
         DataTableIterator it1, it2, kt;

         bool isExact(getTableInterval(sat, ttag, Nhalf, it1, it2, haveVelocity));
         if(isExact && haveVelocity) {
            for(unsigned int i=0; i<Nhalf; i++) ++it1;
            PositionRecord rec(it1->second);
            return rec.Vel;
         }

         // pull data out of the data table
         CommonTime ttag0(it1->first);
         vector<double> times,D[3];       // D will be either Pos or Vel

         kt = it1;
         while(1) {
            times.push_back(kt->first - ttag0);    // sec
            for(i=0; i<3; i++)
               D[i].push_back(haveVelocity ? kt->second.Vel[i] : kt->second.Pos[i]);
            if(kt == it2) break;
            ++kt;
         };

         // interpolate
         Triple Vel;
         double dt(ttag-ttag0), err;
         for(i=0; i<3; i++) {
            if(haveVelocity)
               Vel[i] = LagrangeInterpolation(times,D[i],dt,err);
            else {
               // interpolate positions(km) to get velocity // err is dummy
               LagrangeInterpolation(times,D[i],dt,err,Vel[i]);
               Vel[i] *= 10000.;                                  // km/s -> dm/s
            }
         }

         return Vel;
      }
      catch(InvalidRequest& e) { GPSTK_RETHROW(e); }
   }

   // Return the acceleration for the given satellite at the given time
   // @param[in] sat the SatID of the satellite of interest
   // @param[in] ttag the time (CommonTime) of interest
   // @return Triple containing the acceleration ECEF XYZ meters/second/second
   // @throw InvalidRequest if result cannot be computed, for example because
   //  a) the time t does not lie within the time limits of the data table
   //  b) checkDataGap is true and there is a data gap
   //  c) checkInterval is true and the interval is larger than maxInterval
   //  d) neither velocity nor acceleration data are present
   Triple PositionSatStore::getAcceleration(const SatID& sat, const CommonTime& ttag)
      const throw(InvalidRequest)
   {
      if(!haveVelocity && !haveAcceleration) {
         InvalidRequest e("Neither velocity nor acceleration data are present");
         GPSTK_THROW(e);
      }

      try {
         int i;
         DataTableIterator it1, it2, kt;

         bool isExact(getTableInterval(sat,ttag,Nhalf,it1,it2,haveAcceleration));
         if(isExact && haveAcceleration) {
            // exact match, and have acceleration data
            for(unsigned int i=0; i<Nhalf; i++) ++it1;
            PositionRecord rec(it1->second);
            return rec.Acc;
         }

         // pull data out of the data table
         vector<double> times,D[3];                // D will be either Vel or Acc
         CommonTime ttag0(it1->first);
         kt = it1;
         while(1) {
            times.push_back(kt->first - ttag0);    // sec
            for(i=0; i<3; i++)
               D[i].push_back(haveAcceleration ? kt->second.Acc[i]:kt->second.Vel[i]);
            if(kt == it2) break;
            ++kt;
         };

         // interpolate
         Triple Acc;
         double dt(ttag-ttag0), err;
         for(i=0; i<3; i++) {
            if(haveAcceleration) {
               Acc[i] = LagrangeInterpolation(times,D[i],dt,err);
            }
            else {
               LagrangeInterpolation(times,D[i],dt,err,Acc[i]);   // err is dummy
               Acc[i] *= 0.1;                                     // dm/s/s -> m/s/s
            }
         }

         return Acc;
      }
      catch(InvalidRequest& e) { GPSTK_RETHROW(e); }
   }

   // Add a PositionRecord to the store.
   void PositionSatStore::addPositionRecord(const SatID& sat, const CommonTime& ttag,
                                            const PositionRecord& rec)
      throw(InvalidRequest)
   {
      try {
         checkTimeSystem(ttag.getTimeSystem());

         int i;
         if(!haveVelocity)
            for(i=0; i<3; i++)
               if(rec.Vel[i] != 0.0) { haveVelocity = true; break; }
         if(!haveAcceleration)
            for(i=0; i<3; i++)
               if(rec.Acc[i] != 0.0) { haveAcceleration = true; break; }

         if(tables.find(sat) != tables.end() &&
            tables[sat].find(ttag) != tables[sat].end()) {
                  // record already exists in table
            PositionRecord& oldrec(tables[sat][ttag]);
            oldrec.Pos = rec.Pos;
            oldrec.sigPos = rec.sigPos;
            if(haveVelocity) { oldrec.Vel = rec.Vel; oldrec.sigVel = rec.sigVel; }
            if(haveAcceleration) { oldrec.Acc = rec.Acc; oldrec.sigAcc = rec.sigAcc; }
         }
         else {   // create a new entry in the table
             updateInterval(sat, ttag);
            tables[sat][ttag] = rec;
         }
      }
      catch(InvalidRequest& ir) { GPSTK_RETHROW(ir); }
   }

   // Add position data (only) to the store
   void PositionSatStore::addPositionData(const SatID& sat, const CommonTime& ttag,
                     const Triple& Pos, const Triple& Sig)
      throw(InvalidRequest)
   {
      try {
         checkTimeSystem(ttag.getTimeSystem());

         if(tables.find(sat) != tables.end() &&
            tables[sat].find(ttag) != tables[sat].end()) {
                  // record already exists in table
            PositionRecord& oldrec(tables[sat][ttag]);
            oldrec.Pos = Pos;
            oldrec.sigPos = Sig;
         }
         else {   // create a new entry in the table
            PositionRecord rec;
            rec.Pos = Pos;
            rec.sigPos = Sig;
            rec.Vel = rec.sigVel = rec.Acc = rec.sigAcc = Triple(0,0,0);

            updateInterval(sat, ttag);
            tables[sat][ttag] = rec;
         }
      }
      catch(InvalidRequest& ir) { GPSTK_RETHROW(ir); }
   }

   // Add velocity data (only) to the store
   void PositionSatStore::addVelocityData(const SatID& sat, const CommonTime& ttag,
                        const Triple& Vel, const Triple& Sig)
      throw(InvalidRequest)
   {
      try {
         checkTimeSystem(ttag.getTimeSystem());

         haveVelocity = true;

         if(tables.find(sat) != tables.end() &&
            tables[sat].find(ttag) != tables[sat].end()) {
                  // record already exists in table
            PositionRecord& oldrec(tables[sat][ttag]);
            oldrec.Vel = Vel;
            oldrec.sigVel = Sig;
         }
         else {   // create a new entry in the table
            PositionRecord rec;
            rec.Vel = Vel;
            rec.sigVel = Sig;
            rec.Pos = rec.sigPos = rec.Acc = rec.sigAcc = Triple(0,0,0);

            updateInterval(sat, ttag);
            tables[sat][ttag] = rec;
         }
      }
      catch(InvalidRequest& ir) { GPSTK_RETHROW(ir); }
   }

   // Add acceleration data (only) to the store
   void PositionSatStore::addAccelerationData(const SatID& sat,
                        const CommonTime& ttag, const Triple& Acc, const Triple& Sig)
       throw(InvalidRequest)
   {
      try {
         checkTimeSystem(ttag.getTimeSystem());

         haveAcceleration = true;

         if(tables.find(sat) != tables.end() &&
            tables[sat].find(ttag) != tables[sat].end()) {
                  // record already exists in table
            PositionRecord& oldrec(tables[sat][ttag]);
            oldrec.Acc = Acc;
            oldrec.sigAcc = Sig;
         }
         else {   // create a new entry in the table
            PositionRecord rec;
            rec.Acc = Acc;
            rec.sigAcc = Sig;
            rec.Vel = rec.sigVel = rec.Pos = rec.sigPos = Triple(0,0,0);

            updateInterval(sat, ttag);
            tables[sat][ttag] = rec;
         }
      }
      catch(InvalidRequest& ir) { GPSTK_RETHROW(ir); }
   }

   //@}

}  // End of namespace gpstk
