// This file is part of LatMRG.
//
// Copyright (C) 2012-2023  The LatMRG authors, under the supervision
// of Pierre L'Ecuyer at Universit� de Montr�al.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef LATMRG_CHRONO_H
#define LATMRG_CHRONO_H

#include <string>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;


#undef USE_ANSI_CLOCK

namespace LatMRG {

  /**
   * This class acts as an interface to the system clock to compute the CPU
   * time used by parts of a program. Even though the ANSI/ISO macro
   * <tt>CLOCKS_PER_SEC = 1000000</tt> is the number of clock ticks per second
   * for the value returned by the ANSI-C standard `clock` function (so this
   * function returns the number of microseconds), on some systems where the
   * 32-bit type `long` is used to measure time, the value returned by `clock`
   * wraps around to negative values after about 36 minutes. On some other
   * systems where time is measured using the 32-bit type `unsigned long`, the
   * clock may wrap around to 0 after about 72 minutes. When the macro
   * <tt>USE_ANSI_CLOCK</tt> is undefined, a non-ANSI-C clock is used. On
   * Linux-Unix systems, it calls the POSIX function `times` to get the CPU
   * time used by a program. On a Windows platform (when the macro
   * <tt>HAVE_WINDOWS_H</tt> is defined), the Windows function
   * `GetProcessTimes` will be used to measure the CPU time used by programs.
   * On Linux\f$|\f$Unix platforms, if the macro <tt>USE_ANSI_CLOCK</tt> is
   * defined, the timers will call the ANSI C `clock` function. When
   * <tt>USE_ANSI_CLOCK</tt> is left undefined, class `Chrono` gets the CPU
   * time used by a program via an alternate non-ANSI C timer based on the
   * POSIX (The Portable Operating System Interface) function `times`, assuming
   * this function is available. The POSIX standard is described in the IEEE
   * Std 1003.1-2001 document (see The Open Group web site at
   * [http://www.opengroup.org/onlinepubs/007904975/toc.htm](http://www.opengroup.org/onlinepubs/007904975/toc.htm)).
   *
   * Every object `Chrono` acts as an independent *stopwatch*. Several such
   * stopwatchs can run at any given time. An object of type `Chrono` must be
   * declared for each of them. The method `init` resets the stopwatch to zero,
   * `val` returns its current reading, and `write` writes this reading to the
   * current output. The returned value includes part of the execution time of
   * the functions from class `Chrono`. The `TimeFormat` allows one to choose
   * the kind of time units that are used.
   *
   * Below is an example of how the functions may be used. A stopwatch named
   * `timer` is declared and created. After 2.1 seconds of CPU time have been
   * consumed, the stopwatch is read and reset. Then, after an additional 330
   * seconds (or 5.5 minutes) of CPU time the stopwatch is read again, printed
   * to the output and deleted.
   *
   * <tt>Chrono timer; <br>
   * (...) (<em>suppose 2.1 CPU seconds are used here</em>.)<br>
   * double t = timer.val (Chrono::SEC); // Here, t = 2.1 <br>
   * timer.init(); <br>
   * (...) (<em>suppose 330 CPU seconds are used here</em>.) <br>
   * t = timer.val (Chrono::MIN); // Here, t = 5.5 <br>
   * timer.write (Chrono::HMS); // Prints: 00:05:30.00 </tt>
   */

  class Chrono {
    public:

      /**
       * Types of units in which the time on a `Chrono` can be read or printed: in
       * seconds (<tt>SEC</tt>), minutes (<tt>MIN</tt>), hours (<tt>HOUR</tt>),
       * days (<tt>DAYS</tt>), or in the `HH:MM:SS.xx` format, with hours, minutes,
       * seconds and hundreths of a second (<tt>HMS</tt>).
       */
      enum TimeFormat { SEC, MIN, HOURS, DAYS, HMS };

      /**
       * Constructor for a stopwatch; initializes it to zero. One may
       * reinitialize it later by calling `init`.
       */
      Chrono();

      /**
       * Destructor.
       */
      ~Chrono() {}

      /**
       * (Re)Initializes this stopwatch to zero.
       */
      void init ();

      /**
       * Returns the CPU time measured by this `Chrono` since the last call
       * to `init()`. The parameter `unit` specifies the time unit.
       * Restriction: `unit = HMS` is not allowed here; it will cause an
       * error.
       */
      double val (TimeFormat unit);

      /**
       * Prints, on standard output, the CPU time measured by this `Chrono`
       * since its last call to `init()`. The parameter `unit` specifies the
       * time unit.
       */
      void write (TimeFormat unit);

      /**
       * Returns as a string the CPU time measured by this `Chrono` since its
       * last call to `init()`. The time format used is `HMS`.
       */
      std::string toString();

      /**
       * Returns `true` if this `Chrono` has reached the time `limit` (in
       * seconds), otherwise returns `false`.
       */
      bool timeOver (double limit);

    private:

      /// Microseconds
      unsigned long microsec;

      /// Seconds
      unsigned long second;

      /** Function returning the CPU time used by the program since it was
       * started. This function depends on the operation system and is not
       * intended to be manipulated directly.
       * */
      void tick();
  };

  /**
   * Returns the value of the duration from `timer`.
   */
  std::string toString (Chrono& timer);

};

//============================================================================
// Implementation

namespace {

#ifdef HAVE_WINDOWS_H
#include <windows.h>

  HLANDLE currentProcess = NULL;

  /*
   * A helper function for converting FILETIME to a LONGLONG [safe from memory
   * alignment point of view].
   */
  ULONGLONG fileTimeToInt64 (const FILETIME * time)
  {
    ULARGE_INTEGER _time;
    _time.LowPart = time->dwLowDateTime;
    _time.HighPart = time->dwHighDateTime;
    return _time.QuadPart;
  }
#endif
   
}

//===========================================================================
//Alternative implementation of tick for different operating systems
  
#ifdef HAVE_WINDOWS_H

void Chrono::tick () {
    if (currentProcess == NULL)
      currentProcess = GetCurrentProcess();
    FILETIME creationTime, exitTime, kernelTime, userTime;
    /* Strongly inspired from
     * http://www.javaworld.com/javaworld/javaqa/2002-11/01-qa-1108-cpu.html */
    GetProcessTimes (currentProcess, &creationTime, &exitTime,
        &kernelTime, &userTime);
    ULONGLONG rawTime = fileTimeToInt64 (&kernelTime) +
      fileTimeToInt64 (&userTime);
    /* We have to divide by 10000 to get milliseconds out of
     * the computed time */
    second = static_cast<unsigned long>(rawTime / 10000000);
    microsec = static_cast<unsigned long>((rawTime % 10000000) / 10);
}

#elif defined(USE_ANSI_CLOCK)
  
// ANSI C timer
void Chrono::tick () {
    clock_t t;
    double y;

    t = clock ();
    y = (static_cast<double> (t)) / CLOCKS_PER_SEC;
    second = static_cast<unsigned long>(y);
    microsec = static_cast<unsigned long>((y - second) * 1000000);
}

#else
  // POSIX timer

#include <sys/times.h>
#include <unistd.h>

void Chrono::tick () {
    struct tms us;
    long TICKS, z;

    TICKS = sysconf(_SC_CLK_TCK);
    if (TICKS == -1) {
      cout << "Chrono.cc:   'sysconf(_SC_CLK_TCK)' failed\n";
    }
    z = times (&us);
    if (z == -1) {
      cout << "Chrono.cc:   timer times failed\n";
    }
    /* CPU time = user time + system time */
    microsec = us.tms_utime + us.tms_stime;
    second = microsec / TICKS;
    microsec = (microsec % TICKS) * 1000000 / TICKS;
}

#endif


//===========================================================================

void Chrono::init () {
     tick();
}


Chrono::Chrono() {
    init();
}

//===========================================================================

double Chrono::val (TimeFormat Unit) {
     Chrono now;
     now.tick();
     double temps;                     // Time elapsed, in seconds
     temps = (static_cast<double>(now.microsec) -
         static_cast<double>(microsec)) / 1.E+6 +
       static_cast<double>(now.second) -
       static_cast<double>(second);

     switch (Unit) {
       case SEC:
         return temps;
       case MIN:
         return temps * 1.666666667E-2;
       case HOURS:
         return temps * 2.777777778E-4;
       case DAYS:
         return temps * 1.157407407E-5;
       case HMS:
         cerr << "Chrono.val: HMS is a wrong arg for Time Unit";
     }
     return 0.0;
}

//===========================================================================

void Chrono::write (TimeFormat Form) {
     double temps;
     if (Form != HMS)
       temps = val (Form);
     else
       temps = 0.0;
     switch (Form) {
       case SEC:
         cout << setw(10) << setprecision (2) << temps;
         cout << " seconds";
         break;
       case MIN:
         cout << setw(10) << setprecision (2) << temps;
         cout << " minutes";
         break;
       case HOURS:
         cout << setw(10) << setprecision (2) << temps;
         cout << " hours";
         break;
       case DAYS:
         cout << setw(10) << setprecision (2) << temps;
         cout << " days";
         break;
       case HMS:
         temps = val (SEC);
         long heure = static_cast<long> (temps * 2.777777778E-4);
         if (heure > 0)
           temps -= static_cast<double> (heure) * 3600.0;
         long minute = static_cast<long> (temps * 1.666666667E-2);
         if (minute > 0)
           temps -= static_cast<double> (minute) * 60.0;
         long seconde = static_cast<long> (temps);
         long centieme = static_cast<long>
           (100.0 * (temps - static_cast<double> (seconde)));
         cout << setw(2) << setfill('0') << right << heure << ":" ;
         cout << setw(2) << setfill('0') << right << minute << ":" ;
         cout << setw(2) << setfill('0') << right << seconde << "." ;
         cout << setw(2) << setprecision(2) << centieme;
         break;
     }
}

//===========================================================================

bool Chrono::timeOver (double limit) {
     double temps = val (SEC);
     if (temps >= limit)
       return true;
     else
       return false;
}

//===========================================================================

std::string Chrono::toString () {
     double temps = val (SEC);
     long heure = (long) (temps * 2.777777778E-4);
     if (heure > 0)
       temps -= heure * 3600.0;
     long minute = (long) (temps * 1.666666667E-2);
     if (minute > 0)
       temps -= minute * 60.0;
     long seconde = (long) temps;
     long frac = (long) (100.0 * (temps - seconde));

     std::ostringstream sortie;
     sortie << heure << ":" << minute << ":" << seconde << "." << frac;
     return sortie.str ();
}

#endif
