cristinel.ababei@ndsu.edu
December 2009, Fargo ND


Synopsis
========
This is "reliablenoc" tool: a Branch-and-Bound (BB) mapping algorithm
for regular Networks-on-Chip (NoCs). Mapping is multiobjective: energy-
and reliability-aware. The tool is developed on top of "nocmap 1.2"
tool (from CMU), which is an energy-aware only mapping algorithm.


Installation
============
The latest version of the tool can be downloaded from:
http://venus.ece.ndsu.nodak.edu/~cris/software.html
The tool was developed on a Linux machine running Fedora 8. However,
it should be compile-able on Windows too. On linux, first edit the 
Makefile to reflect the location where you want to compile and link. 
Then, just type:
> make


How to use the tool
===================
Type "nocmap" at the command prompt to see the help menu.


Examples:
---------
nocmap -tilenum 9 -synthesis_routing odd_even -traffic_config tests/mpeg4.config -me -mr -alpha 0.0 -v
nocmap -tilenum 9 -synthesis_routing odd_even -traffic_config tests/mpeg4.config -me -mr -alpha 0.6 -v
nocmap -tilenum 9 -synthesis_routing odd_even -traffic_config tests/mpeg4.config -me -mr -alpha 0.0 -annealing -v
nocmap -tilenum 9 -synthesis_routing odd_even -traffic_config tests/mpeg4.config -me -mr -alpha 0.6 -annealing -v
nocmap -tilenum 16 -synthesis_routing odd_even -traffic_config tests/telecom.config -me -mr -alpha 0.0 -v
nocmap -tilenum 16 -synthesis_routing odd_even -traffic_config tests/telecom.config -me -mr -alpha 0.6 -v
nocmap -tilenum 16 -synthesis_routing odd_even -traffic_config tests/telecom.config -me -mr -alpha 0.0 -annealing -v
nocmap -tilenum 16 -synthesis_routing odd_even -traffic_config tests/telecom.config -me -mr -alpha 0.6 -annealing -v
nocmap -tilenum 25 -synthesis_routing odd_even -traffic_config tests/ami25.config -me -mr -alpha 0.0 -v
nocmap -tilenum 25 -synthesis_routing odd_even -traffic_config tests/ami25.config -me -mr -alpha 0.6 -v
nocmap -tilenum 25 -synthesis_routing odd_even -traffic_config tests/ami25.config -me -mr -alpha 0.0 -annealing -v
nocmap -tilenum 25 -synthesis_routing odd_even -traffic_config tests/ami25.config -me -mr -alpha 0.6 -annealing -v
nocmap -tilenum 49 -synthesis_routing odd_even -traffic_config tests/ami49.config -me -mr -alpha 0.0 -v
nocmap -tilenum 49 -synthesis_routing odd_even -traffic_config tests/ami49.config -me -mr -alpha 0.6 -v
nocmap -tilenum 49 -synthesis_routing odd_even -traffic_config tests/ami49.config -me -mr -alpha 0.0 -annealing -v
nocmap -tilenum 49 -synthesis_routing odd_even -traffic_config tests/ami49.config -me -mr -alpha 0.6 -annealing -v

Notes:
------
1.  When tool is run with options "-mr -alpha 0.0", mapping is done 
essentially with energy only as objective. However, cost calculation
inside the BB algorithm is done with normalization.
2.  To run the mapping tool with energy only as objective and 
without cost normalization do not use "-mr" option. Example:
> nocmap -tilenum 9 -synthesis_routing odd_even -traffic_config tests/mpeg4.config -me -v
3.  There may be differences in the solutions obtained by
running the tool as described in the previous two notes. They
should not be, or they should not be too large. I am still trying
to figure out why these diferences exist.


Credits
=======
Jingcao Hu (while at Carnegie Mellon University) developed "nocmap 1.2"
from which I developed "reliablenoc". Many thanks to SLD group (Radu
Marculescu) at Carnegie Mellon University that make "nocmap 1.2" source
code available at: http://www.ece.cmu.edu/~sld/software/nocmap.php


Copyright
=========
Copyright 2009 by Cristinel Ababei, cristinel.ababei@ndsu.edu
This Copyright notice applies to all files, called hereafter 
"The Software".
Permission to use, copy, and modify this software and its 
documentation is hereby granted only under the following 
terms and conditions.  Both the above copyright notice and 
this permission notice must appear in all copies of the 
software, derivative works or modified versions, and any 
portions thereof, and both notices must appear in supporting 
documentation.  Permission is granted only for non-commercial 
use.  For commercial use, please contact the author.
This software may be distributed (but not offered for sale 
or transferred for compensation) to third parties, provided 
such third parties agree to abide by the terms and conditions
of this notice.
The Software is provided "as is", and the authors, the 
North Dakota State University (NDSU), as well as any and
all previous authors (of portions or modified portions of
the software) disclaim all warranties with regard to this 
software, including all implied warranties of merchantability
and fitness.  In no event shall the authors or NDSU or any and
all previous authors be liable for any special, direct, 
indirect, or consequential damages or any damages whatsoever
resulting from loss of use, data or profits, whether in an
action of contract, negligence or other tortious action,
arising out of or in connection with the use or performance
of this software.


Copyright notice of "nocmap.1.2" from Carnegie Mellon University
================================================================
The author of "nocmap.1.2" is Jingcao Hu (jingcao@ece.cmu.edu)

Copyright (c) 200[X] Carnegie Mellon University.  All rights
reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in
   the documentation and/or other materials provided with the
   distribution.

3. The end-user documentation included with the redistribution,
   if any, must include the following acknowledgment:
      "This product includes software developed by
       Carnegie Mellon University (http://www.cmu.edu/)."
   Alternately, this acknowledgment may appear in the software itself,
   if and wherever such third-party acknowledgments normally appear.

4. The names "Carnegie Mellon University" and "CMU" must
   not be used to endorse or promote products derived from this
   software.

5. Products derived from this software may not be called
    "Carnegie Mellon University" or "CMU",
   nor may "Carnegie Mellon University"  or "CMU" appear in their name.

THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE APACHE SOFTWARE FOUNDATION OR
ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.
