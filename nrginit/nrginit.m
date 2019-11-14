(* 
  "NRG Ljubljana" - Numerical renormalization group for multiple
  impurities and an arbitrary number of channels
  
  Copyright (C) 2005-2016 Rok Zitko
  
     This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.
  
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
  
     You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

   Contact information:
   Rok Zitko
   F5 - Condensed matter physics
   "Jozef Stefan" Institute
   Jamova 39
   SI-1000 Ljubljana
   Slovenia

   rok.zitko@ijs.si

   $Id: nrginit.m 869 2009-11-04 16:54:09Z zitko $
*)

(* Default search path with some convenient defaults. *)
(* Rationale: from more specialized to more generic... *)
PACKAGEPATH = {"~",
               "~/nrg",
               NRGDIR};

(**************************************************************)

Print["NRG Ljubljana -- impurity problem initialization code"];
Print["Dir: ", Directory[]];
Print["Path: ", PACKAGEPATH];
SYMTYPE = "runtime";

(* Load sneg first! *)
Get["sneg.m", Path->PACKAGEPATH];

res=Get["initialparse.m", Path->PACKAGEPATH];
If[res == $Failed,
  Print["Can't load initialparse.m"];
  Exit[1];
];

MyPrint["Parsing parameters"];
res = parse["param"];
If[res == $Failed,
  Print["Failed parsing parameter file."];
  Exit[1];
];
PARSED=True;

Get["initial.m", Path->PACKAGEPATH];

makedata["data"];

Print["Success!"];

