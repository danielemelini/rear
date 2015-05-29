#!/bin/bash
#
#
#  Copyright Daniele Melini & Giorgio Spada, 2014
#
#  This file is part of REAR.
#
#  REAR is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  REAR is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with REAR.  If not, see <http://www.gnu.org/licenses/>.
#
#
#
if [ $# -ne 2 ]; then 
   echo "Usage: compare_map.sh    <FILE1>   <FILE2>"
   echo ""
   exit
fi
#
FILE1=$1
FILE2=$2
#
echo ""
echo ""
echo "Comparing $FILE1 and $FILE2 ..."
echo ""
echo "Min/max ABSOLUTE differences: "
paste $FILE1 $FILE2  |\
    awk '{ print ($3-$9),($4-$10),($5-$11),($6-$12)}' |\
    minmax -H2 -C |\
    awk '{ printf "  - UrDOT:      %12.4e / %12.4e\n  - UthetaDOT:  %12.4e / %12.4e\n  - UphiDOT:    %12.4e / %12.4e\n  - NDOT:       %12.4e / %12.4e\n", $1,$2,$3,$4,$5,$6,$7,$8}'    
#    awk '{ print "  - UrDOT:      "$1" / "$2"\n  - UthetaDOT:  "$3" / "$4"\n  - UphiDOT:    "$5" / "$6"\n  - NDOT:       "$7" / "$8 }'    
#
echo ""
#
echo "Min/max RELATIVE differences: "
paste $FILE1 $FILE2  |\
    awk '{ if ( (NR>2) && ($3 != 0) && ($4 != 0) && ($5 != 0) && ($6 != 0) ) print (($3-$9)/$3),(($4-$10)/$4),(($5-$11)/$5),(($6-$12)/$6)}' |\
    minmax -C | \
    awk '{ printf "  - UrDOT:      %12.4e / %12.4e\n  - UthetaDOT:  %12.4e / %12.4e\n  - UphiDOT:    %12.4e / %12.4e\n  - NDOT:       %12.4e / %12.4e\n", $1,$2,$3,$4,$5,$6,$7,$8}'    
#
echo ""
echo ""


#    awk '{ print "  - UrDOT:      "$1" / "$2"\n  - UthetaDOT:  "$3" / "$4"\n  - UphiDOT:    "$5" / "$6"\n  - NDOT:       "$7" / "$8 }'    
