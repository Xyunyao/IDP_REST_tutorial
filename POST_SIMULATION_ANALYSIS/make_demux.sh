#!/bin/bash

index=$1
restdir=$2

if [ $3 -eq 8 ] ; then

gmx trjcat -f $restdir/0/prod.xtc $restdir/1/prod.xtc $restdir/2/prod.xtc $restdir/3/prod.xtc $restdir/4/prod.xtc $restdir/5/prod.xtc $restdir/6/prod.xtc $restdir/7/prod.xtc -demux $index -o 0.xtc 1.xtc 2.xtc 3.xtc 4.xtc 5.xtc 6.xtc 7.xtc

elif [ $3 -eq 16 ] ; then

gmx trjcat -f $restdir/0/prod.xtc $restdir/1/prod.xtc $restdir/2/prod.xtc $restdir/3/prod.xtc $restdir/4/prod.xtc $restdir/5/prod.xtc $restdir/6/prod.xtc $restdir/7/prod.xtc $restdir/8/prod.xtc $restdir/9/prod.xtc $restdir/10/prod.xtc $restdir/11/prod.xtc $restdir/12/prod.xtc $restdir/13/prod.xtc $restdir/14/prod.xtc $restdir/15/prod.xtc -demux $index -o 0.xtc 1.xtc 2.xtc 3.xtc 4.xtc 5.xtc 6.xtc 7.xtc 8.xtc 9.xtc 10.xtc 11.xtc 12.xtc 13.xtc 14.xtc 15.xtc

elif [ $3 -eq 10 ] ; then

gmx trjcat -f $restdir/0/prod.xtc $restdir/1/prod.xtc $restdir/2/prod.xtc $restdir/3/prod.xtc $restdir/4/prod.xtc $restdir/5/prod.xtc $restdir/6/prod.xtc $restdir/7/prod.xtc $restdir/8/prod.xtc $restdir/9/prod.xtc -demux $index -o 0.xtc 1.xtc 2.xtc 3.xtc 4.xtc 5.xtc 6.xtc 7.xtc 8.xtc 9.xtc

elif [ $3 -eq 20 ] ; then

gmx trjcat -f $restdir/0/prod.xtc $restdir/1/prod.xtc $restdir/2/prod.xtc $restdir/3/prod.xtc $restdir/4/prod.xtc $restdir/5/prod.xtc $restdir/6/prod.xtc $restdir/7/prod.xtc $restdir/8/prod.xtc $restdir/9/prod.xtc $restdir/10/prod.xtc $restdir/11/prod.xtc $restdir/12/prod.xtc $restdir/13/prod.xtc $restdir/14/prod.xtc $restdir/15/prod.xtc $restdir/16/prod.xtc $restdir/17/prod.xtc $restdir/18/prod.xtc $restdir/19/prod.xtc -demux $index -o 0.xtc 1.xtc 2.xtc 3.xtc 4.xtc 5.xtc 6.xtc 7.xtc 8.xtc 9.xtc 10.xtc 11.xtc 12.xtc 13.xtc 14.xtc 15.xtc 16.xtc 17.xtc 18.xtc 19.xtc

fi

