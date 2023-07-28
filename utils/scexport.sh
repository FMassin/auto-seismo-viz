#!/bin/bash
EVENTID=$2
STORE=$3
PREF="/opt/seiscomp/bin/seiscomp exec" 
ENVALIAS="sceewvenv"
#FDSNWS='arclink.ethz.ch:8080'
if [ !  -n "$2" ]; then
    STORE='./'
fi

ssh -o RemoteCommand=none -T $1 <<EOI

IFS="," read -r ID TI LA LO DE ML DESC <<< \$(curl "http://localhost:8080/fdsnws/event/1/query?eventid=${EVENTID}&format=csv")
BEG=\$(date -d "\$TI -1 min" "+%FT%T")
END=\$(date -d "\$TI +3 min" "+%FT%T")

IFS="=" read -r TRASH DBINV <<< \$( grep "^database.inventory" ~/.seiscomp/global.cfg )
IFS="=" read -r TRASH DBLOC <<< \$( grep "^database" ~/.seiscomp/global.cfg |grep -v inventory|grep -v config)
IFS="=" read -r FDSNWS <<< \$( grep "^recordstream" ~/.seiscomp/global.cfg |sed 's/.*;fdsnws\///;s/\/fdsnws.*//')

echo " "
echo "      pwd:" \$(pwd)
echo "     pref:" $PREF
echo "    store:" $STORE 
echo " event id:" $EVENTID
echo " env. al.:" $ENVALIAS
echo "  inv. db:" \$DBINV
echo "   ev. db:" \$DBLOC
echo "   fdsnws:" \$FDSNWS
echo "    begin:" \$BEG
echo "      end:" \$END

mkdir -p $STORE 

# curl "http://\$FDSNWS/fdsnws/dataselect/1/query?starttime=\$BEG&endtime=\$END" |$PREF scmssort -v -E -u |$PREF $ENVALIAS --config.skipDataOlderThan=999999999  -u test --dump --debug -I - > $STORE/${EVENTID//\//_}.mseed 2>$STORE/${EVENTID//\//_}.env.log || tail $STORE/${EVENTID//\//_}.env.log
$PREF scxmldump -f -E ${EVENTID} -d \$DBLOC -PAMFm > $STORE/${EVENTID//\//_}.xml 2>$STORE/${EVENTID//\//_}.dumpev.log || tail $STORE/${EVENTID//\//_}.dumpev.log
# $PREF scxmldump -If -d \$DBINV > $STORE/${EVENTID//\//_}.inv.xml 2>$STORE/${EVENTID//\//_}.dumpinv.log || tail $STORE/${EVENTID//\//_}.dumpinv.log
# curl "http://\$FDSNWS/fdsnws/station/1/query?starttime=\$BEG&endtime=\$END&level=channel&format=stationxml" > $STORE/${EVENTID//\//_}.inv.xml 2>$STORE/${EVENTID//\//_}.dumpinv.log || tail $STORE/${EVENTID//\//_}.dumpinv.log 
exit
EOI

scp $1:$STORE/${EVENTID//\//_}"*xml" ./
