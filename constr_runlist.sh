#!/bin/bash
#usage: ./constr_runlist.sh runlist

# 1st argument: Path to runlistfile

# Currently implemented scheme in main.C: runNr runName pos pdgID energy angle (supports runs 29-75)

# This script creates a runlist supporting runs 19-107 outputting two schemes:
# runs 19-87  -> runNr runName MP pdgID energy angle WCidentifier
# runs 88-107 -> runNr runName MP pdgID energy angle WCidentifier sidePostitionX sidePostitionY
# prints each line into the passed runlist file
# overwrites passed runlist file

# LIST of existing naming schemes:
# runNr 1-18 ->
# "1_18_dist_scan_long"
# runNr 19-28 ->
# scheme 1: runNr_particleNameEnergyvalue_posNr
# runNr 29-75 ->
# scheme 2: runNr_particleNameEnergyvalue_posNr_WOMpos
# runNr 76-87 ->
# scheme 3: runNr_particleNameEnergyvalue_posNr_angleValue_WOMpos
# runNr 88-107 ->
# scheme 4: runNr_particleNameEnergyvalue_scanBD_xPos_yPos_angleValue_WOMpos
# runNr 108 ->
# "108_muon6_posS1_angle90_AB"
# runNr 109 ->
# "109_trigger_test1"

##########################
## PARSING runNr 19-107 ##
##########################

here=`pwd`
data=$here/data

runlist=$1

parse ()
{
	dir_name=$0
	rl_file=$1

	# construct array. contains the elements separated by "_" delimiter in $dir_name
	IFS="_" read -r -a fields <<< "$dir_name"
	
	nfields=${#fields[@]}

	# get run number
	runNr=${fields[0]}

	# get pdgID, beam energy
	runParticle=${fields[1]} 
	case $runParticle in
		muon[6])
			pdgID="-13"
			energy="6"
			;;
		pion[1-6])
			pdgID="211"
			energy=$(echo $runParticle | cut -c 5)
			;;
		e5)
			pdgID="-11"
			energy=5
			;;
		*)
			echo "UNKNOWN particle description in run $runNr"
			;;
	esac

	# get measurement position
	runMP=${fields[2]}
	case $runMP in
		pos[0-9] | pos[1][0-9] ) # 0°/30° measurements, no xy-coordinates
			frontMP=$(echo $runMP | cut -c 4-)
			;;
		scanBD) # 90° WOM scans, hardcoded $sideMP
			side_pos_x=${fields[3]}
			side_pos_y=${fields[4]}
			case $side_pos_y in
				0)
					case $side_pos_x in
						0)sideMP=18;; 6)sideMP=19;; 12)sideMP=20;; 18)sideMP=21;; 24)sideMP=22;;
					esac
					;;
				1)
					case $side_pos_x in
						0)sideMP=23;; 6)sideMP=24;; 12)sideMP=25;; 18)sideMP=26;; 24)sideMP=27;;
					esac
					;;
				2)
					case $side_pos_x in
						0)sideMP=28;; 6)sideMP=29;; 12)sideMP=30;; 18)sideMP=31;; 24)sideMP=32;;
					esac
					;;
				3)
					case $side_pos_x in
						0)sideMP=33;; 6)sideMP=34;; 12)sideMP=35;; 18)sideMP=36;; 24)sideMP=37;;
					esac
					;;
			esac
			;;
		*)
			echo "UNKNOWN position description in run $runNr"
			;;
	esac

	# get angle
	case $nfields in
		3)
			angle="0"
			;;
		4)	
			angle="0"
			;;
		5)	
			angle=$(echo "${fields[3]}" | cut -c 6-)
			;;
		7)	
			angle=$(echo "${fields[5]}" | cut -c 6-)
			;;
		*)
			echo "UNKNOWN angele description in run $runNr"
			;;

	esac
	# get WaveCatcher identifier
	case $nfields in
		3)
			WC="AB"
			;;
		[4-7])
			WC=${fields[@]: -1:1} # accesses last element in ${fields}.
			;;
		*)
			echo "UNKNOWN WaveCatcher identifier in run $runNr"
			;;
	esac

	######################
	## PRINT TO RUNLIST ##
	######################
	# line_to_runlsit
	case $nfields in
		[3-5])
			line_to_runlsit="$runNr $dir_name $frontMP $pdgID $energy $angle $WC"
			;;
		7)
			line_to_runlsit="$runNr $dir_name $sideMP $pdgID $energy $angle $WC $side_pos_x $side_pos_y"
			;;			
	esac

	echo "$line_to_runlsit"
	echo "$line_to_runlsit" >> "$rl_file"


}
export -f parse

# initialize new runlist, first three lines comments
DATE=`date '+%Y-%m-%d %H:%M:%S'`
echo "-> new runlist created: $DATE <-" > $runlist # fist line in runlist. overwrite old runlist
echo "runs 19-87  -> runNr runName MP pdgID energy angle WCidentifier" >> $runlist
echo "runs 88-107 -> runNr runName MP pdgID energy angle WCidentifier sidePostitionX sidePostitionY" >> $runlist


######################
## INITIATE PARSING ##
######################

ls $data | sort -n | xargs -n 1 bash -c "parse $runlist"
