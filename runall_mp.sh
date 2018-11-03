#!/bin/bash
#usage: ./runall_mp.sh case runlist

# Handles WaveCatcher data for processing. Executes processing in multiple threads using BSD xargs
# Uses different cases depending on which format the runlist is available
# 1st argument: case
# 2nd argument: runlist text file
# case == 1: No runlist needed
# 			 Gathers list of all WaveCatcher direcories in ./data and passes it to "read" executable
#			 Expects naming scheme of WaveCatcher run directory as follows:
#			 "RunNr_ParticleName&Energy_Position_WOMName"
#			 e.g. "44_muon6_pos5_AB"
# case == 2: Customized to parse format of TB18 runlist file created by using "constr_runlist.sh"
# 			 Expects white space separeted text file wherein each line corresponds to a single run
# 			 Expected format of each runlist line:
# 			 runs 19-87  -> runNr runName MP pdgID energy angle WCidentifier
# 			 runs 88-107 -> runNr runName MP pdgID energy angle WCidentifier sidePostitionX sidePostitionY


# Verify selected case
case_selector=$1
echo ""; echo "selected case: $case_selector"

# Verify passed runlist file
if [[ $case_selector != 1 ]]
then
	runlist=$2 
	echo "used runlist: $runlist"; echo "";
fi 

#################
### DEBUGGING ###
#################

# Test functions for fast debugging with xargs
test_fn()
{
	echo -n "-> "; for a in "$0"; do echo -n "\"$a\" "; done; echo
	# sleep 1 # show xargs parallel mode 
}
export -f test_fn

test_fn2()
{
	echo $0; echo $1; echo $2; echo $3;
}
export -f test_fn2


#################
# HANDLE CASES ##
#################

case $case_selector in
	1)	
		# No runlist needed. Fetches list of files to be processed by "ls" command from ./data 
		# Expects naming scheme of WaveCatcher run directory as follows:
		# "RunNr_ParticleName&Energy_Position_WOMName"
		# e.g. "44_muon6_pos5_AB"

		work_data()
		{
			here=`pwd`

			if [ ! -d "$here/runs" ]; then
  				mkdir $here/runs
			fi

			runName=$0
			# echo $runName # dummy functionality for debugging

			mkdir $here/runs/$runName
			if [ ! -e $here/runs/$runName/$runName.list ]; then
				ls $here/data/$runName | grep \.bin > $here/runs/$runName/$runName.list
			fi

			# Parse arguments passed to "read" executable 
			inFileList=$here/runs/$runName/$runName.list
			inDataFolder=$here/data/$runName/
			outFile=$here/runs/$runName/out.root

			# reads all characters infront of first "_"-delimiter in 
			runNr=$(echo $runName | cut -d "_" -f 1) 

			time $here/read $inFileList $inDataFolder $outFile $runNr
		}

		# temporary declare bash function "work_data" to PATH variable
		export -f work_data

		# pipe the list of files stored in "./data" into xargs
		# "-n 1" option assures that onyl one argument (here single line of files list) is used
		# "-P 8" option specifies number of parallel threads used
		# xargs utilizes bash that executes function "work_data" 
		ls ./data | xargs -n 1 -P 8 bash -c "work_data"
		;;
	2)	
		work_data()
		{
			here=`pwd`
			rl_line=$0 # passed runlist line

			# Parse run list line.
			# construct array. contains the elements separated by " " delimiter in $rl_line
			IFS=" " read -r -a fields <<< "$rl_line"
			nfields=${#fields[@]} # number of elements in array

			# debugging
			# for element in "${fields[@]}"; do
			# 	echo -n "$element " 
			# done; echo ""

			runNr=${fields[0]}
			runName=${fields[1]}
			MP=${fields[2]}
			pdgID=${fields[3]}
			energy=${fields[4]}
			angle=${fields[5]}
			WC=${fields[6]}
			# xy coordinates in case of 90° WOM scan (runs 88-107)
			case $nfields in
				7) ;;
				9)
					side_pos_x=${fields[7]}
					side_pos_y=${fields[8]}
					;;
				*)
					echo "UNKNOWN runlsit format" ;;
			esac

			# create directories where the analized date will be stored in
			if [[ ! -d "$here/runs" ]]; then mkdir $here/runs; fi
			if [[ ! -d "$here/runs/$runName" ]]; then mkdir $here/runs/$runName; fi
			
			# create list of WaveCatcher files
			if [[ ! -e $here/runs/$runName/$runName.list ]]; then
				ls $here/data/$runName | grep \.bin > $here/runs/$runName/$runName.list
			fi

			inFileList=$here/runs/$runName/$runName.list
			inDataFolder=$here/data/$runName/
			outFile=$here/runs/$runName/out.root

			echo "$runNr $runName $MP $pdgID $energy $angle $WC $side_pos_x $side_pos_y"

			# process data
			case $nfields in
				7)
					time $here/read $inFileList $inDataFolder $outFile $runNr $MP $pdgID $energy $angle $WC	
					;;
				9)
					time $here/read $inFileList $inDataFolder $outFile $runNr $MP $pdgID $energy $angle $WC	$side_pos_x $side_pos_y
					;;					
			esac
		}

		# temporary declare bash function "work_data" to PATH variable
		export -f work_data

		# tail reads runlist starting from 4th line
		# tr translates read EOL to NULL
		# Pipe runlist into xargs that executes child command once per runlist line
		# -0 options splits input around NULL bytes
		# "-n 1" option insures that onyl one command per line is executed
		# "-P 8" option specifies number of parallel threads used, here 8 threads
		# xargs utilizes bash that executes function "work_data"
		tail -n +4 "$runlist" | tr "\n" "\0" | xargs -0 -n 1 -P 8 bash -c "work_data"

		# debugging
		# tail -n +4 "test_runlist.txt" | tr "\n" "\0" | xargs -0 -P 8 -n 1 bash -c "test_fn"

		;;
	*)
		echo "proper case argument nedded. valid arguments: 1,2 "
		;;

esac



