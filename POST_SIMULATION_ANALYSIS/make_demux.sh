#!/bin/bash

#!/bin/bash

# Created by Jaya Krishna
# Edits by Korey Reid





index=""
restdir=""
name=""
log=""
verbose=""

usage() {
	echo "Usage: $0 [Options]"
	echo "Options:"
	echo " -h, --help       Display this help message"
	echo " -v, --verbose    Enable verbosity"
    echo " -i, --index      replica index xvg file"
    echo " -d, --restdir    replica exchange directory"
    echo " -n, --name       trajectory name, i.e. whole.xtc"
	echo " -l, --log        STDIO log File"
}

has_argument() {
	[[ ("$1" == *=* && -n ${1#*=}) || ( ! -z "$2" && "$2" != -*) ]] ;
}

get_argument() {
	echo "${2:-${1#*=}}"
}

handle_flags() {
	while [ $# -gt 0 ]
	do
		case $1 in 
			-h | --help)
				usage
				exit 0
				;;
			-v | --verbose)
				verbose=" -v"
				;;
			-d | --restdir*)
				if ! has_argument $@
				then
					echo "Replica Exchange Directory not specified." >&2
					usage
					exit 1
				fi

				restdir=$(extract_argument $@)

				shift
				;;
            -n | --name*)
				if ! has_argument $@
				then
					echo "Replica xtc Name not specified." >&2
					usage
					exit 1
				fi

				name=$(extract_argument $@)
                if [[ $name == *.xtc ]]
                then
                    temp=$(basename $name .xtc)
                    name=$temp
                fi

				shift
				;;
            -i | --index*)
                if ! has_argument $@
                then
                    echo "Must provide the replica index xvg file with -i or --index"
                    usage
                    exit 1
                fi

                index=$(extract_argument $@)

                shift
                ;;
            -l | --log*)
                if ! has_argument $@
                then
                    echo "When using the log flag you must provide a file name"
                    usage
                    exit 1
                fi

                log=$(extract_argument $@)

                shift
                ;;
			*)
				echo "Invalid option:$1" >&2
				usage
				exit 1
				;;
		esac
		shift
	done
}

demultiplex () {
    if [ ! -n "$name" ]
    then
        echo "must provide a trajectory name."
        usage
        exit 1
    fi

    if [ ! -d "$restdir" ]
    then
        echo "must provide a Replica Exchange Directory that exists."
        if [ -n "$restdir" ]
        then
            echo "Directory provided is $restdir."
        fi
        usage 
        exit 1
    fi

    if [ ! -n "$index" ]
    then
        echo "must provide the number of replicas."
        usage 
        exit 1
    fi

    case ${index} in
        8)
            gmx trjcat -f $restdir/0/$name.xtc $restdir/1/$name.xtc $restdir/2/$name.xtc $restdir/3/$name.xtc $restdir/4/$name.xtc $restdir/5/$name.xtc $restdir/6/$name.xtc $restdir/7/$name.xtc -demux $index -o 0.xtc 1.xtc 2.xtc 3.xtc 4.xtc 5.xtc 6.xtc 7.xtc            
            ;;
        10)
            gmx trjcat -f $restdir/0/$name.xtc $restdir/1/$name.xtc $restdir/2/$name.xtc $restdir/3/$name.xtc $restdir/4/$name.xtc $restdir/5/$name.xtc $restdir/6/$name.xtc $restdir/7/$name.xtc $restdir/8/$name.xtc $restdir/9/$name.xtc -demux $index -o 0.xtc 1.xtc 2.xtc 3.xtc 4.xtc 5.xtc 6.xtc 7.xtc 8.xtc 9.xtc
            ;;
        16)
            gmx trjcat -f $restdir/0/$name.xtc $restdir/1/$name.xtc $restdir/2/$name.xtc $restdir/3/$name.xtc $restdir/4/$name.xtc $restdir/5/$name.xtc $restdir/6/$name.xtc $restdir/7/$name.xtc $restdir/8/$name.xtc $restdir/9/$name.xtc $restdir/10/$name.xtc $restdir/11/$name.xtc $restdir/12/$name.xtc $restdir/13/$name.xtc $restdir/14/$name.xtc $restdir/15/$name.xtc -demux $index -o 0.xtc 1.xtc 2.xtc 3.xtc 4.xtc 5.xtc 6.xtc 7.xtc 8.xtc 9.xtc 10.xtc 11.xtc 12.xtc 13.xtc 14.xtc 15.xtc
            ;;
        20)
            gmx trjcat -f $restdir/0/$name.xtc $restdir/1/$name.xtc $restdir/2/$name.xtc $restdir/3/$name.xtc $restdir/4/$name.xtc $restdir/5/$name.xtc $restdir/6/$name.xtc $restdir/7/$name.xtc $restdir/8/$name.xtc $restdir/9/$name.xtc $restdir/10/$name.xtc $restdir/11/$name.xtc $restdir/12/$name.xtc $restdir/13/$name.xtc $restdir/14/$name.xtc $restdir/15/$name.xtc $restdir/16/$name.xtc $restdir/17/$name.xtc $restdir/18/$name.xtc $restdir/19/$name.xtc -demux $index -o 0.xtc 1.xtc 2.xtc 3.xtc 4.xtc 5.xtc 6.xtc 7.xtc 8.xtc 9.xtc 10.xtc 11.xtc 12.xtc 13.xtc 14.xtc 15.xtc 16.xtc 17.xtc 18.xtc 19.xtc
            ;;
    esac

}


# Main

handle_flags "$@"

if [ "$verbose" = true ]
then
	echo "Verbose mode enable."
    set -x
	set -v
fi

if [ -n "$log" ]
then
	echo "Log file specified: $log"
    exec > >(tee ${log}) 2>&1
fi

demultiplex

echo "Demultiplexing has completed."