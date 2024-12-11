#!/bin/bash

# Default variables

verbose_mode=false
output_file=""


# Function to display script usage

usage() {
	echo "Usage: $0 [Options]"
	echo "Options:"
	echo " -h, --help	Display this help message"
	echo " -v, --verbose	Enable verbosity"
	echo " -o, --ofile	Output FILE"
	echo " -f, --ifile	Input FILE"
}

has_argument() {
	[[ ("$1" == *=* && -n ${1#*=}) || ( ! -z "$2" && "$2" != -*) ]] ;
}

get_argument() {
	echo "${2:-${1#*=}}"
}

# Flag handler
handle_flags() {
	while [ $# -gt 0 ]
	do
		case $1 in 
			-h | --help)
				usage
				exit 0
				;;
			-v | --verbose)
				verbose_mode=true
				;;
			-f | --file*)
				if ! has_argument $@
				then
					echo "File not specified." >&2
					usage
					exit 1
				fi

				output_file=$(extract_argument $@)

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


# Main

handle_flags "$@"


# Perform the desired actions

if [ "$verbose_mode" = true ]
then
	echo "Verbose mode enable."
fi

if [ -n "$output_file" ]
then
	echo "Output file specified: $output_file"
fi

