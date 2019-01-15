#!/bin/bash
while IFS='' read -r line || [[ -n "$line" ]]; do
    	jupyter-nbconvert $line --to markdown
	echo "Text read from file: $line"
done < "$1"
