#!/bin/bash

copy_modify_and_rename(){

  local source_directory="$1"
  local destination_directory="$2"
  local old_string_in_filename="$3"
  local new_string_in_filename="$4"
  local old_string="$5"
  local new_string="$6"
  
  # remove all old files from directory
  rm "$destination_directory"/*.yaml
  
  # Create the destination directory if it doesn't exist
  mkdir -p "$destination_directory"
	
# Loop through each file in the source directory
  for file in "$source_directory"/*
  do
    # Check if it's a file
    if [ -f "$file" ]; then
      # Extract filename and extension
      filename=$(basename -- "$file")
      extension="${filename##*.}"
      filename="${filename%.*}"

      # Replace string in filename
      new_filename="${filename//$old_string_in_filename/$new_string_in_filename}.$extension"

      # Copy the file to the destination directory
      cp -f "$file" "$destination_directory/$new_filename"

      # Change the specific line in the copied file that starts with a specific string
      sed -i "s/$old_string/$new_string/g" "$destination_directory/$new_filename"
    fi
  done
}


copy_modify_and_rename "ljsve" "lj" "ljsve" "lj" "Lennard-Jones SVE" "Lennard-Jones"
copy_modify_and_rename "ljsve" "miesve" "ljsve" "miesve" "Lennard-Jones SVE" "Mie SVE"
copy_modify_and_rename "ljsve" "mie" "ljsve" "mie" "Lennard-Jones SVE" "Mie"
copy_modify_and_rename "ljsve" "miesvefixed" "ljsve" "miesvefixed" "Lennard-Jones SVE" "Mie SVE fixed"
copy_modify_and_rename "ljsve" "miefixed" "ljsve" "miefixed" "Lennard-Jones SVE" "Mie fixed"


