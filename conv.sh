#!/bin/bash

# Directory containing the files (default to the current directory)
TARGET_DIR=${1:-./src/INTERF}

# Check if the directory exists
if [[ ! -d "$TARGET_DIR" ]]; then
    echo "Error: Directory $TARGET_DIR does not exist."
    exit 1
fi

# Iterate over all .f files in the target directory
for file in "$TARGET_DIR"/*.f; do
    # Skip if no .f files are found
    if [[ ! -e "$file" ]]; then
        echo "No .f files found in $TARGET_DIR."
        exit 0
    fi

    # Extract the base name of the file (without extension)
    base_name=$(basename "$file" .f)

    # Define the new file name with .f90 extension
    new_file="$TARGET_DIR/$base_name.f90"

    # Rename the file
    mv "$file" "$new_file"

    echo "Renamed $file to $new_file"
done

echo "All .f files have been converted to .f90."
