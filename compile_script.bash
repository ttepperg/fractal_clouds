#!/bin/bash

# Run this script in order to compile and generate an executable 'fractalClouds'on the project's root directory.
# NOTE: Make sure that 'gcc' is linked to the appropriate C compiler, or replace 'gcc' below by the appropriate C compiler.


# Remove current executable (if present) for a fresh start
echo -e "\nRemoving current executable (if present) for a fresh start\n"
rm -f fractalClouds

# Cd into source directory
echo -e "Changing into source directory\n"
cd src/.

# Compile, link, and create executable
echo -e "Compiling..."
gcc dnrutil.c dlubksb.c dludcmp.c dsavgol.c dfour1.c dfourn.c dran1.c dgasdev.c dgammln.c dgcf.c dgser.c dgammq.c dfit.c dlfit.c dcovsrt.c dgaussj.c alloc_funcs.c fractalClouds.c -o fractalClouds
echo -e "Done.\n"

# Move executable to project's root directory
echo -e "Moving executable 'fractalClouds' to project's root directory\n"
mv fractalClouds ../.

# Cd into project's root directory
echo -e "Returning to project's root directory\n"
cd ../.
echo -e "Done.\n"

echo -e "Run the executable by invoking in the shell:\n\n?>./fractalClouds\n"

exit

