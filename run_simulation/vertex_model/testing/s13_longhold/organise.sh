#!/bin/bash

# Create the backup directory if it doesn't exist
mkdir -p prevSol

# Move solution and data files to prevSol
mv solution* prevSol/ 2>/dev/null || true
mv TimeData.txt prevSol/ 2>/dev/null || true
mv AvgStressStrainData.txt prevSol/ 2>/dev/null || true
mv *.err prevSol/ 2>/dev/null || true
mv *.out prevSol/ 2>/dev/null || true
mv restart_info* prevSol/ 2>/dev/null || true
mv elasticener.txt prevSol/ 2>/dev/null || true
mv Delelasticener.txt prevSol/ 2>/dev/null || true
mv activeener.txt prevSol/ 2>/dev/null || true
mv VMconfig.cfg prevSol/ 2>/dev/null || true

# Copy job script and move the restart script
cp VMrunjob.sh prevSol/ 2>/dev/null || true
mv VMrunjobRestart.sh prevSol/ 2>/dev/null || true

# Display the remaining contents sorted by modification time
echo "Files remaining in the current directory:"
ls -ltr
