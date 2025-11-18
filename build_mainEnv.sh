#!/bin/bash

# Setting the database directory
Main_dir="$( dirname -- "$( readlink -f -- "$0"; )"; )"
aux_dir="${Main_dir}/aux/"

echo "Donwloading MACSE..."
wget -O ${aux_dir}/macse.jar "https://www.agap-ge2pop.org/wp-content/uploads/macse/releases/macse_v2.07.jar"


