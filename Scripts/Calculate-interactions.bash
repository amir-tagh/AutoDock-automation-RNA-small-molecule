#!/bin/bash

here=$(pwd)

for i in {1..7}
 do
  echo "i is $i"	  
   #mkdir -p Interactions_pose${i}
    #cd Interactions_pose${i}
     #cp $here/CBR-pose_state0${i}.pdb .
      #cp $here/*.py .
       cd Docking-G1-active${i}
        now=$(pwd)
	 echo "we are at $now"
	  sleep 2
	   cp $here/extract_and_convert_lowest_model.py .
            python extract_and_convert_lowest_model.py *.dlg output_directory
             cd output_directory
              now_2=$(pwd)
	       echo "we are at $now_2"
                sleep 2
          cp $here/RNA.pdb .
           cp lowest_energy_model.pdb ligand.pdb
            #sed -i '/LIG/d' RNA.pdb
	     sed -i '/CONECT/d' RNA.pdb
              sed -i '/LIG/!d' ligand.pdb
	       obabel -ipdb ligand.pdb -osdf -O ligand.sdf
	       
             
    ~/Downloads/scripts/fingeRNAt.py -r RNA.pdb -l ligand.sdf -detail -print -verbose
    
 cd $here
done
