#bin/bash
for E in {250..5001..250}
do
x=$(( 100000+$E*100 ))
nuwro -i params.txt -p "beam_energy=$E" -p "beam_particle = 14" -p "number_of_test_events = $x" -p "mec_scaling = 0" -o nieves0numu_$E.root
nuwro -i params.txt -p "beam_energy=$E" -p "beam_particle = 14" -p "number_of_test_events = $x" -p "mec_scaling = 1" -o nieves1numu_$E.root
nuwro -i params.txt -p "beam_energy=$E" -p "beam_particle = -14" -p "number_of_test_events = $x" -p "mec_scaling = 0" -o nieves0anumu_$E.root
nuwro -i params.txt -p "beam_energy=$E" -p "beam_particle = -14" -p "number_of_test_events = $x" -p "mec_scaling = 1" -o nieves1anumu_$E.root
done

