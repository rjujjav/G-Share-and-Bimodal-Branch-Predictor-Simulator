# Initial values for ./sim
iterations=14
K=0
M1=0
N=0
M2=5
input_file="gcc_trace.txt"

# Name of the output file to store misprediction rates
output_file="gcc_misprediction_rates.txt"

# Initialize the output file
echo -n > "$output_file"

# Loop to run the command multiple times
for i in $(seq 1 $iterations)
do
  output=$(./sim bimodal $M2 $input_file)
  misprediction_rate=$(echo "$output" | grep -oP 'misprediction rate: \K\d+\.\d+')
  echo "M2=$M2 misprediction rate: $misprediction_rate" >> "$output_file"
  # Increment the 'M2' value by 1 for the next iteration
  M2=$((M2+1))
done