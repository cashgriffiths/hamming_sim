START=$1
END=$2
N=$3
HS=$4

if  [ "$HS" != "hard" ] && [ "$HS" != "soft" ]
then
    echo "Error: Method must be either hard or soft. Got \"$HS\"."
    exit
fi

for SNR in $(seq $START $END)
do
    echo "SNR: $SNR dB"
    time ./a.out $SNR $N $HS | grep ^b | cat > /dev/null
    echo
done
