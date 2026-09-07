#!/bin/bash -e

./run1.bash 0 ""			# avgq=0.8981     avgtc=0.6627

./run1.bash 1 "-m_is 0.024"	# avgq=0.9017     avgtc=0.6704
./run1.bash 2 "-m_is 0.006" # avgq=0.8954     avgtc=0.6539

./run1.bash 3 "-m_il 0.016" # avgq=0.8979     avgtc=0.6621
./run1.bash 4 "-m_il 0.004" # avgq=0.8990     avgtc=0.6630

./run1.bash 5 "-is_is 0.50"	# avgq=0.8991     avgtc=0.6644
./run1.bash 6 "-is_is 0.20" # avgq=0.8983     avgtc=0.6648

./run1.bash 7 "-il_il 0.95" # avgq=0.8980     avgtc=0.6626
./run1.bash 8 "-il_il 0.80" # avgq=0.8984     avgtc=0.6628

./run1.bash 9 "-m_is 0.024 -m_il 0.016 -is_is 0.50 -il_il 0.95" # avgq=0.9018 avgtc=0.6715

./run1.bash 11  "-m_is 0.030" # avgq=0.9022     avgtc=0.6707
./run1.bash 12  "-m_is 0.036" # avgq=0.9030     avgtc=0.6706
./run1.bash 13  "-m_is 0.048" # avgq=0.9042     avgtc=0.6774
./run1.bash 14  "-m_is 0.024 -is_is 0.50" # avgq=0.9015     avgtc=0.6694
./run1.bash 15  "-m_is 0.024 -is_is 0.65" # avgq=0.9019     avgtc=0.6705
./run1.bash 16  "-m_is 0.036 -is_is 0.50" # avgq=0.9036     avgtc=0.6743
./run1.bash 17  "-m_is 0.024 -m_il 0.012" # avgq=0.9021     avgtc=0.6717
./run1.bash 18  "-m_is 0.036 -m_il 0.012 -is_is 0.50" # avgq=0.9047     avgtc=0.6774
./run1.bash 19  "-m_is 0.024 -s_is 0.04" #  avgq=0.9020     avgtc=0.6704
./run1.bash 20  "-m_is 0.036 -m_il 0.012 -is_is 0.50 -il_il 0.95" # avgq=0.9048     avgtc=0.6771

./run1.bash 21  "-m_is 0.060"
./run1.bash 22  "-m_is 0.072"
./run1.bash 23  "-m_is 0.096"
./run1.bash 24  "-m_is 0.048 -m_il 0.012 -is_is 0.50"
./run1.bash 25  "-m_is 0.060 -m_il 0.012 -is_is 0.50"
./run1.bash 26  "-m_is 0.072 -m_il 0.012 -is_is 0.50"
./run1.bash 27  "-m_is 0.048 -m_il 0.016 -is_is 0.50"
./run1.bash 28  "-m_is 0.060 -m_il 0.016 -is_is 0.50"
./run1.bash 29  "-m_is 0.048 -m_il 0.012 -is_is 0.65"
./run1.bash 30  "-m_is 0.060 -m_il 0.012 -is_is 0.50 -il_il 0.95"


./run1.bash 31  "-m_is 0.048 -m_il 0.012 -is_is 0.70"
./run1.bash 32  "-m_is 0.048 -m_il 0.012 -is_is 0.75"
./run1.bash 33  "-m_is 0.048 -m_il 0.012 -is_is 0.80"
./run1.bash 34  "-m_is 0.036 -m_il 0.012 -is_is 0.65"
./run1.bash 35  "-m_is 0.060 -m_il 0.012 -is_is 0.65"
./run1.bash 36  "-m_is 0.048 -m_il 0.008 -is_is 0.65"
./run1.bash 37  "-m_is 0.048 -m_il 0.016 -is_is 0.65"
./run1.bash 38  "-m_is 0.060 -m_il 0.016 -is_is 0.65"
./run1.bash 39  "-m_is 0.048 -m_il 0.012 -is_is 0.65 -il_il 0.95"
./run1.bash 40  "-m_is 0.048 -m_il 0.012 -is_is 0.70 -il_il 0.95"

./run1.bash 31  "-m_is 0.048 -m_il 0.012 -is_is 0.70" # avgq=0.9062     avgtc=0.6816
./run1.bash 32  "-m_is 0.048 -m_il 0.012 -is_is 0.75" # avgq=0.9063     avgtc=0.6817
./run1.bash 33  "-m_is 0.048 -m_il 0.012 -is_is 0.80" # avgq=0.9064     avgtc=0.6823
./run1.bash 34  "-m_is 0.036 -m_il 0.012 -is_is 0.65" # avgq=0.9050     avgtc=0.6772
./run1.bash 35  "-m_is 0.060 -m_il 0.012 -is_is 0.65" # avgq=0.9063     avgtc=0.6818
./run1.bash 36  "-m_is 0.048 -m_il 0.008 -is_is 0.65" # avgq=0.9060     avgtc=0.6793
./run1.bash 37  "-m_is 0.048 -m_il 0.016 -is_is 0.65" # avgq=0.9054     avgtc=0.6805
./run1.bash 38  "-m_is 0.060 -m_il 0.016 -is_is 0.65" # avgq=0.9051     avgtc=0.6791
./run1.bash 39  "-m_is 0.048 -m_il 0.012 -is_is 0.65 -il_il 0.95" # avgq=0.9063     avgtc=0.6801
./run1.bash 40  "-m_is 0.048 -m_il 0.012 -is_is 0.70 -il_il 0.95" # avgq=0.9063     avgtc=0.6799

./run1.bash 41  "-m_is 0.048 -m_il 0.012 -is_is 0.85"
./run1.bash 42  "-m_is 0.048 -m_il 0.012 -is_is 0.90"
./run1.bash 43  "-m_is 0.048 -m_il 0.012 -is_is 0.95"
./run1.bash 44  "-m_is 0.048 -m_il 0.012 -is_is 1.20"
./run1.bash 45  "-m_is 0.060 -m_il 0.012 -is_is 0.80"
./run1.bash 46  "-m_is 0.060 -m_il 0.012 -is_is 0.90"
./run1.bash 47  "-m_is 0.042 -m_il 0.012 -is_is 0.80"
./run1.bash 48  "-m_is 0.054 -m_il 0.012 -is_is 0.80"
./run1.bash 49  "-m_is 0.048 -m_il 0.010 -is_is 0.80"
./run1.bash 50  "-m_is 0.048 -m_il 0.014 -is_is 0.80"

./run1.bash 41  "-m_is 0.048 -m_il 0.012 -is_is 0.85" # avgq=0.9058     avgtc=0.6812
./run1.bash 42  "-m_is 0.048 -m_il 0.012 -is_is 0.90" # avgq=0.9056     avgtc=0.6791
./run1.bash 43  "-m_is 0.048 -m_il 0.012 -is_is 0.95" # avgq=0.9056     avgtc=0.6783
./run1.bash 44  "-m_is 0.048 -m_il 0.012 -is_is 1.20" # avgq=0.9050     avgtc=0.6780
./run1.bash 45  "-m_is 0.060 -m_il 0.012 -is_is 0.80" # avgq=0.9063     avgtc=0.6813
./run1.bash 46  "-m_is 0.060 -m_il 0.012 -is_is 0.90" # avgq=0.9067     avgtc=0.6844
./run1.bash 47  "-m_is 0.042 -m_il 0.012 -is_is 0.80" # avgq=0.9057     avgtc=0.6786
./run1.bash 48  "-m_is 0.054 -m_il 0.012 -is_is 0.80" # avgq=0.9065     avgtc=0.6828
./run1.bash 49  "-m_is 0.048 -m_il 0.010 -is_is 0.80" # avgq=0.9063     avgtc=0.6823
./run1.bash 50  "-m_is 0.048 -m_il 0.014 -is_is 0.80" # avgq=0.9063     avgtc=0.6804

./run1.bash 51  "-m_is 0.060 -m_il 0.012 -is_is 0.85" # avgq=0.9067     avgtc=0.6844
./run1.bash 52  "-m_is 0.060 -m_il 0.012 -is_is 0.88" # avgq=0.9069     avgtc=0.6848
./run1.bash 53  "-m_is 0.060 -m_il 0.012 -is_is 0.92" # avgq=0.9066     avgtc=0.6844
./run1.bash 54  "-m_is 0.060 -m_il 0.012 -is_is 0.95" # avgq=0.9065     avgtc=0.6828
./run1.bash 55  "-m_is 0.054 -m_il 0.012 -is_is 0.90" # avgq=0.9063     avgtc=0.6816
./run1.bash 56  "-m_is 0.066 -m_il 0.012 -is_is 0.90" # avgq=0.9057     avgtc=0.6790
./run1.bash 57  "-m_is 0.072 -m_il 0.012 -is_is 0.90" # avgq=0.9059     avgtc=0.6788
./run1.bash 58  "-m_is 0.060 -m_il 0.010 -is_is 0.90" # avgq=0.9069     avgtc=0.6846
./run1.bash 59  "-m_is 0.060 -m_il 0.014 -is_is 0.90" # avgq=0.9064     avgtc=0.6826
./run1.bash 60  "-m_is 0.066 -m_il 0.012 -is_is 0.95" # avgq=0.9059     avgtc=0.6811

./run1.bash 61  "-m_is 0.060 -m_il 0.012 -is_is 0.86" # avgq=0.9068     avgtc=0.6846
./run1.bash 62  "-m_is 0.060 -m_il 0.012 -is_is 0.87" # avgq=0.9068     avgtc=0.6848
./run1.bash 63  "-m_is 0.060 -m_il 0.012 -is_is 0.89" # avgq=0.9069     avgtc=0.6849 <<
./run1.bash 64  "-m_is 0.060 -m_il 0.010 -is_is 0.88" # avgq=0.9069     avgtc=0.6846
./run1.bash 65  "-m_is 0.060 -m_il 0.011 -is_is 0.88" # avgq=0.9069     avgtc=0.6848
./run1.bash 66  "-m_is 0.057 -m_il 0.012 -is_is 0.88" # avgq=0.9069     avgtc=0.6847
./run1.bash 67  "-m_is 0.063 -m_il 0.012 -is_is 0.88" # avgq=0.9060     avgtc=0.6814
./run1.bash 68  "-m_is 0.057 -m_il 0.010 -is_is 0.88" # avgq=0.9065     avgtc=0.6830
./run1.bash 69  "-m_is 0.060 -m_il 0.010 -is_is 0.87" # avgq=0.9068     avgtc=0.6843
./run1.bash 70  "-m_is 0.060 -m_il 0.010 -is_is 0.89" # avgq=0.9069     avgtc=0.6846

./run1.bash 71  "-m_is 0.060 -m_il 0.012 -is_is 0.89"	# avgq=0.9069     avgtc=0.6849
./run1.bash 72  "-m_is 0.060 -m_il 0.012 -is_is 0.89 -il_il 0.92"	# avgq=0.9070     avgtc=0.6849
./run1.bash 73  "-m_is 0.060 -m_il 0.012 -is_is 0.89 -il_il 0.85"	# avgq=0.9066     avgtc=0.6831
./run1.bash 74  "-m_is 0.060 -m_il 0.012 -is_is 0.89 -s_is 0.03"	# avgq=0.9066     avgtc=0.683
./run1.bash 75  "-m_is 0.060 -m_il 0.012 -is_is 0.89 -s_il 0.12"	# avgq=0.9066     avgtc=0.6840
./run1.bash 76  "-m_is 0.060 -m_il 0.012 -is_is 0.89 -s_il 0.24"	# avgq=0.9071     avgtc=0.6853 <<
./run1.bash 77  "-m_is 0.060 -m_il 0.011 -is_is 0.89"	# avgq=0.9069     avgtc=0.6848
./run1.bash 78  "-m_is 0.060 -m_il 0.012 -is_is 0.885"	# avgq=0.9070     avgtc=0.6849
./run1.bash 79  "-m_is 0.059 -m_il 0.012 -is_is 0.89"	# avgq=0.9070     avgtc=0.6848
./run1.bash 80  "-m_is 0.061 -m_il 0.012 -is_is 0.89"	# avgq=0.9064     avgtc=0.6835
