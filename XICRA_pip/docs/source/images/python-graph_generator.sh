for i in `dir ../../../XICRA/scripts/`; do if [[ $i =~ .*__.* ]]; then continue; else if [[ $i =~ .*pyc.* ]]; then continue; else echo "##############"; echo "##### Create pyan graph for file: "$i; echo "####################"; name=(${i//.py/}); `pyan --dot-rankdir TB --make-svg -f $name".dot" -n -u -c -G -e --focus=$name --only-child ../../../XICRA/*/*py`; fi; fi; done


