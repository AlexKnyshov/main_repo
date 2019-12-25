for f in *.fa*; do echo $f; awk -F'|' '{if ($1 ~ ">") {print ">"$2} else {print $1}}' $f > $f.fas; done
