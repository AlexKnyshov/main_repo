from datetime import datetime
import sys

inf1 = sys.argv[1]

with open(inf1) as infh:
	prevdate = None
	totaltime = 0
	for l in infh:
		datetime_object = datetime.strptime(l.strip(), '%a %b %d %H:%M:%S %Z %Y')
		if prevdate:
			timedif = (datetime_object-prevdate).total_seconds() / 60
			print datetime_object, timedif
			prevdate = datetime_object
		else:
			timedif = 0
			print datetime_object, 0
			prevdate = datetime_object
		totaltime += timedif
print totaltime