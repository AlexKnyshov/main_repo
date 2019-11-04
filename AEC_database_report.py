import sys


dbname = sys.argv[1]
if len(sys.argv) == 3:
	filtername = sys.argv[2]
else:
	filtername = None

# degree_sign= u'\N{DEGREE SIGN}'
degree_sign="&#176"
def add_item(item, d):
	usi = []
	#this part multiplies records by number of specimens (in case several on same pin)
	for r in range(int(item[24])):
		if item[0] == '':	
			sys.stdout.write("error, record without USI\n")
			sys.exit()
			# usi.append("NA")
		else:
			usi.append(item[0])
	#main nested dictionary fill up
	if item[3] in d: #country
		if item[4] in d[item[3]]: #state prov
			if item[5] in d[item[3]][item[4]]: #sec prov
				if item[6] in d[item[3]][item[4]][item[5]]: #locality
					if item[7] in d[item[3]][item[4]][item[5]][item[6]]: #lat
						if item[8] in d[item[3]][item[4]][item[5]][item[6]][item[7]]: #long
							if item[9] in d[item[3]][item[4]][item[5]][item[6]][item[7]][item[8]]: #start date
								if item[10] in d[item[3]][item[4]][item[5]][item[6]][item[7]][item[8]][item[9]]: #end date
									if item[11] in d[item[3]][item[4]][item[5]][item[6]][item[7]][item[8]][item[9]][item[10]]: #collector
										if item[13] in d[item[3]][item[4]][item[5]][item[6]][item[7]][item[8]][item[9]][item[10]][item[11]]: #museum
											if item[12] in d[item[3]][item[4]][item[5]][item[6]][item[7]][item[8]][item[9]][item[10]][item[11]][item[13]]: #sex
												d[item[3]][item[4]][item[5]][item[6]][item[7]][item[8]][item[9]][item[10]][item[11]][item[13]][item[12]] += usi
											else:
												d[item[3]][item[4]][item[5]][item[6]][item[7]][item[8]][item[9]][item[10]][item[11]][item[13]][item[12]] = usi
										else:
											d[item[3]][item[4]][item[5]][item[6]][item[7]][item[8]][item[9]][item[10]][item[11]][item[13]] = {item[12] : usi}
									else:
										d[item[3]][item[4]][item[5]][item[6]][item[7]][item[8]][item[9]][item[10]][item[11]] = {item[13] : {item[12] : usi}}
								else:
									d[item[3]][item[4]][item[5]][item[6]][item[7]][item[8]][item[9]][item[10]] = {item[11] : {item[13] : {item[12] : usi}}}
							else:
								d[item[3]][item[4]][item[5]][item[6]][item[7]][item[8]][item[9]] = {item[10] : {item[11] : {item[13] : {item[12] : usi}}}}
						else:
							d[item[3]][item[4]][item[5]][item[6]][item[7]][item[8]] = {item[9] : {item[10] : {item[11] : {item[13] : {item[12] : usi}}}}}
					else:
						d[item[3]][item[4]][item[5]][item[6]][item[7]] = {item[8] : {item[9] : {item[10] : {item[11] : {item[13] : {item[12] : usi}}}}}}
				else:
					d[item[3]][item[4]][item[5]][item[6]] = {item[7] : {item[8] : {item[9] : {item[10] : {item[11] : {item[13] : {item[12] : usi}}}}}}}
			else:
				d[item[3]][item[4]][item[5]] = {item[6] : {item[7] : {item[8] : {item[9] : {item[10] : {item[11] : {item[13] : {item[12] : usi}}}}}}}}
		else:
			d[item[3]][item[4]] = {item[5] : {item[6] : {item[7] : {item[8] : {item[9] : {item[10] : {item[11] : {item[13] : {item[12] : usi}}}}}}}}}
	else:
		d[item[3]] = {item[4] : {item[5] : {item[6] : {item[7] : {item[8] : {item[9] : {item[10] : {item[11] : {item[13] : {item[12] : usi}}}}}}}}}}
	return d

def print_data(d):
	for country, val1 in sorted(d.items()):
		#print country
		msg1 = "<b>"+country+":</b> "
		sys.stdout.write(msg1)
		for province, val2 in sorted(val1.items()):
			#print province if not none
			if province != "none" and province != "None":
				msg1a = "<b>"+province+":</b> "
				sys.stdout.write(msg1a)
			for province2, val3 in sorted(val2.items()):
				#print secondary province, if not none
				if province2 != "none" and province2 != "None":
					msg2 = "<b><i>"+province2+" Co.:</b></i> "
					sys.stdout.write(msg2)
				for locality, val4 in sorted(val3.items()):
					for lat, val5 in val4.items():
						for lon, val6 in val5.items():
							#transform coordinates to positive and append degree N S E W
							if lat != "" and lon != "":
								if float(lat) >= 0:
									lat1 = lat+degree_sign+"N"
								else:
									lat1 = lat[1:]+degree_sign+"S"
								if float(lon) >= 0:
									lon1 = lon+degree_sign+"E"
								else:
									lon1 = lon[1:]+degree_sign+"W"
								msg3 = locality+", "+lat1+", "+lon1
							else:
								msg3 = locality
							#print locality
							sys.stdout.write(msg3)
							col_counter = 0
							for start_date, val7 in val6.items():
								for end_date, val8 in val7.items():
									for collector, val9 in val8.items():
										#prepare and print data and collector
										if end_date == "":
											total_date = start_date
										else:
											total_date = start_date+" - "+end_date
										if col_counter > 0:
											msg3a = "; "+total_date+", "+collector
										else:
											msg3a = ", "+total_date+", "+collector
										sys.stdout.write(msg3a)
										col_counter += 1
										museum_counter = 0
										for museum, val10 in val9.items():
											for sex, USI in sorted(val10.items(), reverse=True):
												#prepare and print specimen data grouped by museum >> sex >> USI
												if sex == "Adult Male":
													sex1 = ";m"
												elif sex == "Adult Female":
													sex1 = ";f"
												elif sex == "Adult sex unknown":
													sex1 = ";u"
												else:
													sex1 = " "+sex
												final_USI = []
												current_range = []
												sorted_USIs = sorted(range(len(USI)), key=lambda k: USI[k])
												for f in range(len(sorted_USIs)):
													if current_range == []:
														current_range.append(USI[sorted_USIs[f]])
													else:
														if int(USI[sorted_USIs[f]].split(" ")[1])-int(USI[sorted_USIs[f-1]].split(" ")[1]) > 1:
															if len(current_range) > 1:
																final_USI.append(current_range[0]+"-"+current_range[-1])
															else:
																final_USI.append(current_range[0])
															current_range = [USI[sorted_USIs[f]]]
														else:
															current_range.append(USI[sorted_USIs[f]])
												if len(current_range) == 1:
													final_USI.append(current_range[0])
												elif len(current_range) == 2:
													final_USI.append(current_range[0]+", "+current_range[1])
												else:
													final_USI.append(current_range[0]+"-"+current_range[-1])
												if museum_counter > 0:
													msg4 = ", "+str(len(USI))+sex1+" ("+", ".join(final_USI)+")"
												else:
													msg4 = ", "+str(len(USI))+sex1+" ("+", ".join(final_USI)+")"
												sys.stdout.write(msg4)
												museum_counter += 1
											#print museum after all groups output
											sys.stdout.write( " ("+museum+")")
							#print period after the entire locality output
							sys.stdout.write(". ")
	sys.stdout.write("\n")


if filtername:
	list1 = set()
	with open(filtername) as filterhandle:
		for l in filterhandle:
			list1.add(l.rstrip())

with open(dbname) as dbhandle:
	main_dict = {} # [genus_species] = {{holotype}, {paratypes}, {nontypes}}
	c = 0
	for line in dbhandle:
		line = line.rstrip().split("\t")
		if line[0] != "PBIUSI" and line[0] != '': #only parse lines with USI and not the header
			if not filtername or line[0] in list1: #only parse lines if no filter file given or if USI is in filter
				#create new species
				if line[1]+"_"+line[2] not in main_dict:
					main_dict[line[1]+"_"+line[2]] = {"Holotype":{},"Paratype":{},"Nontype":{}}
				#populate holotype data
				if line[15] == "Holotype":
					main_dict[line[1]+"_"+line[2]]["Holotype"] = add_item(line, main_dict[line[1]+"_"+line[2]]["Holotype"])
				#populate paratype data
				elif line[15] == "Paratype":
					main_dict[line[1]+"_"+line[2]]["Paratype"] = add_item(line, main_dict[line[1]+"_"+line[2]]["Paratype"])
				#populate other data
				else:
					main_dict[line[1]+"_"+line[2]]["Nontype"] = add_item(line, main_dict[line[1]+"_"+line[2]]["Nontype"])
				c+=1 
sys.stdout.write("<p>Total records reported: "+str(c)+"\n")
# print output
for species, dat in sorted(main_dict.items()):
	sys.stdout.write("<h3>"+species+"</h3>\n")
	sys.stdout.write("<p>Holotype: ")
	print_data(main_dict[species]["Holotype"])
	sys.stdout.write("<p>Paratypes: ")
	print_data(main_dict[species]["Paratype"])
	sys.stdout.write("<p>Other Specimens Examined: ")
	print_data(main_dict[species]["Nontype"])