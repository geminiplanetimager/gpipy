'''
Generate FoV map of GPI, and plot path of companion as a function of hour angle.
Requires gpi scheduler source code

Robert De Rosa
v0.1 -- Initial version
v0.2 -- Added dashed line for offsets tested in VC2 and VC3
v0.3 -- Rotated to IFS coordinate frame. fov_offx and fov_offy should be the command issued to pnc
v0.4 -- Fixed bug with HA = -24 when LST = RA for Northern targets
v0.5 -- Converted to function. fov.plot(...)
'''

from helperFun import *   #also imports numpy as np
import datetime
import matplotlib.pyplot as pyplot

def plot(targ, targ_ra, targ_de, comp_rho, comp_theta, offset_x, offset_y, hms=False, deg=False, band='H'):

	np.set_printoptions(threshold=np.nan)

	#Gemini South coordinates
	lat = '-30:14:26.700'
	lon = '-70:44:12.096'

	#Convert latitute/longitude to decimal/radians
	lat_deg = dms2dd(lat)
	lon_deg = dms2dd(lon)
	lat_rad = lat_deg * pi/180.0
	lon_rad = lon_deg * pi/180.0

	#Define IFS angle
	ifs_ang = (24.5 * (pi/180.0))

	#Convert from hh:mm:ss.ss to decimal hours and degrees
	if (hms is True):
		targ_ra_hrs = dms2dd(targ_ra)
		targ_ra_deg = targ_ra_hrs * 15.0
		targ_de_deg = dms2dd(targ_de)
		targ_ra_rad = targ_ra_deg * pi/180.0
		targ_de_rad = targ_de_deg * pi/180.0

	if (deg is True):
		targ_ra_hrs = targ_ra / 15.0
		targ_ra_deg = targ_ra
		targ_de_deg = targ_de
		targ_ra_rad = targ_ra_deg * pi/180.0
		targ_de_rad = targ_de_deg * pi/180.0

	if (((hms is False) and (deg is False)) or ((hms is True) and (deg is True))):
		print 'Specify either hms or decimal'
		return

	#Convert companion position angle to radians
	comp_theta_rad = comp_theta * (pi/180.0)

	#Construct hour angle array
	min_ha = -6.0
	max_ha = 6.01
	ha_arr = np.arange(min_ha,max_ha,0.0001,dtype=np.float32)
	ha_arr_rad = (ha_arr * 15.0) * (pi/180.0)

	#Calculate parallactic angle
	pa_arr = (-1.0) * np.arctan2((-np.cos(lat_rad)*np.sin(ha_arr_rad)),(np.sin(lat_rad)*np.cos(targ_de_rad) - np.cos(lat_rad)*np.cos(ha_arr_rad)*np.sin(targ_de_rad)))
	
	#Wrap if declination greater than latitutde
	if targ_de_deg > lat_deg:
		pa_arr = (((pa_arr + (2.0*pi)) % (2.0*pi)))

	#Rotate PA to IFS coordinate frame
	pa_arr += ifs_ang

	#Define GPI FoV
	fov_x = np.array([-1.33,  1.33,  1.33, -1.33, -1.33])
	fov_y = np.array([-1.33, -1.33,  1.33,  1.33, -1.33])

	#Rotate to IFS coordinate frame
	fov_x_rot = (fov_x*np.cos(ifs_ang) - fov_y*np.sin(ifs_ang)) + offset_x
	fov_y_rot = (fov_x*np.sin(ifs_ang) + fov_y*np.cos(ifs_ang)) + offset_y

	#Convert companion rho/theta to x, y.. x ascending to left
	comp_x  = comp_rho * np.cos(((comp_theta_rad - pa_arr) + 1.0*pi) % (2.0*pi))
	comp_y  = comp_rho * np.sin(((comp_theta_rad - pa_arr) + 1.0*pi) % (2.0*pi))
	comp_x *= (-1.0)

	#txt_rho, lin_rho specifies radial distance of text label, and end of line connecting it to curve
	txt_rho = 2.00
	lin_rho = 1.75

	#Generate x,y coordinates for these
	txt_x  = txt_rho * np.cos(((comp_theta_rad - pa_arr) + 1.0*pi) % (2.0*pi))
	txt_y  = txt_rho * np.sin(((comp_theta_rad - pa_arr) + 1.0*pi) % (2.0*pi))
	txt_x *= (-1.0)
	lin_x  = lin_rho * np.cos(((comp_theta_rad - pa_arr) + 1.0*pi) % (2.0*pi))
	lin_y  = lin_rho * np.sin(((comp_theta_rad - pa_arr) + 1.0*pi) % (2.0*pi))
	lin_x *= (-1.0)

	#Define plot
	fig = pyplot.figure(figsize=[10,10])
	pyplot.xlim([-3.0,3.0])
	pyplot.ylim([-3.0,3.0])
	pyplot.xlabel('Offset from star (arc sec)')
	pyplot.ylabel('Offset from star (arc sec)')
	pyplot.plot(0.0, 0.0, '*', markersize=20., color='blue')
	pyplot.plot(offset_x, offset_y, '+', markersize=5, color='black')
	pyplot.plot(fov_x_rot, fov_y_rot, color='black')
	pyplot.plot(comp_x, comp_y, color='red')

	#Box for allowed offsets
	xbox = (np.array([-500,-500,300,300,-500],dtype=np.float32) * 0.00141) + offset_x
	ybox = (np.array([-900,300,300,-900,-900],dtype=np.float32) * 0.00141) + offset_y
	pyplot.plot(xbox,ybox, color='blue', linestyle='dashed')

	#Plot path of satellite spots eminating outwards from 45,135,225,315 - 3.5 deg
	spot_theta = ((np.array([45.0,135.0,225.0,315.0])-3.5)*(np.pi/180.0))
	spot_theta += ifs_ang

	if band == 'Y':
		spot_start = ((20*0.947108e-6)/7.8)*206264.806
		spot_end = ((20*1.14209e-6)/7.8)*206264.806
	if band == 'J':
		spot_start = ((20*1.11407e-6)/7.8)*206264.806
		spot_end = ((20*1.34973e-6)/7.8)*206264.806
	if band == 'H':
		spot_start = ((20*1.49461e-6)/7.8)*206264.806
		spot_end = ((20*1.79739e-6)/7.8)*206264.806
	if band == 'K1':
		spot_start = ((20*1.88609e-6)/7.8)*206264.806
		spot_end = ((20*2.19511e-6)/7.8)*206264.806
	if band == 'K2':
		spot_start = ((20*2.10741e-6)/7.8)*206264.806
		spot_end = ((20*2.39639e-6)/7.8)*206264.806

	spot_x_start = spot_start * np.cos(spot_theta)
	spot_y_start = spot_start * np.sin(spot_theta)

	spot_x_end = spot_end * np.cos(spot_theta)
	spot_y_end = spot_end * np.sin(spot_theta)

	#Spots are at 20lambda/D

	for i in xrange(0,4):
		pyplot.plot([spot_x_start[i],spot_x_end[i]],[spot_y_start[i],spot_y_end[i]],color='black')

	d_ha = 0.25 #loop through every 15 mins
	ha_vals = np.arange(min_ha,max_ha,d_ha,dtype=np.float32)

	for i in xrange(0,len(ha_vals)):
		sym_col = 'red'
		sym_size = 5.0

		if (ha_vals[i] % 1) == 0:
			sym_size = 7.5

		if ha_vals[i] == 0:
			sym_col = 'black'
			sym_size = 10.0

		ind = np.abs (ha_arr - ha_vals[i]).argmin()
		pyplot.plot(comp_x[ind], comp_y[ind], 'o', color=sym_col, markersize=sym_size)

		if ((ha_vals[i] >= -1.0) and (ha_vals[i] <= 1.0)) or (((ha_vals[i] < -1.0) or (ha_vals[i] > 1.0)) and ((ha_vals[i] % 1) == 0)):
			ang = np.arctan2(txt_y[ind],txt_x[ind]) * (180.0/pi)
			pyplot.plot([comp_x[ind],lin_x[ind]],[comp_y[ind],lin_y[ind]],color='grey', linestyle='--')
			pyplot.text(txt_x[ind],txt_y[ind],str(ha_vals[i])+'h', rotation = ang, ha='center', va='center')

	pyplot.arrow(-2.5,2.0,0.5,0.0, head_width=0.1, color='black')
	pyplot.arrow(-2.5,2.0,0.0,0.5, head_width=0.1, color='black')
	pyplot.text(-2.5+0.7,2.0-0.075,'X', fontsize=20)
	pyplot.text(-2.5-0.075,2.0+0.7,'Y', fontsize=20)

	if targ_de_deg < lat_deg:
		xa = -1.30
		ya = 2.55
		arrow_x = 0.5
		arrow_y = 0.0
		xb = (arrow_x * np.cos(ifs_ang) - arrow_y * np.sin(ifs_ang)) + xa
		yb = (arrow_x * np.sin(ifs_ang) + arrow_y * np.cos(ifs_ang)) + ya
		pyplot.arrow(xa,ya,xb-xa,yb-ya, head_width=0.1, color='black')
		pyplot.text(xb+0.20,yb+0.0,'N', fontsize=20)

		arrow_x = 0.0
		arrow_y = -0.5
		xb = (arrow_x * np.cos(ifs_ang) - arrow_y * np.sin(ifs_ang)) + xa
		yb = (arrow_x * np.sin(ifs_ang) + arrow_y * np.cos(ifs_ang)) + ya
		pyplot.arrow(xa,ya,xb-xa,yb-ya, head_width=0.1, color='black')
		pyplot.text(xb+0.10,yb-0.3,'E', fontsize=20)

	else:
		xa = -0.5
		ya = 2.20
		arrow_x = -0.5
		arrow_y = 0.0
		xb = (arrow_x * np.cos(ifs_ang) - arrow_y * np.sin(ifs_ang)) + xa
		yb = (arrow_x * np.sin(ifs_ang) + arrow_y * np.cos(ifs_ang)) + ya
		pyplot.arrow(xa,ya,xb-xa,yb-ya, head_width=0.1, color='black')
		pyplot.text(xb-0.30,yb-0.1,'N', fontsize=20)

		xa = -0.5
		ya = 2.20
		arrow_x = 0.0
		arrow_y = 0.5
		xb = (arrow_x * np.cos(ifs_ang) - arrow_y * np.sin(ifs_ang)) + xa
		yb = (arrow_x * np.sin(ifs_ang) + arrow_y * np.cos(ifs_ang)) + ya
		pyplot.arrow(xa,ya,xb-xa,yb-ya, head_width=0.1, color='black')
		pyplot.text(xb-0.15,yb+0.15,'E', fontsize=20)

	pyplot.text(-1.7,1.6,'at Hour Angle = 0.0',fontsize=16)
	pyplot.text(1.5,3.35,'x offset = '+str(offset_x)+'asec',fontsize=16,horizontalalignment='left')
	pyplot.text(1.5,3.15,'y offset = '+str(offset_y)+'asec',fontsize=16,horizontalalignment='left')
	pyplot.text(0.00,3.25,targ,fontsize=25,horizontalalignment='center')

	fig.savefig("gpi_fov.png",dpi=150)
