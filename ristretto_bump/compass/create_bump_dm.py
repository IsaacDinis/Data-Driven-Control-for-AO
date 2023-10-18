p_geom = wao.supervisor.config.p_geom
p_tel =  wao.supervisor.config.p_tel
p_dm =  wao.supervisor.config.p_dms[1]
pixsize = wao.supervisor.config.p_geom.get_pixsize()
diam = p_tel.diam

xpos = p_dm._xpos-p_dm._n1
ypos = p_dm._ypos-p_dm._n1

# place bump
xpos0 = 240.0 
ypos0 = 185.8

# check bump position
plt.scatter(xpos, ypos, marker='.', color="red")
plt.scatter(xpos0, ypos0, marker='.', color="green")


influ = p_dm._influ
influ0 = p_dm._influ[:,:,0]
influ0 = np.expand_dims(influ0, axis=2)

i10 = xpos0-p_dm._xpos[0]-p_dm._n1 
j10 = ypos0-p_dm._xpos[0]-p_dm._n1

xcenter = p_geom.cent
ycenter = p_geom.cent

xpos0 += p_dm._n1
ypos0 += p_dm._n1
i10 += p_dm._n1
j10 += p_dm._n1

file_name = 'bump.fits'
dm_custom = utils.write_dm_custom_fits(file_name,i10,j10,influ0,xpos0,ypos0,xcenter,ycenter,pixsize,diam)