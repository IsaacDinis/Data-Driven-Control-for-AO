p_geom = wao.supervisor.config.p_geom
p_tel =  wao.supervisor.config.p_tel
p_dm =  wao.supervisor.config.p_dms[1]
pixsize = wao.supervisor.config.p_geom.get_pixsize()
diam = p_tel.diam

xpos = p_dm._xpos-p_dm._n1
ypos = p_dm._ypos-p_dm._n1

i1 = p_dm._i1.copy()
j1 = p_dm._j1.copy()
influ = p_dm._influ.copy()

# place bump
xpos0 = 240.0 
ypos0 = 185.8


influ0 = p_dm._influ[:,:,0]
influ0 = np.expand_dims(influ0, axis=2)

i10 = xpos0-20.5
j10 = ypos0-20.5 

p_dm._i1
xcenter = p_geom.cent
ycenter = p_geom.cent

xpos = np.append(xpos,xpos0)
ypos = np.append(ypos,ypos0)

i1 = np.append(i1,i10)
j1 = np.append(j1,j10)

influ = np.append(influ,influ0,axis = 2)

xpos += p_dm._n1
ypos += p_dm._n1
i1 += p_dm._n1
j1 += p_dm._n1

file_name = 'bump.fits'
dm_custom = utils.write_dm_custom_fits(file_name,i1,j1,influ,xpos,ypos,xcenter,ycenter,pixsize,diam)
