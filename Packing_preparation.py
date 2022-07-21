from yade import pack, plot,qt  # Import standard Yade package, see https://yade-dev.gitlab.io/trunk/

import matplotlib; matplotlib.rc('axes',grid=True) # Import Matlab external package
import pylab
# Parameters
nRead=readParamsFromTable(
	compFricDegree = 25, # contact friction during the confining phase
	key='_sample_genesis_', # simulation's name
	unknownOk=True
)
from yade.params import table
key              = table.key
compFricDegree   = table.compFricDegree 
unbalancedThread = 1e-3
triaxThread      = 1e-5
damp             = 0.2 # damping coefficient
young            = 2.5e9 # contact stiffness, Young's Modulus
seed             = 1
x=z              = 0.3*0.001 # box size
y                = 1.2*0.001
mn,mx            = Vector3(0,0,0), Vector3(x,y,z) # box
folderName       = 'compact_Young_2.5GPa_Fric_25_largeCompact'

O.materials.append(FrictMat(young=young,poisson=0.33,frictionAngle=radians(compFricDegree),density=1400,label='spheres'))
O.materials.append(FrictMat(young=1000*young,poisson=0.30,frictionAngle=0,density=0,label='walls'))
walls=aabbWalls([mn,mx],thickness=0.1*0.001,material='walls',color=(0.3,0.3,0.3))
wallIds=O.bodies.append(walls)

psdSizes=[0.01447*0.001,0.01613*0.001,0.01746*0.001,0.01874*0.001,0.02005*0.001,0.02151*0.001\
    ,0.02328*0.001,0.02568*0.001,0.02973*0.001,0.03386*0.001]
psdCumm=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
sp=pack.SpherePack()
sp.makeCloud(minCorner=mn,maxCorner=mx,psdSizes=psdSizes,psdCumm=psdCumm,distributeMass=True,num=4000,seed=seed)
'''
psdSizes=[0.02005*0.001]
psdCumm=[1.0]
sp=pack.SpherePack()
sp.makeCloud(minCorner=mn,maxCorner=mx,psdSizes=psdSizes,psdCumm=psdCumm,distributeMass=True,num=4000,seed=seed)
'''
'''
pylab.plot(*sp.psd(bins=30,mass=True))
pylab.legend()
pylab.show()
'''
O.bodies.append([sphere(center,rad,material='spheres',color=(0.5,0.5,0.5)) for center,rad in sp])
triax=TriaxialStressController(
	stressMask             = 2,
        wall_bottom_activated  = False,
	internalCompaction     = False,
        updatePorosity         = True
)

O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
		#[Ip2_FrictMat_FrictMat_FrictPhys()],
		#[Law2_ScGeom_FrictPhys_CundallStrack()]
                [Ip2_FrictMat_FrictMat_MindlinPhys()],
		[Law2_ScGeom_MindlinPhys_Mindlin()]
	),

	GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=10,timestepSafetyCoefficient=0.5),
	triax,
	TriaxialStateRecorder(iterPeriod=200,file='WallStresses'+table.key+'0.04',label='triaxRec'),
	NewtonIntegrator(gravity=(0,-9.81,0),damping=damp)
]
triax.goal1 = triax.goal3=0
triax.goal2 = -56.8e6
qt.View()
##########################################################
while 1:
    O.run(100,True)
    unb      = unbalancedForce()
    topWall  = max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
    top      = triax.stress(triax.wall_top_id)[1]
    bottom   = triax.stress(triax.wall_bottom_id)[1]
    meanS    = (triax.stress(triax.wall_right_id)[0]+triax.stress(triax.wall_top_id)[1]+triax.stress(triax.wall_front_id)[2])/3
    kinetic  = 0
    volumeTotal = 0
    for i in O.bodies:
        if isinstance(i.shape,Sphere):
            volume          = 4/3*3.1415926*pow(i.shape.radius,3)
            volumeTotal    += volume
            kinetic        += 0.5*i.mat.density*volume*(pow(i.state.vel[0],2)+pow(i.state.vel[1],2)+pow(i.state.vel[1],2))
    por      = 1-volumeTotal/x/z/topWall
    print ('unbalanced force:',unb,' | y_top', top,' | y_bottom', bottom,'| kinetic',kinetic,'| porosity',por)
    if unb < 10*unbalancedThread and abs((triax.goal2-top)/(triax.goal2))<10*triaxThread and kinetic < 10*1e-4:
        a = open('/home/zhiguo/Yade-pyFile/sintering/real_distribution/compact_Young_2.5GPa_Fric_25/' + folderName +'/sphere_positions'+str(-triax.goal2/1e6)+'.txt', "w+")
        for i in O.bodies:
            if isinstance(i.shape,Sphere):
                a.write(str(i.id)+' '+str(i.state.pos[0])+' '+str(i.state.pos[1])+' '+str(i.state.pos[2])+' '+str(i.shape.radius)+'\n')
        a.close()

        b = open('/home/zhiguo/Yade-pyFile/sintering/real_distribution/compact_Young_2.5GPa_Fric_25/' + folderName +'/interactions'+str(-triax.goal2/1e6)+'.txt', "w+")
        count = 0
        for i in O.interactions:
            if isinstance(O.bodies[i.id1].shape,Sphere) and isinstance(O.bodies[i.id2].shape,Sphere):
                count += 1
                b.write(str(count)+' '+str(i.id1)+' '+str(O.bodies[i.id1].state.pos[0])+' '+str(O.bodies[i.id1].state.pos[1])+' '\
                        +str(O.bodies[i.id1].state.pos[2])+' '+str(O.bodies[i.id1].shape.radius)+' '\
                        +str(i.id2)+' '+str(O.bodies[i.id2].state.pos[0])+' '+str(O.bodies[i.id2].state.pos[1])+' '\
                        +str(O.bodies[i.id2].state.pos[2])+' '+str(O.bodies[i.id2].shape.radius)+'\n')
        b.close()

        c = open('/home/zhiguo/Yade-pyFile/sintering/real_distribution/compact_Young_2.5GPa_Fric_25/' + folderName +'/boun_wall_'+str(-triax.goal2/1e6)+'.txt', "w+")
        for i in O.interactions:
            if i.id1<6:
                c.write(str(count) + ' ' + str(i.id1) + ' ' + str(O.bodies[i.id2].state.pos[0]) + ' ' + str(O.bodies[i.id2].state.pos[1]) + ' ' \
                        + str(O.bodies[i.id2].state.pos[2]) + ' ' + str(O.bodies[i.id2].shape.radius) + '\n')
                count += 1
            if i.id2<6:
                c.write(str(count) + ' ' + str(i.id2) + ' ' + str(O.bodies[i.id1].state.pos[0]) + ' ' + str(O.bodies[i.id1].state.pos[1]) + ' ' \
                        + str(O.bodies[i.id1].state.pos[2]) + ' ' + str(O.bodies[i.id1].shape.radius) + '\n')
                count += 1
        c.close()
        print ("porosity:",por)
        break

#########################################################
for i in range(1,2):
    triax.goal2=-1e4
    while 1:
        O.run(100,True)
        unb=unbalancedForce()
        topWall=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
        top=triax.stress(triax.wall_top_id)[1]
        bottom=triax.stress(triax.wall_bottom_id)[1]
        meanS=(triax.stress(triax.wall_right_id)[0]+triax.stress(triax.wall_top_id)[1]+triax.stress(triax.wall_front_id)[2])/3
        kinetic=0
        volumeTotal=0
        for i in O.bodies:
            if isinstance(i.shape,Sphere):
                volume=4/3*3.1415926*pow(i.shape.radius,3)
                volumeTotal+=volume
                kinetic+=0.5*i.mat.density*volume*(pow(i.state.vel[0],2)+pow(i.state.vel[1],2)+pow(i.state.vel[1],2))
        por=1-volumeTotal/x/z/topWall
        print ('unbalanced force:',unb,' | y_top', top,' | y_bottom', bottom,'| kinetic',kinetic,'| porosity',por)
        if unb < 10*unbalancedThread and abs((triax.goal2-top)/(triax.goal2))<100*triaxThread and kinetic < 10*1e-4:
            a = open('/home/zhiguo/Yade-pyFile/sintering/real_distribution/compact_Young_2.5GPa_Fric_25/' + folderName +'/sphere_positions'+str(-triax.goal2/1e6)+'.txt', "w+")
            for i in O.bodies:
                if isinstance(i.shape,Sphere):
                    a.write(str(i.state.pos[0])+' '+str(i.state.pos[1])+' '+str(i.state.pos[2])+' '+str(i.shape.radius)+'\n')
            a.close()

            b = open('/home/zhiguo/Yade-pyFile/sintering/real_distribution/compact_Young_2.5GPa_Fric_25/' + folderName +'/interactions'+str(-triax.goal2/1e6)+'.txt', "w+")
            count = 0
            for i in O.interactions:
                if isinstance(O.bodies[i.id1].shape, Sphere) and isinstance(O.bodies[i.id2].shape, Sphere):
                    count += 1
                    b.write(str(count) + ' ' + str(i.id1) + ' ' + str(O.bodies[i.id1].state.pos[0]) + ' ' + str(
                        O.bodies[i.id1].state.pos[1]) + ' ' \
                            + str(O.bodies[i.id1].state.pos[2]) + ' ' + str(O.bodies[i.id1].shape.radius) + ' ' \
                            + str(i.id2) + ' ' + str(O.bodies[i.id2].state.pos[0]) + ' ' + str(
                        O.bodies[i.id2].state.pos[1]) + ' ' \
                            + str(O.bodies[i.id2].state.pos[2]) + ' ' + str(O.bodies[i.id2].shape.radius) + '\n')
            b.close()

            c = open('/home/zhiguo/Yade-pyFile/sintering/real_distribution/compact_Young_2.5GPa_Fric_25/' + folderName + '/boun_wall_' + str(-triax.goal2 / 1e6) + '.txt',"w+")
            for i in O.interactions:
                if i.id1 < 6:
                    c.write(str(count) + ' ' + str(i.id1) + ' ' + str(O.bodies[i.id2].state.pos[0]) + ' ' + str(O.bodies[i.id2].state.pos[1]) + ' ' \
                            + str(O.bodies[i.id2].state.pos[2]) + ' ' + str(O.bodies[i.id2].shape.radius) + '\n')
                    count += 1
                if i.id2 < 6:
                    c.write(str(count) + ' ' + str(i.id2) + ' ' + str(O.bodies[i.id1].state.pos[0]) + ' ' + str(O.bodies[i.id1].state.pos[1]) + ' ' \
                            + str(O.bodies[i.id1].state.pos[2]) + ' ' + str(O.bodies[i.id1].shape.radius) + '\n')
                    count += 1
            c.close()
            print ("porosity:",por)
            break

#########################################################
for i in range(1,21):
    triax.goal2 = - 20e6 - 20e6 * (i - 1)

    while 1:
        O.run(100,True)
        unb=unbalancedForce()
        topWall=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
        top=triax.stress(triax.wall_top_id)[1]
        bottom=triax.stress(triax.wall_bottom_id)[1]
        meanS=(triax.stress(triax.wall_right_id)[0]+triax.stress(triax.wall_top_id)[1]+triax.stress(triax.wall_front_id)[2])/3
        kinetic=0
        volumeTotal=0
        for i in O.bodies:
            if isinstance(i.shape,Sphere):
                volume=4/3*3.1415926*pow(i.shape.radius,3)
                volumeTotal+=volume
                kinetic+=0.5*i.mat.density*volume*(pow(i.state.vel[0],2)+pow(i.state.vel[1],2)+pow(i.state.vel[1],2))
        por=1-volumeTotal/x/z/topWall
        print ('unbalanced force:',unb,' | y_top', top,' | y_bottom', bottom,'| kinetic',kinetic,'| porosity',por)
        if unb<unbalancedThread and abs((triax.goal2-top)/(triax.goal2))<triaxThread and kinetic<1e-4:
            a = open('/home/zhiguo/Yade-pyFile/sintering/real_distribution/compact_Young_2.5GPa_Fric_25/' + folderName +'/sphere_positions'+str(-triax.goal2/1e6)+'.txt', "w+")
            for i in O.bodies:
                if isinstance(i.shape,Sphere):
                    a.write(str(i.state.pos[0])+' '+str(i.state.pos[1])+' '+str(i.state.pos[2])+' '+str(i.shape.radius)+'\n')
            a.close()
            b = open('/home/zhiguo/Yade-pyFile/sintering/real_distribution/compact_Young_2.5GPa_Fric_25/' + folderName + '/interactions' + str(-triax.goal2 / 1e6) + '.txt',"w+")
            count = 0
            for i in O.interactions:
                if isinstance(O.bodies[i.id1].shape, Sphere) and isinstance(O.bodies[i.id2].shape, Sphere):
                    count += 1
                    b.write(str(count) + ' ' + str(i.id1) + ' ' + str(O.bodies[i.id1].state.pos[0]) + ' ' + str(O.bodies[i.id1].state.pos[1]) + ' ' \
                            + str(O.bodies[i.id1].state.pos[2]) + ' ' + str(O.bodies[i.id1].shape.radius) + ' ' \
                            + str(i.id2) + ' ' + str(O.bodies[i.id2].state.pos[0]) + ' ' + str(O.bodies[i.id2].state.pos[1]) + ' ' \
                            + str(O.bodies[i.id2].state.pos[2]) + ' ' + str(O.bodies[i.id2].shape.radius) + '\n')
            b.close()

            c = open(
                '/home/zhiguo/Yade-pyFile/sintering/real_distribution/compact_Young_2.5GPa_Fric_25/' + folderName + '/boun_wall_' + str(-triax.goal2 / 1e6) + '.txt',
                "w+")
            for i in O.interactions:
                if i.id1 < 6:
                    c.write(str(count) + ' ' + str(i.id1) + ' ' + str(O.bodies[i.id2].state.pos[0]) + ' ' + str(O.bodies[i.id2].state.pos[1]) + ' ' \
                            + str(O.bodies[i.id2].state.pos[2]) + ' ' + str(O.bodies[i.id2].shape.radius) + '\n')
                    count += 1
                if i.id2 < 6:
                    c.write(str(count) + ' ' + str(i.id2) + ' ' + str(O.bodies[i.id1].state.pos[0]) + ' ' + str(O.bodies[i.id1].state.pos[1]) + ' ' \
                            + str(O.bodies[i.id1].state.pos[2]) + ' ' + str(O.bodies[i.id1].shape.radius) + '\n')
                    count += 1
            c.close()
            print ("porosity:",por)
            break
