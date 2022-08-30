import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
import pylab as pl
import random as rnd
import time
from scipy.integrate import odeint

# Total population, N.
N = 1000
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 1, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 0.4, 0.03
# A grid of time points (in days)
t = np.linspace(0, 160, 160)

# The SIR model differential equations.
def deriv(y, t, N, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I
    dRdt = gamma * I
    return dSdt, dIdt, dRdt

# Initial conditions vector
y0 = S0, I0, R0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma))
S, I, R = ret.T

# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S/1000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I/1000, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R/1000, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number (1000s)')
ax.set_ylim(0,1.2)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()
    
r_container=1000
num=100
r_ball=10
step=0.0006
fs=2.5

t_start=time.time()
time_list=[0]

for vel in [5000]:     
    m=1.67e-29
    k=1.3806488e-23
    def add_container(r):
        patch_container=pl.Circle([0,0],r,fc="None",ec="k")
        ax1.add_patch(patch_container)

    def add_ball(num):
        ax1.add_patch(Ball_patch[num])
    
    def initialize(sim):
        global Ball_data
        global Ball_class
        global speed_generation

        def velocity():
            v_x=rnd.uniform(-sim.v_ball,sim.v_ball)
            sgn=rnd.uniform(-1,1)
            if sgn>=0:
                v_y=(sim.v_ball**2-v_x**2)**0.5
            else:
                v_y=-(sim.v_ball**2-v_x**2)**0.5

            return [v_x,v_y]

        def check_pos(pos):
            r=(pos[0]**2+pos[1]**2)**0.5

            flag=True
            if len(ball_position)>0:
                for i in ball_position:
                    if ((i[0]-pos[0])**2+(i[1]-pos[1])**2)**0.5<1.1*2*sim.r_ball:
                        flag=False
                        break

            if r<=sim.r_container-1.1*sim.r_ball and flag==True:
                return True
            else:
                return False

        def position():
            pos=[rnd.uniform(-sim.r_container,sim.r_container),rnd.uniform(-sim.r_container,sim.r_container)]
            while check_pos(pos)==False:
                pos=[rnd.uniform(-sim.r_container,sim.r_container),rnd.uniform(-sim.r_container,sim.r_container)]

            ball_position.append(pos)
            
        add_container(sim.r_container)

        Ball_data=[None for i in range(0,sim.num)]

        ball_position=[]
        for i in range(0,sim.num):
            position()

        for i in range(0,sim.num):
            Ball_data[i]=[sim.r_ball,ball_position[i],velocity()]

        Ball_class=[None for i in range(0,sim.num)]
        for i in range(0,sim.num):
            Ball_class[i]=Ball(Ball_data[i][0],Ball_data[i][1],Ball_data[i][2])
                
            Ball_patch.append(pl.Circle(Ball_class[i].pos,Ball_class[i].radius,fc="b",ec="b"))
            add_ball(i)

    def next_collision(sim):
        global collision_with_container
        global num_collision
        
        collision_status=[[[i],[[j,None] for j in range(i+1,sim.num)]+[["c",None]]] for i in range(0,sim.num)]
        for i in range(0,sim.num):
            pos=Ball_class[i].pos
            vel=Ball_class[i].velocity

            a=(vel[0]**2+vel[1]**2)
            b=2*(pos[0]*vel[0]+pos[1]*vel[1])
            c=pos[0]**2+pos[1]**2-(sim.r_container-sim.r_ball)**2
            if a!=0:
                t=(-b+(b**2-4*a*c)**0.5)/(2*a)
            else:
                t=1e10
            
            collision_status[i][1][-1][-1]=t

            for j in range(0,len(collision_status[i][1])-1):
                pos_i=Ball_class[i].pos
                vel_i=Ball_class[i].velocity
                pos_j=Ball_class[1+i+j].pos
                vel_j=Ball_class[1+i+j].velocity

                dx=pos_i[0]-pos_j[0]
                dy=pos_i[1]-pos_j[1]
                dvx=vel_i[0]-vel_j[0]
                dvy=vel_i[1]-vel_j[1]

                a=dvx**2+dvy**2
                b=2*(dx*dvx+dy*dvy)
                c=dx**2+dy**2-4*sim.r_ball**2

                if a==0:
                    collision_status[i][1][j][-1]=1e10
                elif round(b**2-4*a*c,9)<0:
                    collision_status[i][1][j][-1]=1e10
                elif round(b**2-4*a*c,9)==0:
                    collision_status[i][1][j][-1]=1e10
                else:
                    t=(-b-(b**2-4*a*c)**0.5)/(2*a)
                    if round(t,9)<0:
                        collision_status[i][1][j][-1]=1e10
                    else:
                        collision_status[i][1][j][-1]=t

        t=1e9
        for i in range(0,sim.num):
            for j in collision_status[i][1]:
                if j[-1]<t:
                    collision=[[i,j[0],j[-1]]]
                    t=j[-1]
                elif j[-1]==t:
                    collision.append([i,j[0],j[-1]])

        num_collision+=1

        for i in collision:
            if i[1]=="c":
                collision_with_container+=1
                break

        return collision

    def ball_par_change(collision_l):
        changed=[]
        for i in collision_l:
            collision=i
            changed.append(collision[0])
            if collision[1]!="c":
                changed.append(collision[1])

            if collision[1]=="c":
                x0=Ball_class[collision[0]].pos[0]
                y0=Ball_class[collision[0]].pos[1]
                vx0=Ball_class[collision[0]].velocity[0]
                vy0=Ball_class[collision[0]].velocity[1]
                t=collision[2]

                x1=x0+t*vx0
                y1=y0+t*vy0

                if y1!=0:
                    Nx=1
                    Ny=-x1/y1
                else:
                    Nx=0
                    Ny=1

                pos1_m=(x1**2+y1**2)**0.5
                upos1_x=x1/pos1_m
                upos1_y=y1/pos1_m
                N_m=(Nx+Ny**2)**0.5
                uN_x=Nx/N_m
                uN_y=Ny/N_m

                vp1=-(vx0*upos1_x+vy0*upos1_y)
                vn1=vx0*uN_x+vy0*uN_y
                vx1=vp1*upos1_x+vn1*uN_x
                vy1=vp1*upos1_y+vn1*uN_y

                Ball_class[collision[0]].pos=[x1,y1]
                Ball_class[collision[0]].velocity=[vx1,vy1]
            else:
                x0_0=Ball_class[collision[0]].pos[0]
                y0_0=Ball_class[collision[0]].pos[1]
                vx0_0=Ball_class[collision[0]].velocity[0]
                vy0_0=Ball_class[collision[0]].velocity[1]

                x1_0=Ball_class[collision[1]].pos[0]
                y1_0=Ball_class[collision[1]].pos[1]
                vx1_0=Ball_class[collision[1]].velocity[0]
                vy1_0=Ball_class[collision[1]].velocity[1]

                t=collision[2]

                x0_1=x0_0+t*vx0_0
                y0_1=y0_0+t*vy0_0

                x1_1=x1_0+t*vx1_0
                y1_1=y1_0+t*vy1_0

                dvx01_0=vx0_0-vx1_0
                dvy01_0=vy0_0-vy1_0
                dx01_1=x0_1-x1_1
                dy01_1=y0_1-y1_1

                if dy01_1!=0:
                    Ndx=1
                    Ndy=-dx01_1/dy01_1
                else:
                    Ndx=0
                    Ndy=1

                drm=(dx01_1**2+dy01_1**2)**0.5
                udx=dx01_1/drm
                udy=dy01_1/drm

                Nm=(Ndx**2+Ndy**2)**0.5
                uNdx=Ndx/Nm
                uNdy=Ndy/Nm

                vh01_1=dvx01_0*uNdx+dvy01_0*uNdy
                vp01_1=0
                vx01_1=vh01_1*uNdx
                vy01_1=vh01_1*uNdy
                vx0_1=vx01_1+vx1_0
                vy0_1=vy01_1+vy1_0

                vh11_1=0
                vp11_1=dvx01_0*udx+dvy01_0*udy
                vx11_1=vp11_1*udx
                vy11_1=vp11_1*udy
                vx1_1=vx11_1+vx1_0
                vy1_1=vy11_1+vy1_0

                Ball_class[collision[0]].pos=[x0_1,y0_1]
                Ball_class[collision[1]].pos=[x1_1,y1_1]
                Ball_class[collision[0]].velocity=[vx0_1,vy0_1]
                Ball_class[collision[1]].velocity=[vx1_1,vy1_1]

        for i in range(0,len(Ball_class)):
            if i not in changed:
                Ball_class[i].pos=[Ball_class[i].pos[0]+t*Ball_class[i].velocity[0],Ball_class[i].pos[1]+t*Ball_class[i].velocity[1]]

    class Ball():
        def __init__(self,r,pos,v):
            self.radius=r
            self.pos=pos
            self.velocity=v

        def __str__(s):
            msg="Radius: "+repr(s.radius)+", Position: "+repr(s.pos)+", velocity: "+repr(s.velocity)
            return msg

    class Simulation():
        def __init__(self,rc,num,rb,vb):
            self.r_container=rc
            self.num=num
            self.r_ball=rb
            self.v_ball=vb

    Ball_patch=[]

    sim1=Simulation(r_container,num,r_ball,vel)
    
    fig1=pl.figure(figsize=(2*fs,2*fs))
    ax1=pl.axes(xlim=(-sim1.r_container, sim1.r_container), ylim=(-sim1.r_container, sim1.r_container),aspect=1)
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_title("simulation")
    initialize(sim1)
    
    time_l=[0]
    Ek_l=[]
    momentum_x=[]
    momentum_y=[]
    num_collision=0
    collision_with_container=0
    nclp=0
    lpct=0
    
    Ek=0
    for i in Ball_class:
        Ek+=0.5*m*(i.velocity[0]**2+i.velocity[1]**2)

    Ek_l.append(Ek)
    
    
    if sim1.num<120:
        n=15
    else:
        n=int(sim1.num/8)+1

    mo_x=0
    mo_y=0
    for i in Ball_class:
        mo_x+=m*i.velocity[0]
        mo_y+=m*i.velocity[1]

    momentum_x.append(mo_x)
    momentum_y.append(mo_y)

    pressure=0
    
    while True:
        collision_l=next_collision(sim1)
        t=collision_l[0][-1]
        if t>step:
            dt=t/(t//step)
        else:
            dt=t
            
        for index in range(0,int(t/dt)):
            for i in range(0,sim1.num):
                Ball_patch[i].center=[Ball_class[i].pos[0]+dt*(index+1)*Ball_class[i].velocity[0],Ball_class[i].pos[1]+dt*(index+1)*Ball_class[i].velocity[1]]
                
            pl.pause(dt)

        ball_par_change(collision_l)

        time_l.append(time_l[-1]+t)
        
        delta_t=int(time.time()-t_start)
        
        if delta_t>time_list[-1]:
            time_list.append(delta_t)
            inf=int(I[delta_t])
            rec=int(R[delta_t])
            sus=1000-inf-rec
            
            for index in range(0,rec//10):
                Ball_patch[index].set_color((0,1,0,1))
            
            if rec%10!=0:
                Ball_patch[rec//10].set_color((1/(1+(rec%10)/(10-rec%10)),1-1/(1+(rec%10)/(10-rec%10)),0,1))
            
            for index in range(rec//10+1,rec//10+inf//10):
                Ball_patch[index].set_color((1,0,0,1))
                
            if sus%10!=0:
                Ball_patch[rec//10+inf//10].set_color((1-1/(1+(10-sus%10)/(sus%10)),0,1/(1+(10-sus%10)/(sus%10)),1))