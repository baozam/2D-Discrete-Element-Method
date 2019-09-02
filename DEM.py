# -*- cording: utf-8 -*-
#Discrete Element Method 2-Dimension edit by kazama 2019/09
import sys
import numpy as np
import math
import time
import random
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.Qt import (QColor, QTimer)
from PyQt5.QtCore import (QLineF, QPointF, QRectF, Qt, QTimer)
from PyQt5.QtGui import (QBrush, QColor, QPainter)
from PyQt5.QtWidgets import (QApplication, QGraphicsView, QGraphicsScene, QGraphicsItem,QGridLayout, QVBoxLayout, QHBoxLayout, QSizePolicy, QLabel, QLineEdit, QPushButton)

class DEM():
    def __init__(self):

        self.particle_num = 0
        self.particle_radius = 0
        self.particle_wall_radius = 0
        self.particle_wall_num = 0
        self.den =0
        self.g = 0
        self.Eras_n = 0
        self.Eras_s = 0
        self.Visc_n = 0
        self.Visc_s = 0
        self.Frac = 0
        self.Eras_wall_n = 0
        self.Eras_wall_s = 0
        self.Visc_wall_n = 0
        self.Visc_wall_s = 0
        self.Frac_wall = 0
        self.m = 0.0
        self.Ir = 0.0
        self.dt = 0.0

        self.iteration = 0
        self.ax = 0.0
        self.ay = 0.0
        self.aphai = 0.0
        self.un = 0.0
        self.us = 0.0
        self.vn = 0.0
        self.vs = 0.0
        self.vx = []
        self.vy = []
        self.vphai = []
        self.dx = []
        self.dy = []
        self.dphai = []
        self.fn = []
        self.fs = []
        self.fx = []
        self.fy = []
        self.fm = []
        self.fEn = []
        self.fEs = []
        self.fVn = []
        self.fVs = []
        self.phai = []

    def main_calc(self,ax,ay,aphai,vx,vy,vphai,dx,dy,dphai,xp,yp,xw,yw,phai,fn,fs,fx,fy,fm,fEn,fEs,fVn,fVs,Ir,dt,m,particle_num,particle_wall_num,particle_radius,particle_wall_radius,Eras_n,Eras_s,Eras_wall_n,Eras_wall_s,Visc_n,Visc_s,Visc_wall_n,Visc_wall_s,Frac,Frac_wall):
        self.ax = ax
        self.ay = ay
        self.aphai = aphai
        self.vx = vx
        self.vy = vy
        self.vphai = vphai
        self.dx = dx
        self.dy = dy
        self.dphai = dphai
        self.x = xp
        self.y = yp
        self.particle_wall_lineX = xw
        self.particle_wall_lineY = yw
        self.phai = phai
        self.fx = fx
        self.fy = fy
        self.fm = fm
        self.fn = fn
        self.fs = fs
        self.fVn = fVn
        self.fVs = fVs
        self.fEs = fEs
        self.fEn = fEn
        self.Ir = Ir
        self.dt = dt
        self.m = m
        self.particle_num = particle_num
        self.particle_wall_num = particle_wall_num
        self.particle_radius = particle_radius
        self.particle_wall_radius = particle_wall_radius
        self.Eras_n = Eras_n
        self.Eras_s = Eras_s
        self.Eras_wall_n = Eras_wall_n
        self.Eras_wall_s = Eras_wall_s
        self.Visc_n = Visc_n
        self.Visc_s = Visc_s
        self.Visc_wall_n = Visc_wall_n
        self.Visc_wall_s = Visc_wall_s
        self.Frac = Frac
        self.Frac_wall = Frac_wall
        self.g = 9.81
        self.fx = [0.0]*self.particle_num
        self.fy = [0.0]*self.particle_num
        self.fm = [0.0]*self.particle_num
        j = 0
        l = 0
        for j in range(self.particle_num):
            k = 1 + j
            while k < self.particle_num:
                lx = self.x[j]-self.x[k]
                ly = self.y[j]-self.y[k]
                ldpp = math.sqrt(lx**2+ly**2)
                if ldpp < self.particle_wall_radius +self.particle_radius:
                    sin = ly/ldpp
                    cos = lx/ldpp
                    self.fEn,self.fEs,self.fVn,self.fVs,self.fn,    self.fs = self.particle_particle_collision(j,   k,lx,ly,sin,cos,ldpp,self.fEn,self.fEs,    self.fVn,self.fVs,self.fn,self.fs,self.dx,  self.dy,self.dphai,self.vx,self.vy,   self.vphai)
                    if self.fn[j] <= 0.0:
                        self.fn[j] = 0.0
                    elif self.fn[k] <= 0.0:
                        self.fn[k] = 0.0
                    elif abs(self.fs[j]) > self.Frac*self.fn[j]:
                        self.fs[j] = self.Frac*abs(self.fn[j])  *self.fs[j]/abs(self.fs[j])
                    elif abs(self.fs[k]) > self.Frac*self.fn[k]:
                        self.fs[k] = self.Frac*abs(self.fn[k])  *self.fs[k]/abs(self.fs[k])
                    self.fx[j] += -self.fn[j]*cos+self.fs[j]*sin
                    self.fy[j] += -self.fn[j]*sin-self.fs[j]*cos
                    self.fm[j] -= self.particle_radius*self.fs[j]
                    self.fx[k] += self.fn[k]*cos-self.fs[k]*sin
                    self.fy[k] += self.fn[k]*sin+self.fs[k]*cos
                    self.fm[k] -= self.particle_radius*self.fs[k]
                else:
                    self.fEn[j] = 0.0
                    self.fEs[j] = 0.0
                    self.fEn[k] = 0.0
                    self.fEs[k] = 0.0
                k += 1
        while l < self.particle_num:
            n0 = 1 + l
            while n0 < self.particle_wall_num:
                lx = self.x[l]-self.particle_wall_lineX[n0]
                ly = self.y[l]-self.particle_wall_lineY[n0]
                ldpw = math.sqrt(lx**2+ly**2)
                if ldpw < self.particle_wall_radius +self.particle_radius:
                    self.fEn,self.fEs,self.fVn,self.fVs,self.fn,    self.fs,self.fx,self.fy,self.fm =   self.particle_wall_collision(l,lx,ly,ldpw,    self.fEn,self.fEs,self.fVn,self.fVs,self.fn,    self.fs,self.dx,self.dy,self.dphai,self.vx, self.vy,self.vphai,self.fx,self.fy,self.fm,  self.Ir)
                else:
                    self.fEn[l] = 0.0
                    self.fEs[l] = 0.0
                n0 += 1
            l += 1
        for n1 in range(self.particle_num):
            self.fy[n1] += -self.g*self.m
        #Parallel(n_jobs=-1)([delayed(update_position)(n,ax,ay, aphai,vx,vy,vphai,dx,dy,dphai,x,y,phai) for n in range   (particle_num)])
        self.vx,self.vy,self.vphai,self.dx,self.dy,self.dphai,  self.x,self.y,self.phai,self.ay = self.update_position    (self.ax,self.ay,self.aphai,self.vx,self.vy,self.vphai, self.dx,self.dy,self.dphai,self.x,self.y,self.phai,  self.fx,self.fy,self.fm,self.Ir)
        print(self.dt,self.y)
        return self.x,self.y,self.particle_wall_lineX,  self.particle_wall_lineY,self.particle_num,   self.particle_wall_num,self.particle_radius,   self.particle_wall_radius

    def update_position(self,ax,ay,aphai,vx,vy,vphai,dx,dy,dphai,x,y,phai,fx,fy,fm,Ir):
        for n in range(self.particle_num):
            self.ax = self.fx[n]/self.m
            self.ay = self.fy[n]/self.m
            self.aphai = self.fm[n]/Ir
            self.vx[n] += self.ax*self.dt
            self.vy[n] += self.ay*self.dt
            self.vphai[n] += self.aphai*self.dt
            self.dx[n] = self.vx[n]*self.dt
            self.dy[n] = self.vy[n]*self.dt
            self.dphai[n] = self.vphai[n]*self.dt
            self.x[n] += self.dx[n]
            self.y[n] += self.dy[n]
            self.phai[n] += self.dphai[n]
        return self.vx,self.vy,self.vphai,self.dx,self.dy,self.dphai,self.x,y,self.phai,self.ay

    def particle_wall_collision(self,l,lx,ly,ldpw,fEn,fEs,fVn,fVs,fn,fs,dx,dy,dphai,vx,vy,vphai,fx,fy,fm,Ir):
        sin = ly/ldpw
        cos = lx/ldpw
        ldx = dx[l]
        ldy = dy[l]
        lvx = vx[l]
        lvy = vy[l]
        self.un = ldx*cos+ldy*sin
        self.us = -ldx*sin+ldy*cos+self.particle_radius*self.dphai[l]
        self.vn = lvx*cos+lvy*sin
        self.vs = -lvx*sin+lvy*cos+self.particle_radius*self.vphai[l]
        self.fEn[l] += self.Eras_wall_n*self.un
        self.fEs[l] += self.Eras_wall_s*self.us
        self.fVn[l] = self.Visc_wall_n*self.vn
        self.fVs[l] = self.Visc_wall_s*self.vs
        self.fn[l] += self.fEn[l]+self.fVn[l]
        self.fs[l] += self.fEs[l]+self.fVs[l]
        if self.fn[l] <= 0.0:
            self.fn[l] = 0.0
        elif abs(self.fs[l]) > self.Frac_wall*self.fn[l]:
            self.fs = self.Frac*abs(self.fn[l])*self.fs[l]/abs(self.fs[l])
        self.fx[l] += self.fn[l]*cos-self.fs[l]*sin
        self.fy[l] += self.fn[l]*sin+self.fs[l]*cos
        self.fm[l] -= self.particle_radius*self.fs[l]
        return self.fEn,self.fEs,self.fVn,self.fVs,self.fn,self.fs,self.fx,self.fy,self.fm

    def particle_particle_collision(self,j,k,lx,ly,sin,cos,ldpp,fEn,fEs,fVn,fVs,fn,fs,dx,dy,dphai,vx,vy,vphai):
        ldx = self.dx[j]-self.dx[k]
        ldy = self.dy[j]-self.dy[k]
        lvx = self.vx[j]-self.vx[k]
        lvy = self.vy[j]-self.vy[k]
        self.un = ldx*cos+ldy*sin
        self.us = -ldx*sin+ldy*cos+(self.particle_radius*dphai[j]+self.particle_radius*self.dphai[k])
        self.vn = lvx*cos+lvy*sin
        self.vs = -lvx*sin+lvy*cos+(self.particle_radius*vphai[j]+self.particle_radius*self.vphai[k])
        self.fEn[j] += self.Eras_n*self.un
        self.fEs[j] += self.Eras_s*self.us
        self.fVn[j] = self.Visc_n*self.vn
        self.fVs[j] = self.Visc_s*self.vs
        self.fn[j] += self.fEn[j]+self.fVn[j]
        self.fs[j] += self.fEs[j]+self.fVs[j]
        self.fEn[k] += self.Eras_n*self.un
        self.fEs[k] += self.Eras_s*self.us
        self.fVn[k] = self.Visc_n*self.vn
        self.fVs[k] = self.Visc_s*self.vs
        self.fn[k] += self.fEn[k]+self.fVn[k]
        self.fs[k] += self.fEs[k]+self.fVs[k]
        return self.fEn,self.fEs,self.fVn,self.fVs,self.fn,self.fs


class MyWindow(FigureCanvas):
    def __init__(self,parent=None,title=None,size=(1000,1000), width=4, height=3, dpi=100):
        #super().__init__()
        self.particle_x = 6
        self.particle_y = 4
        self.particle_num = self.particle_x*self.particle_y
        self.particle_wall_x = 40
        self.particle_wall_y = 40
        self.particle_radius = 1
        self.particle_wall_radius = 2
        self.particle_wall_num = self.particle_wall_x+self.particle_wall_y+2
        self.den = 1
        self.g = 9.81
        self.Eras_n = 2*10
        self.Eras_s = 500
        self.Visc_n = 500
        self.Visc_s = 100
        self.Frac = 100000
        self.Eras_wall_n = 2*10
        self.Eras_wall_s = 500
        self.Visc_wall_n = 100
        self.Visc_wall_s = 100
        self.Frac_wall = 10000
        self.m = 4/3*math.pi*self.den*self.particle_radius**3
        self.Ir = math.pi*self.den*self.particle_radius**4.0/2.0
        self.dt = 2*(self.m/self.Eras_n)**0.5*10**-2
        self.ax = 0.0
        self.ay = 0.0
        self.aphai = 0.0
        self.vx = [0 for i in range(self.particle_num)]
        self.vy = [0 for i in range(self.particle_num)]
        self.vphai = [0 for i in range(self.particle_num)]
        self.dx = [0 for i in range(self.particle_num)]
        self.dy = [0 for i in range(self.particle_num)]
        self.dphai = [0 for i in range(self.particle_num)]
        self.fn = [0 for i in range(self.particle_num)]
        self.fs = [0 for i in range(self.particle_num)]
        self.fx = [0 for i in range(self.particle_num)]
        self.fy = [0 for i in range(self.particle_num)]
        self.fm = [0 for i in range(self.particle_num)]
        self.fEn = [0 for i in range(self.particle_num)]
        self.fEs = [0 for i in range(self.particle_num)]
        self.fVn = [0 for i in range(self.particle_num)]
        self.fVs = [0 for i in range(self.particle_num)]
        self.phai = [0 for i in range(self.particle_num)]

        self.iteration = 0
    
    
        self.xp = []
        for j in range(self.particle_y):
            self.xp.extend([(i+2)*self.particle_radius*2 for i in range(self.particle_x)])
        self.yp = []
        for j in range(self.particle_y):
            self.yp.extend([(j+2)*self.particle_radius*2+self.particle_radius for i in range(self.particle_x)])
    
        self.xw = []
        self.yw = []
        self.xw.extend([ 0 for i in range(self.particle_wall_y+1)])
        self.xw.extend([ (i+1)*self.particle_wall_radius for i in range(self.particle_wall_x)])
        self.xw.extend([ (self.particle_wall_x+1)*self.particle_wall_radius for j in range(self.particle_wall_y+1)])
        self.yw.extend([ i*self.particle_wall_radius for i in range(self.particle_wall_y+1)])
        self.yw.extend([ 0 for i in range(self.particle_wall_x)])
        self.yw.extend([ i*self.particle_wall_radius for i in range(self.particle_wall_y+1)])
    
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)

        super(MyWindow, self).__init__(self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

        self.resize(*size)
        self.setWindowTitle(title)
        self.setStyleSheet("background-color : white;")
        self.show()
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.OnTimer)
        self.timer.start(1000)

    def OnTimer(self):
        self.height()
        self.width()

        dem = DEM()
        self.xp,self.yp,self.xw,self.yw,self.particle_num,self.particle_wall_num,self.particle_radius,self.particle_wall_radius = dem.main_calc(self.ax,self.ay,self.aphai,self.vx,self.vy,self.vphai,self.dx,self.dy,self.dphai,self.xp,self.yp,self.xw,self.yw,self.phai,self.fn,self.fs,self.fx,self.fy,self.fm,self.fEn,self.fEs,self.fVn,self.fVs,self.Ir,self.dt,self.m,self.particle_num,self.particle_wall_num,self.particle_radius,self.particle_wall_radius,self.Eras_n,self.Eras_s,self.Eras_wall_n,self.Eras_wall_s,self.Visc_n,self.Visc_s,self.Visc_wall_n,self.Visc_wall_s,self.Frac,self.Frac_wall)
        self.iteration += 1
        if self.iteration%100 == 0:
            print(self.iteration)
        self.update_figure()

        self.update()

    def update_figure(self):
        self.axes.clear()
        self.axes = self.fig.add_subplot(111)
        self.axes.scatter(self.xp,self.yp,marker='o',color='red',s=self.particle_radius*50)
        self.axes.scatter(self.xw,self.yw,marker='o',color='black',s=self.particle_wall_radius*50)
        self.draw()


def Main():
    dem = DEM()

    #dem.main_calc(ax,ay,aphai,vx,vy,vphai,dx,dy,dphai,xp,yp,xw,yw,phai,fn,fs,fx,fy,fm,fEn,fEs,fVn,fVs,Ir,dt,m,particle_num,particle_wall_num,particle_radius,particle_wall_radius,Eras_n,Eras_s,Eras_wall_n,Eras_wall_s,Visc_n,Visc_s,Visc_wall_n,Visc_wall_s,Frac,Frac_wall)

    app  = QApplication(sys.argv)
    window = MyWindow(title="Discrete Element Method 2D kazama", size=(500, 500), width=4, height=3, dpi=100)
    sys.exit(app.exec_())

if __name__=="__main__":
    Main()