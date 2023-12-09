#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 11:33:42 2019

@author: gibies
"""
import sys
import os
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.path.dirname(CURR_PATH)
OBSLIB=os.environ.get('OBSLIB',PKGHOME+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',PKGHOME+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSNML',PKGHOME+"/nml")
sys.path.append(OBSNML)
import obsmod
import glob
from pathlib import Path
import numpy
import tkinter
import tkFont
from tkinter import Canvas
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

class obsdata:
    def mainloop(self):
        self.mainbox.mainloop()
    
    def fwrdbtn_click(self):
          if self.obstypelist is None: 
            if self.scropt is not None: self.getsrcopt()
          else:
            if self.obstype is None: self.getobstype()
            if self.srcopt is 1: 
                with open(self.infile, "rb") as self.obsfile:
                    if self.subtypelist is None: 
                        self.list_subtype()
                    else:
                        if self.enamlist is None: 
                            self.list_ename()
                        else: 
                            self.addelenam()
            else:
                if self.subtypelist is None: 
                    self.list_subtype()
                else:
                    if self.enamlist is None: 
                        self.list_ename()
                    else: 
                        self.addelenam()
        
    
    def getpath(self,scropt=None):
        if scropt is not None: self.srcopt=scropt
        self.getsrcopt(self.srcopt)
        
        
    def getsrcopt(self,srcopt=None):
        if srcopt is None:
            selection = "You have not made any selection "
        else:
            self.srcopt=srcopt  #.get()
        self.path=self.ebx_path.get()

        if self.srcopt is 0:
            self.path=self.path_datac+"/gl_bufr"
            self.filextn="bufr"
        if self.srcopt is 1:
            self.path=self.path_datac+"/gl_obstore"
            self.filextn="obstore"
        if self.srcopt is 2:
            self.path=self.path_datac+"/gl_odb2"
            self.filextn="odb"
        selection = "You selected the option " +self.filextn
        self.message(selection)  
        self.obstypelist=[Path(obsfile).stem for obsfile in glob.glob(self.path+"/*."+self.filextn)]
        self.message(self.obstypelist)
        self.lbl2=tkinter.Label(self.frm31, text='ObsType')
        self.lbl2.grid(row=1, column=0,sticky="w")
        choices = self.obstypelist
        self.otyp = tkinter.StringVar(self.mainbox)
        # initial value
        self.otyp.set(self.obstypelist[0])
        self.opt_otyp = tkinter.OptionMenu(self.frm31, self.otyp, *self.obstypelist)
        #self.opt_stp = tkinter.OptionMenu(self.frm31, self.ptr_styp, *self.subtypelist)
        #self.opt_enam = tkinter.OptionMenu(self.frm31, self.eselect, *self.enamlist)
        self.opt_otyp.config(width=30)
        self.opt_otyp.grid(row=1, column=1,sticky="w")
        self.btn_opnobs.grid(row=1, column=2)
        ########################
    
    def openobstore(self):
        if self.obstype is None: 
                    self.getobstype()
        if self.srcopt is 1:
            with open(self.infile, "rb") as self.obsfile:
                    if self.subtypelist is None: 
                        self.list_subtype()
                    else:
                        if self.enamlist is None: 
                            self.list_ename()
                        else: 
                            self.addelenam()
        else:
            if self.subtypelist is None: 
                self.list_subtype()
            else:
                if self.enamlist is None: 
                    self.list_ename()
                else: 
                    self.addelenam()
    
    def getobstype(self):
        self.opt_otyp.config(state=tkinter.DISABLED)
        self.message("OBSTYPE") 
        self.obstype=self.otyp.get()
        if self.srcopt is 1:
            self.infile=glob.glob(self.path_datac+"/gl_"+self.filextn+"/"+self.obstype[0].upper() + self.obstype[1:]+"."+self.filextn)[0]
        else:
            self.infile=glob.glob(self.path_datac+"/gl_"+self.filextn+"/"+self.obstype+"."+self.filextn)[0]
        self.message(self.infile)  
            ############
            
    def list_subtype(self):
        if self.srcopt is 1:
            self.subtypelist=obsmod.obstore_list_subtype(self.obsfile).tolist()
        else:
            self.subtypelist=obsmod.odb_list_subtype(self.infile).tolist()
        indx = range(1,len(self.subtypelist)+1)
        self.ptr_indx = tkinter.StringVar(self.mainbox)
        # initial value
        self.ptr_indx.set(0)
        self.opt_indx = tkinter.OptionMenu(self.frm42, self.ptr_indx, *indx, command=lambda:self.ptr_styp.set(self.subtypelist[self.ptr_indx.get()]))
        self.opt_indx.config(width=10)
        self.opt_indx.grid(row=0, column=0,sticky="w")
        
        self.ptr_styp = tkinter.StringVar(self.mainbox)
        # initial value
        self.ptr_styp.set(self.subtypelist[0])
        self.opt_styp = tkinter.OptionMenu(self.frm42, self.ptr_styp, *self.subtypelist, command=lambda:self.ptr_indx.set(0))
        self.opt_styp.config(width=10)
        self.opt_styp.grid(row=0, column=1,sticky="w")
        self.btn_opnobs.grid(row=2, column=2)

    def list_ename(self):
#        if self.srcopt is 1:
#                self.subtype=obsmod.obstore_read_index_subtype(self.obsfile,self.ptr_indx.get())
#                self.stypindx=obsmod.obstore_read_subtype_index(self.obsfile,self.ptr_styp.get())
        if int(self.ptr_indx.get()) is not 0:
            self.stypindx=int(self.ptr_indx.get())
            self.subtype=self.subtypelist[int(self.ptr_indx.get())-1]
            self.ptr_styp.set(self.subtype)
        else:
            self.subtype=self.ptr_styp.get()
        if self.stypindx is None:
            self.stypindx=numpy.where(numpy.array(self.subtypelist) == int(self.subtype))[0][0] + 1
            self.message(self.stypindx)
            self.subtype=self.ptr_styp.get()
            self.ptr_indx.set(self.stypindx)
        self.opt_styp.config(state=tkinter.DISABLED)
        self.opt_indx.config(state=tkinter.DISABLED)
        if self.srcopt is 1:
            self.enamlist=obsmod.getelenams(self.obsfile,self.obs_index_nml,self.subtype).tolist()
        else:
            self.enamlist=obsmod.odb_list_varname(self.infile)   #["lat","lon","varno","obsvalue","obs_error","fg_depar"]  #,"an_depar"
        choices = self.enamlist
        self.display.insert(tkinter.END, choices)
        self.display.insert(tkinter.END, "\n")
        self.eselect = tkinter.StringVar(self.mainbox)
        # initial value
        self.eselect.set(choices[0])
        self.opt_enam = tkinter.OptionMenu(self.frm31, self.eselect, *self.enamlist)
        self.opt_enam.config(width=30)
        self.opt_enam.grid(row=3, column=1,sticky="w")
        #self.btn_opnobs=tkinter.Button(self.frm31, text='>>', width=2, command=lambda:self.fwrdbtn_click() )
        self.btn_opnobs.grid(row=3, column=2)
        
        self.lbox_selected = tkinter.Listbox(self.frm32) 
        self.lbox_selected.grid(row=0, column=0)
        self.lbox_selected.config(height=7)
        
        self.btn_selall=tkinter.Button(self.frm43, text='Select All', width=10, command=lambda:self.selectall() )
        self.btn_selall.grid(row=0, column=1)
        self.btn_select=tkinter.Button(self.frm32, text='Selectlist', width=10, command=lambda:self.select() )
        self.btn_select.grid(row=1, column=0)
        self.btn_usrqry=tkinter.Button(self.frm32, text='>>', width=2, command=lambda:self.filteroption() )
        self.btn_usrqry.grid(row=0, column=1)
        
        self.btn_opnobs.grid(row=3, column=2)
        self.lbox_usrqry = tkinter.Listbox(self.frm33) 
        self.lbox_usrqry.grid(row=4, column=0)
        self.lbox_usrqry.config(width=25,height=5)
        
        self.btn_query=tkinter.Button(self.frm33, text='Query', width=10, command=lambda:self.myjob() )
        self.btn_query.grid(row=5, column=0)
        
    def addelenam(self):
        self.message("ELENAM") 
        self.ename=self.eselect.get()
        self.eselcount=self.eselcount+1
        self.message("SELECTED "+str(self.eselcount)+".) "+self.ename)
        self.lbox_selected.insert(tkinter.END, self.ename) 
        
    def filteroption(self):
        self.frm51 = tkinter.Frame(self.frm33) 
        self.frm51.grid(row=0, column=0)
        fltropt = tkinter.IntVar() 
        self.rb_fltr_sngl=tkinter.Radiobutton(self.frm51, text='single', variable=fltropt, value=1, command=lambda:self.minmax(fltropt.get()) )
        self.rb_fltr_rnge=tkinter.Radiobutton(self.frm51, text='range', variable=fltropt, value=2, command=lambda:self.minmax(fltropt.get()) )
        self.rb_fltr_sngl.grid(row=0, column=0)
        self.rb_fltr_rnge.grid(row=0, column=1)
        
    def minmax(self,fltropt=None):
        if fltropt is not None: 
            self.fltropt=fltropt
        else:
            self.filteroption()
        selection = "You selected the option " +str(self.fltropt)
        self.message(selection)  
        self.elenam=self.lbox_selected.get(self.lbox_selected.curselection())
        if self.srcopt is 2: self.varnolist=obsmod.odb_get_varnolist(self.elenam)
        if self.srcopt is 2 and self.fltropt is 1: self.data=obsmod.odb_filter_varno(self.infile,self.elenam)
        else:
            if self.srcopt is 1:
                with open(self.infile, "rb") as self.obsfile:
                    self.data=obsmod.query_obstore(self.obsfile,self.obs_index_nml,indx=self.stypindx,subtype=self.subtype,selectlist=[self.elenam],userquery=[])
            else:
                self.data=obsmod.query_odb(self.infile,self.odb_index_nml,subtype=self.subtype,selectlist=[self.elenam],varnolist=self.varnolist,userquery=[]) #odb_get_varnolist(varnamlist

        self.frm52 = tkinter.Frame(self.frm33) 
        self.frm52.grid(row=1, column=0)
        self.lbl_elenam=tkinter.Label(self.frm52, text=self.elenam)
        self.ebx_qmin = tkinter.Entry(self.frm52, width=10) 
        self.ebx_qmax = tkinter.Entry(self.frm52, width=10) 
        self.ebx_qmin.grid(row=0, column=0,sticky="w")
        self.lbl_elenam.grid(row=0, column=1,sticky="w")
        self.ebx_qmax.grid(row=0, column=2,sticky="w")
        self.ebx_qmin.insert(tkinter.END, self.data.min().values[0])
        self.ebx_qmax.insert(tkinter.END, self.data.max().values[0])
        self.btn_usrqry=tkinter.Button(self.frm33, text="Filter", width=10, command=lambda:self.usrqry() )
        self.btn_usrqry.grid(row=2, column=0)
        
        
    def usrqry(self):
        self.lbox_usrqry.insert(tkinter.END, self.elenam+">="+str(self.ebx_qmin.get()))
        self.lbox_usrqry.insert(tkinter.END, self.elenam+"<="+str(self.ebx_qmax.get()))
        
    def select(self):
        self.opt_enam.config(state=tkinter.DISABLED)
        self.btn_opnobs.configure(state=tkinter.DISABLED)
        self.selectlist=numpy.array(self.lbox_selected.get(0, tkinter.END)).tolist()
        if self.srcopt is 2: self.varnolist=obsmod.odb_get_varnolist(self.selectlist)
        self.userquery=numpy.array(self.lbox_usrqry.get(0, tkinter.END)).tolist()
        self.message(self.selectlist)
        self.message(self.userquery)
        
    def outfilename_update(self):
        self.outfile=self.ebx_outfile.get()
        self.message(self.outfile)
        
    def savedata(self):
        self.outfilename_update()
        obsmod.obs_frame_ascii(self.outfile,self.data)
    
    def message(self,msg):
        self.display.insert(tkinter.END, msg)
        self.display.insert(tkinter.END, "\n")

    def selectall(self):
        self.lbox_selected.insert(tkinter.END, *self.enamlist)
        
    def clearall(self): 
#            """when clear button is pressed,clears the text area"""
        self.display.delete('1.0',tkinter.END) 
        self.opt_otyp.configure(state=tkinter.NORMAL)
        self.srcopt=None
        self.obstypelist=None
        self.obstype=None
        self.infile=None
        self.selectlist=None
        self.userquery=None
        self.subtypelist=None
        self.stypindx=None
        self.enamlist=None
        self.btn_opnobs.grid(row=0, column=2)
        self.btn_opnobs.configure(state=tkinter.NORMAL)
        self.opt_stp.destroy()
        self.lbox_selected.destroy()
        self.lbox_usrqry.destroy()
        self.lbl_elenam.destroy()
        self.ebx_qmin.destroy()
        self.ebx_qmax.destroy()
        #self.display.insert(END, "") 
        if self.ename is not None: 
                self.eselect=None
                self.ename=None
                self.eselcount=0
                self.opt_enam.grid_forget
                self.opt_enam.destroy()
                self.lbox_selected.grid_forget
                self.lbox_selected.destroy()
        if self.subtype is not None: 
                self.ptr_styp=None
                self.subtype=None
        self.ctime_menu()
   #self.getpath()
    #    self.openobstore()        

    def ctime_menu(self):
        self.srcopt=None
        self.cylctimelist=[Path(ctpath).stem for ctpath in glob.glob(self.path_suite+"/share/cycle/*")]
        self.ptr_ctime = tkinter.StringVar(self.mainbox)
        # initial value
        self.ptr_ctime.set(self.cylctimelist[0])
        self.opt_ctime = tkinter.OptionMenu(self.frm21, self.ptr_ctime, *self.cylctimelist ) 
        self.opt_ctime.grid(row=0, column=2,sticky="w")
        self.opt_ctime.config(width=30)
        self.btn_opnobs=tkinter.Button(self.frm21, text='>>', width=2, command=lambda:self.datac() )
        self.btn_opnobs.grid(row=0, column=3)

    def datac(self, ctime=None):
        if ctime is not None:
                self.cylctime = ctime 
        else:
                self.cylctime = self.ptr_ctime.get()
        self.path_datac=self.path_suite+"/share/cycle/"+str(self.cylctime)
        self.path=self.path_datac
        self.ebx_path = tkinter.Entry(self.frm21, width=100) 
        self.ebx_path.grid(row=0, column=1,sticky="w")
        self.ebx_path.insert(tkinter.END, self.path)
        self.message(self.path)
        self.btn_opnobs=tkinter.Button(self.frm31, text='>>', width=2, command=lambda:self.fwrdbtn_click() )
        self.btn_opnobs.grid(row=0, column=2)
        srcopt = tkinter.IntVar() 
        self.rb1=tkinter.Radiobutton(self.frm41, text='bufr', variable=srcopt, value=0, command=lambda:self.getsrcopt(srcopt.get()))
        self.rb1.grid(row=0, column=0)
        self.rb1=tkinter.Radiobutton(self.frm41, text='obstore', variable=srcopt, value=1, command=lambda:self.getsrcopt(srcopt.get()))
        self.rb1.grid(row=0, column=1)
        self.rb2=tkinter.Radiobutton(self.frm41, text='odb', variable=srcopt, value=2, command=lambda:self.getsrcopt(srcopt.get()))
        self.rb2.grid(row=0, column=2)

    
    def myjob(self):
        self.select()
        if self.srcopt is 1:
            with open(self.infile, "rb") as self.obsfile:
                self.data=obsmod.query_obstore(self.obsfile,self.obs_index_nml,indx=self.stypindx,subtype=self.subtype,selectlist=self.selectlist,userquery=self.userquery)
        else:
            self.message(self.varnolist)
            self.data=obsmod.query_odb(self.infile,self.odb_index_nml,subtype=self.subtype,selectlist=["lat","lon","varno","obsvalue","obs_error","fg_depar"],varnolist=self.varnolist,userquery=self.userquery)
        self.lbl_outpath=tkinter.Label(self.frm23, text='Output File')
        self.lbl_outpath.grid(row=0,sticky="w")
        self.ebx_outfile = tkinter.Entry(self.frm23, width=60) 
        self.ebx_outfile.grid(row=0, column=1,sticky="w")
        self.ebx_outfile.insert(tkinter.END, self.outfile)
        self.btn_savdat=tkinter.Button(self.frm23, text='Save', width=10, command=lambda:self.savedata() )
        self.btn_savdat.grid(row=0, column=3)
        self.message(self.data)
        #self.plotbox = tkinter.Entry(self.frm12)
        #self.plotbox.grid(row=0, column=0)
        self.plot()
        
    def plot (self):
#        x=numpy.array ([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
#        v= numpy.array ([16,16.31925,17.6394,16.003,17.2861,17.3131,19.1259,18.9694,22.0003,22.81226])
#        p= numpy.array ([16.23697,     17.31653,     17.22094,     17.68631,     17.73641 ,    18.6368,
#            19.32125,     19.31756 ,    21.20247  ,   22.41444   ,  22.11718  ,   22.12453])

        fig = Figure(figsize=(6.5,4.5))
        a = fig.add_subplot(111)
        a.scatter(self.data.Longitude,self.data.Latitude,color='red')
        #a.plot(p, range(2 +max(x)),color='blue')
        #a.invert_yaxis()

        a.set_title ("Observation Location", fontsize=16)
        a.set_ylabel("Latitude", fontsize=14)
        a.set_xlabel("Longitude", fontsize=14)

        canvas = FigureCanvasTkAgg(fig, master=self.frm12)
        canvas.get_tk_widget().grid(row=0, column=0)
        canvas.draw()

    def __init__(self,mainbox=None, ROSE_SUITE_DIR=None, OUTFILE=None, obs_index_nml=None, odb_index_nml=None):
        if ROSE_SUITE_DIR is not None: self.path_suite=ROSE_SUITE_DIR
        if OUTFILE is not None: self.outfile=OUTFILE
        if mainbox is not None: 
            self.mainbox=mainbox
        else:
            self.mainbox=tkinter.Tk()
        self.path=self.path_suite
        self.mainbox.title("Monitobs.2019.04")
        self.cylctime="20200520T0000Z"
        self.srcopt=None
        self.obstypelist=None
        self.infile=None
        self.obstype=None
        self.selectlist=None
        self.userquery=None
        self.subtypelist=None
        self.stypindx=None
        self.enamlist=None
        self.ptr_styp=None
        self.subtype=None
        self.eselect=None
        self.ename=None
        self.eselcount=0
        self.obs_index_nml=obs_index_nml
        self.odb_index_nml=odb_index_nml
        self.frm1 = tkinter.Frame(self.mainbox) 
        self.frm1.grid(row=0)
        self.frm11 = tkinter.Frame(self.frm1) 
        self.frm11.grid(row=0, column=0)
        self.frm12 = tkinter.Frame(self.frm1) 
        self.frm12.grid(row=0, column=1)
                
        self.frm20 = tkinter.Frame(self.frm11) 
        self.frm20.grid(row=0, column=0)
        self.frm21 = tkinter.Frame(self.frm11) 
        self.frm21.grid(row=1, column=0)
        self.frm22 = tkinter.Frame(self.frm11) 
        self.frm22.grid(row=2, column=0)
        self.frm23 = tkinter.Frame(self.frm11) 
        self.frm23.grid(row=3, column=0)
        self.frm24 = tkinter.Frame(self.frm11) 
        self.frm24.grid(row=4, column=0)
        
        self.frm_ctime = tkinter.Frame(self.frm21) 
        self.frm_ctime.grid(row=0, column=1)

        self.frm31 = tkinter.Frame(self.frm22) 
        self.frm31.grid(row=0, column=0)
        self.frm32 = tkinter.Frame(self.frm22)
        self.frm32.grid(row=0, column=1)
        self.frm33 = tkinter.Frame(self.frm22) 
        self.frm33.grid(row=0, column=2)
        
        self.frm41 = tkinter.Frame(self.frm31) 
        self.frm41.grid(row=0, column=1)
        self.frm42 = tkinter.Frame(self.frm31) 
        self.frm42.grid(row=2, column=1)
        self.frm43 = tkinter.Frame(self.frm31) 
        self.frm43.grid(row=5, column=1)
        
        self.frm7 = tkinter.Frame(self.mainbox) 
        self.frm7.grid(row=4, column=0)

        
        font1 = tkFont.Font(family="Comic",weight="bold",overstrike=True)
        lbl_header=tkinter.Label(self.frm20, text="MonitObs", fg = "red", font = "font1 50 bold").pack()
        lbl_header=tkinter.Label(self.frm20, text="[  An Obrervation monitoring and database management system developed at NCMRWF, MoES, India  ]", fg = "blue", font = "font1 8 bold italic").pack()
        
        self.lbl_inpath=tkinter.Label(self.frm21, text='Path')
        self.lbl_inpath.grid(row=0,sticky="w")
        self.ebx_path = tkinter.Entry(self.frm21, width=100) 
        self.ebx_path.grid(row=0, column=1,sticky="w")
        self.ebx_path.insert(tkinter.END, self.path)

        
        
#        self.Lb_selectlist = Listbox(self.frm2) 
#        self.Lb.insert(1, 'Surface') 
#        self.Lb.insert(2, 'Sonde') 
#        self.Lb.insert(3, 'Aircraft') 
#        self.Lb.insert(4, 'Satwind') 
        
        
        self.btn_clear=tkinter.Button(self.frm43, text='Clear', width=10, command=lambda:self.clearall() )
        self.btn_clear.grid(row=0, column=0)
        
        
        #self.display = tkinter.Text(self.frm12, height=20, width=80) 
        #self.display.grid(row=0, column=0)
        self.display = tkinter.Text(self.frm7, height=20, width=200) 
        self.display.grid(row=0, column=0)
        #self.display.insert(END, "")
        #self.plot()
        
	self.ctime_menu()
        self.mainloop()
        



