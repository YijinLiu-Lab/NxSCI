'''
Created on 2 May 2012

XMFlex for Python: Toolbox for collecting X-ray transmission image
                       data using Xradia TXM system.

@author:
Mirna Lerotic
2nd Look Consulting
http://www.2ndlookconsulting.com/
'''

from __future__ import division
import wx
import wx.lib.intctrl
import wx.lib.masked.numctrl
import numpy as np

import copy
import os
import time
import datetime
import cPickle as pickle

# from pywinauto.application import Application
from time import gmtime, strftime

import matplotlib as mplot

mplot.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas

import logging

logging.basicConfig(filename='XMFlexPy.log', filemode='w', level=logging.DEBUG)
logging.getLogger(__name__)

import txrm
import scriptgen
import ssrl_sis

__version__ = '1.2.11'

# Verbose
v = False

# Use timer
UseTimer = False

# Debugging
DebugYieldON = True
DebugTextCtrlON = True

# Debug limiting textctrl size while executing a scan
TCMaxLines = 500

SettingsFileName = 'XMF_config.txt'

# Temporary fix for mono movement. Wait time in seconds
DelayMono = False
mono_wait_time = 60  # secs
mono_threshold = 50  # eV

XradiaAxisNameList = ['x', 'y', 'z', 't', 'zpx', 'zpy', 'zpz',
                      'cdx', 'cdy', 'cdz', 'cdpitch', 'cdyaw',
                      'prx', 'pry', 'prz', 'phx', 'phy', 'phz',
                      'energy', 'detx', 'dety', 'detz']


# #----------------------------------------------------------------------
# def calc_R(xc, yc, x, y):
#     """ calculate the distance of each 2D points from the center (xc, yc) """
#     return sqrt((x-xc)**2 + (y-yc)**2)
#
# def f_2(c,x,y):
#     """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
#     xc,yc = c
#     Ri = calc_R(xc, yc, x, y)
#     return Ri - Ri.mean()
#
# #----------------------------------------------------------------------
# #   Fits a circle in x,y plane
# # Result is center point (yc,xc) and radius R.
# def CircleFit(x, y):
#
#     x_m = np.mean(x)
#     y_m = np.mean(y)
#     center_estimate = x_m, y_m
#     center_2, ier = optimize.leastsq(f_2, center_estimate, args=(x,y))
#
#     xc_2, yc_2 = center_2
#     Ri_2       = calc_R(xc_2, yc_2)
#     R_2        = Ri_2.mean()
#     residu_2   = sum((Ri_2 - R_2)**2)

# ----------------------------------------------------------------------
#   Fits a circle in x,y plane
# Result is center point (yc,xc) and radius R.  A is an optional
# output describing the circle's equation:
#   x**2+y**2+a(1)*x+a(2)*y+a(3)=0
def CircleFit(x, y):
    n = x.size
    xx = x * x
    yy = y * y
    xy = x * y
    A = np.array([[np.sum(x), np.sum(y), n],
                  [np.sum(xy), np.sum(yy), np.sum(y)],
                  [np.sum(xx), np.sum(xy), np.sum(x)]])
    B = np.array([[-np.sum(xx + yy)], [-np.sum(xx * y + yy * y)], [-np.sum(xx * x + xy * y)]])
    a = np.linalg.solve(A, B)
    a.resize((3))
    xc = -0.5 * a[0];
    yc = -0.5 * a[1];
    R = np.sqrt((a[0] ** 2 + a[1] ** 2) / 4 - a[2]);

    return xc, yc, R


# ----------------------------------------------------------------------
class ScriptData:
    def __init__(self):
        self.enable_energy = 0
        self.enable_tomo = 0
        self.enable_mosaic = 0
        self.enable_multi = 0

        self.wait_energy = 0
        self.wait_tomo = 0
        self.wait_mosaic = 0
        self.wait_multi = 0

        self.loadedfromscript = 0

        self.NRepeatScan = 1
        self.waitsecs = 0

        # Main loop layers
        self.layer_eng = 1
        self.layer_ang = 2
        self.layer_mos = 3
        self.layer_me = 4

        # Single TOMO image
        self.singletomo = 0

        # Energy Point List
        self.EPlist = []
        self.EPkeys = ['Energy', 'ZPx', 'ZPy', 'ZPz', 'ZPID', 'Detx', 'Dety', 'Detz']

        # Energy Region Params:
        self.RegionParams = np.zeros((1, 7))
        self.EngRegionlist = []
        self.ERegkeys = ['Start', 'Steps', 'EngStep', 'Stop', 'ZPxStart', 'ZPyStart', 'ZPzStart']
        self.DetectorPars = []

        self.defaultregion = dict(zip(self.ERegkeys, [0.0, 1, 0.5, 0.5, 0.0, 0.0, 0.0]))
        self.EngRegionlist.append(dict(zip(self.ERegkeys, [0.0, 1, 0.5, 0.5, 0.0, 0.0, 0.0])))

        self.sortengs = True

        # Multi-exposure
        self.n_exposures = 1
        self.averageotfly = 0

        # Scan Settings:
        self.Exptime = 1.0
        self.useionchamber = 0
        self.Binning = 1
        self.SampleName = 'SampleName'
        self.ScriptName = ' '
        self.OutputPath = ''

        # Reting frame adding
        self.n_retiga = 1

        # Zoneplate Params:
        self.ZP_EngList = []
        self.ZP_xList = []
        self.ZP_yList = []
        self.ZP_zList = []

        # RefImg Params:
        self.RepeatRefExposure = 0
        self.ref4everycertainnumofexp = 0
        self.RefX = 0
        self.RefY = 0
        self.RefZ = 0
        self.useRefZ = False
        self.RefTheta = 0
        self.useRefTheta = False
        self.RefBinning = 1
        self.RefExpTime = 1.0
        self.ABBA = 0
        self.refaverageotfly = 0
        self.refrelativepos = 0

        # Tomo Sample Stage Coordinates
        self.XCorr = 0.0
        self.ZCorr = 0.0
        self.SampleXinZeroDegree = 0.0
        self.SampleYinZeroDegree = 0.0
        self.SampleZinZeroDegree = 0.0
        self.samplexs = []
        self.sampleys = []
        self.samplezs = []
        self.samplets = []

        # Mosaic Params:
        self.mosaic_left = 1
        self.mosaic_right = 1
        self.mosaic_up = 1
        self.mosaic_down = 1
        self.mosaic_overlap = 0.2
        self.mosaic_centraltile = 1

        # XRM file info
        self.pixelsize = 0.014
        self.width = 100
        self.height = 100
        self.defaultpixelsize = 0.014
        self.FOV = 0.0


""" ------------------------------------------------------------------------------------------------"""


class ResumeScanFrame(wx.Frame):

    def __init__(self):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title="Resume Uncompleted Scan", size=(850, 700))

        self.mainframe = wx.GetApp().TopWindow

        self.SetBackgroundColour("White")

        self.dir = ''
        self.scriptlist = []
        self.xrm_files = []

        vboxtop = wx.BoxSizer(wx.VERTICAL)
        vboxtop.Add((0, 20))
        hboxtop = wx.BoxSizer(wx.HORIZONTAL)

        hboxtop.Add((20, 0))

        panel = wx.Panel(self, -1)
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        vbox = wx.BoxSizer(wx.VERTICAL)
        # vbox.Add((0,50))

        button_1 = wx.Button(panel, -1, 'Load Uncompleted Script')
        self.Bind(wx.EVT_BUTTON, self.OnLoadScript, id=button_1.GetId())
        vbox.Add(button_1, 0)
        vbox.Add((0, 20))

        button_2 = wx.Button(panel, -1, 'Select Directory with Images')
        self.Bind(wx.EVT_BUTTON, self.OnSelectDir, id=button_2.GetId())
        vbox.Add(button_2, 0)
        vbox.Add((0, 10))

        self.tc_path1 = wx.TextCtrl(panel, -1, size=((400, -1)), style=wx.TE_READONLY | wx.TE_RICH)
        vbox.Add(self.tc_path1, 0)
        vbox.Add((0, 10))

        self.tc_files = wx.TextCtrl(panel, -1, style=wx.TE_MULTILINE | wx.TE_READONLY |
                                                     wx.TE_RICH | wx.TE_BESTWRAP)
        self.tc_files.SetMinSize((250, 450))
        vbox.Add(self.tc_files, 1, wx.EXPAND)

        vbox.Add((0, 20))
        hbox.Add(vbox, 1, wx.EXPAND)
        panel.SetSizer(hbox)

        hboxtop.Add(panel, -1)
        hboxtop.Add((20, 0))

        # panel 2
        panel2 = wx.Panel(self, -1)
        vbox21 = wx.BoxSizer(wx.VERTICAL)

        self.tc_scriptlist = wx.TextCtrl(panel2, -1, style=wx.TE_MULTILINE | wx.TE_READONLY |
                                                           wx.TE_RICH | wx.VSCROLL | wx.HSCROLL)
        self.tc_scriptlist.SetMinSize((250, 565))
        vbox21.Add(self.tc_scriptlist, 1, wx.EXPAND)

        panel2.SetSizer(vbox21)
        hboxtop.Add(panel2, -1, )

        hboxtop.Add((20, 0))

        vboxtop.Add(hboxtop, 0, wx.EXPAND)

        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        hbox2.Add((20, 0))
        button_resume = wx.Button(self, -1, 'Resume Scan', size=(150, -1))
        self.Bind(wx.EVT_BUTTON, self.OnResume, id=button_resume.GetId())
        hbox2.Add(button_resume, 0)

        hbox2.Add((20, 0))

        button_close = wx.Button(self, -1, 'Cancel', size=(150, -1))
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=button_close.GetId())
        hbox2.Add(button_close, 0)

        vboxtop.Add(hbox2, 0, wx.EXPAND)

        self.SetSizer(vboxtop)

    # ----------------------------------------------------------------------
    def OnSelectDir(self, event):
        """
        Browse for directory with images
        """

        try:
            dialog = wx.DirDialog(None, "Choose a directory with images",
                                  style=wx.DD_DIR_MUST_EXIST | wx.DD_CHANGE_DIR)
            # dialog.SetWildcard(wildcard)
            if dialog.ShowModal() == wx.ID_OK:
                directory = dialog.GetPath()

            self.dir = directory

            self.tc_path1.SetValue(self.dir)

            self.ListFiles()


        except:
            print
            'Error could select a directory.'
            wx.MessageBox("Error could select a directory.")
            import sys;
            print
            sys.exc_info()

        # ----------------------------------------------------------------------

    def OnLoadScript(self, evt):

        wildcard = "TXT files (*.txt)|*.txt"
        dialog = wx.FileDialog(None, "Please select a script file",
                               wildcard=wildcard,
                               style=wx.OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            filepath = dialog.GetPath()
            filedir = dialog.GetDirectory()
        else:
            print
            'ERROR - cannot open file.'
            return

        f = open(str(filepath), 'r')
        f
        scriptlist = []

        for line in f:
            line = line.rstrip('\r\n')

            if line == '':
                continue

            if line.startswith(';'):
                comment = line.replace(';', '')
                script = ';;; %s' % (comment)
                scriptlist.append(script)
                continue

            thiscommand = line.split(' ')

            if thiscommand[0] == 'setexp':
                exposure = float(thiscommand[1])
                script = '%s %.3f' % ('setexp', exposure)
                scriptlist.append(script)
                if v: print
                'set exposure ', exposure

            elif thiscommand[0] == 'setbinning':
                binning = int(thiscommand[1])
                script = '%s %.0f' % ('setbinning', binning)
                scriptlist.append(script)
                if v: print
                'set binning ', binning

            elif thiscommand[0] == 'setfilter':
                filter = int(thiscommand[1])
                script = '%s %.0f' % ('setfilter', filter)
                scriptlist.append(script)
                if v: print
                'set setfilter ', filter


            elif thiscommand[0] == 'sete':
                eng = float(thiscommand[1])
                script = 'sete %.3f' % (eng)
                scriptlist.append(script)
                if v: print
                'set energy ', eng


            elif thiscommand[0] == 'moveto':
                axisname = thiscommand[1]
                position = float(thiscommand[2])
                if axisname == 'E':
                    script = 'sete %.3f' % (position)
                    scriptlist.append(script)
                    if v: print
                    'move energy to %f' % (position)
                else:
                    script = 'moveto %s %.3f' % (axisname, position)
                    scriptlist.append(script)
                    if v: print
                    'move motor %s to %f' % (axisname, position)


            elif thiscommand[0] == 'moveby':
                axisname = thiscommand[1]
                position = float(thiscommand[2])
                script = 'moveby %s %.3f' % (axisname, position)
                scriptlist.append(script)
                if v: print
                'move motor %s by %f' % (axisname, position)

            elif thiscommand[0] == 'collectxrf':
                filename = thiscommand[1]
                script = '%s %s' % ('collectxrf', filename)
                scriptlist.append(script)
                if v: print
                'collectxrf ', filename

            elif thiscommand[0] == 'collect':
                filename = thiscommand[1]
                script = '%s %s' % ('collect', filename)
                scriptlist.append(script)
                if v: print
                'collect ', filename

            elif thiscommand[0] == 'collectaverage':
                filename = thiscommand[1]
                totalnimages = thiscommand[2]
                script = '%s %s %s' % ('collectaverage', filename, str(totalnimages))
                scriptlist.append(script)
                if v: print
                'collectaverage ', filename, totalnimages

            elif thiscommand[0] == 'collectnormalized':
                filename = thiscommand[1]
                script = '%s %s' % ('collectnormalized', filename)
                scriptlist.append(script)
                if v: print
                'collectnormalized ', filename

            elif thiscommand[0] == 'collecttomo':
                filename = thiscommand[1]
                totalnimages = thiscommand[2]
                startangle = thiscommand[3]
                endangle = thiscommand[4]
                script = '%s %s %s %s %s' % ('collecttomo', filename, totalnimages, startangle, endangle)
                scriptlist.append(script)
                if v: print
                'collecttomo ', filename, totalnimages, startangle, endangle

            elif thiscommand[0] == 'collectnormalizedtomo':
                filename = thiscommand[1]
                totalnimages = thiscommand[2]
                startangle = thiscommand[3]
                endangle = thiscommand[4]
                script = '%s %s %s %s %s' % ('collectnormalizedtomo', filename, totalnimages, startangle, endangle)
                scriptlist.append(script)
                if v: print
                'collectnormalizedtomo ', filename, totalnimages, startangle, endangle

            elif thiscommand[0] == 'wait':
                waittimesecs = float(thiscommand[1])
                script = 'wait %f' % (waittimesecs)
                scriptlist.append(script)
                if v: print
                'wait %f' % (waittimesecs)

            elif thiscommand[0] == 'recordioncurrent':
                script = 'recordioncurrent'
                scriptlist.append(script)
                if v: print
                'recordioncurrent'

            elif thiscommand[0] == 'setreadout':
                readout = int(thiscommand[1])
                script = 'setreadout %d' % (readout)
                scriptlist.append(script)
                if v: print
                'setreadout %d' % (readout)


            elif thiscommand[0] == 'shiftobjective':
                nobjective = int(thiscommand[1])
                script = 'shiftobjective %d' % (nobjective)
                scriptlist.append(script)
                if v: print
                'shiftobjective %d' % (nobjective)

            elif thiscommand[0] == 'runshortcut':
                command = ' '.join(thiscommand[1:])
                script = 'runshortcut %s' % (command)
                scriptlist.append(script)
                if v: print
                'runshortcut %d' % (command)

            else:
                print
                'ERROR - unrecognized command in the script.'
                logging.error('ERROR - unrecognized command in the script.')

        f.close()

        self.scriptlist = scriptlist

        self.tc_scriptlist.Clear()

        for i in range(len(self.scriptlist)):
            self.tc_scriptlist.AppendText('{0:s}\n'.format(self.scriptlist[i]))

        self.tc_scriptlist.SetInsertionPoint(0)

        self.SetFocus()

        self.UpdateScript()

    # ----------------------------------------------------------------------
    def ListFiles(self):

        self.xrm_files = [x for x in os.listdir(self.dir) if x.endswith('xrm')]

        self.tc_files.Clear()

        for i in range(len(self.xrm_files)):
            self.tc_files.AppendText('{0:s}\n'.format(self.xrm_files[i]))

        self.tc_files.SetInsertionPoint(0)

        self.SetFocus()

        self.UpdateScript()

    # ----------------------------------------------------------------------
    def UpdateScript(self):

        if (len(self.scriptlist) == 0) or (len(self.xrm_files) == 0):
            return

        # look at the files in the folder and find the one that was taken last
        lastimg_line = 0
        for i in range(len(self.xrm_files)):
            for j in range(len(self.scriptlist)):
                if self.xrm_files[i] in self.scriptlist[j]:
                    thisline = j
                    if thisline > lastimg_line:
                        lastimg_line = thisline

        # Check is sete is present in last 10 lines - if yes repeat from sete
        for i in range(lastimg_line - 1, max(lastimg_line - 10, 0), -1):
            if 'sete ' in self.scriptlist[i]:
                lastimg_line = i
                break

        # Remove all the collect image calls before last line
        added_commands = []
        for i in range(lastimg_line - 1, -1, -1):
            if 'collect' in self.scriptlist[i]:
                self.scriptlist.pop(i)

            else:
                thiscommand = self.scriptlist[i].split(' ')
                if thiscommand[0] == 'moveto':
                    moveaxis = thiscommand[0] + thiscommand[1]
                    if moveaxis not in added_commands:
                        added_commands.append(moveaxis)
                    else:
                        self.scriptlist.pop(i)
                else:
                    if thiscommand[0] not in added_commands:
                        added_commands.append(thiscommand[0])
                    else:
                        self.scriptlist.pop(i)

        self.tc_scriptlist.Clear()

        for i in range(len(self.scriptlist)):
            self.tc_scriptlist.AppendText('{0:s}\n'.format(self.scriptlist[i]))

        self.tc_scriptlist.SetInsertionPoint(0)

    # ----------------------------------------------------------------------
    def OnResume(self, evt):
        self.mainframe.ResumeUncompletedScan(self.scriptlist, self.dir)
        self.Close(True)

    # ----------------------------------------------------------------------
    def OnClose(self, evt):

        self.Close(True)


""" ------------------------------------------------------------------------------------------------"""


class MosaicFrame(wx.Frame):

    def __init__(self, scriptdata):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title="Mosaic Scan Info", size=(530, 400))

        self.scriptdata = scriptdata
        self.mainframe = wx.GetApp().TopWindow

        self.SetBackgroundColour("White")

        vboxtop = wx.BoxSizer(wx.VERTICAL)

        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
        hbox = wx.BoxSizer(wx.HORIZONTAL)

        # Mosaic Image
        self.mosaicimage = np.zeros((300, 300, 3), dtype=np.uint8)

        mimage = wx.ImageFromData(300, 300, self.mosaicimage)
        bitmosaic = wx.BitmapFromImage(mimage)
        self.imageCtrl = wx.StaticBitmap(panel, wx.ID_ANY, bitmosaic)

        # Central Tile

        sizer2 = wx.StaticBoxSizer(wx.StaticBox(panel, -1, 'Central Tile Options'), orient=wx.VERTICAL)
        vbox20 = wx.BoxSizer(wx.VERTICAL)
        vbox20.Add((0, 10))

        self.rb_wctile = wx.RadioButton(panel, -1, 'With Central Tile', style=wx.RB_GROUP)
        self.rb_woctile = wx.RadioButton(panel, -1, 'Without Central Tile')
        self.Bind(wx.EVT_RADIOBUTTON, self.onrb_withctile, id=self.rb_wctile.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.onrb_withctile, id=self.rb_woctile.GetId())

        if (self.scriptdata.mosaic_centraltile == 1):
            self.rb_wctile.SetValue(True)
            self.rb_woctile.SetValue(False)
        else:
            self.rb_wctile.SetValue(False)
            self.rb_woctile.SetValue(True)

        vbox20.Add(self.rb_wctile, 0, wx.EXPAND)
        vbox20.Add((0, 5))
        vbox20.Add(self.rb_woctile, 0, wx.EXPAND)
        vbox20.Add((0, 10))

        sizer2.Add(vbox20, 1, wx.LEFT | wx.RIGHT | wx.EXPAND, 10)

        # Mosaic Design

        sizer3 = wx.StaticBoxSizer(wx.StaticBox(panel, -1, 'Mosaic Design'), orient=wx.VERTICAL)
        vbox31 = wx.BoxSizer(wx.VERTICAL)
        vbox31.Add((0, 10))

        gridtop = wx.FlexGridSizer(cols=2, vgap=5, hgap=10)

        text1 = wx.StaticText(panel, label="Up")
        ebox1 = wx.lib.masked.NumCtrl(panel, size=(70, -1),
                                      style=wx.TE_CENTRE,
                                      integerWidth=7, fractionWidth=0,
                                      min=0, max=100,
                                      value=self.scriptdata.mosaic_up)

        self.Bind(wx.lib.masked.EVT_NUM, self.OnUp, ebox1)

        text2 = wx.StaticText(panel, label="Down")
        ebox2 = wx.lib.masked.NumCtrl(panel, size=(70, -1),
                                      style=wx.TE_CENTRE,
                                      integerWidth=7, fractionWidth=0,
                                      min=0, max=100,
                                      value=self.scriptdata.mosaic_down)

        self.Bind(wx.lib.masked.EVT_NUM, self.OnDown, ebox2)

        text3 = wx.StaticText(panel, label="Left")
        ebox3 = wx.lib.masked.NumCtrl(panel, size=(70, -1),
                                      style=wx.TE_CENTRE,
                                      integerWidth=7, fractionWidth=0,
                                      min=0, max=100,
                                      value=self.scriptdata.mosaic_left)

        self.Bind(wx.lib.masked.EVT_NUM, self.OnLeft, ebox3)

        text4 = wx.StaticText(panel, label="Right")
        ebox4 = wx.lib.masked.NumCtrl(panel, size=(70, -1),
                                      style=wx.TE_CENTRE,
                                      integerWidth=7, fractionWidth=0,
                                      min=0, max=100,
                                      value=self.scriptdata.mosaic_right)

        self.Bind(wx.lib.masked.EVT_NUM, self.OnRight, ebox4)

        text5 = wx.StaticText(panel, label="Overlap")
        ebox5 = wx.lib.masked.numctrl.NumCtrl(panel, -1,
                                              value=self.scriptdata.mosaic_overlap,
                                              integerWidth=3,
                                              fractionWidth=3, min=0, max=1,
                                              limited=True,
                                              allowNone=True)

        self.Bind(wx.lib.masked.EVT_NUM, self.OnOverlap, ebox5)

        gridtop.Add(text1, 0)
        gridtop.Add(ebox1, 1)
        gridtop.Add(text2, 0)
        gridtop.Add(ebox2, 1)
        gridtop.Add(text3, 0)
        gridtop.Add(ebox3, 1)
        gridtop.Add(text4, 0)
        gridtop.Add(ebox4, 1)
        gridtop.Add(text5, 0)
        gridtop.Add(ebox5, 1)

        vbox31.Add(gridtop)

        sizer3.Add(vbox31, 1, wx.LEFT | wx.RIGHT | wx.EXPAND, 10)

        hbox.Add(self.imageCtrl, 1, wx.LEFT | wx.RIGHT, 20)

        vboxr = wx.BoxSizer(wx.VERTICAL)
        vboxr.Add(sizer2, 0, wx.EXPAND)
        vboxr.Add((0, 15))
        vboxr.Add(sizer3, 0, wx.EXPAND)
        hbox.Add(vboxr, 0, wx.EXPAND)

        vbox.Add((0, 20))
        vbox.Add(hbox)

        button_close = wx.Button(panel, -1, 'Submit')
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=button_close.GetId())
        vbox.Add(button_close, 0, wx.RIGHT | wx.ALIGN_RIGHT, 20)

        panel.SetSizer(vbox)

        vboxtop.Add(panel, 1, wx.EXPAND)

        self.SetSizer(vboxtop)

        self.DrawMosaicImage()

    # ----------------------------------------------------------------------
    def onrb_withctile(self, evt):
        state = self.rb_wctile.GetValue()

        if state:
            self.scriptdata.mosaic_centraltile = 1
        else:
            self.scriptdata.mosaic_centraltile = 0

        self.DrawMosaicImage()

    # ----------------------------------------------------------------------
    def OnUp(self, event):
        ctl = event.GetEventObject()
        value = ctl.GetValue()

        if (value < 0) or (value > 100):
            return

        self.scriptdata.mosaic_up = value

        if v: print
        self.scriptdata.mosaic_up

        if self.scriptdata.mosaic_up != None:
            self.DrawMosaicImage()

    # ----------------------------------------------------------------------
    def OnDown(self, event):
        ctl = event.GetEventObject()

        value = ctl.GetValue()

        if (value < 0) or (value > 100):
            return

        self.scriptdata.mosaic_down = value
        if v: print
        self.scriptdata.mosaic_down

        if self.scriptdata.mosaic_down != None:
            self.DrawMosaicImage()

    # ----------------------------------------------------------------------
    def OnLeft(self, event):
        ctl = event.GetEventObject()
        value = ctl.GetValue()

        if (value < 0) or (value > 100):
            return
        self.scriptdata.mosaic_left = value
        if v: print
        self.scriptdata.mosaic_left

        if self.scriptdata.mosaic_left != None:
            self.DrawMosaicImage()

    # ----------------------------------------------------------------------
    def OnRight(self, event):
        ctl = event.GetEventObject()
        value = ctl.GetValue()

        if (value < 0) or (value > 100):
            return
        self.scriptdata.mosaic_right = value
        if v: print
        self.scriptdata.mosaic_right

        if self.scriptdata.mosaic_right != None:
            self.DrawMosaicImage()

    # ----------------------------------------------------------------------
    def OnOverlap(self, event):
        ctl = event.GetEventObject()

        self.scriptdata.mosaic_overlap = ctl.GetValue()

        self.DrawMosaicImage()

    # ----------------------------------------------------------------------
    def CheckMosaic(self):
        # Check if there are any tiles left if the change is made

        havetile = True

        if (self.scriptdata.mosaic_centraltile == 0) and (self.scriptdata.mosaic_up == 0) and \
                (self.scriptdata.mosaic_down == 0) and (self.scriptdata.mosaic_left == 0) and \
                (self.scriptdata.mosaic_right == 0):
            havetile = False

        return havetile

    # ----------------------------------------------------------------------
    def DrawRect(self, x1, y1, x2, y2, color):

        if color == 'g':
            r = 0
            g = 1
            b = 0
        elif color == 'y':
            r = 1
            g = 1
            b = 0
        elif color == 'k':
            r = 0
            g = 0
            b = 0

        for x in range(int(x1), int(x2)):
            for y in range(int(y1), int(y2)):
                self.mosaicimage[y, x, 0] = 255 * r
                self.mosaicimage[y, x, 1] = 255 * g
                self.mosaicimage[y, x, 2] = 255 * b

    # ----------------------------------------------------------------------
    def DrawMosaicImage(self):

        if self.CheckMosaic() == False:
            self.scriptdata.mosaic_centraltile = 1
            self.rb_wctile.SetValue(True)

        self.mosaicimage = np.zeros((300, 300, 3), dtype=np.uint8)
        self.mosaicimage[:, ] = 255
        mimage = wx.ImageFromData(300, 300, self.mosaicimage)
        bmp = mimage.ConvertToBitmap()
        self.imageCtrl.SetBitmap(bmp)

        if (self.scriptdata.mosaic_centraltile == 1):
            NColumns = self.scriptdata.mosaic_left + self.scriptdata.mosaic_right + 1
            NRows = self.scriptdata.mosaic_up + self.scriptdata.mosaic_down + 1
            startingx = (
                                    self.scriptdata.mosaic_left + 0.5) * 10 - self.scriptdata.mosaic_left * self.scriptdata.mosaic_overlap * 10
            startingy = (
                                    self.scriptdata.mosaic_down + 0.5) * 10 - self.scriptdata.mosaic_down * self.scriptdata.mosaic_overlap * 10
        else:
            NColumns = self.scriptdata.mosaic_left + self.scriptdata.mosaic_right
            NRows = self.scriptdata.mosaic_up + self.scriptdata.mosaic_down
            startingx = (
                            self.scriptdata.mosaic_left) * 10 - self.scriptdata.mosaic_left * self.scriptdata.mosaic_overlap * 10 + \
                        10 / 2 * self.scriptdata.mosaic_overlap
            startingy = (
                            self.scriptdata.mosaic_down) * 10 - self.scriptdata.mosaic_down * self.scriptdata.mosaic_overlap * 10 + \
                        10 / 2 * self.scriptdata.mosaic_overlap
            if ((self.scriptdata.mosaic_up == 0) or (self.scriptdata.mosaic_down == 0)):
                if self.scriptdata.mosaic_up == 0:
                    diffff = -1
                else:
                    diffff = 1

                startingy = startingy + diffff * 5 - 10 / 2 * self.scriptdata.mosaic_overlap * diffff

            if ((self.scriptdata.mosaic_left == 0) or (self.scriptdata.mosaic_right == 0)):
                if self.scriptdata.mosaic_left == 0:
                    diffff = 1
                else:
                    diffff = -1

                startingx = startingx + diffff * 5 - 10 / 2 * self.scriptdata.mosaic_overlap * diffff

        moswidth = (NRows - 1) * 10 * (1 - self.scriptdata.mosaic_overlap) + 10
        mosheight = (NColumns - 1) * 10 * (1 - self.scriptdata.mosaic_overlap) + 10
        mossize = max(moswidth, mosheight)
        self.mosaicimage = np.zeros((int(mossize), int(mossize), 3), dtype=np.uint8)
        self.mosaicimage[:, ] = 255

        for x in np.arange(0, ((NColumns - 1) * 10 * (1 - self.scriptdata.mosaic_overlap) + 1),
                           (10 * (1 - self.scriptdata.mosaic_overlap))):
            for y in np.arange(0, ((NRows - 1) * 10 * (1 - self.scriptdata.mosaic_overlap) + 1),
                               (10 * (1 - self.scriptdata.mosaic_overlap))):
                if np.mod((x + y) / (10.0 * (1.0 - self.scriptdata.mosaic_overlap)), 2):
                    rectcolor = 'g'
                else:
                    rectcolor = 'y'
                self.DrawRect(x, y, x + 10, y + 10, rectcolor)

        overlapcolor = 'k'
        for x in np.arange(0, ((NColumns - 1) * 10 * (1 - self.scriptdata.mosaic_overlap)) + 1,
                           (10 * (1 - self.scriptdata.mosaic_overlap))):
            for y in np.arange(0, ((NRows - 1) * 10 * (1 - self.scriptdata.mosaic_overlap)) + 1,
                               (10 * (1 - self.scriptdata.mosaic_overlap))):
                if x > 0:
                    self.DrawRect(x, y, int(x + 10 * (self.scriptdata.mosaic_overlap)), y + 10, overlapcolor)

                if y > 0:
                    self.DrawRect(x, y, x + 10, int(y + 10 * (self.scriptdata.mosaic_overlap)), overlapcolor)

        startingx = int(startingx)
        startingy = int(startingy)
        # Draw red rectangle around central tile
        for x in range(int(startingx - 5), int(startingx + 5)):
            self.mosaicimage[startingy - 5, x, 0] = 255
            self.mosaicimage[startingy - 5, x, 1] = 0
            self.mosaicimage[startingy - 5, x, 2] = 0
            self.mosaicimage[startingy + 4, x, 0] = 255
            self.mosaicimage[startingy + 4, x, 1] = 0
            self.mosaicimage[startingy + 4, x, 2] = 0
        for y in range(int(startingy - 5), int(startingy + 5)):
            self.mosaicimage[y, startingx - 5, 0] = 255
            self.mosaicimage[y, startingx - 5, 1] = 0
            self.mosaicimage[y, startingx - 5, 2] = 0
            self.mosaicimage[y, startingx + 4, 0] = 255
            self.mosaicimage[y, startingx + 4, 1] = 0
            self.mosaicimage[y, startingx + 4, 2] = 0

            # Draw red dot in the center of the mosaic
        self.mosaicimage[startingy, startingx - 1, 0] = 255
        self.mosaicimage[startingy, startingx - 1, 1] = 0
        self.mosaicimage[startingy, startingx - 1, 2] = 0
        self.mosaicimage[startingy, startingx, 0] = 255
        self.mosaicimage[startingy, startingx, 1] = 0
        self.mosaicimage[startingy, startingx, 2] = 0
        self.mosaicimage[startingy - 1, startingx - 1, 0] = 255
        self.mosaicimage[startingy - 1, startingx - 1, 1] = 0
        self.mosaicimage[startingy - 1, startingx - 1, 2] = 0
        self.mosaicimage[startingy - 1, startingx, 0] = 255
        self.mosaicimage[startingy - 1, startingx, 1] = 0
        self.mosaicimage[startingy - 1, startingx, 2] = 0

        mimage = wx.ImageFromData(mossize, mossize, self.mosaicimage)

        neww = 300
        newh = 300

        mimage = mimage.Scale(neww, newh)

        # Flip image vertically because Matlab and wxBitmap have different (0,0) positions
        mimage = mimage.Mirror(horizontally=False)

        bmp = mimage.ConvertToBitmap()

        self.imageCtrl.SetBitmap(bmp)

    # ----------------------------------------------------------------------
    def OnClose(self, evt):
        self.mainframe.ShowEnabledTC()
        self.Close(True)


""" ------------------------------------------------------------------------------------------------"""


class AngleFrame(wx.Frame):

    def __init__(self, scriptdata, microscope, SampleXYZAboveRotation):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title="Tomo Scan Info", size=(430, 200))

        self.scriptdata = scriptdata
        self.microscope = microscope

        self.SampleXYZAboveRotation = SampleXYZAboveRotation

        self.mainframe = wx.GetApp().TopWindow

        self.SetBackgroundColour("White")

        vboxtop = wx.BoxSizer(wx.VERTICAL)

        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)

        # Sample Stage Coordinates Options
        self.usetxrm = True
        sizer = wx.StaticBoxSizer(wx.StaticBox(panel, -1, 'Select the source of sample stage coordinates'),
                                  orient=wx.VERTICAL)
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        hbox1.Add((0, 10))

        self.rb_txrmfile = wx.RadioButton(panel, -1, 'From an existing *.TXRM file', style=wx.RB_GROUP)
        self.rb_scalibration = wx.RadioButton(panel, -1, 'Info from stage calibration')

        self.Bind(wx.EVT_RADIOBUTTON, self.Onrb_sscoordinates, id=self.rb_txrmfile.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.Onrb_sscoordinates, id=self.rb_scalibration.GetId())

        hbox1.Add(self.rb_txrmfile, 0, wx.EXPAND)
        hbox1.Add((15, 0))
        hbox1.Add(self.rb_scalibration, 0, wx.EXPAND)

        sizer.Add(hbox1, 1, wx.LEFT | wx.RIGHT, 10)

        vbox.Add((0, 20))
        vbox.Add(sizer, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, 20)

        button_close = wx.Button(panel, -1, 'Continue')
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=button_close.GetId())
        vbox.Add(button_close, 0, wx.ALL | wx.ALIGN_RIGHT, 20)

        panel.SetSizer(vbox)

        vboxtop.Add(panel, 1, wx.EXPAND)

        self.SetSizer(vboxtop)

    # ----------------------------------------------------------------------
    def Onrb_sscoordinates(self, evt):
        state = self.rb_txrmfile.GetValue()

        if state:
            self.usetxrm = True
        else:
            self.usetxrm = False

    # ----------------------------------------------------------------------
    def OnClose(self, evt):

        if self.usetxrm:

            wildcard = "TXRM files (*.txrm)|*.txrm"
            dialog = wx.FileDialog(None, "Please select a txrm file",
                                   wildcard=wildcard,
                                   style=wx.OPEN)
            if dialog.ShowModal() == wx.ID_OK:
                filepath = dialog.GetPath()
                filename = dialog.GetFilename()
                filedir = dialog.GetDirectory()
            else:
                return

            if self.microscope:
                txrminfo = self.microscope.get_txrm_info(filepath)

                if txrminfo == None:
                    wx.MessageBox("Warning - Could not read .txrm file.")
                    return

                self.scriptdata.samplets = np.array(txrminfo['MotPosT'])
                self.scriptdata.samplexs = np.array(txrminfo['MotPosX'])
                self.scriptdata.sampleys = np.array(txrminfo['MotPosY'])
                self.scriptdata.samplezs = np.array(txrminfo['MotPosZ'])
                self.scriptdata.FOV = txrminfo['FOV']

                self.scriptdata.pixelsize = txrminfo['pixelsize']
                self.scriptdata.width = txrminfo['width']
                self.scriptdata.height = txrminfo['height']

            #                 print self.scriptdata.samplets
            #                 print self.scriptdata.samplexs
            #                 print 'self.scriptdata.pixelsize', self.scriptdata.pixelsize
            #                 print 'self.scriptdata.width', self.scriptdata.width
            #                 print 'self.scriptdata.height', self.scriptdata.height
            #                 print 'self.scriptdata.FOV', self.scriptdata.FOV

            else:

                Ctxrm = txrm.txrm()

                txrminfo = Ctxrm.get_txrm_info(filepath)

                NumofMotors, xxx, yyy, zzz, ttt, zpxxx, zpyyy, zpzzz = Ctxrm.get_motor_index(txrminfo)

                self.scriptdata.samplets = np.array(txrminfo['MotPos'][ttt::NumofMotors])
                self.scriptdata.samplexs = np.array(txrminfo['MotPos'][xxx::NumofMotors])
                self.scriptdata.sampleys = np.array(txrminfo['MotPos'][yyy::NumofMotors])
                self.scriptdata.samplezs = np.array(txrminfo['MotPos'][zzz::NumofMotors])
                self.scriptdata.FOV = txrminfo['FOV']

                self.scriptdata.pixelsize = txrminfo['pixelsize']
                # self.scriptdata.Binning = txrminfo['HBin']
                self.scriptdata.width = txrminfo['width']
                self.scriptdata.height = txrminfo['height']

            #                 print self.scriptdata.samplets
            #                 print self.scriptdata.samplexs
            #                 print 'self.scriptdata.pixelsize', self.scriptdata.pixelsize
            #                 print 'self.scriptdata.width', self.scriptdata.width
            #                 print 'self.scriptdata.height', self.scriptdata.height
            #                 print 'self.scriptdata.FOV', self.scriptdata.FOV

            zeroindex = np.argmin(np.abs(self.scriptdata.samplets))

            self.scriptdata.XCorr = 0.0
            self.scriptdata.ZCorr = 0.0
            self.scriptdata.SampleXinZeroDegree = float(self.scriptdata.samplexs[zeroindex])
            self.scriptdata.SampleYinZeroDegree = float(np.mean(self.scriptdata.sampleys))
            self.scriptdata.SampleZinZeroDegree = float(self.scriptdata.samplezs[zeroindex])

            self.scriptdata.OutputPath = filedir

            self.mainframe.ShowEnabledTC()

            self.mainframe.txrmloaded = True


        else:
            TomoSSCoordinatesFrame(self.scriptdata, self.microscope, self.SampleXYZAboveRotation).Show()
        self.Close(True)


# ------------------------------------------------------------------------------------------------
class TomoSSCoordinatesFrame(wx.Frame):

    def __init__(self, scriptdata, microscope, SampleXYZAboveRotation):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title="Tomo Sample Stage Coordinates", size=(330, 450))

        self.scriptdata = scriptdata
        self.microscope = microscope
        self.mainframe = wx.GetApp().TopWindow

        self.SampleXYZAboveRotation = SampleXYZAboveRotation

        self.SetBackgroundColour("White")

        vboxtop = wx.BoxSizer(wx.VERTICAL)

        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)

        vbox.Add((0, 10))

        grid1 = wx.FlexGridSizer(cols=2, vgap=5, hgap=15)

        text1 = wx.StaticText(panel, label="Sample X Correction:")
        self.ebox1 = wx.lib.masked.numctrl.NumCtrl(panel, style=wx.TE_CENTRE,
                                                   value=0.0, integerWidth=7,
                                                   fractionWidth=4,
                                                   allowNegative=True,
                                                   foregroundColour="Black",
                                                   signedForegroundColour="Black",
                                                   limited=False,
                                                   allowNone=True)
        # self.Bind(wx.lib.masked.EVT_NUM, self.OnEnergy, self.ebox1)
        if self.SampleXYZAboveRotation == 0:
            self.ebox1.Disable()

        text2 = wx.StaticText(panel, label="Sample Z Correction:")
        self.ebox2 = wx.lib.masked.numctrl.NumCtrl(panel, style=wx.TE_CENTRE,
                                                   value=0.0, integerWidth=7,
                                                   fractionWidth=4,
                                                   allowNegative=True,
                                                   foregroundColour="Black",
                                                   signedForegroundColour="Black",
                                                   limited=False,
                                                   allowNone=True)
        # self.Bind(wx.lib.masked.EVT_NUM, self.OnZPx, self.ebox2)
        if self.SampleXYZAboveRotation == 0:
            self.ebox2.Disable()

        text3 = wx.StaticText(panel, label="Sample X in 0 degree:")
        self.ebox3 = wx.lib.masked.numctrl.NumCtrl(panel, style=wx.TE_CENTRE,
                                                   value=0.0, integerWidth=7,
                                                   fractionWidth=4,
                                                   allowNegative=True,
                                                   foregroundColour="Black",
                                                   signedForegroundColour="Black",
                                                   limited=False,
                                                   allowNone=True)
        # self.Bind(wx.lib.masked.EVT_NUM, self.OnZPy, self.ebox3)

        text4 = wx.StaticText(panel, label="Sample Y in 0 degree:")
        self.ebox4 = wx.lib.masked.numctrl.NumCtrl(panel, style=wx.TE_CENTRE,
                                                   value=0.0, integerWidth=7,
                                                   fractionWidth=4,
                                                   allowNegative=True,
                                                   foregroundColour="Black",
                                                   signedForegroundColour="Black",
                                                   limited=False,
                                                   allowNone=True)
        # self.Bind(wx.lib.masked.EVT_NUM, self.OnZPz, self.ebox4)

        text5 = wx.StaticText(panel, label="Sample Z in 0 degree:")
        self.ebox5 = wx.lib.masked.numctrl.NumCtrl(panel, style=wx.TE_CENTRE,
                                                   value=0.0, integerWidth=7,
                                                   fractionWidth=4,
                                                   allowNegative=True,
                                                   foregroundColour="Black",
                                                   signedForegroundColour="Black",
                                                   limited=False,
                                                   allowNone=True)
        # self.Bind(wx.lib.masked.EVT_NUM, self.OnZPz, self.ebox4)

        text6 = wx.StaticText(panel, label="Field of View [um]:")
        self.ebox6 = wx.lib.masked.numctrl.NumCtrl(panel, style=wx.TE_CENTRE,
                                                   value=0.0, integerWidth=7,
                                                   fractionWidth=4,
                                                   allowNegative=True,
                                                   foregroundColour="Black",
                                                   signedForegroundColour="Black",
                                                   limited=False,
                                                   allowNone=True)
        # self.Bind(wx.lib.masked.EVT_NUM, self.OnZPz, self.ebox4)

        grid1.Add(text1, 0)
        grid1.Add(self.ebox1, 1)
        grid1.Add(text2, 0)
        grid1.Add(self.ebox2, 1)
        grid1.Add(text3, 0)
        grid1.Add(self.ebox3, 1)
        grid1.Add(text4, 0)
        grid1.Add(self.ebox4, 1)
        grid1.Add(text5, 0)
        grid1.Add(self.ebox5, 1)
        grid1.Add(text6, 0)
        grid1.Add(self.ebox6, 1)

        vbox.Add(grid1, 0, wx.ALL, 20)

        hbox = wx.BoxSizer(wx.HORIZONTAL)

        st_path = wx.StaticText(panel, label="Path: ")
        hbox.Add(st_path, 0)

        self.tc_path = wx.TextCtrl(panel, -1, size=((100, -1)), style=wx.TE_READONLY | wx.TE_RICH)
        hbox.Add(self.tc_path, -1, wx.LEFT, 5)
        vbox.Add(hbox, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 20)
        vbox.Add((0, 10))

        st = wx.StaticText(panel, label="Desired Angles [in StartDeg:StepDeg:EndDeg]:")
        self.tc_angles = wx.TextCtrl(panel, -1, size=((150, -1)), style=wx.TE_RICH | wx.VSCROLL,
                                     value='-90:1:89')
        vbox.Add(st, 0, wx.LEFT | wx.RIGHT, 20)
        vbox.Add((0, 3))
        vbox.Add(self.tc_angles, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, 20)

        grid2 = wx.FlexGridSizer(cols=2, vgap=5, hgap=25)

        button_infofromtxrm = wx.Button(panel, -1, 'Get Info From Quick *.TXRM', size=((160, -1)))
        self.Bind(wx.EVT_BUTTON, self.OnInfoFromQTxrm, id=button_infofromtxrm.GetId())
        button_infofromxrm = wx.Button(panel, -1, 'Get Info From Quick *.XRM', size=((160, -1)))
        self.Bind(wx.EVT_BUTTON, self.OnInfoFromQXrm, id=button_infofromxrm.GetId())
        if self.SampleXYZAboveRotation == 0:
            button_infofromxrm.Disable()
        button_finercal = wx.Button(panel, -1, 'Finer Stage Calibration', size=((160, -1)))
        self.Bind(wx.EVT_BUTTON, self.OnFineStageCal, id=button_finercal.GetId())
        if self.SampleXYZAboveRotation == 0:
            button_finercal.Disable()

        self.button_close = wx.Button(panel, -1, 'Submit')
        self.Bind(wx.EVT_BUTTON, self.OnSubmit, id=self.button_close.GetId())
        # vbox.Add(button_close, 0, wx.LEFT| wx.RIGHT | wx.ALIGN_RIGHT ,20)
        self.button_close.Disable()

        grid2.Add(button_infofromtxrm, 0)
        grid2.Add(wx.StaticText(panel, label=" "), 1)
        grid2.Add(button_infofromxrm, 0)
        grid2.Add(wx.StaticText(panel, label=" "), 1)
        grid2.Add(button_finercal, 0)
        grid2.Add(self.button_close, 1)

        vbox.Add(grid2, 1, wx.ALL, 20)

        panel.SetSizer(vbox)

        vboxtop.Add(panel, 1, wx.EXPAND)

        self.SetSizer(vboxtop)

    # ----------------------------------------------------------------------
    def OnInfoFromQTxrm(self, event):

        wildcard = "TXRM files (*.txrm)|*.txrm"
        dialog = wx.FileDialog(None, "Please select a txrm file",
                               wildcard=wildcard,
                               style=wx.OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            filepath = dialog.GetPath()
            filename = dialog.GetFilename()
            filedir = dialog.GetDirectory()
        else:
            return

        if self.microscope:
            txrminfo = self.microscope.get_txrm_info(filepath)

            length = txrminfo['ImagesTaken']
            samplets = np.array(txrminfo['MotPosT'])
            samplexs = np.array(txrminfo['MotPosX'])
            sampleys = np.array(txrminfo['MotPosY'])
            samplezs = np.array(txrminfo['MotPosZ'])
            FieldOfView = txrminfo['FOV']

            self.scriptdata.pixelsize = txrminfo['pixelsize']
            self.scriptdata.width = txrminfo['width']
            self.scriptdata.height = txrminfo['height']


        else:

            Ctxrm = txrm.txrm()

            txrminfo = Ctxrm.get_txrm_info(filepath)

            NumofMotors, xxx, yyy, zzz, ttt, zpxxx, zpyyy, zpzzz = Ctxrm.get_motor_index(txrminfo)

            length = txrminfo['ImagesTaken']
            samplets = np.array(txrminfo['MotPos'][ttt::NumofMotors])
            samplexs = np.array(txrminfo['MotPos'][xxx::NumofMotors])
            sampleys = np.array(txrminfo['MotPos'][yyy::NumofMotors])
            samplezs = np.array(txrminfo['MotPos'][zzz::NumofMotors])
            FieldOfView = txrminfo['FOV']

            self.scriptdata.pixelsize = txrminfo['pixelsize']
            # self.scriptdata.Binning = self.txrminfo['HBin']
            self.scriptdata.width = txrminfo['width']
            self.scriptdata.height = txrminfo['height']

        self.tc_path.SetValue(filedir)
        self.ebox6.SetValue(float(FieldOfView))
        zeroindex = np.argmin(np.abs(samplets))
        self.ebox3.SetValue(float(samplexs[zeroindex]))
        self.ebox4.SetValue(float(np.mean(sampleys)))
        self.ebox5.SetValue(float(samplezs[zeroindex]))

        if (self.SampleXYZAboveRotation == 1):
            if (abs(np.amax(samplexs) - np.amin(samplexs)) > 2.0) and \
                    (abs(np.amax(samplezs) - np.amin(samplezs)) > 2.0):
                fittedxc, fittedzc, R = CircleFit(samplexs, samplezs)
                self.ebox1.SetValue(float(samplexs[zeroindex] - fittedxc))
                self.ebox2.SetValue(float(-samplezs[zeroindex] + fittedzc))
                self.scriptdata.XCorr = float(samplexs[zeroindex] - fittedxc)
                self.scriptdata.ZCorr = float(-samplezs[zeroindex] + fittedzc)
            else:
                self.ebox1.SetValue(0.0)
                self.ebox2.SetValue(0.0)
                self.scriptdata.XCorr = 0.0
                self.scriptdata.ZCorr = 0.0
        else:
            self.ebox1.SetValue(0.0)
            self.ebox2.SetValue(0.0)
            self.scriptdata.XCorr = 0.0
            self.scriptdata.ZCorr = 0.0

        self.scriptdata.SampleXinZeroDegree = float(samplexs[zeroindex])
        self.scriptdata.SampleYinZeroDegree = float(np.mean(sampleys))
        self.scriptdata.SampleZinZeroDegree = float(samplezs[zeroindex])
        self.scriptdata.FOV = FieldOfView
        self.scriptdata.OutputPath = filedir

        self.button_close.Enable()

        self.Raise()

    # ----------------------------------------------------------------------
    def OnInfoFromQXrm(self, event):

        wildcard = "XRM files (*.xrm)|*.xrm"
        dialog = wx.FileDialog(None, "Please select a xrm file",
                               wildcard=wildcard,
                               style=wx.OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            filepath = dialog.GetPath()
            filename = dialog.GetFilename()
            filedir = dialog.GetDirectory()
        else:
            return

        if self.microscope:
            xrminfo = self.microscope.get_xrm_info(filepath)

            length = xrminfo['ImagesTaken']
            samplets = np.array([0])
            samplexs = np.array(xrminfo['MotPosX'])
            sampleys = np.array(xrminfo['MotPosY'])
            samplezs = np.array(xrminfo['MotPosZ'])
            FieldOfView = xrminfo['FOV']

            self.scriptdata.pixelsize = xrminfo['pixelsize']
            self.scriptdata.width = xrminfo['width']
            self.scriptdata.height = xrminfo['height']



        else:

            Ctxrm = txrm.txrm()

            xrminfo = Ctxrm.get_xrm_info(filepath)

            NumofMotors, xxx, yyy, zzz, ttt, zpxxx, zpyyy, zpzzz = Ctxrm.get_motor_index(xrminfo)

            length = xrminfo['ImagesTaken']
            samplets = np.array(xrminfo['MotPos'][ttt::NumofMotors])
            samplexs = np.array(xrminfo['MotPos'][xxx::NumofMotors])
            sampleys = np.array(xrminfo['MotPos'][yyy::NumofMotors])
            samplezs = np.array(xrminfo['MotPos'][zzz::NumofMotors])
            FieldOfView = xrminfo['FOV']

            self.scriptdata.pixelsize = xrminfo['pixelsize']
            # self.scriptdata.Binning = self.txrminfo['HBin']
            self.scriptdata.width = xrminfo['width']
            self.scriptdata.height = xrminfo['height']

        self.tc_path.SetValue(filedir)
        self.ebox6.SetValue(float(FieldOfView))
        zeroindex = 0
        self.ebox3.SetValue(float(samplexs))
        self.ebox4.SetValue(float(np.mean(sampleys)))
        self.ebox5.SetValue(float(samplezs))

        self.ebox1.SetValue(0.0)
        self.ebox2.SetValue(0.0)

        self.scriptdata.XCorr = 0.0
        self.scriptdata.ZCorr = 0.0
        self.scriptdata.SampleXinZeroDegree = float(samplexs)
        self.scriptdata.SampleYinZeroDegree = float(np.mean(sampleys))
        self.scriptdata.SampleZinZeroDegree = float(samplezs)
        self.scriptdata.FOV = FieldOfView
        self.scriptdata.OutputPath = filedir

        self.button_close.Enable()

        self.Raise()

    # ----------------------------------------------------------------------
    def OnFineStageCal(self, event):

        FineStageCalFrame(self, self.scriptdata, self.microscope).Show()

    # ----------------------------------------------------------------------
    def ShowFineStageCal(self, FineStageCalibrationInfo):

        self.scriptdata.XCorr = FineStageCalibrationInfo[0]
        self.scriptdata.ZCorr = FineStageCalibrationInfo[1]
        self.scriptdata.SampleXinZeroDegree = FineStageCalibrationInfo[2]
        self.scriptdata.SampleYinZeroDegree = FineStageCalibrationInfo[3]
        self.scriptdata.SampleZinZeroDegree = FineStageCalibrationInfo[4]
        self.scriptdata.FOV = FineStageCalibrationInfo[5]
        self.scriptdata.OutputPath = FineStageCalibrationInfo[6]

        self.ebox1.SetValue(float(self.scriptdata.XCorr))
        self.ebox2.SetValue(float(self.scriptdata.ZCorr))
        self.ebox3.SetValue(float(self.scriptdata.SampleXinZeroDegree))
        self.ebox4.SetValue(float(self.scriptdata.SampleYinZeroDegree))
        self.ebox5.SetValue(float(self.scriptdata.SampleZinZeroDegree))
        self.ebox6.SetValue(self.scriptdata.FOV)
        self.tc_path.SetValue(self.scriptdata.OutputPath)

        self.button_close.Enable()

    # ----------------------------------------------------------------------
    def OnSubmit(self, evt):

        self.scriptdata.XCorr = self.ebox1.GetValue()
        self.scriptdata.ZCorr = self.ebox2.GetValue()
        self.scriptdata.SampleXinZeroDegree = self.ebox3.GetValue()
        self.scriptdata.SampleYinZeroDegree = self.ebox4.GetValue()
        self.scriptdata.SampleZinZeroDegree = self.ebox5.GetValue()
        self.scriptdata.FOV = self.ebox6.GetValue()

        Cx = self.scriptdata.XCorr
        Cz = self.scriptdata.ZCorr
        X0 = self.scriptdata.SampleXinZeroDegree
        Y0 = self.scriptdata.SampleYinZeroDegree
        Z0 = self.scriptdata.SampleZinZeroDegree

        Cx = (-Cx + X0)
        Cz = (Cz + Z0)

        FOV = self.scriptdata.FOV

        # Get Angles from TC
        try:
            anglestext = self.tc_angles.GetValue()
            anglesets = anglestext.split(",")
            for iang in range(len(anglesets)):
                tangles = anglesets[iang]
                langles = tangles.split(":")
                theta1 = float(langles[0])
                theta2 = float(langles[2])
                dtheta = float(langles[1])
                if iang == 0:
                    Theta = np.arange(theta1, theta2 + dtheta, dtheta)
                else:
                    Theta = np.append(Theta, np.arange(theta1, theta2 + dtheta, dtheta))
        except:
            wx.MessageBox(
                "Could not load angles - please input desired angles as StartDeg:StepDeg:EndDeg separated by comma.")
            logging.error(
                "Could not load angles - please input desired angles as StartDeg:StepDeg:EndDeg separated by comma.",
                exc_info=True)
            return

        samplets = np.sort(Theta)
        samplexs = np.zeros((samplets.size))
        sampleys = np.ones((samplets.size)) * Y0
        samplezs = np.zeros((samplets.size))

        # print samplets

        r = ((Cx - X0) ** 2 + (Cz - Z0) ** 2) ** 0.5

        pi = np.math.pi

        try:
            # Disable numpy runtime warning, will check for errors after the calculation is done
            np.seterr(invalid='ignore')
            if X0 != Cx:
                for ii in range(Theta.size):
                    if (np.sin(samplets[ii] / 180. * pi) < 0):
                        samplexs[ii] = np.real((4 * r ** 2 * np.sin((pi * samplets[ii]) / 360) ** 2 +
                                                Cx ** 2 + Cz ** 2 - X0 ** 2 - Z0 ** 2 - r ** 2 -
                                                (2 * Cz * (Cx ** 2 * Cz + Cz * X0 ** 2 + Cx ** 2 * Z0 - Cz * Z0 ** 2 -
                                                           Cz ** 2 * Z0 - Cx * (
                                                                       - Cx ** 4 + 4 * Cx ** 3 * X0 - 2 * Cx ** 2 * Cz ** 2 +
                                                                       4 * Cx ** 2 * Cz * Z0 - 6 * Cx ** 2 * X0 ** 2 - 2 * Cx ** 2 * Z0 ** 2 +
                                                                       8 * Cx ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                       2 * Cx ** 2 * r ** 2 + 4 * Cx * Cz ** 2 * X0 - 8 * Cx * Cz * X0 * Z0 +
                                                                       4 * Cx * X0 ** 3 + 4 * Cx * X0 * Z0 ** 2 -
                                                                       16 * Cx * X0 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 -
                                                                       4 * Cx * X0 * r ** 2 - Cz ** 4 + 4 * Cz ** 3 * Z0 - 2 * Cz ** 2 * X0 ** 2 -
                                                                       6 * Cz ** 2 * Z0 ** 2 +
                                                                       8 * Cz ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                       2 * Cz ** 2 * r ** 2 + 4 * Cz * X0 ** 2 * Z0 + 4 * Cz * Z0 ** 3 -
                                                                       16 * Cz * Z0 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 -
                                                                       4 * Cz * Z0 * r ** 2 - X0 ** 4 - 2 * X0 ** 2 * Z0 ** 2 +
                                                                       8 * X0 ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                       2 * X0 ** 2 * r ** 2 - Z0 ** 4 +
                                                                       8 * Z0 ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                       2 * Z0 ** 2 * r ** 2 - 16 * r ** 4 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 4
                                                                       + 8 * r ** 4 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 - r ** 4) ** 0.5 +
                                                           X0 ** 2 * Z0 + X0 * (
                                                                       - Cx ** 4 + 4 * Cx ** 3 * X0 - 2 * Cx ** 2 * Cz ** 2 +
                                                                       4 * Cx ** 2 * Cz * Z0 - 6 * Cx ** 2 * X0 ** 2 - 2 * Cx ** 2 * Z0 ** 2 +
                                                                       8 * Cx ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                       2 * Cx ** 2 * r ** 2 + 4 * Cx * Cz ** 2 * X0 - 8 * Cx * Cz * X0 * Z0 +
                                                                       4 * Cx * X0 ** 3 + 4 * Cx * X0 * Z0 ** 2 -
                                                                       16 * Cx * X0 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 -
                                                                       4 * Cx * X0 * r ** 2 - Cz ** 4 + 4 * Cz ** 3 * Z0 - 2 * Cz ** 2 * X0 ** 2 -
                                                                       6 * Cz ** 2 * Z0 ** 2 +
                                                                       8 * Cz ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                       2 * Cz ** 2 * r ** 2 + 4 * Cz * X0 ** 2 * Z0 + 4 * Cz * Z0 ** 3 -
                                                                       16 * Cz * Z0 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 -
                                                                       4 * Cz * Z0 * r ** 2 - X0 ** 4 - 2 * X0 ** 2 * Z0 ** 2 +
                                                                       8 * X0 ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                       2 * X0 ** 2 * r ** 2 - Z0 ** 4 +
                                                                       8 * Z0 ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                       2 * Z0 ** 2 * r ** 2 - 16 * r ** 4 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 4
                                                                       + 8 * r ** 4 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 - r ** 4) ** 0.5 -
                                                           Cz * r ** 2 + Z0 * r ** 2 + Cz ** 3 + Z0 ** 3 - 2 * Cx * Cz * X0 -
                                                           2 * Cx * X0 * Z0 + 4 * Cz * r ** 2 * np.sin(
                                                            (pi * samplets[ii]) / 360) ** 2
                                                           - 4 * Z0 * r ** 2 * np.sin(
                                                            (pi * samplets[ii]) / 360) ** 2)) / (2 * Cx ** 2
                                                                                                 - 4 * Cx * X0 + 2 * Cz ** 2 - 4 * Cz * Z0 + 2 * X0 ** 2 + 2 * Z0 ** 2) +
                                                (2 * Z0 * (Cx ** 2 * Cz + Cz * X0 ** 2 + Cx ** 2 * Z0 - Cz * Z0 ** 2 -
                                                           Cz ** 2 * Z0 - Cx * (
                                                                       - Cx ** 4 + 4 * Cx ** 3 * X0 - 2 * Cx ** 2 * Cz ** 2 +
                                                                       4 * Cx ** 2 * Cz * Z0 - 6 * Cx ** 2 * X0 ** 2 - 2 * Cx ** 2 * Z0 ** 2 +
                                                                       8 * Cx ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                       2 * Cx ** 2 * r ** 2 + 4 * Cx * Cz ** 2 * X0 - 8 * Cx * Cz * X0 * Z0 +
                                                                       4 * Cx * X0 ** 3 + 4 * Cx * X0 * Z0 ** 2 -
                                                                       16 * Cx * X0 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 -
                                                                       4 * Cx * X0 * r ** 2 - Cz ** 4 + 4 * Cz ** 3 * Z0 - 2 * Cz ** 2 * X0 ** 2 -
                                                                       6 * Cz ** 2 * Z0 ** 2 +
                                                                       8 * Cz ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                       2 * Cz ** 2 * r ** 2 + 4 * Cz * X0 ** 2 * Z0 + 4 * Cz * Z0 ** 3 -
                                                                       16 * Cz * Z0 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 -
                                                                       4 * Cz * Z0 * r ** 2 - X0 ** 4 - 2 * X0 ** 2 * Z0 ** 2 +
                                                                       8 * X0 ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                       2 * X0 ** 2 * r ** 2 - Z0 ** 4 +
                                                                       8 * Z0 ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                       2 * Z0 ** 2 * r ** 2 - 16 * r ** 4 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 4
                                                                       + 8 * r ** 4 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 - r ** 4) ** 0.5 +
                                                           X0 ** 2 * Z0 + X0 * (
                                                                       - Cx ** 4 + 4 * Cx ** 3 * X0 - 2 * Cx ** 2 * Cz ** 2 +
                                                                       4 * Cx ** 2 * Cz * Z0 - 6 * Cx ** 2 * X0 ** 2 - 2 * Cx ** 2 * Z0 ** 2 +
                                                                       8 * Cx ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                       2 * Cx ** 2 * r ** 2 + 4 * Cx * Cz ** 2 * X0 - 8 * Cx * Cz * X0 * Z0 +
                                                                       4 * Cx * X0 ** 3 + 4 * Cx * X0 * Z0 ** 2 -
                                                                       16 * Cx * X0 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 -
                                                                       4 * Cx * X0 * r ** 2 - Cz ** 4 + 4 * Cz ** 3 * Z0 - 2 * Cz ** 2 * X0 ** 2 -
                                                                       6 * Cz ** 2 * Z0 ** 2 +
                                                                       8 * Cz ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                       2 * Cz ** 2 * r ** 2 + 4 * Cz * X0 ** 2 * Z0 + 4 * Cz * Z0 ** 3 -
                                                                       16 * Cz * Z0 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 -
                                                                       4 * Cz * Z0 * r ** 2 - X0 ** 4 - 2 * X0 ** 2 * Z0 ** 2 +
                                                                       8 * X0 ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                       2 * X0 ** 2 * r ** 2 - Z0 ** 4 +
                                                                       8 * Z0 ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                       2 * Z0 ** 2 * r ** 2 - 16 * r ** 4 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 4
                                                                       + 8 * r ** 4 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 - r ** 4) ** 0.5 -
                                                           Cz * r ** 2 + Z0 * r ** 2 + Cz ** 3 + Z0 ** 3 - 2 * Cx * Cz * X0 -
                                                           2 * Cx * X0 * Z0 + 4 * Cz * r ** 2 * np.sin(
                                                            (pi * samplets[ii]) / 360) ** 2
                                                           - 4 * Z0 * r ** 2 * np.sin(
                                                            (pi * samplets[ii]) / 360) ** 2)) / (2 * Cx ** 2
                                                                                                 - 4 * Cx * X0 + 2 * Cz ** 2 - 4 * Cz * Z0 + 2 * X0 ** 2 +
                                                                                                 2 * Z0 ** 2)) / (
                                                           2 * Cx - 2 * X0))

                        samplezs[ii] = np.real((Cx ** 2 * Cz + Cz * X0 ** 2 + Cx ** 2 * Z0 -
                                                Cz * Z0 ** 2 - Cz ** 2 * Z0 - Cx * (- Cx ** 4 + 4 * Cx ** 3 * X0 -
                                                                                    2 * Cx ** 2 * Cz ** 2 + 4 * Cx ** 2 * Cz * Z0 - 6 * Cx ** 2 * X0 ** 2 -
                                                                                    2 * Cx ** 2 * Z0 ** 2 +
                                                                                    8 * Cx ** 2 * r ** 2 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 2 +
                                                                                    2 * Cx ** 2 * r ** 2 + 4 * Cx * Cz ** 2 * X0 - 8 * Cx * Cz * X0 * Z0 +
                                                                                    4 * Cx * X0 ** 3 + 4 * Cx * X0 * Z0 ** 2 -
                                                                                    16 * Cx * X0 * r ** 2 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 2 -
                                                                                    4 * Cx * X0 * r ** 2 - Cz ** 4 + 4 * Cz ** 3 * Z0 - 2 * Cz ** 2 * X0 ** 2 -
                                                                                    6 * Cz ** 2 * Z0 ** 2 +
                                                                                    8 * Cz ** 2 * r ** 2 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 2 +
                                                                                    2 * Cz ** 2 * r ** 2 + 4 * Cz * X0 ** 2 * Z0 + 4 * Cz * Z0 ** 3 -
                                                                                    16 * Cz * Z0 * r ** 2 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 2 -
                                                                                    4 * Cz * Z0 * r ** 2 - X0 ** 4 - 2 * X0 ** 2 * Z0 ** 2 +
                                                                                    8 * X0 ** 2 * r ** 2 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 2 +
                                                                                    2 * X0 ** 2 * r ** 2 - Z0 ** 4 +
                                                                                    8 * Z0 ** 2 * r ** 2 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 2 +
                                                                                    2 * Z0 ** 2 * r ** 2 - 16 * r ** 4 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 4
                                                                                    + 8 * r ** 4 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 2 - r ** 4) ** 0.5 +
                                                X0 ** 2 * Z0 + X0 * (
                                                            - Cx ** 4 + 4 * Cx ** 3 * X0 - 2 * Cx ** 2 * Cz ** 2 +
                                                            4 * Cx ** 2 * Cz * Z0 - 6 * Cx ** 2 * X0 ** 2 - 2 * Cx ** 2 * Z0 ** 2 +
                                                            8 * Cx ** 2 * r ** 2 * np.sin(
                                                        (pi * samplets[ii]) / 360) ** 2 +
                                                            2 * Cx ** 2 * r ** 2 + 4 * Cx * Cz ** 2 * X0 - 8 * Cx * Cz * X0 * Z0 +
                                                            4 * Cx * X0 ** 3 + 4 * Cx * X0 * Z0 ** 2 -
                                                            16 * Cx * X0 * r ** 2 * np.sin(
                                                        (pi * samplets[ii]) / 360) ** 2 -
                                                            4 * Cx * X0 * r ** 2 - Cz ** 4 + 4 * Cz ** 3 * Z0 - 2 * Cz ** 2 * X0 ** 2 -
                                                            6 * Cz ** 2 * Z0 ** 2 +
                                                            8 * Cz ** 2 * r ** 2 * np.sin(
                                                        (pi * samplets[ii]) / 360) ** 2 +
                                                            2 * Cz ** 2 * r ** 2 + 4 * Cz * X0 ** 2 * Z0 + 4 * Cz * Z0 ** 3 -
                                                            16 * Cz * Z0 * r ** 2 * np.sin(
                                                        (pi * samplets[ii]) / 360) ** 2 -
                                                            4 * Cz * Z0 * r ** 2 - X0 ** 4 - 2 * X0 ** 2 * Z0 ** 2 +
                                                            8 * X0 ** 2 * r ** 2 * np.sin(
                                                        (pi * samplets[ii]) / 360) ** 2 +
                                                            2 * X0 ** 2 * r ** 2 - Z0 ** 4 +
                                                            8 * Z0 ** 2 * r ** 2 * np.sin(
                                                        (pi * samplets[ii]) / 360) ** 2 +
                                                            2 * Z0 ** 2 * r ** 2 - 16 * r ** 4 * np.sin(
                                                        (pi * samplets[ii]) / 360) ** 4
                                                            + 8 * r ** 4 * np.sin(
                                                        (pi * samplets[ii]) / 360) ** 2 - r ** 4) ** 0.5 -
                                                Cz * r ** 2 + Z0 * r ** 2 + Cz ** 3 + Z0 ** 3 - 2 * Cx * Cz * X0 -
                                                2 * Cx * X0 * Z0 + 4 * Cz * r ** 2 * np.sin(
                                    (pi * samplets[ii]) / 360) ** 2
                                                - 4 * Z0 * r ** 2 * np.sin((pi * samplets[ii]) / 360) ** 2) / (
                                                           2 * Cx ** 2 -
                                                           4 * Cx * X0 + 2 * Cz ** 2 - 4 * Cz * Z0 + 2 * X0 ** 2 + 2 * Z0 ** 2))

                    else:
                        samplexs[ii] = np.real((4 * r ** 2 * np.sin((pi * samplets[ii]) / 360) ** 2 + Cx ** 2
                                                + Cz ** 2 - X0 ** 2 - Z0 ** 2 - r ** 2 - (2 * Cz * (Cx ** 2 * Cz +
                                                                                                    Cz * X0 ** 2 + Cx ** 2 * Z0 - Cz * Z0 ** 2 - Cz ** 2 * Z0 + Cx * (
                                                                                                                - Cx ** 4
                                                                                                                + 4 * Cx ** 3 * X0 - 2 * Cx ** 2 * Cz ** 2 + 4 * Cx ** 2 * Cz * Z0 -
                                                                                                                6 * Cx ** 2 * X0 ** 2 - 2 * Cx ** 2 * Z0 ** 2 +
                                                                                                                8 * Cx ** 2 * r ** 2 * np.sin(
                                                                                                            (pi *
                                                                                                             samplets[
                                                                                                                 ii]) / 360) ** 2 +
                                                                                                                2 * Cx ** 2 * r ** 2 + 4 * Cx * Cz ** 2 * X0 - 8 * Cx * Cz * X0 * Z0 +
                                                                                                                4 * Cx * X0 ** 3 + 4 * Cx * X0 * Z0 ** 2 -
                                                                                                                16 * Cx * X0 * r ** 2 * np.sin(
                                                                                                            (pi *
                                                                                                             samplets[
                                                                                                                 ii]) / 360) ** 2 -
                                                                                                                4 * Cx * X0 * r ** 2 - Cz ** 4 + 4 * Cz ** 3 * Z0 - 2 * Cz ** 2 * X0 ** 2 -
                                                                                                                6 * Cz ** 2 * Z0 ** 2 +
                                                                                                                8 * Cz ** 2 * r ** 2 * np.sin(
                                                                                                            (pi *
                                                                                                             samplets[
                                                                                                                 ii]) / 360) ** 2 +
                                                                                                                2 * Cz ** 2 * r ** 2 + 4 * Cz * X0 ** 2 * Z0 + 4 * Cz * Z0 ** 3 -
                                                                                                                16 * Cz * Z0 * r ** 2 * np.sin(
                                                                                                            (pi *
                                                                                                             samplets[
                                                                                                                 ii]) / 360) ** 2 -
                                                                                                                4 * Cz * Z0 * r ** 2 - X0 ** 4 - 2 * X0 ** 2 * Z0 ** 2 +
                                                                                                                8 * X0 ** 2 * r ** 2 * np.sin(
                                                                                                            (pi *
                                                                                                             samplets[
                                                                                                                 ii]) / 360) ** 2 +
                                                                                                                2 * X0 ** 2 * r ** 2 - Z0 ** 4 +
                                                                                                                8 * Z0 ** 2 * r ** 2 * np.sin(
                                                                                                            (pi *
                                                                                                             samplets[
                                                                                                                 ii]) / 360) ** 2 +
                                                                                                                2 * Z0 ** 2 * r ** 2 - 16 * r ** 4 * np.sin(
                                                                                                            (pi *
                                                                                                             samplets[
                                                                                                                 ii]) / 360) ** 4
                                                                                                                + 8 * r ** 4 * np.sin(
                                                                                                            (pi *
                                                                                                             samplets[
                                                                                                                 ii]) / 360) ** 2 - r ** 4) ** 0.5 +
                                                                                                    X0 ** 2 * Z0 - X0 * (
                                                                                                            - Cx ** 4 + 4 * Cx ** 3 * X0 - 2 * Cx ** 2 * Cz ** 2 +
                                                                                                            4 * Cx ** 2 * Cz * Z0 - 6 * Cx ** 2 * X0 ** 2 - 2 * Cx ** 2 * Z0 ** 2 +
                                                                                                            8 * Cx ** 2 * r ** 2 * np.sin(
                                                                                                            (pi *
                                                                                                             samplets[
                                                                                                                 ii]) / 360) ** 2 +
                                                                                                            2 * Cx ** 2 * r ** 2 + 4 * Cx * Cz ** 2 * X0 - 8 * Cx * Cz * X0 * Z0 +
                                                                                                            4 * Cx * X0 ** 3 + 4 * Cx * X0 * Z0 ** 2 -
                                                                                                            16 * Cx * X0 * r ** 2 * np.sin(
                                                                                                            (pi *
                                                                                                             samplets[
                                                                                                                 ii]) / 360) ** 2 -
                                                                                                            4 * Cx * X0 * r ** 2 - Cz ** 4 + 4 * Cz ** 3 * Z0 - 2 * Cz ** 2 * X0 ** 2 -
                                                                                                            6 * Cz ** 2 * Z0 ** 2 +
                                                                                                            8 * Cz ** 2 * r ** 2 * np.sin(
                                                                                                            (pi *
                                                                                                             samplets[
                                                                                                                 ii]) / 360) ** 2 +
                                                                                                            2 * Cz ** 2 * r ** 2 + 4 * Cz * X0 ** 2 * Z0 + 4 * Cz * Z0 ** 3 -
                                                                                                            16 * Cz * Z0 * r ** 2 * np.sin(
                                                                                                            (pi *
                                                                                                             samplets[
                                                                                                                 ii]) / 360) ** 2 -
                                                                                                            4 * Cz * Z0 * r ** 2 - X0 ** 4 - 2 * X0 ** 2 * Z0 ** 2 +
                                                                                                            8 * X0 ** 2 * r ** 2 * np.sin(
                                                                                                            (pi *
                                                                                                             samplets[
                                                                                                                 ii]) / 360) ** 2 +
                                                                                                            2 * X0 ** 2 * r ** 2 - Z0 ** 4 +
                                                                                                            8 * Z0 ** 2 * r ** 2 * np.sin(
                                                                                                            (pi *
                                                                                                             samplets[
                                                                                                                 ii]) / 360) ** 2 +
                                                                                                            2 * Z0 ** 2 * r ** 2 - 16 * r ** 4 * np.sin(
                                                                                                        (pi *
                                                                                                             samplets[
                                                                                                                 ii]) / 360) ** 4
                                                                                                            + 8 * r ** 4 * np.sin(
                                                                                                            (pi *
                                                                                                             samplets[
                                                                                                                 ii]) / 360) ** 2 - r ** 4) ** 0.5 -
                                                                                                    Cz * r ** 2 + Z0 * r ** 2 + Cz ** 3 + Z0 ** 3 - 2 * Cx * Cz * X0 -
                                                                                                    2 * Cx * X0 * Z0 + 4 * Cz * r ** 2 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 2
                                                                                                    - 4 * Z0 * r ** 2 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 2)) / (2 * Cx ** 2
                                                                                 - 4 * Cx * X0 + 2 * Cz ** 2 - 4 * Cz * Z0 + 2 * X0 ** 2 + 2 * Z0 ** 2) +
                                                (2 * Z0 * (Cx ** 2 * Cz + Cz * X0 ** 2 + Cx ** 2 * Z0 - Cz * Z0 ** 2 -
                                                           Cz ** 2 * Z0 + Cx * (
                                                                   - Cx ** 4 + 4 * Cx ** 3 * X0 - 2 * Cx ** 2 * Cz ** 2 +
                                                                   4 * Cx ** 2 * Cz * Z0 - 6 * Cx ** 2 * X0 ** 2 - 2 * Cx ** 2 * Z0 ** 2 +
                                                                   8 * Cx ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                   2 * Cx ** 2 * r ** 2 + 4 * Cx * Cz ** 2 * X0 - 8 * Cx * Cz * X0 * Z0 +
                                                                   4 * Cx * X0 ** 3 + 4 * Cx * X0 * Z0 ** 2 -
                                                                   16 * Cx * X0 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 -
                                                                   4 * Cx * X0 * r ** 2 - Cz ** 4 + 4 * Cz ** 3 * Z0 - 2 * Cz ** 2 * X0 ** 2 -
                                                                   6 * Cz ** 2 * Z0 ** 2 +
                                                                   8 * Cz ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                   2 * Cz ** 2 * r ** 2 + 4 * Cz * X0 ** 2 * Z0 + 4 * Cz * Z0 ** 3 -
                                                                   16 * Cz * Z0 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 -
                                                                   4 * Cz * Z0 * r ** 2 - X0 ** 4 - 2 * X0 ** 2 * Z0 ** 2 +
                                                                   8 * X0 ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                   2 * X0 ** 2 * r ** 2 - Z0 ** 4 +
                                                                   8 * Z0 ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                   2 * Z0 ** 2 * r ** 2 - 16 * r ** 4 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 4
                                                                   + 8 * r ** 4 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 - r ** 4) ** 0.5 +
                                                           X0 ** 2 * Z0 - X0 * (
                                                                   - Cx ** 4 + 4 * Cx ** 3 * X0 - 2 * Cx ** 2 * Cz ** 2 +
                                                                   4 * Cx ** 2 * Cz * Z0 - 6 * Cx ** 2 * X0 ** 2 - 2 * Cx ** 2 * Z0 ** 2 +
                                                                   8 * Cx ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                   2 * Cx ** 2 * r ** 2 + 4 * Cx * Cz ** 2 * X0 - 8 * Cx * Cz * X0 * Z0 +
                                                                   4 * Cx * X0 ** 3 + 4 * Cx * X0 * Z0 ** 2 -
                                                                   16 * Cx * X0 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 -
                                                                   4 * Cx * X0 * r ** 2 - Cz ** 4 + 4 * Cz ** 3 * Z0 - 2 * Cz ** 2 * X0 ** 2 -
                                                                   6 * Cz ** 2 * Z0 ** 2 +
                                                                   8 * Cz ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                   2 * Cz ** 2 * r ** 2 + 4 * Cz * X0 ** 2 * Z0 + 4 * Cz * Z0 ** 3 -
                                                                   16 * Cz * Z0 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 -
                                                                   4 * Cz * Z0 * r ** 2 - X0 ** 4 - 2 * X0 ** 2 * Z0 ** 2 +
                                                                   8 * X0 ** 2 * r ** 2 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 +
                                                                   2 * X0 ** 2 * r ** 2 - Z0 ** 4 +
                                                                   8 * Z0 ** 2 * r ** 2 * np.sin(
                                                               (pi * samplets[ii]) / 360) ** 2 +
                                                                   2 * Z0 ** 2 * r ** 2 - 16 * r ** 4 * np.sin(
                                                               (pi * samplets[ii]) / 360) ** 4
                                                                   + 8 * r ** 4 * np.sin(
                                                                   (pi * samplets[ii]) / 360) ** 2 - r ** 4) ** 0.5 -
                                                           Cz * r ** 2 + Z0 * r ** 2 + Cz ** 3 + Z0 ** 3 - 2 * Cx * Cz * X0 -
                                                           2 * Cx * X0 * Z0 + 4 * Cz * r ** 2 * np.sin(
                                                            (pi * samplets[ii]) / 360) ** 2
                                                           - 4 * Z0 * r ** 2 * np.sin(
                                                            (pi * samplets[ii]) / 360) ** 2)) / (2 * Cx ** 2
                                                                                                 - 4 * Cx * X0 + 2 * Cz ** 2 - 4 * Cz * Z0 + 2 * X0 ** 2 +
                                                                                                 2 * Z0 ** 2)) / (
                                                           2 * Cx - 2 * X0))

                        samplezs[ii] = np.real((Cx ** 2 * Cz + Cz * X0 ** 2 + Cx ** 2 * Z0 -
                                                Cz * Z0 ** 2 - Cz ** 2 * Z0 + Cx * (- Cx ** 4 + 4 * Cx ** 3 * X0 -
                                                                                    2 * Cx ** 2 * Cz ** 2 + 4 * Cx ** 2 * Cz * Z0 - 6 * Cx ** 2 * X0 ** 2 -
                                                                                    2 * Cx ** 2 * Z0 ** 2 +
                                                                                    8 * Cx ** 2 * r ** 2 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 2 +
                                                                                    2 * Cx ** 2 * r ** 2 + 4 * Cx * Cz ** 2 * X0 - 8 * Cx * Cz * X0 * Z0 +
                                                                                    4 * Cx * X0 ** 3 + 4 * Cx * X0 * Z0 ** 2 -
                                                                                    16 * Cx * X0 * r ** 2 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 2 -
                                                                                    4 * Cx * X0 * r ** 2 - Cz ** 4 + 4 * Cz ** 3 * Z0 - 2 * Cz ** 2 * X0 ** 2 -
                                                                                    6 * Cz ** 2 * Z0 ** 2 +
                                                                                    8 * Cz ** 2 * r ** 2 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 2 +
                                                                                    2 * Cz ** 2 * r ** 2 + 4 * Cz * X0 ** 2 * Z0 + 4 * Cz * Z0 ** 3 -
                                                                                    16 * Cz * Z0 * r ** 2 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 2 -
                                                                                    4 * Cz * Z0 * r ** 2 - X0 ** 4 - 2 * X0 ** 2 * Z0 ** 2 +
                                                                                    8 * X0 ** 2 * r ** 2 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 2 +
                                                                                    2 * X0 ** 2 * r ** 2 - Z0 ** 4 +
                                                                                    8 * Z0 ** 2 * r ** 2 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 2 +
                                                                                    2 * Z0 ** 2 * r ** 2 - 16 * r ** 4 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 4
                                                                                    + 8 * r ** 4 * np.sin(
                                            (pi * samplets[ii]) / 360) ** 2 - r ** 4) ** 0.5 +
                                                X0 ** 2 * Z0 - X0 * (
                                                        - Cx ** 4 + 4 * Cx ** 3 * X0 - 2 * Cx ** 2 * Cz ** 2 +
                                                        4 * Cx ** 2 * Cz * Z0 - 6 * Cx ** 2 * X0 ** 2 - 2 * Cx ** 2 * Z0 ** 2 +
                                                        8 * Cx ** 2 * r ** 2 * np.sin(
                                                        (pi * samplets[ii]) / 360) ** 2 +
                                                        2 * Cx ** 2 * r ** 2 + 4 * Cx * Cz ** 2 * X0 - 8 * Cx * Cz * X0 * Z0 +
                                                        4 * Cx * X0 ** 3 + 4 * Cx * X0 * Z0 ** 2 -
                                                        16 * Cx * X0 * r ** 2 * np.sin(
                                                        (pi * samplets[ii]) / 360) ** 2 -
                                                        4 * Cx * X0 * r ** 2 - Cz ** 4 + 4 * Cz ** 3 * Z0 - 2 * Cz ** 2 * X0 ** 2 -
                                                        6 * Cz ** 2 * Z0 ** 2 +
                                                        8 * Cz ** 2 * r ** 2 * np.sin(
                                                        (pi * samplets[ii]) / 360) ** 2 +
                                                        2 * Cz ** 2 * r ** 2 + 4 * Cz * X0 ** 2 * Z0 + 4 * Cz * Z0 ** 3 -
                                                        16 * Cz * Z0 * r ** 2 * np.sin(
                                                    (pi * samplets[ii]) / 360) ** 2 -
                                                        4 * Cz * Z0 * r ** 2 - X0 ** 4 - 2 * X0 ** 2 * Z0 ** 2 +
                                                        8 * X0 ** 2 * r ** 2 * np.sin(
                                                        (pi * samplets[ii]) / 360) ** 2 +
                                                        2 * X0 ** 2 * r ** 2 - Z0 ** 4 +
                                                        8 * Z0 ** 2 * r ** 2 * np.sin(
                                                        (pi * samplets[ii]) / 360) ** 2 +
                                                        2 * Z0 ** 2 * r ** 2 - 16 * r ** 4 * np.sin(
                                                    (pi * samplets[ii]) / 360) ** 4
                                                        + 8 * r ** 4 * np.sin(
                                                        (pi * samplets[ii]) / 360) ** 2 - r ** 4) ** 0.5 -
                                                Cz * r ** 2 + Z0 * r ** 2 + Cz ** 3 + Z0 ** 3 - 2 * Cx * Cz * X0 -
                                                2 * Cx * X0 * Z0 + 4 * Cz * r ** 2 * np.sin(
                                    (pi * samplets[ii]) / 360) ** 2
                                                - 4 * Z0 * r ** 2 * np.sin((pi * samplets[ii]) / 360) ** 2) / (
                                                       2 * Cx ** 2 -
                                                       4 * Cx * X0 + 2 * Cz ** 2 - 4 * Cz * Z0 + 2 * X0 ** 2 + 2 * Z0 ** 2))

            else:
                if Z0 != Cz:
                    for ii in range(Theta.size):
                        samplezs[ii] = ((r ** 2 - (2 * r * np.sin(samplets[ii] / 2 / 180 * pi)) ** 2) / (
                                    Z0 - Cz) + Cz + Z0) / 2
                        samplexs[ii] = np.sign(np.sin(samplets[ii] / 180 / pi)) * (
                                    r ** 2 - (samplezs[ii] - Cz) ** 2) ** 0.5 + Cx
                else:
                    samplexs = np.ones((samplexs.size)) * Cx
                    samplezs = np.ones((samplezs.size)) * Cz
        except:
            wx.MessageBox("Tomo position array calculation is incorrect. Please check sample stage calibration values.")
            logging.error("Tomo position array calculation is incorrect. Please check sample stage calibration values.",
                          exc_info=True)
            return

        # Check if the fit was OK and if position arrays have Nans!
        if (np.isnan(np.sum(samplexs)) or np.isnan(np.sum(samplezs))):
            Nnans = 0

            # Interpolate Nan's if found
            nans, x = self.nan_helper(samplexs)
            Nnans += np.sum(nans)
            if nans.any():
                samplexs[nans] = np.interp(x(nans), x(~nans), samplexs[~nans])
                sampleys[nans] = np.interp(x(nans), x(~nans), sampleys[~nans])
                samplezs[nans] = np.interp(x(nans), x(~nans), samplezs[~nans])
                samplets[nans] = np.interp(x(nans), x(~nans), samplets[~nans])
                print
                'Interpolated values', samplets[nans]
            #                 samplexs = samplexs[~nans]
            #                 sampleys = sampleys[~nans]
            #                 samplezs = samplezs[~nans]
            #                 samplets = samplets[~nans]

            nans, x = self.nan_helper(samplezs)
            Nnans += np.sum(nans)
            if nans.any():
                samplexs[nans] = np.interp(x(nans), x(~nans), samplexs[~nans])
                sampleys[nans] = np.interp(x(nans), x(~nans), sampleys[~nans])
                samplezs[nans] = np.interp(x(nans), x(~nans), samplezs[~nans])
                samplets[nans] = np.interp(x(nans), x(~nans), samplets[~nans])
                print
                'Interpolated values', samplets[nans]
            #                 samplexs = samplexs[~nans]
            #                 sampleys = sampleys[~nans]
            #                 samplezs = samplezs[~nans]
            #                 samplets = samplets[~nans]

            message = "Warning - {0:2d} NaNs found in Tomo position arrays - interpolating the values.".format(Nnans)
            print
            message
            logging.warning(message)

        self.scriptdata.samplexs = samplexs
        self.scriptdata.sampleys = sampleys
        self.scriptdata.samplezs = samplezs
        self.scriptdata.samplets = samplets

        self.Close(True)

        self.mainframe.ShowEnabledTC()
        self.mainframe.txrmloaded = True

    # ----------------------------------------------------------------------
    def nan_helper(self, y):
        """Helper to handle indices and logical indices of NaNs.

        Input:
            - y, 1d numpy array with possible NaNs
        Output:
            - nans, logical indices of NaNs
            - index, a function, with signature indices= index(logical_indices),
              to convert logical indices of NaNs to 'equivalent' indices
        Example:
            >>> # linear interpolation of NaNs
            >>> nans, x= nan_helper(y)
            >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
        """

        return np.isnan(y), lambda z: z.nonzero()[0]


""" ------------------------------------------------------------------------------------------------"""


class FineStageCalFrame(wx.Frame):

    def __init__(self, parent, scriptdata, microscope):
        wx.Frame.__init__(self, parent, title="Fine Stage Calibration", size=(875, 600))

        self.scriptdata = scriptdata
        self.microscope = microscope
        self.parent = parent

        self.SetBackgroundColour("White")

        self.i_img = 0

        self.imgsize = 500

        vboxtop = wx.BoxSizer(wx.VERTICAL)

        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
        hbox = wx.BoxSizer(wx.HORIZONTAL)

        # Data Image

        self.dataimage = np.zeros((self.imgsize * self.imgsize * 3), dtype=np.int8)

        image = wx.ImageFromData(self.imgsize, self.imgsize, self.dataimage)
        bitdata = wx.BitmapFromImage(image)

        from wx.lib.statbmp import GenStaticBitmap
        self.imageCtrl = GenStaticBitmap(panel, wx.ID_ANY, bitdata)
        wx.EVT_LEFT_DOWN(self.imageCtrl, self.OnImageMouseDown)
        self.imageCtrl.Disable()

        # Path
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        self.tc_path = wx.TextCtrl(panel, -1, size=((self.imgsize, -1)), style=wx.TE_READONLY | wx.TE_RICH)
        hbox1.Add(self.tc_path, 0, wx.RIGHT, 20)

        button_browse = wx.Button(panel, -1, 'Browse 4 TXRM', size=((150, -1)))
        self.Bind(wx.EVT_BUTTON, self.OnGetInfoFromTxrm, id=button_browse.GetId())
        hbox1.Add(button_browse, 0, wx.RIGHT, 5)

        self.button_submit = wx.Button(panel, -1, 'Submit Stage Calibration', size=((150, -1)))
        self.button_submit.Disable()
        self.Bind(wx.EVT_BUTTON, self.OnSubmit, id=self.button_submit.GetId())
        hbox1.Add(self.button_submit, 0)

        # Sample XYZ info
        grid1 = wx.FlexGridSizer(cols=3, vgap=3, hgap=2)

        textx = wx.StaticText(panel, label="Sample X")
        texty = wx.StaticText(panel, label="Sample Y")
        textz = wx.StaticText(panel, label="Sample Z")

        self.tc_x = wx.ListCtrl(panel, -1, size=(100, 150), style=wx.LC_REPORT | wx.LC_NO_HEADER | wx.SIMPLE_BORDER)
        self.tc_x.InsertColumn(0, 'x')
        self.tc_x.SetColumnWidth(0, 95)

        self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnXListClick, self.tc_x)

        self.tc_y = wx.ListCtrl(panel, -1, size=(100, 150), style=wx.LC_REPORT | wx.LC_NO_HEADER | wx.SIMPLE_BORDER)
        self.tc_y.InsertColumn(0, 'y')
        self.tc_y.SetColumnWidth(0, 95)
        self.tc_y.Disable()

        self.tc_z = wx.ListCtrl(panel, -1, size=(100, 150), style=wx.LC_REPORT | wx.LC_NO_HEADER | wx.SIMPLE_BORDER)
        self.tc_z.InsertColumn(0, 'z')
        self.tc_z.SetColumnWidth(0, 95)
        self.tc_z.Disable()

        grid1.Add(textx, 0)
        grid1.Add(texty, 1)
        grid1.Add(textz, 2)
        grid1.Add(self.tc_x, 0)
        grid1.Add(self.tc_y, 1)
        grid1.Add(self.tc_z, 2)

        # Info
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        vbox1.Add((0, 10))

        grid2 = wx.FlexGridSizer(cols=4, vgap=5, hgap=10)

        textps = wx.StaticText(panel, label="PixelSize:")
        self.eboxps = wx.lib.masked.numctrl.NumCtrl(panel, style=wx.TE_CENTRE | wx.TE_READONLY,
                                                    value=0.0, integerWidth=5,
                                                    fractionWidth=5, limited=False,
                                                    allowNone=True)

        textpsu = wx.StaticText(panel, label="um")

        textis = wx.StaticText(panel, label="ImgSize")
        self.eboxisx = wx.lib.intctrl.IntCtrl(panel, size=(100, -1),
                                              style=wx.TE_CENTRE | wx.TE_READONLY,
                                              value=0, limited=True)
        self.eboxisy = wx.lib.intctrl.IntCtrl(panel, size=(100, -1),
                                              style=wx.TE_CENTRE | wx.TE_READONLY,
                                              value=0, limited=True)

        textfc = wx.StaticText(panel, label="Fitted Center:")
        self.eboxfc1 = wx.lib.masked.numctrl.NumCtrl(panel, style=wx.TE_CENTRE | wx.TE_READONLY,
                                                     value=0.0, integerWidth=5,
                                                     fractionWidth=5,
                                                     allowNegative=True,
                                                     foregroundColour="Black",
                                                     signedForegroundColour="Black",
                                                     limited=False,
                                                     allowNone=True)

        self.eboxfc2 = wx.lib.masked.numctrl.NumCtrl(panel, style=wx.TE_CENTRE | wx.TE_READONLY,
                                                     value=0.0, integerWidth=5,
                                                     fractionWidth=5,
                                                     allowNegative=True,
                                                     foregroundColour="Black",
                                                     signedForegroundColour="Black",
                                                     limited=False,
                                                     allowNone=True)

        grid2.Add(textps, 0)
        grid2.Add(self.eboxps, 1)
        grid2.Add(wx.StaticText(panel, label=" "), 2)
        grid2.Add(textpsu, 3)
        grid2.Add(textis, 0)
        grid2.Add(self.eboxisx, 1)
        grid2.Add(wx.StaticText(panel, label="X"), 2)
        grid2.Add(self.eboxisy, 3)
        grid2.Add(textfc, 0)
        grid2.Add(self.eboxfc1, 1)
        grid2.Add(wx.StaticText(panel, label=" "), 2)
        grid2.Add(self.eboxfc2, 3)

        # Data Plot
        self.fig = Figure(figsize=(3.3, 3.3), facecolor='w')
        self.canvas = FigCanvas(panel, -1, self.fig)
        self.fig.add_axes()
        self.axes = self.fig.gca()
        mplot.rcParams['font.size'] = 8

        hbox.Add(self.imageCtrl, 1, wx.LEFT | wx.RIGHT, 20)

        vboxr = wx.BoxSizer(wx.VERTICAL)
        vboxr.Add((0, 10))
        vboxr.Add(grid1, 0)
        vboxr.Add((0, 10))
        vboxr.Add(grid2, 0, wx.EXPAND)
        # vboxr.Add((0,0))
        vboxr.Add(self.canvas, 0, wx.LEFT, 10)

        hbox.Add(vboxr, 0, wx.EXPAND)

        vbox.Add((0, 20))
        vbox.Add(hbox1, 0, wx.EXPAND | wx.LEFT, 20)
        vbox.Add(hbox)

        panel.SetSizer(vbox)

        vboxtop.Add(panel, 1, wx.EXPAND)

        self.SetSizer(vboxtop)

    # ----------------------------------------------------------------------
    def OnGetInfoFromTxrm(self, event):

        wildcard = "TXRM files (*.txrm)|*.txrm"
        dialog = wx.FileDialog(self, "Choose txrm file",
                               wildcard=wildcard,
                               style=wx.OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            self.filepath = dialog.GetPath()
            filename = dialog.GetFilename()
            filedir = dialog.GetDirectory()
        else:
            return

        if self.microscope:
            txrminfo = self.microscope.get_txrm_info(self.filepath)

            length = txrminfo['ImagesTaken']
            samplets = np.array(txrminfo['MotPosT'])
            samplexs = np.array(txrminfo['MotPosX'])
            sampleys = np.array(txrminfo['MotPosY'])
            samplezs = np.array(txrminfo['MotPosZ'])
            FieldOfView = txrminfo['FOV']

            self.scriptdata.pixelsize = txrminfo['pixelsize']
            self.scriptdata.width = txrminfo['width']
            self.scriptdata.height = txrminfo['height']

        else:

            self.Ctxrm = txrm.txrm()

            txrminfo = self.Ctxrm.get_txrm_info(self.filepath)

            NumofMotors, xxx, yyy, zzz, ttt, zpxxx, zpyyy, zpzzz = self.Ctxrm.get_motor_index(txrminfo)

            length = txrminfo['ImagesTaken']
            samplets = np.array(txrminfo['MotPos'][ttt::NumofMotors])
            samplexs = np.array(txrminfo['MotPos'][xxx::NumofMotors])
            sampleys = np.array(txrminfo['MotPos'][yyy::NumofMotors])
            samplezs = np.array(txrminfo['MotPos'][zzz::NumofMotors])
            FieldOfView = txrminfo['FOV']

            self.scriptdata.pixelsize = txrminfo['pixelsize']
            self.scriptdata.width = txrminfo['width']
            self.scriptdata.height = txrminfo['height']

        self.tc_path.SetValue(filedir)
        self.scriptdata.OutputPath = filedir

        self.fittedxc, self.fittedzc, R = CircleFit(samplexs, samplezs)

        self.eboxps.SetValue(float(self.scriptdata.pixelsize))
        self.eboxisx.SetValue(int(self.scriptdata.width))
        self.eboxisy.SetValue(int(self.scriptdata.height))
        self.eboxfc1.SetValue(float(self.fittedxc))
        self.eboxfc2.SetValue(float(self.fittedzc))

        # FineStageCalibrationInfo
        self.FineStageCalibrationInfo = []
        zeroindx = np.argmin(np.abs(samplets))
        self.FineStageCalibrationInfo.append(samplexs[zeroindx] - self.fittedxc)
        self.FineStageCalibrationInfo.append(-samplezs[zeroindx] + self.fittedzc)
        self.FineStageCalibrationInfo.append(samplexs[zeroindx])
        self.FineStageCalibrationInfo.append(np.mean(sampleys))
        self.FineStageCalibrationInfo.append(samplezs[zeroindx])
        self.FineStageCalibrationInfo.append(FieldOfView)
        self.FineStageCalibrationInfo.append(filedir)

        self.samplex = samplexs
        self.sampley = sampleys
        self.samplez = samplezs

        self.newsamplex = samplexs.copy()
        self.newsampley = sampleys.copy()
        self.newsamplez = samplezs.copy()

        self.i_img = 0

        self.ShowSamplexyz()
        self.DrawPlot()
        self.DisplayImage(self.i_img + 1)

        self.imageCtrl.Enable()
        self.button_submit.Enable()

    # ----------------------------------------------------------------------
    def OnXListClick(self, event):

        self.i_img = event.m_itemIndex

        event.Skip()

        self.ShowSamplexyz()
        self.DrawPlot()
        self.DisplayImage(self.i_img + 1)

    # ----------------------------------------------------------------------
    def DisplayImage(self, i):

        if self.microscope:
            img = self.microscope.getTXRMimage(self.filepath, i - 1)
        else:
            img = self.Ctxrm.getTXRMimage(self.filepath, i)

        limits = [np.min(img), np.max(img)]
        if (limits[0] != limits[1]):
            delta = 1.0 / (limits[1] - limits[0]);
            img = (img.copy() - limits[0]) * delta * 255

        rgbimg = np.zeros((img.shape[0], img.shape[1], 3), dtype=np.uint8)

        rgbimg[:, :, 0] = img
        rgbimg[:, :, 1] = img
        rgbimg[:, :, 2] = img

        # Draw cursor:
        pixelsize = self.scriptdata.pixelsize

        CofImgx = self.samplex[self.i_img]
        CofImgy = self.sampley[self.i_img]

        markerx = self.newsamplex[self.i_img]
        markery = self.newsampley[self.i_img]

        markerxindx = ((markerx - CofImgx) / pixelsize + img.shape[0] / 2)
        markeryindx = ((markery - CofImgy) / pixelsize + img.shape[1] / 2)
        linelen = min(img.shape[0], img.shape[1]) / 20;

        image = wx.EmptyImage(img.shape[0], img.shape[1])

        for x in range(max(int(markerxindx - linelen), 0), min(int(markerxindx + linelen), img.shape[0] - 1)):
            rgbimg[markeryindx, x, 0] = 255
            rgbimg[markeryindx, x, 1] = 0
            rgbimg[markeryindx, x, 2] = 0

        for y in range(max(int(markeryindx - linelen), 0), min(int(markeryindx + linelen), img.shape[1] - 1)):
            rgbimg[y, markerxindx, 0] = 255
            rgbimg[y, markerxindx, 1] = 0
            rgbimg[y, markerxindx, 2] = 0

        image.SetData(rgbimg)
        # Scale the image to fit bitmap panel:
        neww = self.imgsize
        newh = self.imgsize
        if (img.shape[0] > img.shape[1]):
            newh = img.shape[1] / img.shape[0] * self.imgsize
        elif (img.shape[1] > img.shape[0]):
            neww = img.shape[0] / img.shape[1] * self.imgsize
        imagescaled = image.Scale(neww, newh, wx.IMAGE_QUALITY_HIGH)

        bmp = imagescaled.ConvertToBitmap()
        self.imageCtrl.SetBitmap(bmp)

    # ----------------------------------------------------------------------
    def DrawPlot(self):

        self.fittedxc, self.fittedzc, R = CircleFit(self.newsamplex, self.newsamplez)
        self.axes.clear()
        samplearr = np.array([self.newsamplex, self.newsamplez])

        addbound = (np.max(samplearr) - np.min(samplearr)) / 10

        bounds = [np.min(samplearr) - addbound, np.max(samplearr) + addbound]

        plot = self.axes.plot(self.newsamplex, self.newsamplez, 'c*', self.fittedxc, self.fittedzc, 'b*', markersize=8)

        self.axes.set_xlim(bounds)
        self.axes.set_ylim(bounds)
        self.canvas.draw()

        self.eboxfc1.SetValue(float(self.fittedxc))
        self.eboxfc2.SetValue(float(self.fittedzc))

    # ----------------------------------------------------------------------
    def ShowSamplexyz(self):

        self.tc_x.DeleteAllItems()
        self.tc_y.DeleteAllItems()
        self.tc_z.DeleteAllItems()

        for i in range(len(self.newsamplex)):
            self.tc_x.InsertStringItem(i, '{0:08.2f}'.format(self.newsamplex[i]))
            self.tc_y.InsertStringItem(i, '{0:08.2f}'.format(self.newsampley[i]))
            self.tc_z.InsertStringItem(i, '{0:08.2f}'.format(self.newsamplez[i]))

        self.tc_x.SetItemBackgroundColour(self.i_img, 'light blue')
        self.tc_y.SetItemBackgroundColour(self.i_img, 'light blue')
        self.tc_z.SetItemBackgroundColour(self.i_img, 'light blue')

    # ----------------------------------------------------------------------
    def OnImageMouseDown(self, event):

        # get the position tuple
        pt = event.GetPosition()

        newmarkerX = pt[0]
        newmarkerY = pt[1]

        dispsz1 = self.imgsize
        dispsz2 = self.imgsize
        imgsz1 = self.scriptdata.width
        imgsz2 = self.scriptdata.height
        newmarkerX = float(newmarkerX) / dispsz1 * imgsz1
        newmarkerY = float(newmarkerY) / dispsz2 * imgsz2
        tilecenterX = self.samplex[self.i_img]
        tilecenterY = self.sampley[self.i_img]

        pixelsize = self.scriptdata.pixelsize

        newmarkerYcoordinate = (newmarkerY - imgsz1 / 2) * pixelsize + tilecenterY
        newmarkerXcoordinate = (newmarkerX - imgsz2 / 2) * pixelsize + tilecenterX

        self.newsamplex[self.i_img] = newmarkerXcoordinate
        self.newsampley[self.i_img] = newmarkerYcoordinate

        if self.i_img < len(self.newsamplex) - 1:
            self.i_img = self.i_img + 1

        self.ShowSamplexyz()
        self.DrawPlot()
        self.DisplayImage(self.i_img + 1)

    # ----------------------------------------------------------------------
    def OnSubmit(self, evt):

        self.FineStageCalibrationInfo[1] = -self.FineStageCalibrationInfo[4] + self.fittedzc
        self.FineStageCalibrationInfo[0] = self.FineStageCalibrationInfo[2] - self.fittedxc

        self.parent.ShowFineStageCal(self.FineStageCalibrationInfo)
        self.Close(True)

    # ----------------------------------------------------------------------
    def OnClose(self, evt):

        self.Close(True)


""" ------------------------------------------------------------------------------------------------"""


class EnergyFrame(wx.Frame):

    def __init__(self, scriptdata, microscope, MotorizedDetector=0):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title="Energy Scan Info", size=(850, 640))

        self.scriptdata = scriptdata
        self.microscope = microscope
        self.mainframe = wx.GetApp().TopWindow

        self.MotorizedDetector = MotorizedDetector

        self.i_current_EP = 0
        self.i_region = 0

        self.maxengregions = 10
        self.maxengsteps = 1000

        self.sortengs = True

        self.SetBackgroundColour("White")

        vboxtop = wx.BoxSizer(wx.VERTICAL)

        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)

        # Energy Points
        hbox212 = wx.BoxSizer(wx.HORIZONTAL)
        vbox212 = wx.BoxSizer(wx.VERTICAL)

        self.button_ep1 = wx.Button(panel, -1, 'Add E P')
        self.Bind(wx.EVT_BUTTON, self.OnAddEPInfo, id=self.button_ep1.GetId())
        vbox212.Add(self.button_ep1, 0, wx.EXPAND)
        vbox212.Add((0, 3))

        self.button_ep2 = wx.Button(panel, -1, 'Edit E P')
        self.Bind(wx.EVT_BUTTON, self.OnEditEPInfo, id=self.button_ep2.GetId())
        vbox212.Add(self.button_ep2, 0, wx.EXPAND)
        vbox212.Add((0, 3))

        self.button_ep3 = wx.Button(panel, -1, 'Remove E P')
        self.Bind(wx.EVT_BUTTON, self.OnRemoveEP, id=self.button_ep3.GetId())
        vbox212.Add(self.button_ep3, 0, wx.EXPAND)
        vbox212.Add((0, 3))

        self.button_ep4 = wx.Button(panel, -1, 'Sort E P')
        self.Bind(wx.EVT_BUTTON, self.OnSortEP, id=self.button_ep4.GetId())
        vbox212.Add(self.button_ep4, 0, wx.EXPAND)
        vbox212.Add((0, 3))

        self.button_ep5 = wx.Button(panel, -1, 'Show E Ps')
        self.Bind(wx.EVT_BUTTON, self.OnShowEPs, id=self.button_ep5.GetId())
        vbox212.Add(self.button_ep5, 0, wx.EXPAND)
        vbox212.Add((0, 3))

        hbox212.Add((23, 0))

        hbox212.Add(vbox212, 1, wx.EXPAND)

        hbox212.Add((23, 0))

        self.tc_ep1 = wx.ListCtrl(panel, -1, size=(560, -1), style=wx.LC_REPORT | wx.LC_NO_HEADER)
        self.tc_ep1.InsertColumn(0, 'EPs')
        self.tc_ep1.SetColumnWidth(0, 530)

        self.Bind(wx.EVT_LIST_ITEM_FOCUSED, self.OnEPListClick, self.tc_ep1)
        hbox212.Add(self.tc_ep1, 4, wx.EXPAND | wx.LEFT, 20)

        vbox.Add(hbox212, 0, wx.LEFT | wx.RIGHT, 15)

        vbox.Add((0, 20))

        ''' Unused        
        #Zoneplate parameters

        sizer2 = wx.StaticBoxSizer(wx.StaticBox(panel, -1, 'Zoneplate parameters'), orient=wx.VERTICAL)
        vbox21 = wx.BoxSizer(wx.VERTICAL)
        vbox21.Add((0,10))          

        grid2 = wx.FlexGridSizer(4, 5, vgap=5, hgap=35)

        t_seng = wx.StaticText(panel, label="Start Energy [eV]")
        self.eb_seng = wx.lib.masked.numctrl.NumCtrl(panel, -1, 
                                              value = str(self.scriptdata.ZP_startE),  
                                              integerWidth = 5,
                                              fractionWidth = 2, min = 0, max = 2000,
                                              limited = False)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnZPStartEng, self.eb_seng)

        t_eeng = wx.StaticText(panel, label="End Energy [eV]")
        self.eb_eeng = wx.lib.masked.numctrl.NumCtrl(panel, -1, 
                                              value = str(self.scriptdata.ZP_stopE),  integerWidth = 5,
                                              fractionWidth = 2, min = 0, max = 2000,
                                              limited = False)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnZPEndEng, self.eb_eeng)

        t_ueng = wx.StaticText(panel, label="delta / eV")


        t_szpx = wx.StaticText(panel, label="ZP-x Startvalue")
        self.eb_szpx = wx.lib.masked.numctrl.NumCtrl(panel, -1, 
                                              value = str(self.scriptdata.ZP_startZPx),  
                                              integerWidth = 5,
                                              fractionWidth = 2, min = 0, max = 2000,
                                              limited = False)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnStartZPx, self.eb_szpx)

        t_ezpx = wx.StaticText(panel, label="ZP-x Endvalue")
        self.eb_ezpx = wx.lib.masked.numctrl.NumCtrl(panel, -1, 
                                              value = str(self.scriptdata.ZP_stopZPx),  
                                              integerWidth = 5,
                                              fractionWidth = 2, min = 0, max = 2000,
                                              limited = False)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnEndZPx, self.eb_ezpx)

        self.t_deltazpx = wx.TextCtrl(panel, -1, value=str(self.scriptdata.ZPx_delta), style = wx.TE_READONLY|wx.NO_BORDER)

        t_szpy = wx.StaticText(panel, label="ZP-y Startvalue")
        self.eb_szpy = wx.lib.masked.numctrl.NumCtrl(panel, -1, 
                                              value = str(self.scriptdata.ZP_startZPy),  
                                              integerWidth = 5,
                                              fractionWidth = 2, min = 0, max = 2000,
                                              limited = False)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnStartZPy, self.eb_szpy)

        t_ezpy = wx.StaticText(panel, label="ZP-y Endvalue")
        self.eb_ezpy = wx.lib.masked.numctrl.NumCtrl(panel, -1, 
                                              value = str(self.scriptdata.ZP_stopZPy),  
                                              integerWidth = 5,
                                              fractionWidth = 2, min = 0, max = 2000,
                                              limited = False)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnEndZPy, self.eb_ezpy)

        self.t_deltazpy = wx.TextCtrl(panel, -1, value=str(self.scriptdata.ZPy_delta), style = wx.TE_READONLY|wx.NO_BORDER)

        t_szpz = wx.StaticText(panel, label="ZP-z Startvalue")
        self.eb_szpz = wx.lib.masked.numctrl.NumCtrl(panel, -1, 
                                              value = str(self.scriptdata.ZP_startZPz),  
                                              integerWidth = 5,
                                              fractionWidth = 2, min = 0, max = 2000,
                                              limited = False)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnStartZPz, self.eb_szpz)

        t_ezpz = wx.StaticText(panel, label="ZP-z Endvalue")
        self.eb_ezpz = wx.lib.masked.numctrl.NumCtrl(panel, -1, 
                                              value = str(self.scriptdata.ZP_stopZPz),  
                                              integerWidth = 5,
                                              fractionWidth = 2, min = 0, max = 2000,
                                              limited = False)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnEndZPz, self.eb_ezpz)

        self.t_deltazpz = wx.TextCtrl(panel, -1, value=str(self.scriptdata.ZPz_delta), style = wx.TE_READONLY|wx.NO_BORDER)


        grid2.Add(t_seng,0)
        grid2.Add(self.eb_seng,1)
        grid2.Add(t_eeng,2)
        grid2.Add(self.eb_eeng,3)
        grid2.Add(t_ueng,4)
        grid2.Add(t_szpx,0)
        grid2.Add(self.eb_szpx,1)
        grid2.Add(t_ezpx,2)
        grid2.Add(self.eb_ezpx,3)
        grid2.Add(self.t_deltazpx,4)
        grid2.Add(t_szpy,0)
        grid2.Add(self.eb_szpy,1)
        grid2.Add(t_ezpy,2)
        grid2.Add(self.eb_ezpy,3)
        grid2.Add(self.t_deltazpy,4)
        grid2.Add(t_szpz,0)
        grid2.Add(self.eb_szpz,1)
        grid2.Add(t_ezpz,2)
        grid2.Add(self.eb_ezpz,3)
        grid2.Add(self.t_deltazpz,4)

        vbox21.Add(grid2, 0, wx.LEFT|wx.RIGHT, 10)

        sizer2.Add(vbox21,-1, wx.EXPAND)

        vbox.Add(sizer2, 0, wx.EXPAND|wx.RIGHT, 30)

        vbox.Add((0,20))

        '''

        # Number of regions
        self.n_eng_regions = len(self.scriptdata.EngRegionlist)

        sizer22 = wx.StaticBoxSizer(wx.StaticBox(panel, -1, 'Energy Regions'), orient=wx.VERTICAL)
        vbox22 = wx.BoxSizer(wx.VERTICAL)
        vbox22.Add((0, 10))

        grid2 = wx.FlexGridSizer(cols=2, vgap=10)
        t_nregions = wx.StaticText(panel, label="Number of energy regions", size=(160, -1))
        self.eb_nregions = wx.lib.intctrl.IntCtrl(panel, size=(80, -1), style=wx.TE_CENTRE,
                                                  value=self.n_eng_regions, limited=True,
                                                  allow_none=True)
        self.eb_nregions.SetMin(1)
        self.eb_nregions.SetMax(self.maxengregions)
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnNEngRegions, self.eb_nregions)

        grid2.Add(t_nregions, 0)
        grid2.Add(self.eb_nregions, 1)

        t_nepoints = wx.StaticText(panel, label="Total number of energy points", size=(200, -1))
        self.eb_nepoints = wx.lib.intctrl.IntCtrl(panel, size=(80, -1), style=wx.TE_CENTRE | wx.TE_READONLY,
                                                  value=0, limited=False)

        grid2.Add(t_nepoints, 2)
        grid2.Add(self.eb_nepoints, 3)

        vbox22.Add(grid2, 0, wx.LEFT | wx.RIGHT, 10)
        vbox22.Add((0, 3))

        self.cb_engsort = wx.CheckBox(panel, -1, '  Sort energy points', size=((142, 20)))
        self.cb_engsort.SetValue(self.scriptdata.sortengs)
        self.Bind(wx.EVT_CHECKBOX, self.OnSortEnergy, self.cb_engsort)
        vbox22.Add(self.cb_engsort, 0, wx.LEFT | wx.RIGHT, 10)
        vbox22.Add((0, 5))

        sizer22.Add(vbox22, 1, wx.EXPAND)

        vbox.Add(sizer22, 0, wx.EXPAND | wx.RIGHT, 30)

        vbox.Add((0, 20))

        # Region file
        self.i_region = 0
        hbox3 = wx.BoxSizer(wx.HORIZONTAL)

        sizer3 = wx.StaticBoxSizer(wx.StaticBox(panel, -1, 'Xanes Region Parameters'), orient=wx.VERTICAL)
        vbox31 = wx.BoxSizer(wx.VERTICAL)
        vbox31.Add((0, 10))

        gridtop = wx.FlexGridSizer(cols=2, vgap=5, hgap=10)

        text1 = wx.StaticText(panel, label="Start energy [eV]")
        self.regbox1 = wx.lib.masked.numctrl.NumCtrl(panel, -1,
                                                     value=float(self.scriptdata.EngRegionlist[0]['Start']),
                                                     integerWidth=5,
                                                     fractionWidth=2, min=0, max=20000,
                                                     limited=False,
                                                     allowNone=True)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnRegStartEng, self.regbox1)

        text2 = wx.StaticText(panel, label="Number of steps")

        self.regbox2 = wx.lib.intctrl.IntCtrl(panel, size=(72, -1), style=wx.TE_RIGHT,
                                              value=int(self.scriptdata.EngRegionlist[0]['Steps']), limited=True)
        self.regbox2.SetMin(0)
        self.regbox2.SetMax(self.maxengsteps)
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnRegNSteps, self.regbox2)

        text3 = wx.StaticText(panel, label="Energy step [eV]")
        self.regbox3 = wx.lib.masked.numctrl.NumCtrl(panel, -1,
                                                     value=float(self.scriptdata.EngRegionlist[0]['EngStep']),
                                                     integerWidth=5,
                                                     fractionWidth=2, min=0, max=2000,
                                                     limited=False,
                                                     allowNone=True)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnRegStepEng, self.regbox3)

        text4 = wx.StaticText(panel, label="Stop energy [eV]")
        self.regbox4 = wx.lib.masked.numctrl.NumCtrl(panel, -1, style=wx.TE_READONLY | wx.NO_BORDER,
                                                     value=float(self.scriptdata.EngRegionlist[0]['Stop']),
                                                     integerWidth=5,
                                                     fractionWidth=2, min=0, max=20000,
                                                     limited=False,
                                                     allowNone=True)

        text5 = wx.StaticText(panel, label="ZPx Start Value")
        self.regbox5 = wx.lib.masked.numctrl.NumCtrl(panel, -1, style=wx.TE_READONLY | wx.NO_BORDER,
                                                     value=float(self.scriptdata.EngRegionlist[0]['ZPxStart']),
                                                     integerWidth=5,
                                                     fractionWidth=2,
                                                     allowNegative=True,
                                                     foregroundColour="Black",
                                                     signedForegroundColour="Black",
                                                     limited=False,
                                                     allowNone=True)

        text6 = wx.StaticText(panel, label="ZPy Start Value")
        self.regbox6 = wx.lib.masked.numctrl.NumCtrl(panel, -1, style=wx.TE_READONLY | wx.NO_BORDER,
                                                     value=float(self.scriptdata.EngRegionlist[0]['ZPyStart']),
                                                     integerWidth=5,
                                                     fractionWidth=2,
                                                     allowNegative=True,
                                                     foregroundColour="Black",
                                                     signedForegroundColour="Black",
                                                     limited=False,
                                                     allowNone=True)

        text7 = wx.StaticText(panel, label="ZPz Start Value")
        self.regbox7 = wx.lib.masked.numctrl.NumCtrl(panel, -1, style=wx.TE_READONLY | wx.NO_BORDER,
                                                     value=float(self.scriptdata.EngRegionlist[0]['ZPzStart']),
                                                     integerWidth=5,
                                                     fractionWidth=2,
                                                     allowNegative=True,
                                                     foregroundColour="Black",
                                                     signedForegroundColour="Black",
                                                     limited=False,
                                                     allowNone=True)

        gridtop.Add(text1, 0)
        gridtop.Add(self.regbox1, 1)
        gridtop.Add(text2, 0)
        gridtop.Add(self.regbox2, 1)
        gridtop.Add(text3, 0)
        gridtop.Add(self.regbox3, 1)
        gridtop.Add(text4, 0)
        gridtop.Add(self.regbox4, 1)
        gridtop.Add(text5, 0)
        gridtop.Add(self.regbox5, 1)
        gridtop.Add(text6, 0)
        gridtop.Add(self.regbox6, 1)
        gridtop.Add(text7, 0)
        gridtop.Add(self.regbox7, 1)

        vbox31.Add(gridtop)

        sizer3.Add(vbox31, 1, wx.LEFT | wx.RIGHT | wx.EXPAND, 10)

        hbox3.Add(sizer3, 0, wx.EXPAND)

        self.tc_rp1 = wx.ListCtrl(panel, -1, size=(560, 0), style=wx.LC_REPORT | wx.LC_NO_HEADER)
        self.tc_rp1.InsertColumn(0, 'RPs')
        self.tc_rp1.SetColumnWidth(0, 540)
        font = wx.SystemSettings_GetFont(wx.SYS_DEFAULT_GUI_FONT)
        font.SetPointSize(8)
        self.tc_rp1.SetFont(font)

        self.Bind(wx.EVT_LIST_ITEM_FOCUSED, self.OnRegionListClick, self.tc_rp1)

        hbox3.Add(self.tc_rp1, 1, wx.EXPAND | wx.LEFT, 20)

        vbox.Add(hbox3, 0)

        hbox4 = wx.BoxSizer(wx.HORIZONTAL)

        button_save = wx.Button(panel, -1, 'Save Eng Regions')
        self.Bind(wx.EVT_BUTTON, self.OnSave, id=button_save.GetId())
        hbox4.Add(button_save, 0, wx.RIGHT, 20)

        button_load = wx.Button(panel, -1, 'Load Eng Regions')
        self.Bind(wx.EVT_BUTTON, self.OnLoad, id=button_load.GetId())
        hbox4.Add(button_load, 0, wx.RIGHT, 20)

        button_close = wx.Button(panel, -1, 'Submit')
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=button_close.GetId())
        hbox4.Add(button_close, 0, wx.RIGHT, 30)

        vbox.Add((0, 20))
        vbox.Add(hbox4, 0, wx.ALIGN_RIGHT)

        panel.SetSizer(vbox)

        vboxtop.Add(panel, 1, wx.EXPAND | wx.LEFT | wx.TOP, 20)

        self.SetSizer(vboxtop)

        self.ShowEPList()
        self.ShowRegionList()

    # Energy Point
    # ----------------------------------------------------------------------
    def OnAddEPInfo(self, event):

        self.thisEP = dict(zip(self.scriptdata.EPkeys, [0, 0, 0, 0, 0, 0, 0, 0]))

        EPFrame(self, self.microscope, self.thisEP, MotorizedDetector=self.MotorizedDetector).Show()

    # ----------------------------------------------------------------------
    def AddEPInfo(self):

        self.scriptdata.EPlist.append(self.thisEP.copy())

        self.ShowEPList()

        self.SetRegionParamsAll()
        self.ShowRegionList()

    # ----------------------------------------------------------------------
    def OnEditEPInfo(self, event):
        if (self.i_current_EP < len(self.scriptdata.EPlist)):
            # get the data from current EP
            currEP = self.scriptdata.EPlist[self.i_current_EP]
            self.thisEP = currEP.copy()

            EPFrame(self, self.microscope, self.thisEP, MotorizedDetector=self.MotorizedDetector, edit=True).Show()

    # ----------------------------------------------------------------------
    def EditEPInfo(self):

        self.scriptdata.EPlist[self.i_current_EP] = self.thisEP.copy()

        self.ShowEPList()

        self.SetRegionParamsAll()
        self.ShowRegionList()

    # ----------------------------------------------------------------------
    def OnRemoveEP(self, event):
        if (self.i_current_EP < len(self.scriptdata.EPlist)):
            self.scriptdata.EPlist.pop(self.i_current_EP)

            self.ShowEPList()

            self.SetRegionParamsAll()
            self.ShowRegionList()

    # ----------------------------------------------------------------------
    def OnSortEP(self, event):
        if (len(self.scriptdata.EPlist) > 1):
            self.scriptdata.EPlist.sort(key=lambda k: k['Energy'])

            self.ShowEPList()

            self.SetRegionParamsAll()
            self.ShowRegionList()

    # ----------------------------------------------------------------------
    def ShowEPList(self):

        self.tc_ep1.DeleteAllItems()

        if self.MotorizedDetector == 0:
            for i in range(len(self.scriptdata.EPlist)):
                self.tc_ep1.InsertStringItem(i,
                                             'Energy: {0:07.2f}   ZPx: {1:07.2f}   ZPy: {2:07.2f}   ZPz: {3:07.2f}  ZPID: {4:2.0f}'.format(
                                                 self.scriptdata.EPlist[i]['Energy'],
                                                 self.scriptdata.EPlist[i]['ZPx'],
                                                 self.scriptdata.EPlist[i]['ZPy'],
                                                 self.scriptdata.EPlist[i]['ZPz'],
                                                 self.scriptdata.EPlist[i]['ZPID']))
        else:
            for i in range(len(self.scriptdata.EPlist)):
                self.tc_ep1.InsertStringItem(i,
                                             'Energy: {0:07.2f}   ZPx: {1:07.2f}   ZPy: {2:07.2f}   ZPz: {3:07.2f}  ZPID: {4:2.0f}  Detx: {5:07.2f}   Dety: {6:07.2f}   Detz: {7:07.2f}'.format(
                                                 self.scriptdata.EPlist[i]['Energy'],
                                                 self.scriptdata.EPlist[i]['ZPx'],
                                                 self.scriptdata.EPlist[i]['ZPy'],
                                                 self.scriptdata.EPlist[i]['ZPz'],
                                                 self.scriptdata.EPlist[i]['ZPID'],
                                                 self.scriptdata.EPlist[i]['Detx'],
                                                 self.scriptdata.EPlist[i]['Dety'],
                                                 self.scriptdata.EPlist[i]['Detz']))

            # ----------------------------------------------------------------------

    def OnEPListClick(self, event):

        self.i_current_EP = event.m_itemIndex

    # ----------------------------------------------------------------------
    def OnShowEPs(self, event):

        if len(self.scriptdata.EPlist) >= 2:
            eng = self.scriptdata.ZP_EngList
            zpx = self.scriptdata.ZP_xList
            zpy = self.scriptdata.ZP_yList
            zpz = self.scriptdata.ZP_zList

            import matplotlib.pyplot as plt
            plt.figure(1, figsize=(5, 9))
            plt.subplot(311)
            plt.plot(eng, zpx, 'bo')
            for i in range(len(self.zpxfp)):
                [a, b] = self.zpxfp[i]
                plt.plot(eng, a * eng + b, '-')
            plt.title('ZP X-axis')

            plt.subplot(312)
            plt.plot(eng, zpy, 'bo')
            for i in range(len(self.zpyfp)):
                [a, b] = self.zpyfp[i]
                plt.plot(eng, a * eng + b, '-')
            plt.title('ZP Y-axis')

            plt.subplot(313)
            plt.plot(eng, zpz, 'bo')
            for i in range(len(self.zpzfp)):
                [a, b] = self.zpzfp[i]
                plt.plot(eng, a * eng + b, '-')
            plt.title('ZP Z-axis')
            plt.show()

    '''  ----- Unused -----
#Zone plate Parameters 
#----------------------------------------------------------------------    
    def OnZPStartEng(self, event):   
        self.scriptdata.ZP_startE = self.eb_seng.GetValue()

        self.SetZPParams()
        self.SetRegionParamsAll()

#----------------------------------------------------------------------    
    def OnZPEndEng(self, event):   
        self.scriptdata.ZP_stopE = self.eb_eeng.GetValue()    

        self.SetZPParams()
        self.SetRegionParamsAll()

#----------------------------------------------------------------------    
    def OnStartZPx(self, event):   
        self.scriptdata.ZP_startZPx = self.eb_szpx.GetValue()

        self.SetZPParams()
        self.SetRegionParamsAll()

#----------------------------------------------------------------------    
    def OnEndZPx(self, event):   
        self.scriptdata.ZP_stopZPx = self.eb_ezpx.GetValue()   

        self.SetZPParams()
        self.SetRegionParamsAll()

#----------------------------------------------------------------------    
    def OnStartZPy(self, event):   
        self.scriptdata.ZP_startZPy = self.eb_szpy.GetValue()

        self.SetZPParams()
        self.SetRegionParamsAll()

#----------------------------------------------------------------------    
    def OnEndZPy(self, event):   
        self.scriptdata.ZP_stopZPy = self.eb_ezpy.GetValue()  

        self.SetZPParams()
        self.SetRegionParamsAll()

#----------------------------------------------------------------------    
    def OnStartZPz(self, event):   
        self.scriptdata.ZP_startZPz = self.eb_szpz.GetValue()

        self.SetZPParams()
        self.SetRegionParamsAll()

#----------------------------------------------------------------------    
    def OnEndZPz(self, event):   
        self.scriptdata.ZP_stopZPz = self.eb_ezpz.GetValue()  

        self.SetZPParams()
        self.SetRegionParamsAll()

#----------------------------------------------------------------------           
    def SetZPParams(self): 
        ZP_deltaE = self.scriptdata.ZP_stopE - self.scriptdata.ZP_startE
        if ZP_deltaE == 0:
            ZP_deltaE=0.0001;
            print 'The starting E and the ending E should not be the same....'

        ZPx_delta = self.scriptdata.ZP_stopZPx - self.scriptdata.ZP_startZPx
        ZPy_delta = self.scriptdata.ZP_stopZPy - self.scriptdata.ZP_startZPy
        ZPz_delta = self.scriptdata.ZP_stopZPz - self.scriptdata.ZP_startZPz
        ZPx_delta_per_eV = ZPx_delta/ZP_deltaE
        ZPy_delta_per_eV = ZPy_delta/ZP_deltaE
        ZPz_delta_per_eV = ZPz_delta/ZP_deltaE
        self.scriptdata.ZPx_delta = ZPx_delta_per_eV
        self.scriptdata.ZPy_delta = ZPy_delta_per_eV
        self.scriptdata.ZPz_delta = ZPz_delta_per_eV

        #refresh text on the GUI
        self.t_deltazpx.SetValue(str(self.scriptdata.ZPx_delta))
        self.t_deltazpy.SetValue(str(self.scriptdata.ZPy_delta))
        self.t_deltazpz.SetValue(str(self.scriptdata.ZPz_delta)) '''

    # Region parameters
    # ----------------------------------------------------------------------
    def OnNEngRegions(self, event):
        ctl = event.GetEventObject()
        value = ctl.GetValue()

        if value < 1:
            return

        self.n_eng_regions = value

        oldn_eregions = len(self.scriptdata.EngRegionlist)

        if (self.n_eng_regions > oldn_eregions):
            # Add new default regions
            for i in range(self.n_eng_regions - oldn_eregions):
                self.scriptdata.EngRegionlist.append(
                    dict(zip(self.scriptdata.ERegkeys, [0.0, 1, 0.5, 0.5, 0.0, 0.0, 0.0])))
        elif (oldn_eregions > self.n_eng_regions):
            # Remove regions
            for i in range(oldn_eregions - self.n_eng_regions):
                self.scriptdata.EngRegionlist.pop()

        self.i_region = 0
        self.SetRegionParamsAll()
        self.ShowRegionList()
        self.ShowValsRegionListClick()

    # ----------------------------------------------------------------------
    def OnRegionListClick(self, event):

        self.i_region = event.m_itemIndex
        self.ShowValsRegionListClick()
        self.ShowRegionList()

    # ----------------------------------------------------------------------
    def ShowValsRegionListClick(self):

        # Show the region values
        self.regbox1.SetValue(float(self.scriptdata.EngRegionlist[self.i_region]['Start']))
        self.regbox2.SetValue(int(self.scriptdata.EngRegionlist[self.i_region]['Steps']))
        self.regbox3.SetValue(float(self.scriptdata.EngRegionlist[self.i_region]['EngStep']))
        self.regbox4.SetValue(float(self.scriptdata.EngRegionlist[self.i_region]['Stop']))
        self.regbox5.SetValue(float(self.scriptdata.EngRegionlist[self.i_region]['ZPxStart']))
        self.regbox6.SetValue(float(self.scriptdata.EngRegionlist[self.i_region]['ZPyStart']))
        self.regbox7.SetValue(float(self.scriptdata.EngRegionlist[self.i_region]['ZPzStart']))

    # ----------------------------------------------------------------------
    def OnRegStartEng(self, event):
        self.scriptdata.EngRegionlist[self.i_region]['Start'] = self.regbox1.GetValue()
        self.SetRegionParams(self.i_region)
        self.ShowRegionList()

    # ----------------------------------------------------------------------
    def OnRegNSteps(self, event):
        ctl = event.GetEventObject()
        self.scriptdata.EngRegionlist[self.i_region]['Steps'] = ctl.GetValue()
        self.SetRegionParams(self.i_region)
        self.ShowRegionList()

    # ----------------------------------------------------------------------
    def OnRegStepEng(self, event):
        self.scriptdata.EngRegionlist[self.i_region]['EngStep'] = self.regbox3.GetValue()
        self.SetRegionParams(self.i_region)
        self.ShowRegionList()

    # ----------------------------------------------------------------------
    def OnSortEnergy(self, event):
        if self.cb_engsort.GetValue():
            self.scriptdata.sortengs = True
        else:
            self.scriptdata.sortengs = False

    # ----------------------------------------------------------------------
    def ShowRegionList(self):

        self.tc_rp1.DeleteAllItems()

        for i in range(len(self.scriptdata.EngRegionlist)):
            self.tc_rp1.InsertStringItem(i,
                                         'Region:{0:02d}  Start:{1:7.2f}  Steps:{2:3d}  Energy-step:{3:7.2f}  Stop:{4:7.2f}  ZPx:{5:7.2f}  ZPy:{6:7.2f}  ZPz:{7:7.2f}'.format(
                                             i + 1,
                                             self.scriptdata.EngRegionlist[i]['Start'],
                                             int(self.scriptdata.EngRegionlist[i]['Steps']),
                                             self.scriptdata.EngRegionlist[i]['EngStep'],
                                             self.scriptdata.EngRegionlist[i]['Stop'],
                                             self.scriptdata.EngRegionlist[i]['ZPxStart'],
                                             self.scriptdata.EngRegionlist[i]['ZPyStart'],
                                             self.scriptdata.EngRegionlist[i]['ZPzStart']))

        self.regbox4.SetValue(float(self.scriptdata.EngRegionlist[self.i_region]['Stop']))
        self.regbox5.SetValue(float(self.scriptdata.EngRegionlist[self.i_region]['ZPxStart']))
        self.regbox6.SetValue(float(self.scriptdata.EngRegionlist[self.i_region]['ZPyStart']))
        self.regbox7.SetValue(float(self.scriptdata.EngRegionlist[self.i_region]['ZPzStart']))

        totalEpoints = self.GetTotalNumOfEnergyPoints()

        self.eb_nepoints.SetValue(totalEpoints)

        self.tc_rp1.SetItemBackgroundColour(self.i_region, 'light blue')

    # ----------------------------------------------------------------------
    def GetTotalNumOfEnergyPoints(self):

        # Calculate total number of energy steps
        totalEpoints = 0
        for ii in range(self.scriptdata.RegionParams.shape[0]):
            totalEpoints = totalEpoints + self.scriptdata.RegionParams[ii, 1] + 1

        return int(totalEpoints)

    # ----------------------------------------------------------------------
    def SetRegionParamsAll(self):

        self.scriptdata.RegionParams = np.zeros((self.n_eng_regions, 7))
        for i in range(self.n_eng_regions):
            self.SetRegionParams(i)

        self.ShowRegionList()

    # ----------------------------------------------------------------------
    def SetRegionParams(self, RegionNbr):

        i = RegionNbr
        startE = float(self.scriptdata.EngRegionlist[i]['Start'])
        self.scriptdata.RegionParams[RegionNbr, 0] = startE
        self.scriptdata.RegionParams[RegionNbr, 1] = float(self.scriptdata.EngRegionlist[i]['Steps'])
        self.scriptdata.RegionParams[RegionNbr, 2] = float(self.scriptdata.EngRegionlist[i]['EngStep'])
        endEnergy = self.scriptdata.RegionParams[RegionNbr, 0] + self.scriptdata.RegionParams[RegionNbr, 1] * \
                    self.scriptdata.RegionParams[RegionNbr, 2]

        self.scriptdata.RegionParams[RegionNbr, 3] = endEnergy

        self.scriptdata.ZP_EngList = []
        self.scriptdata.ZP_xList = []
        self.scriptdata.ZP_yList = []
        self.scriptdata.ZP_zList = []

        if (len(self.scriptdata.EPlist) > 1):
            sortedEPlist = list(self.scriptdata.EPlist)
            sortedEPlist.sort(key=lambda k: k['Energy'])

            # Calculate interpolation functions:
            for ie in range(len(sortedEPlist)):
                self.scriptdata.ZP_EngList.append(sortedEPlist[ie]['Energy'])
                self.scriptdata.ZP_xList.append(sortedEPlist[ie]['ZPx'])
                self.scriptdata.ZP_yList.append(sortedEPlist[ie]['ZPy'])
                self.scriptdata.ZP_zList.append(sortedEPlist[ie]['ZPz'])

            self.scriptdata.ZP_EngList = np.array(self.scriptdata.ZP_EngList)
            self.scriptdata.ZP_xList = np.array(self.scriptdata.ZP_xList)
            self.scriptdata.ZP_yList = np.array(self.scriptdata.ZP_yList)
            self.scriptdata.ZP_zList = np.array(self.scriptdata.ZP_zList)

            zpstartx, zpstarty, zpstartz = self.CalcZPPosition(startE)

            self.scriptdata.RegionParams[RegionNbr, 4] = zpstartx[0]
            self.scriptdata.RegionParams[RegionNbr, 5] = zpstarty[0]
            self.scriptdata.RegionParams[RegionNbr, 6] = zpstartz[0]



        else:
            self.scriptdata.RegionParams[RegionNbr, 4] = 0.
            self.scriptdata.RegionParams[RegionNbr, 5] = 0.
            self.scriptdata.RegionParams[RegionNbr, 6] = 0.

        self.scriptdata.EngRegionlist[i]['Stop'] = endEnergy
        self.scriptdata.EngRegionlist[i]['ZPxStart'] = self.scriptdata.RegionParams[RegionNbr, 4]
        self.scriptdata.EngRegionlist[i]['ZPyStart'] = self.scriptdata.RegionParams[RegionNbr, 5]
        self.scriptdata.EngRegionlist[i]['ZPzStart'] = self.scriptdata.RegionParams[RegionNbr, 6]

    # ----------------------------------------------------------------------
    def CalcZPPosition(self, zpeng):

        sg = scriptgen.scriptgen()

        zpbounds, self.zpxfp, self.zpyfp, self.zpzfp = sg.CalcXanesFitParams(self.scriptdata.EPlist)

        zpelist = [zpeng]

        zpx, zpy, zpz = sg.GetZPPositions(zpelist)

        return zpx, zpy, zpz

    # ----------------------------------------------------------------------
    def extrap(self, x, xp, yp):
        """np.interp function with linear extrapolation"""
        y = np.interp(x, xp, yp)
        y = np.where(x < xp[0], yp[0] + (x - xp[0]) * (yp[0] - yp[1]) / (xp[0] - xp[1]), y)
        y = np.where(x > xp[-1], yp[-1] + (x - xp[-1]) * (yp[-1] - yp[-2]) / (xp[-1] - xp[-2]), y)
        return y

    # ----------------------------------------------------------------------
    def OnSave(self, evt):

        wildcard = "PTWZ files (*.ptwz)|*.ptwz"
        dialog = wx.FileDialog(None, "Please select a ptwz file",
                               wildcard=wildcard,
                               style=wx.SAVE)
        if dialog.ShowModal() == wx.ID_OK:
            filepath = dialog.GetPath()
            filename = dialog.GetFilename()
            filedir = dialog.GetDirectory()

            pickleFile = open(filepath, 'wb')
            fversion = 1.0
            pickle.dump(fversion, pickleFile, pickle.HIGHEST_PROTOCOL)

            pickle.dump(self.scriptdata.RegionParams, pickleFile, pickle.HIGHEST_PROTOCOL)
            neng = len(self.scriptdata.EngRegionlist)
            pickle.dump(neng, pickleFile, pickle.HIGHEST_PROTOCOL)
            for i in range(neng):
                pickle.dump(self.scriptdata.EngRegionlist[i], pickleFile, pickle.HIGHEST_PROTOCOL)

            nzp = len(self.scriptdata.EPlist)
            pickle.dump(nzp, pickleFile, pickle.HIGHEST_PROTOCOL)
            for i in range(nzp):
                pickle.dump(self.scriptdata.EPlist[i], pickleFile, pickle.HIGHEST_PROTOCOL)

            pickleFile.close()

        self.Raise()

    # ----------------------------------------------------------------------
    def OnLoad(self, evt):

        wildcard = "PTWZ files (*.ptwz)|*.ptwz"
        dialog = wx.FileDialog(None, "Please select a ptwz file",
                               wildcard=wildcard,
                               style=wx.OPEN)

        if dialog.ShowModal() == wx.ID_OK:
            filepath = dialog.GetPath()
        else:
            print
            'Warning - Could not read file.'
            return

        pickleFile = open(filepath, 'rb')
        fversion = pickle.load(pickleFile)

        self.scriptdata.RegionParams = pickle.load(pickleFile)

        self.scriptdata.EngRegionlist = []
        self.scriptdata.EPlist = []
        neng = pickle.load(pickleFile)
        for i in range(neng):
            self.scriptdata.EngRegionlist.append(pickle.load(pickleFile))

        nzp = pickle.load(pickleFile)
        for i in range(nzp):
            self.scriptdata.EPlist.append(pickle.load(pickleFile))

        pickleFile.close()

        self.n_eng_regions = neng
        self.eb_nregions.SetValue(self.n_eng_regions)

        self.ShowEPList()

        self.SetRegionParamsAll()
        self.ShowRegionList()

        self.Raise()

    # ----------------------------------------------------------------------
    def OnClose(self, evt):
        self.mainframe.ShowEnabledTC()
        self.Close(True)


# ------------------------------------------------------------------------------------------------
class EPFrame(wx.Frame):

    def __init__(self, parent, microscope, thisEP, MotorizedDetector=0, edit=False):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title="Energy Point Info", size=(300, 420))

        self.microscope = microscope
        self.thisEP = thisEP
        self.parent = parent
        self.edit = edit

        self.MotorizedDetector = MotorizedDetector

        self.energy = self.thisEP['Energy']
        self.ZPx = self.thisEP['ZPx']
        self.ZPy = self.thisEP['ZPy']
        self.ZPz = self.thisEP['ZPz']
        self.ZPID = self.thisEP['ZPID']

        self.detx = self.thisEP['Detx']
        self.dety = self.thisEP['Dety']
        self.detz = self.thisEP['Detz']

        self.SetBackgroundColour("White")

        vboxtop = wx.BoxSizer(wx.VERTICAL)

        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)

        vbox.Add((0, 10))

        text1 = wx.StaticText(panel, label="X-ray E [eV]")
        self.ebox1 = wx.lib.masked.numctrl.NumCtrl(panel, -1,
                                                   value=self.thisEP['Energy'], integerWidth=5,
                                                   fractionWidth=2,
                                                   limited=False,
                                                   allowNone=True)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnEnergy, self.ebox1)

        text2 = wx.StaticText(panel, label="ZPx [microns]")
        self.ebox2 = wx.lib.masked.numctrl.NumCtrl(panel, -1,
                                                   value=self.thisEP['ZPx'], integerWidth=5,
                                                   fractionWidth=2,
                                                   allowNegative=True,
                                                   foregroundColour="Black",
                                                   signedForegroundColour="Black",
                                                   limited=False,
                                                   allowNone=True)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnZPx, self.ebox2)

        text3 = wx.StaticText(panel, label="ZPy [microns]")
        self.ebox3 = wx.lib.masked.numctrl.NumCtrl(panel, -1,
                                                   value=self.thisEP['ZPy'], integerWidth=5,
                                                   fractionWidth=2,
                                                   allowNegative=True,
                                                   foregroundColour="Black",
                                                   signedForegroundColour="Black",
                                                   limited=False,
                                                   allowNone=True)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnZPy, self.ebox3)

        text4 = wx.StaticText(panel, label="ZPz [microns]")
        self.ebox4 = wx.lib.masked.numctrl.NumCtrl(panel, -1,
                                                   value=self.thisEP['ZPz'], integerWidth=5,
                                                   fractionWidth=2,
                                                   allowNegative=True,
                                                   foregroundColour="Black",
                                                   signedForegroundColour="Black",
                                                   limited=False,
                                                   allowNone=True)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnZPz, self.ebox4)

        text5 = wx.StaticText(panel, label="ZP ID")
        self.ebox5 = wx.lib.masked.numctrl.NumCtrl(panel, -1,
                                                   value=self.thisEP['ZPID'], integerWidth=8,
                                                   fractionWidth=0,
                                                   allowNegative=False,
                                                   foregroundColour="Black",
                                                   signedForegroundColour="Black",
                                                   limited=True,
                                                   allowNone=True)
        self.ebox5.SetMin(0)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnZPID, self.ebox5)

        text6 = wx.StaticText(panel, label="Det X [microns]")
        self.ebox6 = wx.lib.masked.numctrl.NumCtrl(panel, -1,
                                                   value=self.thisEP['Detx'], integerWidth=5,
                                                   fractionWidth=2,
                                                   allowNegative=True,
                                                   foregroundColour="Black",
                                                   signedForegroundColour="Black",
                                                   limited=False,
                                                   allowNone=True)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnDetx, self.ebox6)

        text7 = wx.StaticText(panel, label="Det Y [microns]")
        self.ebox7 = wx.lib.masked.numctrl.NumCtrl(panel, -1,
                                                   value=self.thisEP['Dety'], integerWidth=5,
                                                   fractionWidth=2,
                                                   allowNegative=True,
                                                   foregroundColour="Black",
                                                   signedForegroundColour="Black",
                                                   limited=False,
                                                   allowNone=True)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnDety, self.ebox7)

        text8 = wx.StaticText(panel, label="Det Z [microns]")
        self.ebox8 = wx.lib.masked.numctrl.NumCtrl(panel, -1,
                                                   value=self.thisEP['Detz'], integerWidth=5,
                                                   fractionWidth=2,
                                                   allowNegative=True,
                                                   foregroundColour="Black",
                                                   signedForegroundColour="Black",
                                                   limited=False,
                                                   allowNone=True)
        self.Bind(wx.lib.masked.EVT_NUM, self.OnDetz, self.ebox8)

        if self.MotorizedDetector == 0:
            text6.Disable()
            self.ebox6.Disable()
            text7.Disable()
            self.ebox7.Disable()
            text8.Disable()
            self.ebox8.Disable()

        gridtop = wx.FlexGridSizer(cols=2, vgap=5, hgap=10)
        gridtop.Add(text1, 0)
        gridtop.Add(self.ebox1, 1)
        gridtop.Add(wx.StaticText(panel, label=' '), 0)
        gridtop.Add(wx.StaticText(panel, label=' '), 1)

        gridtop.Add(text2, 0)
        gridtop.Add(self.ebox2, 1)
        gridtop.Add(text3, 0)
        gridtop.Add(self.ebox3, 1)
        gridtop.Add(text4, 0)
        gridtop.Add(self.ebox4, 1)
        gridtop.Add(text5, 0)
        gridtop.Add(self.ebox5, 1)
        gridtop.Add(wx.StaticText(panel, label=' '), 0)
        gridtop.Add(wx.StaticText(panel, label=' '), 1)

        gridtop.Add(text6, 0)
        gridtop.Add(self.ebox6, 1)
        gridtop.Add(text7, 0)
        gridtop.Add(self.ebox7, 1)
        gridtop.Add(text8, 0)
        gridtop.Add(self.ebox8, 1)

        vbox.Add(gridtop, 1, wx.ALL | wx.EXPAND, 20)

        button_load = wx.Button(panel, -1, 'Load Positions from .xrm File')
        self.Bind(wx.EVT_BUTTON, self.OnLoadPositions, id=button_load.GetId())
        vbox.Add(button_load, -1, wx.LEFT, 20)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add((20, 0))
        button_ok = wx.Button(panel, -1, 'Accept')
        self.Bind(wx.EVT_BUTTON, self.OnOK, id=button_ok.GetId())
        hbox.Add(button_ok, -1, wx.TOP, 10)

        button_close = wx.Button(panel, -1, 'Discard')
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=button_close.GetId())
        hbox.Add(button_close, -1, wx.TOP, 10)
        hbox.Add((20, 0))
        vbox.Add(hbox, 1, wx.EXPAND)
        vbox.Add((20, 0))

        panel.SetSizer(vbox)

        vboxtop.Add(panel, 1, wx.EXPAND)

        self.SetSizer(vboxtop)

    # ----------------------------------------------------------------------
    def OnEnergy(self, event):
        self.energy = self.ebox1.GetValue()

    # ----------------------------------------------------------------------
    def OnZPx(self, event):
        self.ZPx = self.ebox2.GetValue()

    # ----------------------------------------------------------------------
    def OnZPy(self, event):
        self.ZPy = self.ebox3.GetValue()

    # ----------------------------------------------------------------------
    def OnZPz(self, event):
        self.ZPz = self.ebox4.GetValue()

    # ----------------------------------------------------------------------
    def OnZPID(self, event):
        self.ZPID = self.ebox5.GetValue()

    # ----------------------------------------------------------------------
    def OnDetx(self, event):
        self.detx = self.ebox6.GetValue()

    # ----------------------------------------------------------------------
    def OnDety(self, event):
        self.dety = self.ebox7.GetValue()

    # ----------------------------------------------------------------------
    def OnDetz(self, event):
        self.detz = self.ebox8.GetValue()

    # ----------------------------------------------------------------------
    def OnLoadPositions(self, evt):

        wildcard = "XRM files (*.xrm)|*.xrm"
        dialog = wx.FileDialog(None, "Please select a xrm file",
                               wildcard=wildcard,
                               style=wx.OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            filepath = dialog.GetPath()
            filename = dialog.GetFilename()
            filedir = dialog.GetDirectory()
        else:
            return

        zpx = []
        zpy = []
        zpz = []
        detx = []
        dety = []
        detz = []

        if self.microscope:
            zpx, zpy, zpz, detx, dety, detz = self.microscope.get_xrm_zpinfo(filepath)

        else:

            Cxrm = txrm.txrm()

            xrminfo = Cxrm.get_xrm_info(filepath)

            NumofMotors, zpxxx, zpyyy, zpzzz, detxxx, detyyy, detzzz = Cxrm.get_zpmotor_index(xrminfo)

            zpx = xrminfo['MotPos'][zpxxx::zpxxx + 1]
            zpy = xrminfo['MotPos'][zpyyy::zpyyy + 1]
            zpz = xrminfo['MotPos'][zpzzz::zpzzz + 1]

            if detxxx != -1:
                detx = xrminfo['MotPos'][detxxx::detxxx + 1]
            if detyyy != -1:
                dety = xrminfo['MotPos'][detyyy::detyyy + 1]
            if detzzz != -1:
                detz = xrminfo['MotPos'][detzzz::detzzz + 1]

        if len(zpx) > 0:
            self.ZPx = zpx[0]
            self.ebox2.SetValue(self.ZPx)
            self.ZPy = zpy[0]
            self.ebox3.SetValue(self.ZPy)
            self.ZPz = zpz[0]
            self.ebox4.SetValue(self.ZPz)

        if len(detx) > 0:
            self.detx = detx[0]
            self.ebox6.SetValue(self.detx)
        if len(dety) > 0:
            self.dety = dety[0]
            self.ebox7.SetValue(self.dety)
        if len(detz) > 0:
            self.detz = detz[0]
            self.ebox8.SetValue(self.detz)

        self.Raise()

    # ----------------------------------------------------------------------
    def OnOK(self, evt):

        if (self.edit == False):
            for i in range(len(self.parent.scriptdata.EPlist)):
                if self.energy == self.parent.scriptdata.EPlist[i]['Energy']:
                    wx.MessageBox("Energy value is already used. Please input different energy point.")
                    return

        self.thisEP['Energy'] = self.energy
        self.thisEP['ZPx'] = self.ZPx
        self.thisEP['ZPy'] = self.ZPy
        self.thisEP['ZPz'] = self.ZPz
        self.thisEP['ZPID'] = self.ZPID

        if self.MotorizedDetector == 1:
            self.thisEP['Detx'] = self.detx
            self.thisEP['Dety'] = self.dety
            self.thisEP['Detz'] = self.detz

        if (self.edit == False):
            self.parent.AddEPInfo()
        else:
            self.parent.EditEPInfo()

        self.Close(True)

        self.parent.Raise()

    # ----------------------------------------------------------------------
    def OnClose(self, evt):

        self.Destroy()

        self.parent.Raise()


""" ------------------------------------------------------------------------------------------------"""


class MainFrame(wx.Frame):
    def __init__(self, parent, id, title, Winsizex, Winsizey):
        wx.Frame.__init__(self, parent, id, title, size=(Winsizex, Winsizey))

        self.epframe = None

        self.scroll = wx.ScrolledWindow(self, -1, style=wx.VSCROLL)
        self.scroll.SetScrollbars(1, 1, 1, 1)

        self.scriptdata = ScriptData()
        self.scriptlist = []

        self.txrmloaded = False

        self.multiscript = []
        self.multiscriptlist = []
        self.multiscan_foldernum = 1
        self.multiscan_namelist = []

        self.multiscan_userepeat = 0
        self.multiscan_repeat = 1
        self.multiscan_wait = 0

        self.scanrunning = False
        self.pausescan = False
        self.resumescan = False
        self.resumeuncompleted = False
        self.stopscan = False
        self.refreshscriptlist = False
        self.iscanstopped = 0

        self.TurnAgitatorOff = True

        self.IonCurrentReading = 0

        self.HaveRetiga = False

        # Read setting from a TW_settings.txt file
        self.GetSettingsFromFile()

        if UseTimer:
            print
            "Timer Enabled."

        try:
            self.microscope = None
            if self.HaveXradiaAPI:
                import xradia_helper
                self.microscope = xradia_helper.Microscope(motorwait=self.XradiaMotorWait, verbose=v, timer=UseTimer)

        except:
            self.microscope = None
            error_message = ('ERROR - Could not initialize XRADIA API.')
            logging.error(error_message, exc_info=True)
            wx.MessageBox(error_message, 'INIT ERROR',
                          parent=self, style=wx.OK | wx.ICON_ERROR)

        if self.HaveXradiaAPI:
            try:
                # Check if Retiga camera averaging is available
                self.HaveRetiga = self.microscope.GetRetigaCameraAveraging()
            except:
                error_message = ('ERROR - Could not get camera type. Retiga averaging will be disabled')
                logging.error(error_message, exc_info=True)
                wx.MessageBox(error_message, 'INIT ERROR',
                              parent=self, style=wx.OK | wx.ICON_ERROR)

        if self.HaveMono == True:
            try:
                # Check to see if monochromator is available
                self.mono = ssrl_sis.Csis(verbose=v)
                if self.mono.soc is None:
                    error_message = (
                        'ERROR - Could not initialize SSRL SIS connection for monochromator control, continuing without monochromator control.')
                    logging.error(error_message)
                    wx.MessageBox(error_message, 'INIT ERROR', parent=self, style=wx.OK | wx.ICON_ERROR)
                    self.mono = None

                if self.mono:
                    self.mono.SetMonoSoftLimits(self.mono_use_limits, self.mono_l, self.mono_u)
            except:
                self.mono = None
                error_message = (
                    'ERROR - Could not initialize SSRL SIS connection for monochromator control, continuing without monochromator control.')
                logging.error(error_message, exc_info=True)
                wx.MessageBox(error_message, 'INIT ERROR', parent=self, style=wx.OK | wx.ICON_ERROR)

        else:
            self.mono = None

        self.SetBackgroundColour("White")

        vboxT = wx.BoxSizer(wx.VERTICAL)

        # panel 1
        panel1 = wx.Panel(self.scroll, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)

        hbox11 = wx.BoxSizer(wx.HORIZONTAL)

        self.cb_1 = wx.CheckBox(panel1, -1, '  Repeat the defined scan')
        self.Bind(wx.EVT_CHECKBOX, self.OnCBRepeatScan, self.cb_1)

        self.cb_exclusiveControl = wx.CheckBox(panel1, -1, '  Enable Exclusive control of the instrument', (600, 5))
        self.Bind(wx.EVT_CHECKBOX, self.OnCBEnableDisableExclusiveControl, self.cb_exclusiveControl)

        self.ntc1 = wx.lib.intctrl.IntCtrl(panel1, size=(60, -1), style=wx.TE_CENTRE, value=self.scriptdata.NRepeatScan,
                                           limited=True, allow_none=True)
        self.ntc1.SetMin(1)
        self.ntc1.SetMax(1000)
        self.ntc1.Disable()
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnRepeatScanN, self.ntc1)

        tc1 = wx.StaticText(panel1, label='times.')

        hbox11.Add((15, 0))
        hbox11.Add(self.cb_1, 0, wx.EXPAND)
        hbox11.Add((5, 0))
        hbox11.Add(self.ntc1, 0, wx.TOP, 2)
        hbox11.Add((5, 0))
        hbox11.Add(tc1, 0, wx.EXPAND | wx.TOP, 7)

        tc12 = wx.StaticText(panel1, label='Wait between scans ')
        self.ntc12 = wx.lib.intctrl.IntCtrl(panel1, size=(60, -1), style=wx.TE_CENTRE, value=self.scriptdata.waitsecs,
                                            limited=True)
        self.ntc12.SetMin(0)
        self.ntc12.SetMax(9000)
        self.ntc12.Disable()
        tc13 = wx.StaticText(panel1, label='seconds.')

        hbox11.Add((15, 0))
        hbox11.Add(tc12, 0, wx.EXPAND | wx.TOP, 7)
        hbox11.Add((5, 0))
        hbox11.Add(self.ntc12, 0, wx.TOP, 2)
        hbox11.Add((5, 0))
        hbox11.Add(tc13, 0, wx.EXPAND | wx.TOP, 7)

        vbox1.Add(hbox11, 0, wx.EXPAND)

        panel1.SetSizer(vbox1)

        # panel 2
        panel2 = wx.Panel(self.scroll, -1)

        sbox2 = wx.StaticBox(panel2, -1, 'Scan info')
        sizer2 = wx.StaticBoxSizer(sbox2, orient=wx.VERTICAL)
        vbox21 = wx.BoxSizer(wx.VERTICAL)

        # Energy scan
        hbox211 = wx.BoxSizer(wx.HORIZONTAL)
        vbox211 = wx.BoxSizer(wx.VERTICAL)

        self.cb_eng = wx.CheckBox(panel2, -1, '  Enable Energy Scan', size=((142, 20)))
        self.Bind(wx.EVT_CHECKBOX, self.OnEnableEnergy, self.cb_eng)
        vbox211.Add(self.cb_eng, 0, wx.EXPAND)

        self.button_eng = wx.Button(panel2, -1, 'Edit Energy Info')
        self.Bind(wx.EVT_BUTTON, self.OnEditEngInfo, id=self.button_eng.GetId())
        vbox211.Add(self.button_eng, 0, wx.EXPAND)
        vbox211.Add((0, 3))

        hbox_layer = wx.BoxSizer(wx.HORIZONTAL)
        tl = wx.StaticText(panel2, label="Layer:")
        hbox_layer.Add(tl, 0, wx.RIGHT | wx.TOP, 3)
        self.combo_eng = wx.ComboBox(panel2, choices=['1', '2'], style=wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.OnSelectLayerEng, self.combo_eng)
        self.combo_eng.SetToolTip(wx.ToolTip('Layer'))
        self.combo_eng.SetStringSelection(str(self.scriptdata.layer_eng))
        hbox_layer.Add(self.combo_eng, 1, wx.EXPAND)
        vbox211.Add(hbox_layer, 0, wx.EXPAND)

        hbox_w1 = wx.BoxSizer(wx.HORIZONTAL)
        tc1 = wx.StaticText(panel2, label='Wait [s] ')
        self.ntc_ms_waiteng = wx.lib.intctrl.IntCtrl(panel2, size=(70, -1), style=wx.TE_CENTRE, value=0, limited=True)
        self.ntc_ms_waiteng.SetMin(0)
        self.ntc_ms_waiteng.SetMax(9000)
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnSetWaittimeEng, self.ntc_ms_waiteng)

        hbox_w1.Add(tc1, 0, wx.EXPAND | wx.TOP, 7)
        hbox_w1.Add((5, 0))
        hbox_w1.Add(self.ntc_ms_waiteng, 1, wx.EXPAND | wx.TOP, 2)
        vbox211.Add(hbox_w1, 0, wx.EXPAND)

        hbox211.Add(vbox211, 1, wx.EXPAND)

        self.tc_e = wx.TextCtrl(panel2, -1, size=((-1, 80)), style=wx.TE_MULTILINE | wx.TE_READONLY |
                                                                   wx.TE_RICH | wx.VSCROLL)

        hbox211.Add(self.tc_e, 4, wx.EXPAND | wx.LEFT, 20)

        vbox21.Add(hbox211, 1, wx.EXPAND)
        vbox21.Add((0, 10))

        # Angle scan
        hbox212 = wx.BoxSizer(wx.HORIZONTAL)
        vbox212 = wx.BoxSizer(wx.VERTICAL)

        self.cb_ang = wx.CheckBox(panel2, -1, '  Enable Angle Scan  ', size=((142, 20)))
        self.Bind(wx.EVT_CHECKBOX, self.OnEnableTomo, self.cb_ang)
        vbox212.Add(self.cb_ang, 0, wx.EXPAND)

        self.button_ang = wx.Button(panel2, -1, 'Edit Tomo Info')
        self.Bind(wx.EVT_BUTTON, self.OnEditTomoInfo, id=self.button_ang.GetId())
        vbox212.Add(self.button_ang, 0, wx.EXPAND)
        vbox212.Add((0, 3))

        hbox_layer2 = wx.BoxSizer(wx.HORIZONTAL)
        tl2 = wx.StaticText(panel2, label="Layer:")
        hbox_layer2.Add(tl2, 0, wx.RIGHT | wx.TOP, 3)
        self.combo_ang = wx.ComboBox(panel2, choices=['1', '2'], style=wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.OnSelectLayerAng, self.combo_ang)
        self.combo_ang.SetToolTip(wx.ToolTip('Layer'))
        self.combo_ang.SetStringSelection(str(self.scriptdata.layer_ang))
        hbox_layer2.Add(self.combo_ang, 1, wx.EXPAND)
        vbox212.Add(hbox_layer2, 0, wx.EXPAND)

        hbox_w2 = wx.BoxSizer(wx.HORIZONTAL)
        tc1 = wx.StaticText(panel2, label='Wait [s] ')
        self.ntc_ms_waittomo = wx.lib.intctrl.IntCtrl(panel2, size=(70, -1), style=wx.TE_CENTRE, value=0, limited=True)
        self.ntc_ms_waittomo.SetMin(0)
        self.ntc_ms_waittomo.SetMax(9000)
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnSetWaittimeTomo, self.ntc_ms_waittomo)

        hbox_w2.Add(tc1, 0, wx.EXPAND | wx.TOP, 7)
        hbox_w2.Add((5, 0))
        hbox_w2.Add(self.ntc_ms_waittomo, 1, wx.EXPAND | wx.TOP, 2)
        vbox212.Add(hbox_w2, 0, wx.EXPAND)

        self.cb_singletomo = wx.CheckBox(panel2, -1, '  Single Tomo Image', size=((142, 20)))
        self.Bind(wx.EVT_CHECKBOX, self.OnEnableTomo, self.cb_singletomo)
        self.cb_singletomo.Disable()
        vbox212.Add(self.cb_singletomo, 0, wx.EXPAND)

        hbox212.Add(vbox212, 1, wx.EXPAND)

        self.tc_an = wx.TextCtrl(panel2, 1, size=((-1, 80)), style=wx.TE_MULTILINE | wx.TE_READONLY |
                                                                   wx.TE_RICH | wx.VSCROLL)
        hbox212.Add(self.tc_an, 4, wx.EXPAND | wx.LEFT, 20)

        vbox21.Add(hbox212, 1, wx.EXPAND)

        vbox21.Add((0, 10))

        # Mosaic scan
        hbox213 = wx.BoxSizer(wx.HORIZONTAL)
        vbox213 = wx.BoxSizer(wx.VERTICAL)

        self.cb_mos = wx.CheckBox(panel2, -1, '  Enable Mosaic Scan  ', size=((142, 20)))
        self.Bind(wx.EVT_CHECKBOX, self.OnEnableMosaic, self.cb_mos)
        vbox213.Add(self.cb_mos, 0, wx.EXPAND)

        self.button_mos = wx.Button(panel2, -1, 'Edit Mosaic Info')
        self.Bind(wx.EVT_BUTTON, self.OnEditMosaicInfo, id=self.button_mos.GetId())
        vbox213.Add(self.button_mos, 0, wx.EXPAND)
        vbox213.Add((0, 3))

        hbox_layer3 = wx.BoxSizer(wx.HORIZONTAL)
        tl3 = wx.StaticText(panel2, label="Layer:")
        hbox_layer3.Add(tl3, 0, wx.RIGHT | wx.TOP, 3)
        self.combo_mos = wx.ComboBox(panel2, choices=['3'], style=wx.CB_READONLY)
        self.combo_mos.SetToolTip(wx.ToolTip('Layer'))
        self.combo_mos.SetStringSelection(str(self.scriptdata.layer_mos))
        hbox_layer3.Add(self.combo_mos, -1, wx.EXPAND)
        vbox213.Add(hbox_layer3, 0, wx.EXPAND)

        hbox_w3 = wx.BoxSizer(wx.HORIZONTAL)
        tc1 = wx.StaticText(panel2, label='Wait [s] ')
        self.ntc_ms_waitmos = wx.lib.intctrl.IntCtrl(panel2, size=(70, -1), style=wx.TE_CENTRE, value=0, limited=True)
        self.ntc_ms_waitmos.SetMin(0)
        self.ntc_ms_waitmos.SetMax(9000)
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnSetWaittimeMosaic, self.ntc_ms_waitmos)

        hbox_w3.Add(tc1, 0, wx.EXPAND | wx.TOP, 7)
        hbox_w3.Add((5, 0))
        hbox_w3.Add(self.ntc_ms_waitmos, 1, wx.EXPAND | wx.TOP, 2)
        vbox213.Add(hbox_w3, 0, wx.EXPAND)

        hbox213.Add(vbox213, 1, wx.EXPAND)

        self.tc_mos = wx.TextCtrl(panel2, -1, size=((-1, 80)), style=wx.TE_MULTILINE | wx.TE_READONLY |
                                                                     wx.TE_RICH | wx.VSCROLL)
        hbox213.Add(self.tc_mos, 4, wx.EXPAND | wx.LEFT, 20)

        vbox21.Add(hbox213, 1, wx.EXPAND)

        vbox21.Add((0, 10))

        # Multiexposure scan
        hbox214 = wx.BoxSizer(wx.HORIZONTAL)
        vbox214 = wx.BoxSizer(wx.VERTICAL)

        self.cb_me = wx.CheckBox(panel2, -1, '  Enable MultiExposures ', size=((142, 20)))
        self.Bind(wx.EVT_CHECKBOX, self.OnEnableMulti, self.cb_me)
        vbox214.Add(self.cb_me, 0, wx.EXPAND)

        hbox2ne = wx.BoxSizer(wx.HORIZONTAL)
        t_nexp = wx.StaticText(panel2, label="NExposures:")
        self.tc_nexp = wx.lib.intctrl.IntCtrl(panel2, size=(70, -1), style=wx.TE_CENTRE,
                                              value=self.scriptdata.n_exposures,
                                              limited=True, allow_none=True)
        self.tc_nexp.SetMin(1)
        self.tc_nexp.SetMax(10000)
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnNExposures, self.tc_nexp)
        self.tc_nexp.Disable()

        hbox2ne.Add(t_nexp, 0, flag=wx.ALIGN_CENTER_VERTICAL)
        hbox2ne.Add(self.tc_nexp, 0, wx.LEFT, 5)
        vbox214.Add((0, 3))
        vbox214.Add(hbox2ne, 0, wx.EXPAND)
        vbox214.Add((0, 3))

        self.cb_avfly = wx.CheckBox(panel2, -1, '  Average on the fly', size=((142, 20)))
        self.cb_avfly.Disable()
        vbox214.Add(self.cb_avfly, 0, wx.EXPAND)

        hbox_layer4 = wx.BoxSizer(wx.HORIZONTAL)
        tl4 = wx.StaticText(panel2, label="Layer:")
        hbox_layer4.Add(tl4, 0, wx.RIGHT | wx.TOP, 3)
        self.combo_me = wx.ComboBox(panel2, choices=['4'], style=wx.CB_READONLY)
        self.combo_me.SetToolTip(wx.ToolTip('Layer'))
        self.combo_me.SetStringSelection(str(self.scriptdata.layer_me))
        hbox_layer4.Add(self.combo_me, 1, wx.EXPAND)
        vbox214.Add(hbox_layer4, 0, wx.EXPAND)

        hbox_w4 = wx.BoxSizer(wx.HORIZONTAL)
        tc1 = wx.StaticText(panel2, label='Wait [s] ')
        self.ntc_ms_waitmulti = wx.lib.intctrl.IntCtrl(panel2, size=(70, -1), style=wx.TE_CENTRE, value=0, limited=True)
        self.ntc_ms_waitmulti.SetMin(0)
        self.ntc_ms_waitmulti.SetMax(9000)
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnSetWaittimeMulti, self.ntc_ms_waitmulti)

        hbox_w4.Add(tc1, 0, wx.EXPAND | wx.TOP, 7)
        hbox_w4.Add((5, 0))
        hbox_w4.Add(self.ntc_ms_waitmulti, 1, wx.EXPAND | wx.TOP, 2)
        vbox214.Add(hbox_w4, 0, wx.EXPAND)

        hbox214.Add(vbox214, 1, wx.EXPAND)

        self.tc_me = wx.TextCtrl(panel2, -1, size=((-1, 80)), style=wx.TE_MULTILINE | wx.TE_READONLY |
                                                                    wx.TE_RICH | wx.VSCROLL)
        hbox214.Add(self.tc_me, 4, wx.EXPAND | wx.LEFT, 20)

        vbox21.Add(hbox214, 1, wx.EXPAND)
        vbox21.Add((0, 10))

        sizer2.Add(vbox21, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, 10)
        panel2.SetSizer(sizer2)

        # panel 3
        panel3 = wx.Panel(self.scroll, -1)

        sbox3 = wx.StaticBox(panel3, -1, 'Reference Collection')
        sizer3 = wx.StaticBoxSizer(sbox3, orient=wx.VERTICAL)
        vbox31 = wx.BoxSizer(wx.VERTICAL)
        vbox31.Add((0, 2))

        hbox31 = wx.BoxSizer(wx.HORIZONTAL)

        trefx = wx.StaticText(panel3, label="RefX:")
        self.ntc_refx = wx.lib.masked.numctrl.NumCtrl(panel3, -1,
                                                      value=0, integerWidth=5,
                                                      fractionWidth=2,
                                                      allowNegative=True,
                                                      foregroundColour="Black",
                                                      signedForegroundColour="Black",
                                                      limited=False,
                                                      allowNone=True)

        trefy = wx.StaticText(panel3, label="RefY:")
        self.ntc_refy = wx.lib.masked.numctrl.NumCtrl(panel3, -1,
                                                      value=0, integerWidth=5,
                                                      fractionWidth=2,
                                                      allowNegative=True,
                                                      foregroundColour="Black",
                                                      signedForegroundColour="Black",
                                                      limited=False,
                                                      allowNone=True)

        trefz = wx.StaticText(panel3, label="RefZ:")
        self.ntc_refz = wx.lib.masked.numctrl.NumCtrl(panel3, -1,
                                                      value=0, integerWidth=5,
                                                      fractionWidth=2,
                                                      allowNegative=True,
                                                      foregroundColour="Black",
                                                      signedForegroundColour="Black",
                                                      limited=False,
                                                      allowNone=True)
        self.ntc_refz.Disable()

        self.cb_refz = wx.CheckBox(panel3, -1, '')
        self.Bind(wx.EVT_CHECKBOX, self.OnRefZ, self.cb_refz)

        treft = wx.StaticText(panel3, label="RefTheta:")
        self.ntc_reft = wx.lib.masked.numctrl.NumCtrl(panel3, -1,
                                                      value=0, integerWidth=5,
                                                      fractionWidth=2,
                                                      allowNegative=True,
                                                      foregroundColour="Black",
                                                      signedForegroundColour="Black",
                                                      limited=False,
                                                      allowNone=True)
        self.ntc_reft.Disable()

        self.cb_reft = wx.CheckBox(panel3, -1, '')
        self.Bind(wx.EVT_CHECKBOX, self.OnRefTheta, self.cb_reft)

        self.button_refpos = wx.Button(panel3, -1, 'Load Ref Position')
        self.Bind(wx.EVT_BUTTON, self.OnLoadRefPosition, id=self.button_refpos.GetId())

        hbox31.Add(trefx, 0, wx.RIGHT | wx.TOP, 3)
        hbox31.Add(self.ntc_refx, 0, wx.RIGHT, 20)
        hbox31.Add(trefy, 0, wx.RIGHT | wx.TOP, 3)
        hbox31.Add(self.ntc_refy, 0, wx.RIGHT, 20)
        hbox31.Add(trefz, 0, wx.RIGHT | wx.TOP, 3)
        hbox31.Add(self.cb_refz, 0, wx.RIGHT | wx.TOP, 3)
        hbox31.Add(self.ntc_refz, 0, wx.RIGHT, 20)
        hbox31.Add(treft, 0, wx.RIGHT | wx.TOP, 3)
        hbox31.Add(self.cb_reft, 0, wx.RIGHT | wx.TOP, 3)
        hbox31.Add(self.ntc_reft, 0, wx.RIGHT, 20)
        hbox31.Add(self.button_refpos, 0)
        vbox31.Add(hbox31)

        vbox31.Add((0, 5))

        hbox32 = wx.BoxSizer(wx.HORIZONTAL)

        tccol = wx.StaticText(panel3, label="Collect")
        self.ntc_ncol = wx.lib.intctrl.IntCtrl(panel3, size=(70, -1), style=wx.TE_CENTRE, value=0, limited=True)
        self.ntc_ncol.SetMin(0)
        self.ntc_ncol.SetMax(100)
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnCollectNRefImgs, self.ntc_ncol)

        tcrimg = wx.StaticText(panel3, label="Reference Images for every")
        self.ntc_nexp = wx.lib.intctrl.IntCtrl(panel3, size=(70, -1), style=wx.TE_CENTRE, value=0, limited=True)
        self.ntc_nexp.SetMin(0)
        self.ntc_nexp.SetMax(9000)
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnCollectRefImgsEveryNExp, self.ntc_nexp)

        tcrexp = wx.StaticText(panel3, label="exposures of the sample.")

        self.cb_abba = wx.CheckBox(panel3, -1, '  ABBA mode')
        # self.Bind(wx.EVT_CHECKBOX, self.OnEnableTomo, self.cb_ang)

        hbox32.Add(tccol, 0, wx.RIGHT | wx.TOP, 3)
        hbox32.Add(self.ntc_ncol, 0, wx.RIGHT, 10)
        hbox32.Add(tcrimg, 0, wx.RIGHT | wx.TOP, 3)
        hbox32.Add(self.ntc_nexp, 0, wx.RIGHT, 10)
        hbox32.Add(tcrexp, 0, wx.RIGHT | wx.TOP, 3)
        hbox32.Add((10, 0))
        hbox32.Add(self.cb_abba, 0, wx.RIGHT | wx.TOP, 3)
        vbox31.Add(hbox32)
        vbox31.Add((0, 5))

        hbox33 = wx.BoxSizer(wx.HORIZONTAL)

        st33 = wx.StaticText(panel3, label="Reference Exposure Time: ")
        self.ntc_refexptime = wx.lib.masked.numctrl.NumCtrl(panel3, -1,
                                                            value=1, integerWidth=5,
                                                            fractionWidth=2,
                                                            allowNegative=False,
                                                            foregroundColour="Black",
                                                            signedForegroundColour="Black",
                                                            limited=False,
                                                            allowNone=True)

        st34 = wx.StaticText(panel3, label="Reference Binning: ")
        self.combo_RefBinning = wx.ComboBox(panel3, size=(40, -1), choices=['1', '2', '4', '8'], style=wx.CB_READONLY)
        self.combo_RefBinning.SetToolTip(wx.ToolTip('Reference Binning'))
        self.combo_RefBinning.SetStringSelection(str(self.scriptdata.RefBinning))

        hbox33.Add(st33, 0, wx.RIGHT | wx.TOP, 3)
        hbox33.Add(self.ntc_refexptime, 0, wx.RIGHT, 10)
        hbox33.Add(st34, 0, wx.RIGHT | wx.TOP, 3)
        hbox33.Add(self.combo_RefBinning, 0, wx.RIGHT, 10)

        self.cb_refrelpos = wx.CheckBox(panel3, -1, '  Use Relative Positions', size=((142, 20)))
        hbox33.Add(self.cb_refrelpos, 0, wx.RIGHT, 20)

        self.cb_refavfly = wx.CheckBox(panel3, -1, '  Average on the fly', size=((142, 20)))
        self.cb_refavfly.Disable()
        hbox33.Add(self.cb_refavfly, 0, wx.RIGHT, 20)

        vbox31.Add(hbox33)
        vbox31.Add((0, 10))

        sizer3.Add(vbox31, 1, wx.LEFT | wx.RIGHT | wx.EXPAND, 10)
        panel3.SetSizer(sizer3)

        # panel 4
        panel4 = wx.Panel(self.scroll, -1)
        vbox41 = wx.BoxSizer(wx.VERTICAL)

        self.tc_scriptlist = wx.TextCtrl(panel4, -1, style=wx.TE_MULTILINE | wx.TE_READONLY |
                                                           wx.TE_RICH | wx.VSCROLL | wx.HSCROLL)
        self.tc_scriptlist.SetMinSize((250, -1))
        vbox41.Add(self.tc_scriptlist, 5, wx.EXPAND | wx.TOP, 6)

        self.tc_status = wx.TextCtrl(panel4, -1, style=wx.TE_MULTILINE | wx.TE_READONLY |
                                                       wx.TE_RICH | wx.TE_BESTWRAP)
        self.tc_status.SetMinSize((250, -1))
        vbox41.Add(self.tc_status, 1, wx.EXPAND | wx.TOP, 6)

        panel4.SetSizer(vbox41)

        # panel 5
        panel5 = wx.Panel(self.scroll, -1)
        vbox5 = wx.BoxSizer(wx.VERTICAL)

        fgs5 = wx.FlexGridSizer(cols=6, vgap=10, hgap=20)
        fgs5.AddGrowableCol(0)
        fgs5.AddGrowableCol(1)
        fgs5.AddGrowableCol(2)
        fgs5.AddGrowableCol(3)
        fgs5.AddGrowableCol(4)
        fgs5.AddGrowableCol(5)

        self.button_1 = wx.Button(panel5, -1, 'Generate Scan')
        self.Bind(wx.EVT_BUTTON, self.OnGenerateScan, id=self.button_1.GetId())
        fgs5.Add(self.button_1, 0, wx.EXPAND)

        self.button_3 = wx.Button(panel5, -1, 'Load Recipe')
        self.Bind(wx.EVT_BUTTON, self.OnLoadRecipe, id=self.button_3.GetId())
        fgs5.Add(self.button_3, 1, wx.EXPAND)

        self.button_2 = wx.Button(panel5, -1, 'Save Recipe')
        self.Bind(wx.EVT_BUTTON, self.OnSaveRecipe, id=self.button_2.GetId())
        fgs5.Add(self.button_2, 2, wx.EXPAND)

        self.button_4 = wx.Button(panel5, -1, 'Load Script')
        self.Bind(wx.EVT_BUTTON, self.OnLoadScript, id=self.button_4.GetId())
        fgs5.Add(self.button_4, 3, wx.EXPAND)

        self.button_11 = wx.Button(panel5, -1, 'Save Script')
        self.button_11.Disable()
        self.Bind(wx.EVT_BUTTON, self.OnSaveScript, id=self.button_11.GetId())
        fgs5.Add(self.button_11, 4, wx.EXPAND)

        self.button_5 = wx.Button(panel5, -1, 'Add to MultiScan')
        self.button_5.Disable()
        self.Bind(wx.EVT_BUTTON, self.OnAddMultiScan, id=self.button_5.GetId())
        fgs5.Add(self.button_5, 5, wx.EXPAND)

        self.button_6 = wx.Button(panel5, -1, 'Execute Scan')
        self.button_6.Disable()
        self.Bind(wx.EVT_BUTTON, self.OnExecuteScan, id=self.button_6.GetId())
        fgs5.Add(self.button_6, 0, wx.EXPAND)

        self.button_7 = wx.Button(panel5, -1, 'Pause Scan')
        self.button_7.Disable()
        self.Bind(wx.EVT_BUTTON, self.OnPauseScan, id=self.button_7.GetId())
        fgs5.Add(self.button_7, 1, wx.EXPAND)

        self.button_8 = wx.Button(panel5, -1, 'Resume Paused')
        self.button_8.Disable()
        self.Bind(wx.EVT_BUTTON, self.OnResumeScan, id=self.button_8.GetId())
        fgs5.Add(self.button_8, 2, wx.EXPAND)

        self.button_12 = wx.Button(panel5, -1, 'Stop Scan')
        self.button_12.Disable()
        self.Bind(wx.EVT_BUTTON, self.OnStopScan, id=self.button_12.GetId())
        fgs5.Add(self.button_12, 3, wx.EXPAND)

        self.button_14 = wx.Button(panel5, -1, 'Resume Uncompleted')
        # self.button_14.Disable()
        self.Bind(wx.EVT_BUTTON, self.OnResumeUncompletedScan, id=self.button_14.GetId())
        fgs5.Add(self.button_14, 4, wx.EXPAND)

        self.button_13 = wx.ToggleButton(panel5, -1, 'Agitator OFF on Exit')
        self.button_13.SetValue(True)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnTAgitator, id=self.button_13.GetId())
        fgs5.Add(self.button_13, 5, wx.EXPAND)

        vbox5.Add(fgs5, 1, wx.EXPAND)

        panel5.SetSizer(vbox5)

        # panel 6
        panel6 = wx.Panel(self.scroll, -1)
        hbox60 = wx.BoxSizer(wx.HORIZONTAL)

        sbox6 = wx.StaticBox(panel6, -1, 'Scan Settings')
        sizer6 = wx.StaticBoxSizer(sbox6, orient=wx.VERTICAL)
        vbox61 = wx.BoxSizer(wx.VERTICAL)
        vbox61.Add((0, 5))

        hbox61 = wx.BoxSizer(wx.HORIZONTAL)

        st61 = wx.StaticText(panel6, label="Sample Name: ")

        self.tc_samplename = wx.TextCtrl(panel6, -1, size=((200, -1)), style=wx.TE_RICH | wx.VSCROLL,
                                         value=self.scriptdata.SampleName)
        hbox61.Add(st61, 0, flag=wx.ALIGN_CENTER_VERTICAL)
        hbox61.Add(self.tc_samplename, 1, wx.EXPAND | wx.LEFT, 2)

        hbox62 = wx.BoxSizer(wx.HORIZONTAL)

        st62 = wx.StaticText(panel6, label="Exposure Time: ")
        self.ntc_exptime = wx.lib.masked.numctrl.NumCtrl(panel6, -1,
                                                         value=1.0, integerWidth=3,
                                                         fractionWidth=2,
                                                         allowNegative=False,
                                                         foregroundColour="Black",
                                                         signedForegroundColour="Black",
                                                         limited=False,
                                                         allowNone=True)

        self.cb_ioncham = wx.CheckBox(panel6, -1, '  Use Ion Chamber', size=((142, 20)))

        st64 = wx.StaticText(panel6, label="Binning: ")
        self.combo_Binning = wx.ComboBox(panel6, size=(60, -1), choices=['1', '2', '4', '8'], style=wx.CB_READONLY)
        self.combo_Binning.SetToolTip(wx.ToolTip('Binning'))
        self.combo_Binning.SetStringSelection(str(self.scriptdata.Binning))

        st65 = wx.StaticText(panel6, label="Pixel Size [Binning 1]: ")
        self.ntc_pixelsize = wx.lib.masked.numctrl.NumCtrl(panel6, -1,
                                                           value=self.scriptdata.pixelsize, integerWidth=2,
                                                           fractionWidth=4,
                                                           allowNegative=False,
                                                           foregroundColour="Black",
                                                           signedForegroundColour="Black",
                                                           limited=False,
                                                           allowNone=True)

        hbox62.Add(st64, 0, wx.TOP, 3)
        hbox62.Add(self.combo_Binning, 0, wx.RIGHT, 15)
        hbox62.Add(st65, 0, wx.TOP, 3)
        hbox62.Add(self.ntc_pixelsize, 0)

        hbox62a = wx.BoxSizer(wx.HORIZONTAL)
        hbox62a.Add(st62, 0, wx.TOP, 3)
        hbox62a.Add(self.ntc_exptime, 0, wx.RIGHT, 8)
        hbox62a.Add(self.cb_ioncham, 0, wx.LEFT, 20)

        hbox63 = wx.BoxSizer(wx.HORIZONTAL)

        st63 = wx.StaticText(panel6, label="Path: ")

        self.tc_outpath = wx.TextCtrl(panel6, -1, style=wx.TE_RICH | wx.VSCROLL | wx.TE_READONLY,
                                      value=self.scriptdata.OutputPath)
        hbox63.Add(st63, 0, flag=wx.ALIGN_CENTER_VERTICAL)
        hbox63.Add(self.tc_outpath, 1, wx.EXPAND | wx.LEFT, 2)

        hbox64 = wx.BoxSizer(wx.HORIZONTAL)
        st64 = wx.StaticText(panel6, label="Reting Frame Adding - N frames:")
        self.ntc_n_retiga = wx.lib.intctrl.IntCtrl(panel6, size=(70, -1), style=wx.TE_CENTRE,
                                                   value=self.scriptdata.n_retiga,
                                                   limited=True, allow_none=True)
        self.ntc_n_retiga.SetMin(1)
        self.ntc_n_retiga.SetMax(1000)
        # self.Bind(wx.lib.intctrl.EVT_INT, self.OnNExposures, self.ntc_n_retiga)
        if not self.HaveRetiga: self.ntc_n_retiga.Disable()

        hbox64.Add(st64, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 8)
        hbox64.Add(self.ntc_n_retiga, 0, wx.LEFT, 0)

        vbox61.Add(hbox61, 0, wx.EXPAND)
        vbox61.Add((0, 5))
        vbox61.Add(hbox63, 0, wx.EXPAND)
        vbox61.Add((0, 5))
        vbox61.Add(hbox62a)
        vbox61.Add((0, 5))
        vbox61.Add(hbox62)
        vbox61.Add((0, 5))
        vbox61.Add(hbox64)

        sizer6.Add(vbox61, 1, wx.LEFT | wx.RIGHT | wx.EXPAND, 10)

        spanel6 = wx.StaticBox(panel6, -1, 'MultiScan Settings')
        sizer62 = wx.StaticBoxSizer(spanel6, orient=wx.VERTICAL)

        hbox6t = wx.BoxSizer(wx.HORIZONTAL)

        self.tc_multiscan = wx.ListCtrl(panel6, -1, size=(120, -1), style=wx.LC_REPORT | wx.LC_NO_HEADER)
        self.tc_multiscan.InsertColumn(0, '')
        self.tc_multiscan.SetColumnWidth(0, 155)

        self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnMultiScanListClick, self.tc_multiscan)

        hbox6t.Add(self.tc_multiscan, 1, wx.EXPAND)

        vbox62 = wx.BoxSizer(wx.VERTICAL)

        self.cb_multiscript = wx.CheckBox(panel6, -1, '  Show MultiScan List')
        self.Bind(wx.EVT_CHECKBOX, self.OnShowMultiScanList, self.cb_multiscript)

        vbox62.Add(self.cb_multiscript, 0, wx.RIGHT | wx.TOP, 3)

        hbox66 = wx.BoxSizer(wx.HORIZONTAL)

        self.button_9 = wx.Button(panel6, -1, 'Edit', size=(60, -1))
        self.button_9.Disable()
        self.Bind(wx.EVT_BUTTON, self.OnEditMultiscan, id=self.button_9.GetId())
        hbox66.Add(self.button_9, 0)

        self.button_14 = wx.Button(panel6, -1, 'Delete', size=(60, -1))
        self.button_14.Disable()
        self.Bind(wx.EVT_BUTTON, self.OnDeleteMultiscan, id=self.button_14.GetId())
        hbox66.Add(self.button_14, 0)

        self.button_10 = wx.Button(panel6, -1, 'Save', size=(60, -1))
        self.button_10.Disable()
        self.Bind(wx.EVT_BUTTON, self.OnSaveEditMultiscan, id=self.button_10.GetId())
        hbox66.Add(self.button_10, 0)

        hbox69 = wx.BoxSizer(wx.HORIZONTAL)
        self.button_15 = wx.Button(panel6, -1, 'Move Up', size=(75, -1))
        self.button_15.Disable()
        self.Bind(wx.EVT_BUTTON, self.OnMoveUpMultiscan, id=self.button_15.GetId())
        hbox69.Add(self.button_15, 0)

        self.button_16 = wx.Button(panel6, -1, 'Move Down', size=(75, -1))
        self.button_16.Disable()
        self.Bind(wx.EVT_BUTTON, self.OnMoveDownMultiscan, id=self.button_16.GetId())
        hbox69.Add(self.button_16, 0)

        hbox67 = wx.BoxSizer(wx.HORIZONTAL)

        self.cb_ms_repeat = wx.CheckBox(panel6, -1, ' Repeat')
        self.Bind(wx.EVT_CHECKBOX, self.OnCBRepeatMultiScan, self.cb_ms_repeat)

        self.ntc_ms_repeat = wx.lib.intctrl.IntCtrl(panel6, size=(70, -1), style=wx.TE_CENTRE, value=1,
                                                    limited=True, allow_none=True)

        self.ntc_ms_repeat.SetMin(0)
        self.ntc_ms_repeat.SetMax(1000)
        self.ntc_ms_repeat.Disable()
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnRepeatMultiscanSettings, self.ntc_ms_repeat)
        tc13 = wx.StaticText(panel6, label='times')

        hbox67.Add(self.cb_ms_repeat, 0, wx.EXPAND)
        hbox67.Add((5, 0))
        hbox67.Add(self.ntc_ms_repeat, 0, wx.TOP, 2)
        hbox67.Add((5, 0))
        hbox67.Add(tc13, 0, wx.TOP, 5)

        hbox68 = wx.BoxSizer(wx.HORIZONTAL)

        tc12 = wx.StaticText(panel6, label='Wait [s] ')
        self.ntc_ms_wait = wx.lib.intctrl.IntCtrl(panel6, size=(70, -1), style=wx.TE_CENTRE, value=0, limited=True)
        self.ntc_ms_wait.SetMin(0)
        self.ntc_ms_wait.SetMax(9000)
        self.ntc_ms_wait.Disable()
        self.Bind(wx.lib.intctrl.EVT_INT, self.OnRepeatMultiscanSettings, self.ntc_ms_wait)

        hbox68.Add((18, 0))
        hbox68.Add(tc12, 0, wx.EXPAND | wx.TOP, 7)
        hbox68.Add((5, 0))
        hbox68.Add(self.ntc_ms_wait, 0, wx.TOP, 2)

        vbox62.Add((0, 2))
        vbox62.Add(hbox66)
        vbox62.Add(hbox69)
        vbox62.Add((0, 5))
        vbox62.Add(hbox67)
        vbox62.Add((0, 1))
        vbox62.Add(hbox68)

        hbox6t.Add((10, 0))
        hbox6t.Add(vbox62)

        sizer62.Add(hbox6t, 1, wx.EXPAND, 10)

        hbox60.Add(sizer6, 1, wx.EXPAND)
        hbox60.Add((8, 0))
        hbox60.Add(sizer62, 1, wx.EXPAND)

        panel6.SetSizer(hbox60)

        self.Bind(wx.EVT_CLOSE, self.OnClose)

        randomId = wx.NewId()
        self.Bind(wx.EVT_MENU, self.OnClose, id=randomId)
        accel_tbl = wx.AcceleratorTable([(wx.ACCEL_CTRL, ord('Q'), randomId)])
        self.SetAcceleratorTable(accel_tbl)

        # ----------------------------

        hboxT = wx.BoxSizer(wx.HORIZONTAL)

        vboxT2 = wx.BoxSizer(wx.VERTICAL)

        vboxT2.Add(panel2, 1, wx.EXPAND | wx.RIGHT, 5)
        vboxT2.Add((0, 5))
        vboxT2.Add(panel3, 0, wx.EXPAND | wx.RIGHT, 5)
        vboxT2.Add((0, 5))
        vboxT2.Add(panel6, 0, wx.EXPAND | wx.RIGHT, 5)

        hboxT.Add(vboxT2, 3, wx.EXPAND)
        hboxT.Add(panel4, 1, wx.EXPAND)

        vboxT.Add((0, 10))
        vboxT.Add(panel1, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 18)
        vboxT.Add((0, 2))
        vboxT.Add(hboxT, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, 18)
        vboxT.Add((0, 18))
        vboxT.Add(panel5, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 18)
        vboxT.Add((0, 10))

        self.scroll.SetSizer(vboxT)

        self.Centre()
        self.Show(True)

        # self.Maximize()

    # ----------------------------------------------------------------------
    def OnEditEngInfo(self, event):

        if self.epframe:
            self.epframe.Raise()
        else:
            self.epframe = EnergyFrame(self.scriptdata, self.microscope, MotorizedDetector=self.MotorizedDetector)
            self.epframe.Show()

    # ----------------------------------------------------------------------
    def OnEditTomoInfo(self, event):

        AngleFrame(self.scriptdata, self.microscope, self.SampleXYZAboveRotation).Show()

    # ----------------------------------------------------------------------
    def OnEditMosaicInfo(self, event):

        MosaicFrame(self.scriptdata).Show()

    # ----------------------------------------------------------------------
    def OnNExposures(self, event):
        ctl = event.GetEventObject()
        value = ctl.GetValue()
        if value == None:
            value = 1
        self.scriptdata.n_exposures = value

        self.ShowEnabledTC()

    # ----------------------------------------------------------------------
    def OnRepeatScanN(self, event):
        ctl = event.GetEventObject()
        value = ctl.GetValue()
        if value == None:
            value = 1
        self.scriptdata.NRepeatScan = value
        if v: print
        'self.scriptdata.NRepeatScan=', self.scriptdata.NRepeatScan

    # ----------------------------------------------------------------------
    def OnCBEnableDisableExclusiveControl(self, event):
        if self.cb_exclusiveControl.GetValue():
            self.microscope.GetExclusiveControl()
        else:
            self.microscope.GiveUpExclusiveControl()

        # ----------------------------------------------------------------------

    def OnEnableEnergy(self, event):
        if self.cb_eng.GetValue():
            self.scriptdata.enable_energy = 1
        else:
            self.scriptdata.enable_energy = 0

        self.ShowEnabledTC()

    # ----------------------------------------------------------------------
    def OnEnableTomo(self, event):
        if self.cb_ang.GetValue():
            self.scriptdata.enable_tomo = 1
        else:
            self.scriptdata.enable_tomo = 0

        self.ShowEnabledTC()

    # ----------------------------------------------------------------------
    def OnEnableMosaic(self, event):
        if self.cb_mos.GetValue():
            self.scriptdata.enable_mosaic = 1
        else:
            self.scriptdata.enable_mosaic = 0

        self.ShowEnabledTC()

    # ----------------------------------------------------------------------
    def OnEnableMulti(self, event):
        if self.cb_me.GetValue():
            self.scriptdata.enable_multi = 1
            self.tc_nexp.Enable()
        else:
            self.scriptdata.enable_multi = 0
            self.scriptdata.n_exposures = 1
            self.tc_nexp.Disable()

        self.tc_nexp.SetValue(self.scriptdata.n_exposures)

        self.ShowEnabledTC()

    # ----------------------------------------------------------------------
    def OnSetWaittimeEng(self, event):

        self.scriptdata.wait_energy = self.ntc_ms_waiteng.GetValue()

    # ----------------------------------------------------------------------
    def OnSetWaittimeTomo(self, event):
        self.scriptdata.wait_tomo = self.ntc_ms_waittomo.GetValue()

    # ----------------------------------------------------------------------
    def OnSetWaittimeMosaic(self, event):
        self.scriptdata.wait_mosaic = self.ntc_ms_waitmos.GetValue()

    # ----------------------------------------------------------------------
    def OnSetWaittimeMulti(self, event):
        self.scriptdata.wait_multi = self.ntc_ms_waitmulti.GetValue()

    # ----------------------------------------------------------------------
    def OnSelectLayerEng(self, event):
        item = event.GetSelection()
        self.scriptdata.layer_eng = int(item) + 1

        if self.scriptdata.layer_eng == 1:
            self.scriptdata.layer_ang = 2
        else:
            self.scriptdata.layer_ang = 1

        self.combo_eng.SetStringSelection(str(self.scriptdata.layer_eng))
        self.combo_ang.SetStringSelection(str(self.scriptdata.layer_ang))

        self.CheckTomoInerLoop()

    # ----------------------------------------------------------------------
    def OnSelectLayerAng(self, event):
        item = event.GetSelection()
        self.scriptdata.layer_ang = int(item) + 1

        if self.scriptdata.layer_ang == 1:
            self.scriptdata.layer_eng = 2
        else:
            self.scriptdata.layer_eng = 1

        self.combo_eng.SetStringSelection(str(self.scriptdata.layer_eng))
        self.combo_ang.SetStringSelection(str(self.scriptdata.layer_ang))

        self.CheckTomoInerLoop()

    # ----------------------------------------------------------------------
    def OnShowMultiScanList(self, event):
        self.ShowScriptList()

    # ----------------------------------------------------------------------
    def CheckTomoInerLoop(self):

        if (self.scriptdata.enable_tomo == 1) and (self.scriptdata.layer_ang == 2) and (
                self.scriptdata.enable_mosaic == 0) and (self.SampleXYZAboveRotation == 1):
            self.cb_singletomo.Enable()
        else:
            self.cb_singletomo.Disable()
            self.cb_singletomo.SetValue(False)

        if self.cb_singletomo.IsChecked():
            self.scriptdata.singletomo = 1
        else:
            self.scriptdata.singletomo = 0

    # ----------------------------------------------------------------------
    def OnRefZ(self, event):
        if self.cb_refz.GetValue():
            self.scriptdata.useRefZ = True
            self.ntc_refz.Enable()

        else:
            self.scriptdata.useRefZ = False
            self.scriptdata.RefZ = 0.0
            self.ntc_refz.SetValue(self.scriptdata.RefZ)
            self.ntc_refz.Disable()

    # ----------------------------------------------------------------------
    def OnRefTheta(self, event):
        if self.cb_reft.GetValue():
            self.scriptdata.useRefTheta = True
            self.ntc_reft.Enable()

        else:
            self.scriptdata.useRefTheta = False
            self.scriptdata.RefTheta = 0.0
            self.ntc_reft.SetValue(self.scriptdata.RefTheta)
            self.ntc_reft.Disable()

        # ----------------------------------------------------------------------

    def OnLoadRefPosition(self, event):

        wildcard = "XRM files (*.xrm)|*.xrm"
        dialog = wx.FileDialog(None, "Please select a reference xrm file",
                               wildcard=wildcard,
                               style=wx.OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            filepath = dialog.GetPath()
            filename = dialog.GetFilename()
            filedir = dialog.GetDirectory()
        else:
            return

        if not self.microscope:
            xrminfo = self.microscope.get_xrm_info(filepath)

            samplets = xrminfo['MotPosT']
            samplexs = xrminfo['MotPosX']
            sampleys = xrminfo['MotPosY']
            samplezs = xrminfo['MotPosZ']


        else:

            Ctxrm = txrm.txrm()

            xrminfo = Ctxrm.get_xrm_info(filepath)

            NumofMotors, xxx, yyy, zzz, ttt, zpxxx, zpyyy, zpzzz = Ctxrm.get_motor_index(xrminfo)

            samplets = xrminfo['MotPos'][ttt]
            samplexs = xrminfo['MotPos'][xxx]
            sampleys = xrminfo['MotPos'][yyy]
            samplezs = xrminfo['MotPos'][zzz]

        self.cb_refz.SetValue(True)
        self.scriptdata.useRefZ = True
        self.ntc_refz.Enable()

        self.cb_reft.SetValue(True)
        self.scriptdata.useRefTheta = True
        self.ntc_reft.Enable()

        self.ntc_refx.SetValue(samplexs)
        self.ntc_refy.SetValue(sampleys)
        self.ntc_refz.SetValue(samplezs)
        self.ntc_reft.SetValue(samplets)

    # ----------------------------------------------------------------------
    def ShowEnabledTC(self):

        self.CheckTomoInerLoop()

        self.tc_e.Clear()
        self.tc_an.Clear()
        self.tc_mos.Clear()
        self.tc_me.Clear()

        if self.scriptdata.enable_energy:
            self.tc_e.SetValue('Energy Regions:')
            for i in range(len(self.scriptdata.EngRegionlist)):
                self.tc_e.AppendText(
                    '\nRegion:{0:02d} Start:{1:7.2f} Steps:{2:3d} Energy-step:{3:7.2f} Stop:{4:7.2f} ZPx:{5:7.2f} ZPy:{6:7.2f} ZPz:{7:7.2f}'.format(
                        i + 1,
                        self.scriptdata.EngRegionlist[i]['Start'],
                        int(self.scriptdata.EngRegionlist[i]['Steps']),
                        self.scriptdata.EngRegionlist[i]['EngStep'],
                        self.scriptdata.EngRegionlist[i]['Stop'],
                        self.scriptdata.EngRegionlist[i]['ZPxStart'],
                        self.scriptdata.EngRegionlist[i]['ZPyStart'],
                        self.scriptdata.EngRegionlist[i]['ZPzStart']))
                self.tc_e.SetInsertionPoint(0)

            self.cb_eng.SetValue(True)
        else:
            self.cb_eng.SetValue(False)

        if self.scriptdata.enable_tomo:
            self.tc_an.SetValue('Angle Scan Info: \n')
            self.tc_an.AppendText(
                'Xcorr:{0:4.3f} \tZCorr:{1:4.3f} \nSampleXinZeroDegree:{2:4.3f}  SampleYinZeroDegree:{3:4.3f}  SampleZinZeroDegree:{4:4.3f} \nFOV:{5:4.3f}'.format(
                    self.scriptdata.XCorr,
                    self.scriptdata.ZCorr,
                    self.scriptdata.SampleXinZeroDegree,
                    self.scriptdata.SampleYinZeroDegree,
                    self.scriptdata.SampleZinZeroDegree,
                    self.scriptdata.FOV))
            if len(self.scriptdata.samplets) > 0:
                self.tc_an.AppendText('\nNumber of Angles: {0:02d}'.format(len(self.scriptdata.samplets)))

            if (self.scriptdata.singletomo == 1):
                self.tc_an.AppendText('\nAcquire Single Tomo Image (.txrm)')

            self.cb_ang.SetValue(True)
        else:
            self.cb_ang.SetValue(False)

        if self.scriptdata.enable_mosaic:
            self.tc_mos.SetValue('Mosaic Info: \n')
            self.tc_mos.AppendText(
                'Up:{0:3d}  \tDown:{1:3d}  \nLeft:{2:3d}  \tRight:{3:3d}  \nOverlap:{4:7.2f}  \nCentralTile:{5:3d}'.format(
                    self.scriptdata.mosaic_up,
                    self.scriptdata.mosaic_down,
                    self.scriptdata.mosaic_left,
                    self.scriptdata.mosaic_right,
                    self.scriptdata.mosaic_overlap,
                    self.scriptdata.mosaic_centraltile))

            self.cb_mos.SetValue(True)
        else:
            self.cb_mos.SetValue(False)

        if self.scriptdata.enable_multi:
            self.tc_me.SetValue('MultiExposure Enabled:\n')
            self.tc_me.AppendText('Number of multiple exposures: {0:3d} '.format(self.scriptdata.n_exposures))
            self.tc_nexp.Enable()
            self.cb_me.SetValue(True)
        else:
            self.cb_me.SetValue(False)
            self.tc_nexp.Disable()

        if self.scriptdata.NRepeatScan > 1:
            self.cb_1.SetValue(True)
            self.ntc1.SetValue(self.scriptdata.NRepeatScan)
            self.ntc12.SetValue(self.scriptdata.waitsecs)
        else:
            self.cb_1.SetValue(False)
            self.ntc1.SetValue(self.scriptdata.NRepeatScan)
            self.ntc12.SetValue(self.scriptdata.waitsecs)

        self.ntc_ms_waiteng.SetValue(self.scriptdata.wait_energy)
        self.ntc_ms_waittomo.SetValue(self.scriptdata.wait_tomo)
        self.ntc_ms_waitmos.SetValue(self.scriptdata.wait_mosaic)
        self.ntc_ms_waitmulti.SetValue(self.scriptdata.wait_multi)

    # ----------------------------------------------------------------------
    def OnCollectNRefImgs(self, event):
        ctl = event.GetEventObject()
        self.scriptdata.RepeatRefExposure = ctl.GetValue()

    # ----------------------------------------------------------------------
    def OnCollectRefImgsEveryNExp(self, event):
        ctl = event.GetEventObject()
        self.scriptdata.ref4everycertainnumofexp = ctl.GetValue()

    # ----------------------------------------------------------------------
    def OnCBRepeatScan(self, event):

        if self.cb_1.GetValue():
            self.ntc1.Enable()
            self.ntc12.Enable()
        else:
            self.scriptdata.NRepeatScan = 1
            self.scriptdata.waitsecs = 0
            self.ntc1.SetValue(1)
            self.ntc1.Disable()
            self.ntc12.SetValue(0)
            self.ntc12.Disable()

    # ----------------------------------------------------------------------
    def OnCBRepeatMultiScan(self, event):

        if self.cb_ms_repeat.GetValue():
            self.ntc_ms_repeat.Enable()
            self.ntc_ms_wait.Enable()
            self.multiscan_userepeat = 1
        else:
            self.multiscan_userepeat = 0
            self.multiscan_repeat = 1
            self.multiscan_wait = 0
            self.ntc_ms_repeat.SetValue(1)
            self.ntc_ms_repeat.Disable()
            self.ntc_ms_wait.SetValue(0)
            self.ntc_ms_wait.Disable()

    # ----------------------------------------------------------------------
    def OnRepeatMultiscanSettings(self, event):

        self.ShowScriptList()

    # ----------------------------------------------------------------------
    def OnMultiScanListClick(self, event):

        self.ieditmulti = event.m_itemIndex

        self.ShowMultiscanList()

        lineposition = 0
        for i in range(self.ieditmulti):
            lineposition = lineposition + len(self.multiscriptlist[i])

        if self.cb_multiscript.GetValue():
            p1 = self.tc_scriptlist.XYToPosition(0, lineposition)
            self.tc_scriptlist.ShowPosition(p1)

        self.button_10.Disable()

        self.scriptlist = []
        self.scan.scriptlist = []
        self.scriptdata = ScriptData()

        text = self.tc_status.GetValue()

        self.ResetWidgets()

        self.tc_status.SetValue(text)

    # ----------------------------------------------------------------------
    def OnEditMultiscan(self, event):

        self.scriptdata = copy.deepcopy(self.multiscript[self.ieditmulti])

        # self.ShowEnabledTC()
        self.UpdateWidgets()
        self.tc_nexp.SetValue(self.scriptdata.n_exposures)

        self.button_10.Enable()

    # ----------------------------------------------------------------------
    def OnDeleteMultiscan(self, event):

        del self.multiscript[self.ieditmulti]
        del self.multiscriptlist[self.ieditmulti]
        del self.multiscan_namelist[self.ieditmulti]

        for i in range(len(self.multiscriptlist)):
            newmsname = ';;; Multi Scan %d ;;;' % (i + 1)
            self.multiscriptlist[i][0] = newmsname

        # New scan numbering after delete
        for i in range(len(self.multiscan_namelist)):
            items = self.multiscan_namelist[i].split(' ')
            scanname = ' '.join(items[1:])
            self.multiscan_namelist[i] = '{0} {1}'.format(i + 1, scanname)

        self.ieditmulti = 0

        self.ShowMultiscanList()

        self.scriptlist = []
        self.scan.scriptlist = []
        self.scriptdata = ScriptData()
        self.txrmloaded = 0

        self.ResetWidgets()
        self.ShowScriptList()

        if len(self.multiscriptlist) == 0:
            self.button_9.Disable()
            self.button_14.Disable()
            self.button_15.Disable()
            self.button_16.Disable()

        self.UpdateButtonWidgets()

    # ----------------------------------------------------------------------
    def OnSaveEditMultiscan(self, event):

        self.scriptdata.RefX = self.ntc_refx.GetValue()
        self.scriptdata.RefY = self.ntc_refy.GetValue()
        self.scriptdata.RefZ = self.ntc_refz.GetValue()
        self.scriptdata.RefTheta = self.ntc_reft.GetValue()
        if self.cb_abba.GetValue():
            self.scriptdata.ABBA = 1
        else:
            self.scriptdata.ABBA = 0

        if self.cb_avfly.GetValue():
            self.scriptdata.averageotfly = 1
        else:
            self.scriptdata.averageotfly = 0

        if self.cb_refavfly.GetValue():
            self.scriptdata.refaverageotfly = 1
        else:
            self.scriptdata.refaverageotfly = 0

        if self.cb_refrelpos.GetValue():
            self.scriptdata.refrelativepos = 1
        else:
            self.scriptdata.refrelativepos = 0

        self.scriptdata.useRefZ = self.cb_refz.GetValue()
        self.scriptdata.useRefTheta = self.cb_reft.GetValue()

        self.CheckTomoInerLoop()

        self.scriptdata.wait_energy = self.ntc_ms_waiteng.GetValue()
        self.scriptdata.wait_tomo = self.ntc_ms_waittomo.GetValue()
        self.scriptdata.wait_mosaic = self.ntc_ms_waitmos.GetValue()
        self.scriptdata.wait_multi = self.ntc_ms_waitmulti.GetValue()

        # Scan Settings:
        self.scriptdata.Exptime = self.ntc_exptime.GetValue()
        self.scriptdata.Binning = 2 ** (int(self.combo_Binning.GetSelection()))
        self.scriptdata.SampleName = self.tc_samplename.GetValue()
        self.scriptdata.SampleName = self.scriptdata.SampleName.replace(" ", "")
        self.scriptdata.RefBinning = 2 ** (int(self.combo_RefBinning.GetSelection()))
        self.scriptdata.RefExpTime = self.ntc_refexptime.GetValue()
        self.scriptdata.waitsecs = self.ntc12.GetValue()

        self.multiscript[self.ieditmulti] = copy.deepcopy(self.scriptdata)

        multiscriptnum = self.ieditmulti + 1
        self.scan.generate_script(self.scriptdata, self.scriptdata.enable_energy, self.scriptdata.enable_tomo,
                                  self.scriptdata.enable_mosaic, self.scriptdata.enable_multi,
                                  self.scriptdata.layer_eng, self.scriptdata.layer_ang, multiscriptnum=multiscriptnum,
                                  UseRefZ=self.scriptdata.useRefZ, UseRefTheta=self.scriptdata.useRefTheta,
                                  SampleXYZAboveRotation=self.SampleXYZAboveRotation,
                                  TakeSingleTomo=self.scriptdata.singletomo, MotorizedDetector=self.MotorizedDetector)

        newmultilist = []
        newmultilist.append(';;; Multi Scan %d ;;;' % (multiscriptnum))
        for i in range(len(self.scan.scriptlist)):
            newmultilist.append(copy.deepcopy(self.scan.scriptlist[i]))

        self.multiscriptlist[self.ieditmulti] = copy.deepcopy(newmultilist)

        scanname = self.GetScanName(self.scriptdata.enable_energy, self.scriptdata.enable_tomo,
                                    self.scriptdata.enable_mosaic, self.scriptdata.layer_eng, self.scriptdata.layer_ang)

        self.multiscan_namelist[self.ieditmulti] = '{0} {1}'.format(multiscriptnum, scanname)
        self.ShowMultiscanList()

        self.scriptlist = []
        self.scan.scriptlist = []
        self.scriptdata = ScriptData()
        self.txrmloaded = 0

        self.ResetWidgets()
        self.ShowScriptList()

        self.button_10.Disable()

    # ----------------------------------------------------------------------
    def OnMoveUpMultiscan(self, event):

        if self.ieditmulti == 0:
            return

        tempms1 = copy.deepcopy(self.multiscript[self.ieditmulti - 1])
        tempms2 = copy.deepcopy(self.multiscriptlist[self.ieditmulti - 1])
        tempms3 = copy.deepcopy(self.multiscan_namelist[self.ieditmulti - 1])

        self.multiscript[self.ieditmulti - 1] = copy.deepcopy(self.multiscript[self.ieditmulti])
        self.multiscriptlist[self.ieditmulti - 1] = copy.deepcopy(self.multiscriptlist[self.ieditmulti])
        self.multiscan_namelist[self.ieditmulti - 1] = copy.deepcopy(self.multiscan_namelist[self.ieditmulti])

        self.multiscript[self.ieditmulti] = copy.deepcopy(tempms1)
        self.multiscriptlist[self.ieditmulti] = copy.deepcopy(tempms2)
        self.multiscan_namelist[self.ieditmulti] = copy.deepcopy(tempms3)

        for i in range(len(self.multiscriptlist)):
            newmsname = ';;; Multi Scan %d ;;;' % (i + 1)
            self.multiscriptlist[i][0] = newmsname

        # New scan numbering after delete
        for i in range(len(self.multiscan_namelist)):
            items = self.multiscan_namelist[i].split(' ')
            scanname = ' '.join(items[1:])
            self.multiscan_namelist[i] = '{0} {1}'.format(i + 1, scanname)

        self.ieditmulti = self.ieditmulti - 1
        self.ShowMultiscanList()

        self.scriptlist = []
        self.scan.scriptlist = []
        self.scriptdata = ScriptData()
        self.txrmloaded = 0

        self.ResetWidgets()
        self.ShowScriptList()

    # ----------------------------------------------------------------------
    def OnMoveDownMultiscan(self, event):

        if self.ieditmulti == len(self.multiscript) - 1:
            return

        tempms1 = copy.deepcopy(self.multiscript[self.ieditmulti + 1])
        tempms2 = copy.deepcopy(self.multiscriptlist[self.ieditmulti + 1])
        tempms3 = copy.deepcopy(self.multiscan_namelist[self.ieditmulti + 1])

        self.multiscript[self.ieditmulti + 1] = copy.deepcopy(self.multiscript[self.ieditmulti])
        self.multiscriptlist[self.ieditmulti + 1] = copy.deepcopy(self.multiscriptlist[self.ieditmulti])
        self.multiscan_namelist[self.ieditmulti + 1] = copy.deepcopy(self.multiscan_namelist[self.ieditmulti])

        self.multiscript[self.ieditmulti] = copy.deepcopy(tempms1)
        self.multiscriptlist[self.ieditmulti] = copy.deepcopy(tempms2)
        self.multiscan_namelist[self.ieditmulti] = copy.deepcopy(tempms3)

        for i in range(len(self.multiscriptlist)):
            newmsname = ';;; Multi Scan %d ;;;' % (i + 1)
            self.multiscriptlist[i][0] = newmsname

        # New scan numbering after delete
        for i in range(len(self.multiscan_namelist)):
            items = self.multiscan_namelist[i].split(' ')
            scanname = ' '.join(items[1:])
            self.multiscan_namelist[i] = '{0} {1}'.format(i + 1, scanname)

        self.ieditmulti = self.ieditmulti + 1
        self.ShowMultiscanList()

        self.scriptlist = []
        self.scan.scriptlist = []
        self.scriptdata = ScriptData()
        self.txrmloaded = 0

        self.ResetWidgets()
        self.ShowScriptList()

    # ----------------------------------------------------------------------
    def OnGenerateScan(self, evt):

        if (self.scriptdata.enable_energy == 0) and (self.scriptdata.enable_tomo == 0) and (
                self.scriptdata.enable_mosaic == 0):
            warning_message = ('Please select either Energy, Angle or Mosaic before generating a scan.')
            wx.MessageBox(warning_message, 'Warning',
                          parent=self, style=wx.OK | wx.ICON_WARNING)
            return

        if v: print
        'Generating Scan.'

        self.scan = scriptgen.scriptgen()

        self.scriptdata.RefX = self.ntc_refx.GetValue()
        self.scriptdata.RefY = self.ntc_refy.GetValue()
        self.scriptdata.RefZ = self.ntc_refz.GetValue()
        self.scriptdata.RefTheta = self.ntc_reft.GetValue()
        if self.cb_abba.GetValue():
            self.scriptdata.ABBA = 1
        else:
            self.scriptdata.ABBA = 0

        if self.cb_avfly.GetValue():
            self.scriptdata.averageotfly = 1
        else:
            self.scriptdata.averageotfly = 0

        if self.cb_refavfly.GetValue():
            self.scriptdata.refaverageotfly = 1
        else:
            self.scriptdata.refaverageotfly = 0

        if self.cb_refrelpos.GetValue():
            self.scriptdata.refrelativepos = 1
        else:
            self.scriptdata.refrelativepos = 0

        if self.cb_ioncham.GetValue():
            self.scriptdata.useionchamber = 1
        else:
            self.scriptdata.useionchamber = 0

        self.CheckTomoInerLoop()

        self.scriptdata.wait_energy = self.ntc_ms_waiteng.GetValue()
        self.scriptdata.wait_tomo = self.ntc_ms_waittomo.GetValue()
        self.scriptdata.wait_mosaic = self.ntc_ms_waitmos.GetValue()
        self.scriptdata.wait_multi = self.ntc_ms_waitmulti.GetValue()

        # Scan Settings:
        self.scriptdata.Exptime = self.ntc_exptime.GetValue()
        self.scriptdata.Binning = 2 ** (int(self.combo_Binning.GetSelection()))
        self.scriptdata.SampleName = self.tc_samplename.GetValue()
        self.scriptdata.SampleName = self.scriptdata.SampleName.replace(" ", "")
        self.scriptdata.RefBinning = 2 ** (int(self.combo_RefBinning.GetSelection()))
        self.scriptdata.RefExpTime = self.ntc_refexptime.GetValue()
        self.scriptdata.waitsecs = self.ntc12.GetValue()

        # If position data was not taken from Tomo load the info from a .XRM file
        if (self.scriptdata.enable_tomo == 0):
            wildcard = "XRM files (*.xrm)|*.xrm"
            # Allow loading of multiple files for mosaic scans. Path, binning and exp time will be taken from
            # the first file.
            if (self.scriptdata.enable_mosaic == 1):
                dialog = wx.FileDialog(None, "Please select an xrm file",
                                       wildcard=wildcard,
                                       style=wx.OPEN | wx.FD_MULTIPLE)
            else:
                dialog = wx.FileDialog(None, "Please select an xrm file",
                                       wildcard=wildcard,
                                       style=wx.OPEN)
            if dialog.ShowModal() == wx.ID_OK:
                filepath = dialog.GetPath()
                filename = dialog.GetFilename()
                filedir = dialog.GetDirectory()
                # Get all the info from the first file and positions from the other ones
                filepaths = dialog.GetPaths()

            else:
                return

            if self.microscope:
                xrminfo = self.microscope.get_xrm_info(filepath)

                if xrminfo:
                    self.scriptdata.samplets = []
                    self.scriptdata.samplexs = []
                    self.scriptdata.sampleys = []
                    self.scriptdata.samplezs = []

                    self.scriptdata.samplets.append(xrminfo['MotPosT'])
                    self.scriptdata.samplexs.append(xrminfo['MotPosX'])
                    self.scriptdata.sampleys.append(xrminfo['MotPosY'])
                    self.scriptdata.samplezs.append(xrminfo['MotPosZ'])

                    self.scriptdata.pixelsize = xrminfo['pixelsize']
                    self.scriptdata.width = xrminfo['width']
                    self.scriptdata.height = xrminfo['height']
                    self.scriptdata.FOV = xrminfo['FOV']

                    if len(filepaths) > 1:
                        for i in range(1, len(filepaths)):
                            xrminfo = self.microscope.get_xrm_info(filepaths[i])

                            self.scriptdata.samplets.append(xrminfo['MotPosT'])
                            self.scriptdata.samplexs.append(xrminfo['MotPosX'])
                            self.scriptdata.sampleys.append(xrminfo['MotPosY'])
                            self.scriptdata.samplezs.append(xrminfo['MotPosZ'])
                else:
                    errormsg = 'ERROR - Could not open .xrm file.'
                    wx.MessageBox(errormsg, 'ERROR',
                                  parent=self, style=wx.OK | wx.ICON_ERROR)
                    return

            else:

                Ctxrm = txrm.txrm()

                xrminfo = Ctxrm.get_xrm_info(filepath)

                NumofMotors, xxx, yyy, zzz, ttt, zpxxx, zpyyy, zpzzz = Ctxrm.get_motor_index(xrminfo)

                # print ' NumofMotors, xxx, yyy, zzz, ttt, zpxxx, zpyyy, zpzzz',  NumofMotors, xxx, yyy, zzz, ttt, zpxxx, zpyyy, zpzzz

                self.scriptdata.samplets = []
                self.scriptdata.samplexs = []
                self.scriptdata.sampleys = []
                self.scriptdata.samplezs = []

                self.scriptdata.samplets.append(xrminfo['MotPos'][ttt])
                self.scriptdata.samplexs.append(xrminfo['MotPos'][xxx])
                self.scriptdata.sampleys.append(xrminfo['MotPos'][yyy])
                self.scriptdata.samplezs.append(xrminfo['MotPos'][zzz])

                self.scriptdata.pixelsize = xrminfo['pixelsize']
                self.scriptdata.width = xrminfo['width']
                self.scriptdata.height = xrminfo['height']
                self.scriptdata.FOV = xrminfo['FOV']

                if len(filepaths) > 1:
                    for i in range(1, len(filepaths)):
                        xrminfo = Ctxrm.get_xrm_info(filepaths[i])

                        NumofMotors, xxx, yyy, zzz, ttt, zpxxx, zpyyy, zpzzz = Ctxrm.get_motor_index(xrminfo)

                        self.scriptdata.samplets.append(xrminfo['MotPos'][ttt])
                        self.scriptdata.samplexs.append(xrminfo['MotPos'][xxx])
                        self.scriptdata.sampleys.append(xrminfo['MotPos'][yyy])
                        self.scriptdata.samplezs.append(xrminfo['MotPos'][zzz])

            self.scriptdata.OutputPath = filedir

            self.txrmloaded = False

        else:
            if self.txrmloaded == False:
                error_message = ('ERROR - Please load Info from a .txrm file before generating a scan.')
                wx.MessageBox(error_message, 'ERROR',
                              parent=self, style=wx.OK | wx.ICON_ERROR)
                logging.error(error_message)
                return

                # If we have XANES scan check if we have all the energy info
        if (self.scriptdata.enable_energy > 0):
            totalEpoints = 0
            for ii in range(self.scriptdata.RegionParams.shape[0]):
                totalEpoints = totalEpoints + self.scriptdata.RegionParams[ii, 1] + 1
                # If we have more then 1 energy we need to interpolate - check if we have zoneplate params
            if (totalEpoints > 1) and (len(self.scriptdata.EPlist) == 0):
                error_message = ('ERROR - Please input zoneplate parameters before generating a scan.')
                wx.MessageBox(error_message, 'ERROR',
                              parent=self, style=wx.OK | wx.ICON_ERROR)
                logging.error(error_message, exc_info=True)
                return

        self.scriptdata.loadedfromscript = 0

        self.scan.generate_script(self.scriptdata, self.scriptdata.enable_energy, self.scriptdata.enable_tomo,
                                  self.scriptdata.enable_mosaic, self.scriptdata.enable_multi,
                                  self.scriptdata.layer_eng, self.scriptdata.layer_ang,
                                  UseRefZ=self.scriptdata.useRefZ, UseRefTheta=self.scriptdata.useRefTheta,
                                  SampleXYZAboveRotation=self.SampleXYZAboveRotation,
                                  TakeSingleTomo=self.scriptdata.singletomo, MotorizedDetector=self.MotorizedDetector)

        # Check the number of exposures in the inner loop and make sure it's
        # divisible by reference image frequency
        if self.scriptdata.ref4everycertainnumofexp > 1:
            showwarning = 0
            innerloop = 0
            # find inner loop
            if (self.scriptdata.enable_mosaic == 1):
                innerloop = 3
            elif (self.scriptdata.enable_energy == 1) and (self.scriptdata.enable_tomo == 1):
                if (self.scriptdata.layer_eng == 1):
                    innerloop = 2
                else:
                    innerloop = 1
            elif (self.scriptdata.enable_energy == 1) and (self.scriptdata.enable_tomo == 0):
                innerloop = 1
            elif (self.scriptdata.enable_energy == 0) and (self.scriptdata.enable_tomo == 1):
                innerloop = 2

            # Energy
            if innerloop == 1:
                if (float(totalEpoints * self.scriptdata.n_exposures) % self.scriptdata.ref4everycertainnumofexp != 0):
                    showwarning = 1
            # Tomo
            if innerloop == 2:
                if (float(len(
                        self.scriptdata.samplets) * self.scriptdata.n_exposures) % self.scriptdata.ref4everycertainnumofexp != 0):
                    showwarning = 1
            # Mosaic
            if innerloop == 3:
                if (self.scriptdata.mosaic_centraltile == 1):
                    NColumns = self.scriptdata.mosaic_left + self.scriptdata.mosaic_right + 1
                    NRows = self.scriptdata.mosaic_up + self.scriptdata.mosaic_down + 1
                else:
                    NColumns = self.scriptdata.mosaic_left + self.scriptdata.mosaic_right
                    NRows = self.scriptdata.mosaic_up + self.scriptdata.mosaic_down
                # Multipy by number of positions for multiple files
                nsampleexp = NColumns * NRows * len(self.scriptdata.samplexs)
                if (float(nsampleexp * self.scriptdata.n_exposures) % self.scriptdata.ref4everycertainnumofexp != 0):
                    showwarning = 1

            if showwarning == 1:
                wx.MessageBox("Number of exposures is not divisible by reference taking frequency. \nPlease check.")

        self.tc_outpath.SetValue(self.scriptdata.OutputPath)

        self.CheckLocalDrive(self.scriptdata.OutputPath)

        self.scriptlist = self.scan.scriptlist
        self.ShowScriptList()
        self.UpdateButtonWidgets()
        self.UpdateWidgets()

    # ----------------------------------------------------------------------
    def OnAddMultiScan(self, evt):

        self.multiscript.append(copy.deepcopy(self.scriptdata))
        multiscriptnum = len(self.multiscript)
        self.scan.generate_script(self.scriptdata, self.scriptdata.enable_energy, self.scriptdata.enable_tomo,
                                  self.scriptdata.enable_mosaic, self.scriptdata.enable_multi,
                                  self.scriptdata.layer_eng, self.scriptdata.layer_ang, multiscriptnum=multiscriptnum,
                                  UseRefZ=self.scriptdata.useRefZ, UseRefTheta=self.scriptdata.useRefTheta,
                                  SampleXYZAboveRotation=self.SampleXYZAboveRotation,
                                  TakeSingleTomo=self.scriptdata.singletomo, MotorizedDetector=self.MotorizedDetector)

        newmultilist = []
        newmultilist.append(';;; Multi Scan %d ;;;' % (multiscriptnum))
        for i in range(len(self.scan.scriptlist)):
            newmultilist.append(copy.deepcopy(self.scan.scriptlist[i]))

        self.multiscriptlist.append(list(newmultilist))

        scanname = self.GetScanName(self.scriptdata.enable_energy, self.scriptdata.enable_tomo,
                                    self.scriptdata.enable_mosaic, self.scriptdata.layer_eng, self.scriptdata.layer_ang)

        self.multiscan_namelist.append('{0} {1}'.format(multiscriptnum, scanname))
        self.ieditmulti = multiscriptnum - 1
        self.ShowMultiscanList()

        self.scriptlist = []
        self.scan.scriptlist = []
        self.scriptdata = ScriptData()
        self.txrmloaded = 0

        self.ResetWidgets()
        self.ShowScriptList()

        self.button_9.Enable()
        self.button_14.Enable()
        self.button_15.Enable()
        self.button_16.Enable()

        self.UpdateButtonWidgets()

    # ----------------------------------------------------------------------
    def ShowMultiscanList(self):

        self.tc_multiscan.DeleteAllItems()

        for i in range(len(self.multiscan_namelist)):
            self.tc_multiscan.InsertStringItem(i, self.multiscan_namelist[i])

        if len(self.multiscan_namelist) > 0:
            self.tc_multiscan.SetItemBackgroundColour(self.ieditmulti, 'light blue')

    # ----------------------------------------------------------------------
    def OnTAgitator(self, evt):

        value = self.button_13.GetValue()

        if not value:
            self.button_13.SetLabel('Agitator ON on Exit')
            self.TurnAgitatorOff = False
        else:
            self.button_13.SetLabel('Agitator OFF on Exit')
            self.TurnAgitatorOff = True

    # ----------------------------------------------------------------------
    def OnExecuteScan(self, evt):

        self.scanrunning = True
        self.UpdateButtonWidgets()

        self.beamstatus = False

        if v: print
        'Executing scan.'

        if UseTimer:
            print
            'Starting timer on execute scan:'
            StartTime = time.time()

        fg = wx.Colour(156, 156, 156)
        at = wx.TextAttr(fg)

        # If we have Multiscan build command list
        if (self.resumescan == False):
            self.imultiscan = 0
            self.progress = 0

        scriptlist = []
        self.scriptdata.SampleName = self.tc_samplename.GetValue()
        self.scriptdata.SampleName = self.scriptdata.SampleName.replace(" ", "")
        if (len(self.multiscriptlist) > 0):
            if not self.resumescan:
                self.cb_multiscript.SetValue(True)
                self.ShowScriptList()
            scriptlist = self.GetMultiscanScriptList()

        else:
            scriptlist = self.scriptlist

        if self.resumescan == False:
            if (self.resumeuncompleted == False):
                self.CreateDirAndScanInfo(self.scriptdata, scriptlist)

        if self.resumeuncompleted == True:
            self.NewOutputPath = self.scriptdata.OutputPath

        self.status_warning = ''
        if (self.mono == None):
            self.UpdateStatusPanel('Warning - monochoromator control not initialized.')
            self.status_warning = 'Warning - monochoromator control not initialized.\n'
            logging.warning('Warning - monochoromator control not initialized. Cannot use monochromator.')
            if v: print
            'Warning - monochoromator control not initialized. Cannot use monochromator.'

        if (self.microscope == None):
            self.UpdateStatusPanel('ERROR - XradiaAPI not initialized. Cannot run a script.')
            print
            'ERROR - XradiaAPI not initialized. Cannot run a script.'
            logging.error('ERROR - XradiaAPI not initialized. Cannot run a script.')
            self.scanrunning = False
            self.UpdateButtonWidgets()
            return

        if self.refreshscriptlist:
            self.ShowScriptList()
            self.refreshscriptlist = False

        # Turn on agitator before the scan
        self.microscope.TurnOnAgitator()

        # Scan is executed for the first time, not paused and resumed
        if (self.resumescan == False):
            # Set pixel size
            self.scriptdata.pixelsize = self.ntc_pixelsize.GetValue()
            if (self.scriptdata.pixelsize == self.scriptdata.defaultpixelsize):
                self.status_warning = self.status_warning + 'Warning - using default pixel size: ' + str(
                    self.scriptdata.pixelsize) + '\n'
                self.UpdateStatusPanel(self.status_warning)
            self.microscope.SetPixelSize(self.scriptdata.pixelsize)
            if v: print
            "Setting pixel size: ", self.scriptdata.pixelsize

            self.GetStartingMotorPositions(scriptlist)

            self.imultiscan = 0

            # Get Retiga averaging number of frames
            if self.HaveRetiga:
                self.scriptdata.n_retiga = self.ntc_n_retiga.GetValue()

        self.UpdateStatusPanel(self.status_warning + "Scan is running...")
        self.resumescan = False

        if UseTimer:
            ElapsedTime = time.time() - StartTime
            print
            '\nScan startup. Elapsed Time: ', ElapsedTime
            print
            'Command list:'

        for i in range(len(scriptlist)):

            if i % TCMaxLines == 0:
                self.ShowScriptList(startline=i)

            if UseTimer:
                print
                '\nRunning ', scriptlist[i]
                StartTime2 = time.time()

                # Skip the commands if we paused the scan previously and resuming now
            if (i < self.iscanstopped):
                continue

            if DebugYieldON == True:
                # Check if the scan was paused
                wx.SafeYield(win=self)
                if (self.pausescan):
                    self.iscanstopped = i
                    self.UpdateStatusPanel(self.status_warning + "Scan paused.")
                    return
                if (self.stopscan):
                    self.StopScan()
                    return

            # MultiScan settings
            if ('Multi Scan' in scriptlist[i]) and (self.scriptdata.loadedfromscript == 0):
                self.NewOutputPath = self.multifolderlist[self.imultiscan]
                self.imultiscan = self.imultiscan + 1

            if ';;' in scriptlist[i]:
                continue

            if UseTimer:
                ElapsedTime = time.time() - StartTime2
                print
                'Checking for user input. Elapsed Time: ', ElapsedTime
                StartTime2 = time.time()

            thiscommand = scriptlist[i].split(' ')

            if thiscommand[0] == 'setexp':
                exposure = float(thiscommand[1])
                if v: print
                'set exposure ', exposure
                self.scriptdata.Exptime = exposure


            elif thiscommand[0] == 'setbinning':
                binning = int(thiscommand[1])
                if v: print
                'set binning ', binning
                self.scriptdata.Binning = binning


            elif thiscommand[0] == 'setfilter':
                ifilter = int(thiscommand[1])
                if v: print
                'set filter ', ifilter
                if self.HaveXradiaAPI:
                    if (ifilter < self.MaxNumFilters):
                        status = self.microscope.SetFilter(ifilter)
                        if not status:
                            print
                            'ERROR - could not change filter.'
                            logging.error('ERROR - could not change filter.')
                    else:
                        print
                        'ERROR - Max number of filters is ', self.MaxNumFilters, ' trying to set ', ifilter

            elif thiscommand[0] == 'sete':
                eng = float(thiscommand[1])
                if self.EnergyAxis == 0:
                    if self.mono:
                        if v: print
                        'set energy ', eng
                        status = self.mono.MoveMono(eng)
                        if DelayMono: self.MonoDelay(eng)
                        if not status:
                            error_message = ('ERROR in moving monochromator.')
                            logging.error('ERROR in moving monochromator.')
                            wx.MessageBox(error_message, 'INIT ERROR',
                                          parent=self, style=wx.OK | wx.ICON_ERROR)
                            self.scanrunning = False
                            self.UpdateButtonWidgets()
                            return
                else:
                    axisname = 'energy'
                    if v: print
                    'set energy ', eng
                    if self.HaveXradiaAPI:
                        status = self.microscope.MoveAxis(axisname, eng)
                        if not status:
                            error_message = ('ERROR in moving ' + axisname)
                            wx.MessageBox(error_message, 'INIT ERROR',
                                          parent=self, style=wx.OK | wx.ICON_ERROR)
                            logging.error(error_message)
                            self.scanrunning = False
                            self.UpdateButtonWidgets()
                            return


            elif thiscommand[0] == 'moveto':
                axis = thiscommand[1]
                position = float(thiscommand[2])
                if (axis == 'E') or (axis == 'e'):
                    if self.EnergyAxis == 0:
                        if self.mono:
                            status = self.mono.MoveMono(position)
                            if DelayMono: self.MonoDelay(position)
                            if not status:
                                error_message = ('ERROR in moving monochromator.')
                                wx.MessageBox(error_message, 'INIT ERROR',
                                              parent=self, style=wx.OK | wx.ICON_ERROR)
                                logging.error('ERROR in moving monochromator.')
                                self.scanrunning = False
                                self.UpdateButtonWidgets()
                                return
                    else:
                        axisname = 'energy'
                        if v: print
                        'move motor %s to %f' % (axisname, position)
                        if self.HaveXradiaAPI:
                            status = self.microscope.MoveAxis(axisname, position)
                            if not status:
                                error_message = ('ERROR in moving ' + axisname)
                                wx.MessageBox(error_message, 'INIT ERROR',
                                              parent=self, style=wx.OK | wx.ICON_ERROR)
                                logging.error(error_message)
                                self.scanrunning = False
                                self.UpdateButtonWidgets()
                                return
                    if v: print
                    'move energy to %f' % (position)
                elif axis in XradiaAxisNameList:
                    axisname = thiscommand[1]
                    if v: print
                    'move motor %s to %f' % (axisname, position)
                    if self.HaveXradiaAPI:
                        status = self.microscope.MoveAxis(axisname, position)
                        if not status:
                            error_message = ('ERROR in moving ' + axisname)
                            wx.MessageBox(error_message, 'INIT ERROR',
                                          parent=self, style=wx.OK | wx.ICON_ERROR)
                            logging.error(error_message)
                            self.scanrunning = False
                            self.UpdateButtonWidgets()
                            return
                else:
                    # Move other motors using SSRL SIS
                    if self.mono:
                        status = self.mono.MoveMotor(axis, position)
                        if not status:
                            error_message = ('ERROR in moving motor %s.' % axis)
                            wx.MessageBox(error_message, 'INIT ERROR',
                                          parent=self, style=wx.OK | wx.ICON_ERROR)
                            logging.error('ERROR in moving motor %s.' % axis)
                            self.scanrunning = False
                            self.UpdateButtonWidgets()
                            return
                        if v: print
                        'move motor %s to %f' % (axis, position)


            elif thiscommand[0] == 'moveby':
                axis = thiscommand[1].strip()
                position = float(thiscommand[2])
                if axis == 'E':
                    if self.mono:
                        status = self.mono.MoveMonoBy(position)
                        if DelayMono:
                            if position > mono_threshold: time.sleep(mono_wait_time)
                        if not status:
                            error_message = ('ERROR in moving monochromator.')
                            logging.error('ERROR in moving monochromator.')
                            wx.MessageBox(error_message, 'INIT ERROR',
                                          parent=self, style=wx.OK | wx.ICON_ERROR)
                            self.scanrunning = False
                            self.UpdateButtonWidgets()
                            return
                        if v: print
                        'move energy by %f' % (position)
                elif axis in XradiaAxisNameList:
                    axisname = thiscommand[1]
                    if v: print
                    'move motor %s by %f' % (axisname, position)
                    if self.HaveXradiaAPI:
                        status = self.microscope.MoveAxisBy(axisname, position)
                        if not status:
                            error_message = ('ERROR in moving ' + axisname)
                            wx.MessageBox(error_message, 'INIT ERROR',
                                          parent=self, style=wx.OK | wx.ICON_ERROR)
                            logging.error(error_message)
                            self.scanrunning = False
                            self.UpdateButtonWidgets()
                            return

                else:
                    # Move other motors using SSRL SIS

                    print
                    'move motor %s by %f' % (axis, position)
                    if self.mono:
                        print
                        'LOC-1'
                        status = self.mono.MoveMotorBy(axis, position)
                        if not status:
                            error_message = ('ERROR in moving motor %s.' % axis)
                            logging.error('ERROR in moving motor %s.' % axis)
                            wx.MessageBox(error_message, 'INIT ERROR',
                                          parent=self, style=wx.OK | wx.ICON_ERROR)
                            self.scanrunning = False
                            self.UpdateButtonWidgets()
                            return
                        if v: print
                        'move motor %s by %f' % (axis, position)

            elif thiscommand[0] == 'collectxrf':
                filename = thiscommand[1]
                print
                filename
                if v: print
                'collectxrf ', filename
                self.beamstatus = self.CheckBeamStatus(previousstatus=self.beamstatus)
                if UseTimer:
                    ElapsedTime = time.time() - StartTime2
                    print
                    '\tCollect: Checking beam status. Elapsed Time: ', ElapsedTime
                if not self.beamstatus:
                    if (self.pausescan):
                        # Scan was paused while waiting for the beam
                        self.iscanstopped = i
                        self.UpdateStatusPanel(self.status_warning + "Scan paused.")
                        return
                    elif (self.stopscan):
                        self.StopScan()
                        returna

                app = Application(backend="uia").connect(path=r'C:\Program Files (x86)\XIA\ProSpect 1.1\\ProSpect.exe')

                parent_dir = 'C:\\TXMDATA\\202102\\20210303_Jizhou\\'
                directory = 'XRFAuto\\'
                path = os.path.join(parent_dir, directory)
                if not os.path.exists(path):
                    os.mkdir(path)
                    print("Directory '% s' created" % parent_dir + directory)

                print('open xrf signal...')
                app['ProSpect - v1.1.12'].Start_Run.click()

                time.sleep(0.1)

                app['ProSpect - v1.1.12'].Stop_Run.click()

                # ct = os.path.join(path,strftime("%Y-%m-%d_%H-%M-%S", gmtime()))
                ct = filename

                win = app['ProSpect - v1.1.12']
                submenu = app['ProSpect - v1.1.12']['']  # Dropdown submenu is a top-level window
                sub = win['File']
                savemcaMenuItem = (sub.children()[11])
                try:
                    savemcaMenuItem.select()
                    savemcaMenuItem.select()
                except:
                    print('1st error, no matter')

                print('saving signal...')
                Sub = app['ProSpect - v1.1.12'].child_window(title_re="Save MCA data File As...", class_name="#32770")

                Sub.FileNameEdit.type_keys(str(ct))
                Sub.Save.click()

                print('save finished ' + str(ct))


            elif thiscommand[0] == 'collect':
                # Check if the next command is also collect if it is leave the shutter opened and skip beam
                # status check on the next command
                leaveshutteropened = False
                nextcommand = ['']
                for j in range(i + 1, len(scriptlist)):
                    if ';;' in scriptlist[j]:
                        continue
                    nextcommand = scriptlist[j].split(' ')
                    break
                if 'collect' in nextcommand[0]:
                    leaveshutteropened = True

                filename = thiscommand[1]
                if v: print
                'collect ', filename
                self.beamstatus = self.CheckBeamStatus(previousstatus=self.beamstatus)
                if UseTimer:
                    ElapsedTime = time.time() - StartTime2
                    print
                    '\tCollect: Checking beam status. Elapsed Time: ', ElapsedTime
                if not self.beamstatus:
                    if (self.pausescan):
                        # Scan was paused while waiting for the beam
                        self.iscanstopped = i
                        self.UpdateStatusPanel(self.status_warning + "Scan paused.")
                        return
                    elif (self.stopscan):
                        self.StopScan()
                        return
                if self.HaveXradiaAPI:
                    #                     #If crashed scan was resumed, delete image files before acquisition
                    #                     if self.resumeuncompleted == True:
                    #                         if os.path.isfile(os.path.join(self.NewOutputPath, filename)):
                    #                             os.remove(os.path.join(self.NewOutputPath, filename))
                    #                             print 'Deleted file: ', os.path.join(self.NewOutputPath, filename)

                    if not self.microscope.StartSingleAcquisition(str(os.path.join(self.NewOutputPath, filename)),
                                                                  self.scriptdata.Binning,
                                                                  self.scriptdata.Exptime,
                                                                  leaveshutteropened=leaveshutteropened,
                                                                  useretigaaveraging=self.HaveRetiga,
                                                                  n_frames=self.scriptdata.n_retiga):
                        error_message = ('ERROR - Could not collect image.')
                        wx.MessageBox(error_message, 'INIT ERROR',
                                      parent=self, style=wx.OK | wx.ICON_ERROR)
                        logging.error(error_message)
                        self.UpdateStatusPanel("ERROR - Could not collect image.")
                        self.iscanstopped = 0
                        self.scanrunning = False
                        self.UpdateButtonWidgets()
                        return


            elif thiscommand[0] == 'collectaverage':
                # Check if the next command is also collect if it is leave the sutter opened and skip beam
                # status check on the next command
                leaveshutteropened = False
                nextcommand = ['']
                for j in range(i + 1, len(scriptlist)):
                    if ';;' in scriptlist[j]:
                        continue
                    nextcommand = scriptlist[j].split(' ')
                    break
                if 'collect' in nextcommand[0]:
                    leaveshutteropened = True

                filename = thiscommand[1]
                if v: print
                'collectaverage ', filename
                ntotalimages = int(thiscommand[2])
                if v: print
                'Average n images: ', ntotalimages
                self.beamstatus = self.CheckBeamStatus(previousstatus=self.beamstatus)
                if UseTimer:
                    ElapsedTime = time.time() - StartTime2
                    print
                    '\tCollectaverage: Checking beam status. Elapsed Time: ', ElapsedTime
                if not self.beamstatus:
                    if (self.pausescan):
                        # Scan was paused while waiting for the beam
                        self.iscanstopped = i
                        self.UpdateStatusPanel(self.status_warning + "Scan paused.")
                        return
                    elif (self.stopscan):
                        self.StopScan()
                        return
                if self.HaveXradiaAPI:
                    if not self.microscope.StartAveragedAcquisition(str(os.path.join(self.NewOutputPath, filename)),
                                                                    self.scriptdata.Binning,
                                                                    self.scriptdata.Exptime,
                                                                    ntotalimages,
                                                                    leaveshutteropened=leaveshutteropened,
                                                                    useretigaaveraging=self.HaveRetiga,
                                                                    n_frames=self.scriptdata.n_retiga):
                        error_message = ('ERROR - Could not collect average image.')
                        wx.MessageBox(error_message, 'INIT ERROR',
                                      parent=self, style=wx.OK | wx.ICON_ERROR)
                        logging.error(error_message)
                        self.UpdateStatusPanel("ERROR - Could not collect average image.")
                        self.iscanstopped = 0
                        self.scanrunning = False
                        self.UpdateButtonWidgets()
                        return


            elif thiscommand[0] == 'collectnormalized':
                if (self.IonCurrentReading == 0):
                    error_message = (
                        'ERROR - Could not collect normalized image. Recorded ion current reading is 0. \nPlease use recordioncurrent before calling collectnormalized.')
                    wx.MessageBox(error_message, 'INIT ERROR',
                                  parent=self, style=wx.OK | wx.ICON_ERROR)
                    self.StopScan()
                    return
                # Check if the next command is also collect if it is leave the sutter opened and skip beam
                # status check on the next command
                leaveshutteropened = False
                nextcommand = ['']
                for j in range(i + 1, len(scriptlist)):
                    if ';;' in scriptlist[j]:
                        continue
                    nextcommand = scriptlist[j].split(' ')
                    break
                if 'collect' in nextcommand[0]:
                    leaveshutteropened = True

                filename = thiscommand[1]
                if v: print
                'collectnormalized', filename
                self.beamstatus = self.CheckBeamStatus(previousstatus=self.beamstatus)
                if UseTimer:
                    ElapsedTime = time.time() - StartTime2
                    print
                    '\tCollect: Checking beam status. Elapsed Time: ', ElapsedTime
                if not self.beamstatus:
                    if (self.pausescan):
                        # Scan was paused while waiting for the beam
                        self.iscanstopped = i
                        self.UpdateStatusPanel(self.status_warning + "Scan paused.")
                        return
                    elif (self.stopscan):
                        self.StopScan()
                        return

                # This collects an image using the recorded ion chamber current
                # value to scale the exposure time using the equation :
                # exposure time *  Recorded Ion Chamber Reading /current IonChamber Reading
                NewIonCurrentReading = self.microscope.GetIonChamberReading()
                if NewIonCurrentReading != 0:
                    NormExpTime = float(self.scriptdata.Exptime) * self.IonCurrentReading / NewIonCurrentReading
                else:
                    error_message = ('Ion chamber reading is equal to 0. Exposure time will not be scaled.')
                    wx.MessageBox(error_message, 'INIT ERROR',
                                  parent=self, style=wx.OK | wx.ICON_ERROR)
                    NormExpTime = float(self.scriptdata.Exptime)

                if self.HaveXradiaAPI:
                    if not self.microscope.StartSingleAcquisition(str(os.path.join(self.NewOutputPath, filename)),
                                                                  self.scriptdata.Binning,
                                                                  NormExpTime,
                                                                  leaveshutteropened=leaveshutteropened,
                                                                  useretigaaveraging=self.HaveRetiga,
                                                                  n_frames=self.scriptdata.n_retiga):
                        error_message = ('ERROR - Could not collect normalized image.')
                        wx.MessageBox(error_message, 'INIT ERROR',
                                      parent=self, style=wx.OK | wx.ICON_ERROR)
                        logging.error(error_message)
                        self.UpdateStatusPanel("ERROR - Could not collect normalized image.")
                        self.iscanstopped = 0
                        self.scanrunning = False
                        self.UpdateButtonWidgets()
                        return


            elif thiscommand[0] == 'collecttomo':
                # Check if the next command is also collect if it is leave the sutter opened and skip beam
                # status check on the next command
                leaveshutteropened = False
                nextcommand = ['']
                for j in range(i + 1, len(scriptlist)):
                    if ';;' in scriptlist[j]:
                        continue
                    nextcommand = scriptlist[j].split(' ')
                    break
                if 'collect' in nextcommand[0]:
                    leaveshutteropened = True

                filename = thiscommand[1]
                totalimages = int(thiscommand[2])
                startangle = float(thiscommand[3])
                endangle = float(thiscommand[4])
                if v: print
                'collecttomo ', filename
                self.beamstatus = self.CheckBeamStatus(previousstatus=self.beamstatus)
                if UseTimer:
                    ElapsedTime = time.time() - StartTime2
                    print
                    '\tCollect: Checking beam status. Elapsed Time: ', ElapsedTime
                if not self.beamstatus:
                    if (self.pausescan):
                        # Scan was paused while waiting for the beam
                        self.iscanstopped = i
                        self.UpdateStatusPanel(self.status_warning + "Scan paused.")
                        return
                    elif (self.stopscan):
                        self.StopScan()
                        return
                if self.HaveXradiaAPI:
                    if not self.microscope.StartTomoAcquisition(str(os.path.join(self.NewOutputPath, filename)),
                                                                self.scriptdata.Binning,
                                                                self.scriptdata.Exptime,
                                                                totalimages,
                                                                startangle,
                                                                endangle,
                                                                leaveshutteropened=leaveshutteropened,
                                                                useretigaaveraging=self.HaveRetiga,
                                                                n_frames=self.scriptdata.n_retiga):
                        error_message = ('ERROR - Could not collect image.')
                        wx.MessageBox(error_message, 'INIT ERROR',
                                      parent=self, style=wx.OK | wx.ICON_ERROR)
                        logging.error(error_message)
                        self.UpdateStatusPanel("ERROR - Could not collect image.")
                        self.iscanstopped = 0
                        self.scanrunning = False
                        self.UpdateButtonWidgets()
                        return


            elif thiscommand[0] == 'collectnormalizedtomo':
                # Check if the next command is also collect if it is leave the sutter opened and skip beam
                # status check on the next command
                leaveshutteropened = False
                nextcommand = ['']
                for j in range(i + 1, len(scriptlist)):
                    if ';;' in scriptlist[j]:
                        continue
                    nextcommand = scriptlist[j].split(' ')
                    break
                if 'collect' in nextcommand[0]:
                    leaveshutteropened = True

                filename = thiscommand[1]
                totalimages = int(thiscommand[2])
                startangle = float(thiscommand[3])
                endangle = float(thiscommand[4])
                if v: print
                'collectnormalizedtomo ', filename
                self.beamstatus = self.CheckBeamStatus(previousstatus=self.beamstatus)
                if UseTimer:
                    ElapsedTime = time.time() - StartTime2
                    print
                    '\tCollect: Checking beam status. Elapsed Time: ', ElapsedTime
                if not self.beamstatus:
                    if (self.pausescan):
                        # Scan was paused while waiting for the beam
                        self.iscanstopped = i
                        self.UpdateStatusPanel(self.status_warning + "Scan paused.")
                        return
                    elif (self.stopscan):
                        self.StopScan()
                        return

                # This collects an image using the recorded ion chamber current
                # value to scale the exposure time using the equation :
                # exposure time *  Recorded Ion Chamber Reading /current IonChamber Reading
                NewIonCurrentReading = self.microscope.GetIonChamberReading()
                if NewIonCurrentReading != 0:
                    NormExpTime = float(self.scriptdata.Exptime) * self.IonCurrentReading / NewIonCurrentReading
                else:
                    error_message = ('Ion chamber reading is equal to 0. Exposure time will not be scaled.')
                    wx.MessageBox(error_message, 'INIT ERROR',
                                  parent=self, style=wx.OK | wx.ICON_ERROR)
                    NormExpTime = float(self.scriptdata.Exptime)

                if self.HaveXradiaAPI:
                    if not self.microscope.StartTomoAcquisition(str(os.path.join(self.NewOutputPath, filename)),
                                                                self.scriptdata.Binning,
                                                                NormExpTime,
                                                                totalimages,
                                                                startangle,
                                                                endangle,
                                                                normalized=True,
                                                                leaveshutteropened=leaveshutteropened,
                                                                useretigaaveraging=self.HaveRetiga,
                                                                n_frames=self.scriptdata.n_retiga):
                        error_message = ('ERROR - Could not collect image.')
                        wx.MessageBox(error_message, 'INIT ERROR',
                                      parent=self, style=wx.OK | wx.ICON_ERROR)
                        logging.error(error_message)
                        self.UpdateStatusPanel("ERROR - Could not collect image.")
                        self.iscanstopped = 0
                        self.scanrunning = False
                        self.UpdateButtonWidgets()
                        return


            elif thiscommand[0] == 'wait':
                waittimesecs = float(thiscommand[1])
                if v: print
                'wait %f' % (waittimesecs)
                time.sleep(waittimesecs)

            elif thiscommand[0] == 'recordioncurrent':
                self.IonCurrentReading = self.microscope.GetIonChamberReading()
                if v: print
                'Ion Current Reading %f' % (self.IonCurrentReading)

            elif thiscommand[0] == 'setreadout':
                readout = int(thiscommand[1])
                if v: print
                'readout %d' % (readout)
                if readout in [1, 2]:
                    self.microscope.SetReadoutIndex(readout)
                else:
                    wx.MessageBox(error_message, 'Error - setreadout. Readout Index allowed values are 1 and 2. ',
                                  parent=self, style=wx.OK | wx.ICON_ERROR)

            elif thiscommand[0] == 'shiftobjective':
                lensindex = int(thiscommand[1])
                if v: print
                'shiftobjective %f' % (lensindex)
                if not self.microscope.SetCurrentLens(lensindex):
                    wx.MessageBox(error_message, 'Error - shiftobjective. Could not set lens index.',
                                  parent=self, style=wx.OK | wx.ICON_ERROR)

            elif thiscommand[0] == 'runshortcut':
                command = ' '.join(thiscommand[1:])
                if v: print
                'runshortcut %s' % (command)
                os.startfile(command)

            else:
                print
                'ERROR - unrecognized command in the script.'
                logging.error('ERROR - unrecognized command in the script.')
                self.scanrunning = False
                self.UpdateButtonWidgets()
                return

            if DebugTextCtrlON == True:
                itcline = i % TCMaxLines
                p1 = self.tc_scriptlist.XYToPosition(0, itcline)
                p2 = self.tc_scriptlist.XYToPosition(self.tc_scriptlist.GetLineLength(itcline), itcline)
                self.tc_scriptlist.ShowPosition(p1)
                self.tc_scriptlist.SetStyle(p1, p2, at)

            if UseTimer:
                ElapsedTime = time.time() - StartTime2
                print
                '\tTotal Elapsed Time: ', ElapsedTime

            if i % 5 == 0:
                self.progress = int(float(i) / len(scriptlist) * 100)
                self.UpdateStatusPanel(
                    self.status_warning + "Scan is running...\nProgress: " + str(self.progress) + '%')

        # Leave agitation on after the scan Sep-2013
        # Turn off agitator after the scan
        # self.microscope.TurnOffAgitator()

        if UseTimer:
            ElapsedTime = time.time() - StartTime
            print
            '\nTotal Scan Run Time: ', ElapsedTime
            print
            '\n'

        self.StopScan()
        self.UpdateButtonWidgets()

        return

    # ----------------------------------------------------------------------
    def OnPauseScan(self, evt):

        # Shutter could be opened if scan was paused between mulitple image collections
        self.microscope.CloseXrayShutter()

        self.pausescan = True
        self.iscanstopped = 0

        self.button_7.Disable()
        self.button_8.Enable()

    # ----------------------------------------------------------------------
    def OnResumeScan(self, evt):

        self.pausescan = False
        self.resumescan = True

        self.button_7.Enable()
        self.button_8.Disable()

        self.OnExecuteScan(evt)

    # ----------------------------------------------------------------------
    def OnStopScan(self, evt):

        self.stopscan = True

    # ----------------------------------------------------------------------
    def OnResumeUncompletedScan(self, evt):

        ResumeScanFrame().Show()

    # ----------------------------------------------------------------------
    def ResumeUncompletedScan(self, scriptlist, imgdir):

        self.resumeuncompleted = True

        self.scriptlist = scriptlist

        self.scriptdata.OutputPath = imgdir
        self.scriptdata.loadedfromscript = 0

        self.tc_outpath.SetValue(self.scriptdata.OutputPath)

        self.CheckLocalDrive(self.scriptdata.OutputPath)

        self.ShowScriptList()
        self.UpdateButtonWidgets()

    # ----------------------------------------------------------------------
    def StopScan(self):

        if self.HaveXradiaAPI:
            self.microscope.SetFilter(0)
            # Shutter could be opened if scan was paused between mulitple image \ions
            self.microscope.CloseXrayShutter()

        self.scanrunning = False
        self.pausescan = False
        self.resumescan = False
        self.resumeuncompleted = False
        self.iscanstopped = 0
        self.UpdateButtonWidgets()

        self.UpdateStatusPanel(self.status_warning + "Returning motors to starting positions.")
        self.ReturnMotorsToStartingPositions()

        if self.stopscan:
            self.UpdateStatusPanel(self.status_warning + "Scan stopped.")
        else:
            self.UpdateStatusPanel(self.status_warning + "Scan is finished.")

        self.multiscan_foldernum = 1

        self.stopscan = False
        self.refreshscriptlist = True

    # ----------------------------------------------------------------------
    def GetMultiscanScriptList(self):

        if (self.multiscan_userepeat == 0):

            scriptlist = []

            for ii in range(len(self.multiscriptlist)):
                for jj in range(len(self.multiscriptlist[ii])):
                    scriptlist.append(self.multiscriptlist[ii][jj])

        else:

            self.multiscan_repeat = self.ntc_ms_repeat.GetValue()
            if self.multiscan_repeat == None:
                self.multiscan_repeat = 1
            self.multiscan_wait = self.ntc_ms_wait.GetValue()

            scriptlist = []
            scriptlist_single = []

            for ii in range(len(self.multiscriptlist)):
                for jj in range(len(self.multiscriptlist[ii])):
                    scriptlist_single.append(self.multiscriptlist[ii][jj])

            for i in range(self.multiscan_repeat):
                for ii in range(len(scriptlist_single)):
                    scriptlist.append(scriptlist_single[ii])

                if (self.multiscan_wait > 0) and (i < self.multiscan_repeat - 1):
                    script = 'wait %.0f' % (self.multiscan_wait)
                    scriptlist.append(script)

        return scriptlist

    # ----------------------------------------------------------------------
    def CheckBeamStatus(self, previousstatus=False):

        DelaySecs = 1.0
        waiting = True

        updatestatuspanel = True

        while waiting == True:

            if self.HaveXradiaAPI: beamstatus = self.microscope.CheckBeamStatus(self.IonChamber_nobeam)
            if v: print
            'Beam Status =', beamstatus
            if beamstatus:
                if v: print
                'Beam is ON.'
                if previousstatus == False:
                    self.UpdateStatusPanel(
                        self.status_warning + "Scan is running...\nProgress: " + str(self.progress) + '%')
                waiting = False
            else:
                if updatestatuspanel == True:
                    self.UpdateStatusPanel(self.status_warning + "Scan paused, waiting for the beam.")
                    updatestatuspanel = False
                    if DebugYieldON == True:
                        wx.SafeYield(win=self)
                if v: print
                'Beam is OFF.'
                # wait for DelaySecs seconds
                time.sleep(DelaySecs)

            if DebugYieldON == True:
                # Check if a user paused the scan
                wx.SafeYield(win=self)
                if (self.pausescan):
                    return False
                if (self.stopscan):
                    return False

        return True

    # ----------------------------------------------------------------------
    def MonoDelay(self, newposition):

        status, oldposition = self.mono.ShowMonoStatus()
        oldstatus = self.tc_status.GetValue()
        newstatus = oldstatus + '\nMono is moving, waiting ' + str(mono_wait_time) + ' seconds.\n'
        self.UpdateStatusPanel(newstatus)
        wx.SafeYield(win=self)
        moveby = np.abs(oldposition - newposition)
        if moveby > mono_threshold: time.sleep(mono_wait_time)
        self.UpdateStatusPanel(oldstatus)

        return

    # ----------------------------------------------------------------------
    def GetStartingMotorPositions(self, scriptlist):

        if v: print
        "Getting starting motor positions."

        # Motor starting positions 'motor name':[is the motor moved, starting position]
        self.StartingMotorPositions = {'mono': [0, 0.0], 'x': [0, 0.0], 'y': [0, 0.0], 'z': [0, 0.0], 't': [0, 0.0],
                                       'zpx': [0, 0.0], 'zpy': [0, 0.0], 'zpz': [0, 0.0],
                                       'cdx': [0, 0.0], 'cdy': [0, 0.0], 'cdz': [0, 0.0],
                                       'cdtip': [0, 0.0], 'cdtil': [0, 0.0],
                                       'prx': [0, 0.0], 'pry': [0, 0.0], 'prz': [0, 0.0],
                                       'phx': [0, 0.0], 'phy': [0, 0.0], 'phz': [0, 0.0]}

        # Check which motor will be moved by the script
        for i in range(len(scriptlist)):

            if ';;' in scriptlist[i]:
                continue

            thiscommand = scriptlist[i].split(' ')

            if thiscommand[0] == 'sete':
                self.StartingMotorPositions['mono'][0] = 1


            elif thiscommand[0] == 'moveto':
                axis = thiscommand[1]
                if axis == 'E':
                    self.StartingMotorPositions['mono'][0] = 1
                elif axis in self.StartingMotorPositions:
                    self.StartingMotorPositions[thiscommand[1]][0] = 1


            elif thiscommand[0] == 'moveby':
                axis = thiscommand[1]
                if axis == 'E':
                    self.StartingMotorPositions['mono'][0] = 1
                elif axis in self.StartingMotorPositions:
                    self.StartingMotorPositions[thiscommand[1]][0] = 1

                    # Get starting motor positions
        for item in self.StartingMotorPositions:
            if self.StartingMotorPositions[item][0] == 1:
                if (item == 'mono'):
                    if not self.mono: continue
                    status, position = self.mono.ShowMonoStatus()
                    self.StartingMotorPositions[item][1] = position
                else:
                    axisname = item
                    position = self.microscope.GetAxisPosition(axisname)
                    self.StartingMotorPositions[item][1] = position

        return

    # ----------------------------------------------------------------------
    def ReturnMotorsToStartingPositions(self):

        if v: print
        "Returning motors to starting positions."

        for item in self.StartingMotorPositions:
            if self.StartingMotorPositions[item][0] == 1:
                if (item == 'mono'):
                    if not self.mono: continue
                    position = self.StartingMotorPositions[item][1]
                    status = self.mono.MoveMono(position)
                    if not status:
                        error_message = ('ERROR in moving monochromator.')
                        logging.error(error_message)
                        wx.MessageBox(error_message, 'INIT ERROR',
                                      parent=self, style=wx.OK | wx.ICON_ERROR)
                else:
                    axisname = item
                    position = self.StartingMotorPositions[item][1]
                    if not self.HaveXradiaAPI: continue
                    status = self.microscope.MoveAxis(axisname, position)
                    if not status:
                        error_message = ('ERROR in moving ' + axisname)
                        logging.error(error_message)
                        wx.MessageBox(error_message, 'INIT ERROR',
                                      parent=self, style=wx.OK | wx.ICON_ERROR)

                    # ----------------------------------------------------------------------

    def ShowScriptList(self, startline=0):

        scriptlist = []
        if (len(self.multiscriptlist) > 0) and (self.cb_multiscript.GetValue()):
            scriptlist = self.GetMultiscanScriptList()
        else:
            scriptlist = self.scriptlist

        self.tc_scriptlist.Clear()

        if self.scanrunning == False:
            for i in range(len(scriptlist)):
                self.tc_scriptlist.AppendText('{0:s}\n'.format(scriptlist[i]))
        else:
            # Show only TCMaxLines
            for i in range(startline, min(len(scriptlist), startline + TCMaxLines)):
                self.tc_scriptlist.AppendText('{0:s}\n'.format(scriptlist[i]))

        self.tc_scriptlist.SetInsertionPoint(0)

        nexposures, exposuretime = self.GetNExposuresandExpTime(scriptlist)

        if (len(scriptlist) > 0) and (nexposures > 0):
            text = 'Number of Exposures:\t%d\nTotal Exposure Time:\t%s' % (nexposures,
                                                                           str(datetime.timedelta(
                                                                               seconds=exposuretime)))
        else:
            text = ''

        self.UpdateStatusPanel(text)

    # ----------------------------------------------------------------------
    def GetNExposuresandExpTime(self, scriptlist):

        # Calculate the number of exposures
        nexposures = 0
        totalexposuretime = 0
        exposure = 1
        for i in range(len(scriptlist)):
            if ';' in scriptlist[i]:
                continue
            thiscommand = scriptlist[i].split(' ')

            if thiscommand[0] == 'setexp':
                exposure = float(thiscommand[1])

            elif thiscommand[0] == 'collect':
                nexposures += 1
                totalexposuretime += exposure

            elif thiscommand[0] == 'collectaverage':
                avnexposures = int(thiscommand[2])
                nexposures += avnexposures
                totalexposuretime += exposure * avnexposures

            elif thiscommand[0] == 'collecttomo':
                totalimgs = int(thiscommand[2])
                nexposures += totalimgs
                totalexposuretime += exposure * totalimgs

            elif thiscommand[0] == 'collectnormalizedtomo':
                totalimgs = int(thiscommand[2])
                nexposures += totalimgs
                totalexposuretime += exposure * totalimgs

        return nexposures, totalexposuretime

    # ----------------------------------------------------------------------
    def UpdateStatusPanel(self, status):

        self.tc_status.Clear()

        self.tc_status.SetValue(status)
        self.tc_status.Refresh()

    # ----------------------------------------------------------------------
    def UpdateWidgets(self):

        self.UpdateButtonWidgets()

        self.cb_singletomo.SetValue(self.scriptdata.singletomo)

        self.ShowEnabledTC()

        self.combo_eng.SetStringSelection(str(self.scriptdata.layer_eng))
        self.combo_ang.SetStringSelection(str(self.scriptdata.layer_ang))

        self.tc_nexp.SetValue(self.scriptdata.n_exposures)

        self.cb_avfly.SetValue(self.scriptdata.averageotfly)
        self.cb_refavfly.SetValue(self.scriptdata.refaverageotfly)

        self.ntc_refx.SetValue(self.scriptdata.RefX)
        self.ntc_refy.SetValue(self.scriptdata.RefY)
        self.ntc_refz.SetValue(self.scriptdata.RefZ)
        self.ntc_reft.SetValue(self.scriptdata.RefTheta)
        self.ntc_ncol.SetValue(self.scriptdata.RepeatRefExposure)
        self.ntc_nexp.SetValue(self.scriptdata.ref4everycertainnumofexp)
        self.cb_abba.SetValue(self.scriptdata.ABBA)
        self.ntc_refexptime.SetValue(self.scriptdata.RefExpTime)
        self.combo_RefBinning.SetStringSelection(str(self.scriptdata.RefBinning))
        self.cb_refz.SetValue(self.scriptdata.useRefZ)
        self.cb_reft.SetValue(self.scriptdata.useRefTheta)
        self.cb_refrelpos.SetValue(self.scriptdata.refrelativepos)

        self.cb_avfly.SetValue(self.scriptdata.averageotfly)
        self.cb_refavfly.SetValue(self.scriptdata.refaverageotfly)

        if self.scriptdata.useRefZ:
            self.ntc_refz.Enable()
        else:
            self.ntc_refz.Disable()

        if self.scriptdata.useRefTheta:
            self.ntc_reft.Enable()
        else:
            self.ntc_reft.Disable()

        self.tc_samplename.SetValue(self.scriptdata.SampleName)
        self.tc_outpath.SetValue(self.scriptdata.OutputPath)
        self.ntc_exptime.SetValue(self.scriptdata.Exptime)
        self.combo_Binning.SetStringSelection(str(self.scriptdata.Binning))
        self.ntc_pixelsize.SetValue(self.scriptdata.pixelsize)

        self.scriptdata.RefX = self.ntc_refx.GetValue()
        self.scriptdata.RefY = self.ntc_refy.GetValue()
        self.scriptdata.RefZ = self.ntc_refz.GetValue()
        self.scriptdata.RefTheta = self.ntc_reft.GetValue()

        self.cb_ioncham.SetValue(self.scriptdata.useionchamber)

    # ----------------------------------------------------------------------
    def UpdateButtonWidgets(self):

        if self.scanrunning == False:
            self.button_1.Enable()
            self.button_2.Enable()
            self.button_3.Enable()
            self.button_4.Enable()
            self.button_5.Enable()
            self.button_11.Enable()

            self.button_7.Disable()
            self.button_8.Disable()
            self.button_12.Disable()

            self.button_eng.Enable()
            self.button_ang.Enable()
            self.button_mos.Enable()
            self.button_refpos.Enable()
            self.button_13.Enable()
            self.combo_eng.Enable()
            self.combo_ang.Enable()
            self.combo_me.Enable()
            self.combo_mos.Enable()

            if (len(self.scriptlist) > 0) or (len(self.multiscriptlist) > 0):
                self.button_6.Enable()
                self.button_11.Enable()

                self.pausescan = False
                self.resumescan = False
                self.iscanstopped = 0

            else:
                self.button_6.Disable()
                self.button_11.Disable()

            if (len(self.scriptlist) > 0):
                self.button_5.Enable()
            else:
                self.button_5.Disable()

        else:
            self.button_1.Disable()
            self.button_2.Disable()
            self.button_3.Disable()
            self.button_4.Disable()
            self.button_5.Disable()
            self.button_6.Disable()
            self.button_11.Disable()

            self.button_eng.Disable()
            self.button_ang.Disable()
            self.button_mos.Disable()
            self.button_refpos.Disable()
            self.button_13.Disable()
            self.combo_eng.Disable()
            self.combo_ang.Disable()
            self.combo_me.Disable()
            self.combo_mos.Disable()

            self.button_7.Enable()
            self.button_12.Enable()

        # ----------------------------------------------------------------------

    def ResetWidgets(self):
        self.scriptdata.enable_energy = 0
        self.scriptdata.enable_tomo = 0
        self.scriptdata.enable_mosaic = 0
        self.scriptdata.enable_multi = 0

        self.cb_eng.SetValue(False)
        self.cb_ang.SetValue(False)
        self.cb_mos.SetValue(False)
        self.cb_me.SetValue(False)

        self.tc_e.Clear()
        self.tc_an.Clear()
        self.tc_mos.Clear()
        self.tc_me.Clear()

        self.cb_1.SetValue(False)
        self.ntc1.SetValue(self.scriptdata.NRepeatScan)
        self.ntc1.Disable()
        self.ntc12.SetValue(self.scriptdata.waitsecs)
        self.ntc12.Disable()

        self.combo_eng.SetStringSelection(str(self.scriptdata.layer_eng))
        self.combo_ang.SetStringSelection(str(self.scriptdata.layer_ang))

        self.tc_nexp.SetValue(self.scriptdata.n_exposures)

        self.scriptdata.wait_energy = 0
        self.ntc_ms_waiteng.SetValue(self.scriptdata.wait_energy)
        self.scriptdata.wait_tomo = 0
        self.ntc_ms_waittomo.SetValue(self.scriptdata.wait_tomo)
        self.scriptdata.wait_mosaic = 0
        self.ntc_ms_waitmos.SetValue(self.scriptdata.wait_mosaic)
        self.scriptdata.wait_multi = 0
        self.ntc_ms_waitmulti.SetValue(self.scriptdata.wait_multi)

        self.scriptdata.singletomo = 0
        self.cb_singletomo.SetValue(self.scriptdata.singletomo)

        self.ntc_refx.SetValue(0.0)
        self.ntc_refy.SetValue(0.0)
        self.ntc_refz.SetValue(0.0)
        self.ntc_reft.SetValue(0.0)
        self.ntc_ncol.SetValue(0)
        self.ntc_nexp.SetValue(1)
        self.cb_abba.SetValue(False)
        self.ntc_refexptime.SetValue(1.0)
        self.combo_RefBinning.SetStringSelection(str(1))
        self.cb_refz.SetValue(self.scriptdata.useRefZ)
        self.cb_reft.SetValue(self.scriptdata.useRefTheta)
        self.cb_refrelpos.SetValue(self.scriptdata.refrelativepos)

        self.cb_refavfly.SetValue(self.scriptdata.refaverageotfly)
        self.cb_avfly.SetValue(self.scriptdata.averageotfly)

        self.tc_samplename.SetValue(self.scriptdata.SampleName)
        self.tc_outpath.SetValue(self.scriptdata.OutputPath)
        self.ntc_exptime.SetValue(1.0)
        self.combo_Binning.SetStringSelection(str(1))

        self.UpdateStatusPanel('')

    # ----------------------------------------------------------------------
    # Read settings from TW_settings file
    def GetSettingsFromFile(self):

        self.HaveXradiaAPI = 1
        self.HaveMono = 1
        self.IonChamber_nobeam = -0.1
        self.XradiaMotorWait = 0
        self.MaxNumFilters = 16

        # Monochromator soft limits
        self.mono_use_limits = 0
        self.mono_l = 0
        self.mono_u = 10000

        # Display drive warning
        self.DisplayDriveWarning = 0
        self.LocalDrives = []

        # Sample XYZ above Rotation (Sample XYZ does move during a tomography)
        self.SampleXYZAboveRotation = 1

        # Refocusing with constant magnification for systems with a motorized detector stage.
        self.MotorizedDetector = 0

        # Energy axis (0 - SSRL SIS, 1 - Xradia generic axis)
        self.EnergyAxis = 0

        try:
            f = open(SettingsFileName, 'rt')
            for line in f:

                if 'USE_XRADIAAPI' in line:
                    slist = line.split('=')
                    value = ''.join(slist[1:])
                    self.HaveXradiaAPI = int(value)

                if 'USE_MONO' in line:
                    slist = line.split('=')
                    value = ''.join(slist[1:])
                    self.HaveMono = int(value)

                if 'MOTOR_WAIT' in line:
                    slist = line.split('=')
                    value = ''.join(slist[1:])
                    self.XradiaMotorWait = int(value)

                if 'IC_NOBEAM' in line:
                    slist = line.split('=')
                    value = ''.join(slist[1:])
                    self.IonChamber_nobeam = float(value)

                if 'MAX_NUM_FILTERS' in line:
                    slist = line.split('=')
                    value = ''.join(slist[1:])
                    self.MaxNumFilters = int(value)

                if 'MONO_USE_LIMITS' in line:
                    slist = line.split('=')
                    value = ''.join(slist[1:])
                    self.mono_use_limits = int(value)

                if 'MONO_L' in line:
                    slist = line.split('=')
                    value = ''.join(slist[1:])
                    self.mono_l = float(value)

                if 'MONO_U' in line:
                    slist = line.split('=')
                    value = ''.join(slist[1:])
                    self.mono_u = float(value)

                if 'DISPLAY_DRIVE_WARNING' in line:
                    slist = line.split('=')
                    value = ''.join(slist[1:])
                    self.DisplayDriveWarning = int(value)

                if 'LOCAL_DRIVES' in line:
                    slist = line.split('=')
                    value = ''.join(slist[1:])
                    for item in value.split(','):
                        self.LocalDrives.append(item.strip())

                if 'SAMPLEXYZ_ABOVE_ROTATION' in line:
                    slist = line.split('=')
                    value = ''.join(slist[1:])
                    self.SampleXYZAboveRotation = int(value)

                if 'MOTORIZED_DETECTOR' in line:
                    slist = line.split('=')
                    value = ''.join(slist[1:])
                    self.MotorizedDetector = int(value)

                if 'ENERGY_AXIS' in line:
                    slist = line.split('=')
                    value = ''.join(slist[1:])
                    self.EnergyAxis = int(value)

            f.close()
        except:
            print
            'WARNING - Could not read settings from TW_config.txt using defaults.'
            logging.warning('WARNING - Could not read settings from TW_config.txt using defaults.')

        if v:
            print
            'self.HaveXradiaAPI', self.HaveXradiaAPI
            print
            'self.HaveMono', self.HaveMono
            print
            'self.IonChamber_nobeam', self.IonChamber_nobeam
            print
            'XradiaMotorWait', self.XradiaMotorWait
        return

    # ----------------------------------------------------------------------
    def OutputScanInfoFile(self, filename, scriptdata, scriptlist):

        # If the script was loaded from a script file do not output a Info File
        if self.scriptdata.loadedfromscript == 1:
            print
            'Warning - Script was loaded from a script file, not able to output info file.'
            return

        path = scriptdata.OutputPath
        filepath = str(os.path.join(path, filename))

        f = open(filepath, 'w')

        f.write('VERSION 1\n')

        if (scriptdata.enable_energy == 0):
            f.write('ENERGY 0\n')
        else:
            f.write('ENERGY %1.0f\n' % (scriptdata.layer_eng))

        if (scriptdata.enable_tomo == 0):
            f.write('TOMO 0\n')
        else:
            f.write('TOMO %1.0f\n' % (scriptdata.layer_ang))

        if (scriptdata.enable_mosaic == 0):
            f.write('MOSAIC 0\n')
        else:
            f.write('MOSAIC %1.0f\n' % (scriptdata.layer_mos))

        if (scriptdata.enable_multi == 0):
            f.write('MULTIEXPOSURE 0\n')
        else:
            f.write('MULTIEXPOSURE %1.0f\n' % (scriptdata.layer_me))

        f.write('NREPEATSCAN %3.0f\n' % (scriptdata.NRepeatScan))
        f.write('WAITNSECS %3.0f\n' % (scriptdata.waitsecs))

        f.write('NEXPOSURES %3.0f\n' % (scriptdata.n_exposures))
        f.write('AVERAGEONTHEFLY %3.0f\n' % (scriptdata.averageotfly))

        f.write('REFNEXPOSURES %3.0f\n' % (scriptdata.RepeatRefExposure))
        f.write('REF4EVERYEXPOSURES %3.0f\n' % (scriptdata.ref4everycertainnumofexp))
        f.write('REFABBA %1.0f\n' % (scriptdata.ABBA))
        f.write('REFAVERAGEONTHEFLY %1.0f\n' % (scriptdata.refaverageotfly))

        f.write('MOSAICUP %3.0f\n' % (scriptdata.mosaic_up))
        f.write('MOSAICDOWN %3.0f\n' % (scriptdata.mosaic_down))
        f.write('MOSAICLEFT %3.0f\n' % (scriptdata.mosaic_left))
        f.write('MOSAICRIGHT %3.0f\n' % (scriptdata.mosaic_right))
        f.write('MOSAICOVERLAP %1.2f\n' % (scriptdata.mosaic_overlap))
        f.write('MOSAICCENTRALTILE %3.0f\n' % (scriptdata.mosaic_centraltile))

        f.write('FILES\n')

        for i in range(len(scriptlist)):

            if ';;' in scriptlist[i]:
                continue

            thiscommand = scriptlist[i].split(' ')

            if ((thiscommand[0] == 'collect') or (thiscommand[0] == 'collectaverage')
                    or (thiscommand[0] == 'collectnomalized') or (thiscommand[0] == 'collecttomo')
                    or (thiscommand[0] == 'collectnomalizedtomo')):
                filename = thiscommand[1]
                f.write(str(filename + '\n'))

        f.close()
        return

    # ----------------------------------------------------------------------
    def CreateDirAndScanInfo(self, scriptdata, scriptlist):
        # Create new directory to store the scan if it's a single scan
        # Get scan type name
        scanname = self.GetScanName(scriptdata.enable_energy, scriptdata.enable_tomo,
                                    scriptdata.enable_mosaic, scriptdata.layer_eng,
                                    scriptdata.layer_ang, scriptdata.loadedfromscript)

        date = datetime.datetime.now()
        thistime = str(datetime.datetime.now()).split(' ')
        hrmin = thistime[1].split(':')
        timestamp = '_' + date.strftime("%y%m%d") + '_' + hrmin[0] + hrmin[1]

        if (len(self.multiscriptlist) == 0):
            self.NewOutputPath = str(os.path.join(scriptdata.OutputPath, str(scriptdata.SampleName + '_' + scanname)))

            self.NewOutputPath = self.NewOutputPath + str(timestamp)
            if not os.path.exists(self.NewOutputPath):
                os.makedirs(self.NewOutputPath)
                if not os.path.exists(self.NewOutputPath):
                    print
                    'ERROR: did not find the output directory, and could not create a new output directory.'
                    logging.error(
                        'ERROR: did not find the output directory, and could not create a new output directory.')

            # Save scan info in a text file
            scaninfofilename = 'ScanInfo_%s%s.txt' % (scriptdata.SampleName, timestamp)

            self.OutputScanInfoFile(str(os.path.join(self.NewOutputPath, scaninfofilename)),
                                    scriptdata, scriptlist)

            # Save script into a text file
            scrdir, scrfilename = os.path.split(self.NewOutputPath)
            scrfilename = scrfilename[0:-12]
            scriptfilepath = str(os.path.join(self.NewOutputPath, str(scrfilename + '.txt')))
            self.SaveScript(scriptfilepath)

        else:
            # Create all directories for a multiscan
            self.multifolderlist = []

            if self.multiscan_repeat == None:
                self.multiscan_repeat == 1

            for ir in range(self.multiscan_repeat):

                for ims in range(len(self.multiscript)):

                    scriptdata = self.multiscript[ims]
                    scriptlist = self.multiscriptlist[ims]

                    scanname = self.GetScanName(scriptdata.enable_energy, scriptdata.enable_tomo,
                                                scriptdata.enable_mosaic, scriptdata.layer_eng,
                                                scriptdata.layer_ang, scriptdata.loadedfromscript)

                    scanname = '{0:02d}'.format(self.multiscan_foldernum) + scanname
                    self.multiscan_foldernum += 1

                    self.NewOutputPath = str(
                        os.path.join(scriptdata.OutputPath, str(scriptdata.SampleName + '_' + scanname)))

                    self.NewOutputPath = self.NewOutputPath + str(timestamp)
                    if not os.path.exists(self.NewOutputPath):
                        os.makedirs(self.NewOutputPath)
                        if not os.path.exists(self.NewOutputPath):
                            print
                            'ERROR: did not find the output directory, and could not create a new output directory.'
                            logging.error(
                                'ERROR: did not find the output directory, and could not create a new output directory.')

                    self.multifolderlist.append(self.NewOutputPath)

                    # Save scan info in a text file
                    scaninfofilename = 'ScanInfo_%s%s.txt' % (scriptdata.SampleName, timestamp)

                    self.OutputScanInfoFile(str(os.path.join(self.NewOutputPath, scaninfofilename)),
                                            scriptdata, scriptlist)

                    # Save individual scripts into text files
                    scrdir, scrfilename = os.path.split(self.NewOutputPath)
                    scrfilename = scrfilename[0:-12]
                    scriptfilepath = str(os.path.join(self.NewOutputPath, str(scrfilename + '.txt')))
                    f = open(str(scriptfilepath), 'w')
                    for i in range(len(scriptlist)):
                        f.write('{0:s}\n'.format(scriptlist[i]))
                    f.close()

    # ----------------------------------------------------------------------
    def CheckLocalDrive(self, drivepath):

        if self.DisplayDriveWarning == 1:
            if len(self.LocalDrives) > 0:
                driveused = os.path.dirname(drivepath)
                if driveused[0] not in self.LocalDrives:
                    wx.MessageBox("Warning - Data not saved to a local drive: " + str(self.LocalDrives))

    # ----------------------------------------------------------------------
    def GetScanName(self, enable_energy, enable_tomo, enable_mosaic, l1, l2, loadedfromscript=0):

        scanname = ''

        # 2D XANES
        if (enable_energy == 1) and (enable_tomo == 0) and (enable_mosaic == 0):
            scanname = 'XANES'
        # TOMO
        if (enable_energy == 0) and (enable_tomo == 1) and (enable_mosaic == 0):
            scanname = 'TOMO'
        # Mosaic
        if (enable_energy == 0) and (enable_tomo == 0) and (enable_mosaic == 1):
            scanname = 'MOSAIC'

            # TOMO XANES
        if (enable_energy == 1) and (enable_tomo == 1) and (enable_mosaic == 0) and (l1 == 1):
            scanname = 'TOMO-XANES'

        # XANES TOMO
        if (enable_energy == 1) and (enable_tomo == 1) and (enable_mosaic == 0) and (l2 == 1):
            scanname = 'XANES-TOMO'

        # Mosaic XANES
        if (enable_energy == 1) and (enable_tomo == 0) and (enable_mosaic == 1):
            scanname = 'MOSAIC-XANES'

            # Mosaic TOMO
        if (enable_energy == 0) and (enable_tomo == 1) and (enable_mosaic == 1):
            scanname = 'MOSAIC-TOMO'

        # Mosaic TOMO XANES
        if (enable_energy == 1) and (enable_tomo == 1) and (enable_mosaic == 1) and (l1 == 1):
            scanname = 'MOSAIC-TOMO-XANES'

        # Mosaic XANES TOMO
        if (enable_energy == 1) and (enable_tomo == 1) and (enable_mosaic == 1) and (l2 == 1):
            scanname = 'MOSAIC-XANES-TOMO'

        if (loadedfromscript == 1):
            scanname = 'SCRIPT'

        return scanname

    # ----------------------------------------------------------------------
    def OnSaveRecipe(self, evt):

        self.scriptdata.RefX = self.ntc_refx.GetValue()
        self.scriptdata.RefY = self.ntc_refy.GetValue()
        self.scriptdata.RefZ = self.ntc_refz.GetValue()
        self.scriptdata.useRefZ = self.cb_refz.GetValue()
        self.scriptdata.RefTheta = self.ntc_reft.GetValue()
        self.scriptdata.useRefTheta = self.cb_reft.GetValue()
        if self.cb_abba.GetValue():
            self.scriptdata.ABBA = 1
        else:
            self.scriptdata.ABBA = 0
        self.scriptdata.refrelativepos = self.cb_refrelpos.GetValue()

        # Scan Settings:
        self.scriptdata.Exptime = self.ntc_exptime.GetValue()
        self.scriptdata.Binning = 2 ** (int(self.combo_Binning.GetSelection()))
        self.scriptdata.SampleName = self.tc_samplename.GetValue()
        self.scriptdata.SampleName = self.scriptdata.SampleName.replace(" ", "")
        self.scriptdata.RefBinning = 2 ** (int(self.combo_RefBinning.GetSelection()))
        self.scriptdata.RefExpTime = self.ntc_refexptime.GetValue()
        self.scriptdata.waitsecs = self.ntc12.GetValue()
        self.scriptdata.pixelsize = self.ntc_pixelsize.GetValue()

        if self.cb_avfly.GetValue():
            self.scriptdata.averageotfly = 1
        else:
            self.scriptdata.averageotfly = 0

        if self.cb_refavfly.GetValue():
            self.scriptdata.refaverageotfly = 1
        else:
            self.scriptdata.refaverageotfly = 0

        if self.cb_ioncham.GetValue():
            self.scriptdata.useionchamber = 1
        else:
            self.scriptdata.useionchamber = 0

        self.CheckTomoInerLoop()

        self.scriptdata.wait_energy = self.ntc_ms_waiteng.GetValue()
        self.scriptdata.wait_tomo = self.ntc_ms_waittomo.GetValue()
        self.scriptdata.wait_mosaic = self.ntc_ms_waitmos.GetValue()
        self.scriptdata.wait_multi = self.ntc_ms_waitmulti.GetValue()

        wildcard = "PTWZ files (*.ptwz)|*.ptwz"
        dialog = wx.FileDialog(None, "Please select a ptwz file",
                               wildcard=wildcard,
                               style=wx.SAVE)
        if dialog.ShowModal() == wx.ID_OK:
            filepath = dialog.GetPath()
            filename = dialog.GetFilename()
            filedir = dialog.GetDirectory()

            pickleFile = open(filepath, 'wb')
            fversion = 1.0
            pickle.dump(fversion, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.enable_energy, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.enable_tomo, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.enable_mosaic, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.enable_multi, pickleFile, pickle.HIGHEST_PROTOCOL)

            pickle.dump(self.scriptdata.NRepeatScan, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.waitsecs, pickleFile, pickle.HIGHEST_PROTOCOL)

            pickle.dump(self.scriptdata.layer_eng, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.layer_ang, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.layer_mos, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.layer_me, pickleFile, pickle.HIGHEST_PROTOCOL)

            pickle.dump(self.scriptdata.singletomo, pickleFile, pickle.HIGHEST_PROTOCOL)

            pickle.dump(self.scriptdata.wait_energy, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.wait_tomo, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.wait_mosaic, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.wait_multi, pickleFile, pickle.HIGHEST_PROTOCOL)

            pickle.dump(self.scriptdata.RegionParams, pickleFile, pickle.HIGHEST_PROTOCOL)
            neng = len(self.scriptdata.EngRegionlist)
            pickle.dump(neng, pickleFile, pickle.HIGHEST_PROTOCOL)
            for i in range(neng):
                pickle.dump(self.scriptdata.EngRegionlist[i], pickleFile, pickle.HIGHEST_PROTOCOL)
            nzp = len(self.scriptdata.EPlist)
            pickle.dump(nzp, pickleFile, pickle.HIGHEST_PROTOCOL)
            for i in range(nzp):
                pickle.dump(self.scriptdata.EPlist[i], pickleFile, pickle.HIGHEST_PROTOCOL)

            pickle.dump(self.scriptdata.n_exposures, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.averageotfly, pickleFile, pickle.HIGHEST_PROTOCOL)

            pickle.dump(self.scriptdata.Exptime, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.useionchamber, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.Binning, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.SampleName, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.ScriptName, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.OutputPath, pickleFile, pickle.HIGHEST_PROTOCOL)

            pickle.dump(self.scriptdata.ZP_EngList, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.ZP_xList, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.ZP_yList, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.ZP_zList, pickleFile, pickle.HIGHEST_PROTOCOL)

            pickle.dump(self.scriptdata.RepeatRefExposure, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.ref4everycertainnumofexp, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.RefX, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.RefY, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.RefZ, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.useRefZ, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.RefTheta, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.useRefTheta, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.RefBinning, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.RefExpTime, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.ABBA, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.refrelativepos, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.refaverageotfly, pickleFile, pickle.HIGHEST_PROTOCOL)

            pickle.dump(self.scriptdata.XCorr, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.ZCorr, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.SampleXinZeroDegree, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.SampleYinZeroDegree, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.SampleZinZeroDegree, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.FOV, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.samplexs, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.sampleys, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.samplezs, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.samplets, pickleFile, pickle.HIGHEST_PROTOCOL)

            pickle.dump(self.scriptdata.mosaic_left, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.mosaic_right, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.mosaic_up, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.mosaic_down, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.mosaic_overlap, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.mosaic_centraltile, pickleFile, pickle.HIGHEST_PROTOCOL)

            pickle.dump(self.scriptdata.pixelsize, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.width, pickleFile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(self.scriptdata.height, pickleFile, pickle.HIGHEST_PROTOCOL)

            pickleFile.close()

    # ----------------------------------------------------------------------
    def OnLoadRecipe(self, evt):

        wildcard = "PTWZ files (*.ptwz)|*.ptwz"
        dialog = wx.FileDialog(None, "Please select a ptwz file",
                               wildcard=wildcard,
                               style=wx.OPEN)

        if dialog.ShowModal() == wx.ID_OK:
            filepath = dialog.GetPath()
        else:
            print
            'Warning - Could not read recipe file.'
            logging.warning('Warning - Could not read recipe file.')
            return

        pickleFile = open(filepath, 'rb')
        fversion = pickle.load(pickleFile)

        self.scriptdata.enable_energy = pickle.load(pickleFile)
        self.scriptdata.enable_tomo = pickle.load(pickleFile)
        self.scriptdata.enable_mosaic = pickle.load(pickleFile)
        self.scriptdata.enable_multi = pickle.load(pickleFile)

        self.scriptdata.NRepeatScan = pickle.load(pickleFile)
        self.scriptdata.waitsecs = pickle.load(pickleFile)

        self.scriptdata.layer_eng = pickle.load(pickleFile)
        self.scriptdata.layer_ang = pickle.load(pickleFile)
        self.scriptdata.layer_mos = pickle.load(pickleFile)
        self.scriptdata.layer_me = pickle.load(pickleFile)

        self.scriptdata.singletomo = pickle.load(pickleFile)

        self.scriptdata.wait_energy = pickle.load(pickleFile)
        self.scriptdata.wait_tomo = pickle.load(pickleFile)
        self.scriptdata.wait_mosaic = pickle.load(pickleFile)
        self.scriptdata.wait_multi = pickle.load(pickleFile)

        self.scriptdata.RegionParams = pickle.load(pickleFile)

        self.scriptdata.EngRegionlist = []
        neng = pickle.load(pickleFile)
        for i in range(neng):
            self.scriptdata.EngRegionlist.append(pickle.load(pickleFile))

        nzp = pickle.load(pickleFile)
        for i in range(nzp):
            self.scriptdata.EPlist.append(pickle.load(pickleFile))

        self.scriptdata.n_exposures = pickle.load(pickleFile)
        self.scriptdata.averageotfly = pickle.load(pickleFile)

        self.scriptdata.Exptime = pickle.load(pickleFile)
        self.scriptdata.useionchamber = pickle.load(pickleFile)
        self.scriptdata.Binning = pickle.load(pickleFile)
        self.scriptdata.SampleName = pickle.load(pickleFile)
        self.scriptdata.ScriptName = pickle.load(pickleFile)
        self.scriptdata.OutputPath = pickle.load(pickleFile)

        self.scriptdata.ZP_EngList = pickle.load(pickleFile)
        self.scriptdata.ZP_xList = pickle.load(pickleFile)
        self.scriptdata.ZP_yList = pickle.load(pickleFile)
        self.scriptdata.ZP_zList = pickle.load(pickleFile)

        self.scriptdata.RepeatRefExposure = pickle.load(pickleFile)
        self.scriptdata.ref4everycertainnumofexp = pickle.load(pickleFile)
        self.scriptdata.RefX = pickle.load(pickleFile)
        self.scriptdata.RefY = pickle.load(pickleFile)
        self.scriptdata.RefZ = pickle.load(pickleFile)
        self.scriptdata.useRefZ = pickle.load(pickleFile)
        self.scriptdata.RefTheta = pickle.load(pickleFile)
        self.scriptdata.useRefTheta = pickle.load(pickleFile)
        self.scriptdata.RefBinning = pickle.load(pickleFile)
        self.scriptdata.RefExpTime = pickle.load(pickleFile)
        self.scriptdata.ABBA = pickle.load(pickleFile)
        self.scriptdata.refrelativepos = pickle.load(pickleFile)
        self.scriptdata.refaverageotfly = pickle.load(pickleFile)

        self.scriptdata.XCorr = pickle.load(pickleFile)
        self.scriptdata.ZCorr = pickle.load(pickleFile)
        self.scriptdata.SampleXinZeroDegree = pickle.load(pickleFile)
        self.scriptdata.SampleYinZeroDegree = pickle.load(pickleFile)
        self.scriptdata.SampleZinZeroDegree = pickle.load(pickleFile)
        self.scriptdata.FOV = pickle.load(pickleFile)
        self.scriptdata.samplexs = pickle.load(pickleFile)
        self.scriptdata.sampleys = pickle.load(pickleFile)
        self.scriptdata.samplezs = pickle.load(pickleFile)
        self.scriptdata.samplets = pickle.load(pickleFile)

        self.scriptdata.mosaic_left = pickle.load(pickleFile)
        self.scriptdata.mosaic_right = pickle.load(pickleFile)
        self.scriptdata.mosaic_up = pickle.load(pickleFile)
        self.scriptdata.mosaic_down = pickle.load(pickleFile)
        self.scriptdata.mosaic_overlap = pickle.load(pickleFile)
        self.scriptdata.mosaic_centraltile = pickle.load(pickleFile)

        self.scriptdata.pixelsize = pickle.load(pickleFile)
        self.scriptdata.width = pickle.load(pickleFile)
        self.scriptdata.height = pickle.load(pickleFile)

        pickleFile.close()

        self.UpdateWidgets()
        self.ShowEnabledTC()

    # ----------------------------------------------------------------------
    def OnLoadScript(self, evt):

        wildcard = "TXT files (*.txt)|*.txt"
        dialog = wx.FileDialog(None, "Please select a script file",
                               wildcard=wildcard,
                               style=wx.OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            filepath = dialog.GetPath()
            filedir = dialog.GetDirectory()
        else:
            print
            'ERROR - cannot open file.'
            return

        if v: print
        'Loading script from a file.'

        f = open(str(filepath), 'r')

        scriptlist = []

        for line in f:
            line = line.rstrip('\r\n')

            if line == '':
                continue

            if line.startswith(';'):
                comment = line.replace(';', '')
                script = ';;; %s' % (comment)
                scriptlist.append(script)
                continue

            thiscommand = line.split(' ')

            if thiscommand[0] == 'setexp':
                exposure = float(thiscommand[1])
                script = '%s %.3f' % ('setexp', exposure)
                scriptlist.append(script)
                if v: print
                'set exposure ', exposure

            elif thiscommand[0] == 'setbinning':
                binning = int(thiscommand[1])
                script = '%s %.0f' % ('setbinning', binning)
                scriptlist.append(script)
                if v: print
                'set binning ', binning

            elif thiscommand[0] == 'setfilter':
                filter = int(thiscommand[1])
                script = '%s %.0f' % ('setfilter', filter)
                scriptlist.append(script)
                if v: print
                'set setfilter ', filter


            elif thiscommand[0] == 'sete':
                eng = float(thiscommand[1])
                script = 'sete %.3f' % (eng)
                scriptlist.append(script)
                if v: print
                'set energy ', eng


            elif thiscommand[0] == 'moveto':
                axisname = thiscommand[1]
                position = float(thiscommand[2])
                if axisname == 'E':
                    script = 'sete %.3f' % (position)
                    scriptlist.append(script)
                    if v: print
                    'move energy to %f' % (position)
                else:
                    script = 'moveto %s %.3f' % (axisname, position)
                    scriptlist.append(script)
                    if v: print
                    'move motor %s to %f' % (axisname, position)


            elif thiscommand[0] == 'moveby':
                axisname = thiscommand[1]
                position = float(thiscommand[2])
                script = 'moveby %s %.3f' % (axisname, position)
                scriptlist.append(script)
                if v: print
                'move motor %s by %f' % (axisname, position)


            elif thiscommand[0] == 'collect':
                filename = thiscommand[1]
                script = '%s %s' % ('collect', filename)
                scriptlist.append(script)
                if v: print
                'collect ', filename

            elif thiscommand[0] == 'collectxrf':
                filename = thiscommand[1]
                script = '%s %s' % ('collectxrf', filename)
                scriptlist.append(script)
                if v: print
                'collectxrf ', filename

            elif thiscommand[0] == 'collectaverage':
                filename = thiscommand[1]
                totalnimages = thiscommand[2]
                script = '%s %s %s' % ('collectaverage', filename, str(totalnimages))
                scriptlist.append(script)
                if v: print
                'collectaverage ', filename, totalnimages

            elif thiscommand[0] == 'collectnormalized':
                filename = thiscommand[1]
                script = '%s %s' % ('collectnormalized', filename)
                scriptlist.append(script)
                if v: print
                'collectnormalized ', filename

            elif thiscommand[0] == 'collecttomo':
                filename = thiscommand[1]
                totalnimages = thiscommand[2]
                startangle = thiscommand[3]
                endangle = thiscommand[4]
                script = '%s %s %s %s %s' % ('collecttomo', filename, totalnimages, startangle, endangle)
                scriptlist.append(script)
                if v: print
                'collecttomo ', filename, totalnimages, startangle, endangle

            elif thiscommand[0] == 'collectnormalizedtomo':
                filename = thiscommand[1]
                totalnimages = thiscommand[2]
                startangle = thiscommand[3]
                endangle = thiscommand[4]
                script = '%s %s %s %s %s' % ('collectnormalizedtomo', filename, totalnimages, startangle, endangle)
                scriptlist.append(script)
                if v: print
                'collectnormalizedtomo ', filename, totalnimages, startangle, endangle

            elif thiscommand[0] == 'wait':
                waittimesecs = float(thiscommand[1])
                script = 'wait %f' % (waittimesecs)
                scriptlist.append(script)
                if v: print
                'wait %f' % (waittimesecs)

            elif thiscommand[0] == 'recordioncurrent':
                script = 'recordioncurrent'
                scriptlist.append(script)
                if v: print
                'recordioncurrent'

            elif thiscommand[0] == 'setreadout':
                readout = int(thiscommand[1])
                script = 'setreadout %d' % (readout)
                scriptlist.append(script)
                if v: print
                'setreadout %d' % (readout)


            elif thiscommand[0] == 'shiftobjective':
                nobjective = int(thiscommand[1])
                script = 'shiftobjective %d' % (nobjective)
                scriptlist.append(script)
                if v: print
                'shiftobjective %d' % (nobjective)

            elif thiscommand[0] == 'runshortcut':
                command = ' '.join(thiscommand[1:])
                script = 'runshortcut %s' % (command)
                scriptlist.append(script)
                if v: print
                'runshortcut %d' % (command)

            else:
                print
                'ERROR - unrecognized command in the script.'
                logging.error('ERROR - unrecognized command in the script.')

        f.close()
        self.scriptlist = scriptlist

        self.scriptdata.OutputPath = filedir
        self.scriptdata.loadedfromscript = 1

        self.tc_outpath.SetValue(self.scriptdata.OutputPath)

        self.CheckLocalDrive(self.scriptdata.OutputPath)

        self.ShowScriptList()
        self.UpdateButtonWidgets()

    # ----------------------------------------------------------------------
    def OnSaveScript(self, evt):

        wildcard = "TXT files (*.txt)|*.txt"
        dialog = wx.FileDialog(None, "Save as ...",
                               wildcard=wildcard,
                               style=wx.SAVE)
        if dialog.ShowModal() == wx.ID_OK:
            filepath = dialog.GetPath()
            filedir = dialog.GetDirectory()
        else:
            print
            'ERROR - could not save file.'
            return

        if v: print
        'Save script to a file.'

        self.SaveScript(filepath)

    # ----------------------------------------------------------------------
    def SaveScript(self, filepath):

        scriptlist = []
        self.scriptdata.SampleName = self.tc_samplename.GetValue()
        self.scriptdata.SampleName = self.scriptdata.SampleName.replace(" ", "")
        if (len(self.multiscriptlist) > 0):
            self.cb_multiscript.SetValue(True)
            self.ShowScriptList()
            for ii in range(len(self.multiscriptlist)):
                for jj in range(len(self.multiscriptlist[ii])):
                    scriptlist.append(self.multiscriptlist[ii][jj])
        else:
            scriptlist = self.scriptlist

        f = open(str(filepath), 'w')

        for i in range(len(scriptlist)):
            f.write('{0:s}\n'.format(scriptlist[i]))

        f.close()

    # ----------------------------------------------------------------------
    def OnClose(self, evt):

        if self.mono: self.mono.Close()
        if self.microscope:

            if self.TurnAgitatorOff:
                self.microscope.TurnOffAgitator()

            self.microscope.Shutdown()

        self.Destroy()


""" ------------------------------------------------------------------------------------------------"""


def main():
    app = wx.App(False)

    (ssx, ssy) = wx.DisplaySize()

    Winsizex = 1200
    Winsizey = 1010

    if (Winsizex > ssx):
        Winsizex = ssx

    if (Winsizey > ssy):
        Winsizey = ssy - 70

    title = 'XMFlex  v:%s' % (__version__)
    frame = MainFrame(None, -1, title, Winsizex, Winsizey)

    if (ssy <= 1024):
        frame.SetPosition(wx.Point(-1, 0))

    frame.Show()

    app.MainLoop()


# ----------------------------------------------------------------------
if __name__ == '__main__':
    main()



