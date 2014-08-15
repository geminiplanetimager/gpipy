#!/usr/bin/env python
#/opt/python2.7.2/bin/python
import os
import sys

import wx
from ObjectListView import ObjectListView, GroupListView, ColumnDefn

import gpi_pipeline
import numpy as np

import logging
_log = logging.getLogger('gpicaldb')


 
########################################################################
class CalDBViewer(wx.Frame):
    """ Graphical view tool for examining the contents of a GPI
    calibration files directory

    """
 
    #----------------------------------------------------------------------
    def __init__(self, parent=None, id=-1, size=(900,600)):
        """ Viewer window for GPI calibration DB objects 
        Based on objectlistview's simpleExample1.py
        
        """
        wx.Frame.__init__(self, parent=parent, id=id, size=size, title="GPI Calibration DB Viewer")

        self.InitModel()
        self.InitWidgets()
        self.InitObjectListView()

    def InitModel(self):
        """ Set up the database connection and model objects
        """

        _log.info("Querying Calibration DB for files...")
        self.db = gpi_pipeline.GPICalDB()
        _log.info("Reading in FITS headers of calibration files...")
        self.files = self.db.select_calfiles(return_objects=True)
        self.selected = np.arange(len(self.files))

    def InitWidgets(self):
        """
        create widgets 
        """
        # Add a panel so it looks the correct on all platforms
        topPanel = wx.Panel(self, wx.ID_ANY)
        topSizer = wx.BoxSizer(wx.VERTICAL)

        #--- top panel with info and controls
        panel1 = wx.Panel(topPanel, wx.ID_ANY)
        sizer1 = wx.BoxSizer(wx.VERTICAL)
        text = wx.StaticText(panel1, label='GPI Calibration DB: '+self.db.caldir)
        sizer1.Add(text, 1, wx.ALL|wx.EXPAND, border=5)

        panel1a = wx.Panel(panel1)
        sizer1a = wx.BoxSizer(wx.HORIZONTAL)

        text1 = wx.StaticText(panel1a, label='File Type: ')
        choices = list(self.db.unique_caltypes)
        choices.insert(0,'--Any--')
        self.cb_type = wx.ComboBox(panel1a, -1,value=choices[0], choices=choices, style=wx.CB_DROPDOWN|wx.CB_READONLY)
        self.cb_type.Bind(wx.EVT_COMBOBOX, self.OnSelect)

        text2 = wx.StaticText(panel1a, label='     Prism: ')
        choices = list(self.db.unique_dispersers)
        choices.insert(0,'--Any--')
        self.cb_prism = wx.ComboBox(panel1a, -1,value=choices[0], choices=choices, style=wx.CB_DROPDOWN|wx.CB_READONLY)
        self.cb_prism.Bind(wx.EVT_COMBOBOX, self.OnSelect)

        text3 = wx.StaticText(panel1a, label='     IFS Filter: ')
        choices = ['--Any--','Y','J','H','K1','K2']
        self.cb_filter = wx.ComboBox(panel1a, -1,value=choices[0], choices=choices, style=wx.CB_DROPDOWN|wx.CB_READONLY)
        self.cb_filter.Bind(wx.EVT_COMBOBOX, self.OnSelect)

        text4 = wx.StaticText(panel1a, label='     ITime: ')
        choices = list(set(self.db.table['ITIME']))
        choices.sort()
        choices =  ['%.1f'% val for val in choices]
        choices.insert(0,'--Any--')
        self.cb_itime = wx.ComboBox(panel1a, -1,value=choices[0], choices=choices, style=wx.CB_DROPDOWN|wx.CB_READONLY)
        self.cb_itime.Bind(wx.EVT_COMBOBOX, self.OnSelect)



        for item in [text1,self.cb_type,text2,self.cb_prism,text3,self.cb_filter,text4,self.cb_itime]:
            sizer1a.Add(item,1,wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER_HORIZONTAL, border=0)
 
        sizer1.Add(panel1a, 1, wx.ALL|wx.EXPAND, border=5)
        panel1a.SetSizerAndFit(sizer1a)



        panel1.SetSizerAndFit(sizer1)


        #-- second panel with ObjectListView
        panel2 = wx.Panel(topPanel, wx.ID_ANY)
        sizer2 = wx.BoxSizer(wx.VERTICAL)

        #self.myOlv = ObjectListView(panel2, -1, style=wx.LC_REPORT|wx.SUNKEN_BORDER)
        self.myOlv = GroupListView(panel2, -1, style=wx.LC_REPORT|wx.SUNKEN_BORDER)
        sizer2.Add(self.myOlv, 1, wx.ALL|wx.EXPAND, 5)
        panel2.SetSizerAndFit(sizer2)

        #-- third panel, button bar
        panel3 = wx.Panel(topPanel, wx.ID_ANY)
        sizer3 = wx.BoxSizer(wx.HORIZONTAL)
        btn = wx.Button(panel3, label="Refresh")
        btn.Bind(wx.EVT_BUTTON, self.OnRefresh)
        btnToggleGrouping = wx.Button(panel3, label="Toggle Grouping")
        btnToggleGrouping.Bind(wx.EVT_BUTTON, self.OnToggleGrouping)
        #btnGroupTrue = wx.Button(panel3, label="Toggle Grouping")
        #btnGroupTrue.Bind(wx.EVT_BUTTON, self.OnGroupTrue)
  
        sizer3.Add(btn, 0, wx.ALL|wx.CENTER, border=10)
        sizer3.Add(btnToggleGrouping, 0, wx.ALL|wx.CENTER, border=10)
        #sizer3.Add(btnGroupTrue, 0, wx.ALL|wx.CENTER, border=10)
        panel3.SetSizerAndFit(sizer3)

        topSizer.Add(panel1, 0, wx.ALL|wx.EXPAND, border=5)
        topSizer.Add(panel2, 1, wx.ALL|wx.EXPAND, border=5)
        topSizer.Add(panel3, 0, wx.ALL|wx.EXPAND, border=10)

        topPanel.SetSizerAndFit(topSizer)


    def InitObjectListView(self):
        """ Populate the Object List View widget
        """

        self.myOlv.useExpansionColumn=True # add extra initial column for expand/contract group widgets
        self.myOlv.SetColumns([
            ColumnDefn("Filename", "left", 220, "name", stringConverter=os.path.basename),
            ColumnDefn("Type", "left", 130, "filetype"),
            ColumnDefn("Prism", "left", 50, "prism"),
            ColumnDefn("IFS Filter", "left", 50, "filter"),
            ColumnDefn("Lyot Mask", "left", 70, "lyot"),
            ColumnDefn("Apodizer", "left", 70, "apodizer"),
            ColumnDefn("ITime", "right", 70, "itime",stringConverter="%.1f" ),
            ColumnDefn("Readmode", "right", 70, "readmode" ),
            #ColumnDefn("Type", "left", 100, "filetype", stringConverter="%.1f"),
            #ColumnDefn("Last Played", "left", 100, "lastPlayed", stringConverter="%d-%m-%Y"),
            ColumnDefn("Date", "right", 100, "dateobs"),
            ColumnDefn("# combined", "right", 70, "ncombined"),
            ColumnDefn("DRP version", "right", 70, "version")
        ], repopulate=False)
        self.myOlv.SetObjects(self.files)
        self.myOlv.SetSortColumn(self.myOlv.columns[2], resortNow=True)
        # hack work-around for keys not being set correctly:
        self.myOlv.SetShowGroups(False)
        self.myOlv.SetShowGroups(True)



    #----------------------------------------------------------------------
    def OnSelect(self, event=None):
        """ Update which files are displayed, in response to the combo box selections """
        caltype = self.cb_type.GetValue()
        if caltype == '--Any--': caltype=None
        prism = self.cb_prism.GetValue()
        if prism == '--Any--': prism=None
        filter = self.cb_filter.GetValue()
        if filter == '--Any--': filter=None
        itime = self.cb_itime.GetValue()
        if itime == '--Any--': itime=None


        #selection = self.db.select_calfiles(return_indices=True, caltype='Dark File')
        # for debugging:
        self.selection= self.db.select_calfiles(caltype=caltype, ifsfilt=filter, itime=itime, disperser=prism, 
                return_indices=True)
        self.myOlv.SetObjects([self.files[i] for i in self.selection])
        self.myOlv.RepopulateList()
        #self.myOlv.SetShowGroups( not self.myOlv.GetShowGroups())
        #self.myOlv.SetShowGroups( not self.myOlv.GetShowGroups())

    def OnRefresh(self, event):
        """ Update all calibration DB information from disk """
        self.InitModel()
        self.OnSelect()

    def OnToggleGrouping(self, event):
        self.myOlv.SetShowGroups( not self.myOlv.GetShowGroups())

    def OnGroupTrue(self, event):
         self.myOlv.SetShowGroups( True)


#----------------------------------------------------------------------
# Run the program
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    app = wx.App(False)
    frame = CalDBViewer()
    frame.Show()
    app.MainLoop()
