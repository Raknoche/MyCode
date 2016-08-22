from tkinter import *
from tkinter import messagebox
import pymysql as mdb
import pandas as pd
import sys
import os
import re
import datetime
import time
from time import sleep
import numpy as np
import collections
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


class PlotInterface(object):
    def __init__(self, df, EvT_color_set, EvT_size_set):
        self.df = df
        self.ax = plt.gca()
        self.rect = Rectangle((0.1,0.1), 0, 0, fill=False)
        self.rect.set_visible(False)
        self.ax.add_patch(self.rect)
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.press = False
        self.sub_df = None
        self.EvT_color_set = EvT_color_set
        self.EvT_size_set = EvT_size_set
        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax.figure.canvas.mpl_connect('motion_notify_event',self.on_motion)
        self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        
    def on_press(self, event):
        if event.inaxes:
            
            self.press = True
            self.x0 = event.xdata
            self.y0 = event.ydata
            self.rect.set_width(0)
            self.rect.set_height(0)
            self.rect.set_xy((self.x0,self.y0))
            self.rect.set_visible(True)
            #Save the plot background
            self.background = self.ax.figure.canvas.copy_from_bbox(self.ax.bbox)
            


    def on_release(self, event):
        if event.inaxes:
            self.press = False
            self.x1 = event.xdata
            self.y1 = event.ydata
            self.rect.set_visible(False)
            self.ax.figure.canvas.draw()
            self.min_x = min(self.x0,self.x1)
            self.max_x = max(self.x0,self.x1)
            self.min_y = min(self.y0,self.y1)
            self.max_y = max(self.y0,self.y1)
            
            if self.sub_df is not None:
                self.sub_df = self.sub_df.append(self.df[ (self.df['timestamp'] > self.min_x) & (self.df['timestamp'] < self.max_x) & \
                                 (self.df['amplitude'] > self.min_y) & (self.df['amplitude'] < self.max_y)])
            else:
                self.sub_df = self.df[ (self.df['timestamp'] > self.min_x) & (self.df['timestamp'] < self.max_x) & \
                                 (self.df['amplitude'] > self.min_y) & (self.df['amplitude'] < self.max_y)]            
            
            plt.scatter(self.sub_df['timestamp'],self.sub_df['amplitude'],color='blue',s=self.EvT_size_set)  
            self.ax.figure.canvas.draw()

            print('Selected (%d < timestamp < %d) and (%d < amplitude < %d)' % (self.min_x, self.max_x, self.min_y, self.max_y))
            
    def on_motion(self, event):
        if self.press is False: return
        if event.inaxes != self.rect.axes: return
        self.x1 = event.xdata
        self.y1 = event.ydata
        self.rect.set_width(self.x1 - self.x0)
        self.rect.set_height(self.y1 - self.y0)
        self.rect.set_xy((self.x0, self.y0))
        self.ax.figure.canvas.restore_region(self.background)
        self.ax.draw_artist(self.rect)
        self.ax.figure.canvas.blit(self.rect.clipbox)
        
        
        
class App(object):
    def __init__(self, window, ip, user, pswd, db):
        self.window = window
        self.ip = ip
        self.user = user
        self.pswd = pswd
        self.db = db
        
        window.wm_title("Radon Database Interface")
        
        self.current_row=0
        
        '''Query Panel'''
        #Run tag field
        self.runtag_label = Label (window, text= "Data Run Tag: ")
        self.runtag_label.grid(row=self.current_row,column=0)
        self.runtag_text = StringVar()
        Entry(window, textvariable=self.runtag_text).grid(row=self.current_row,column=1)
        self.current_row += 1
        
        #Start time field
        self.start_t_label = Label (window, text= "Start Time (YYYY-MM-DD HH:MM:SS): ")
        self.start_t_label.grid(row=self.current_row,column=0)          
        self.start_t_text = StringVar()
        Entry(window, textvariable=self.start_t_text).grid(row=self.current_row,column=1)
        self.current_row += 1
        
        #End time field
        self.end_t_label = Label (window, text= "End Time (YYYY-MM-DD HH:MM:SS): ")
        self.end_t_label.grid(row=self.current_row, column=0)          
        self.end_t_text = StringVar()
        Entry(window, textvariable=self.end_t_text).grid(row=self.current_row,column=1)
        self.current_row += 1
        
        #Limit field
        self.limit_label = Label (window, text= "Query Limit (Integer)")
        self.limit_label.grid(row=self.current_row, column=0)          
        self.limit_text = StringVar()
        Entry(window, textvariable=self.limit_text).grid(row=self.current_row,column=1)
        self.current_row += 1
        
        #Send Query button
        self.query_button = StringVar()
        self.query_button.set("Send Query")        
        Button(window, textvariable=self.query_button, command=self.sendQuery).grid(row=self.current_row,column=0,columnspan=2)
        self.current_row += 1
        
        #Query status
        self.querystatus = Label(window, text='')
        self.querystatus.grid(row=self.current_row,column=0,columnspan=2)
        self.current_row += 1
    
        '''Line to separate segments'''
        canvas = Canvas(master=window, width=500, height=40)
        canvas.create_line(0, 20, 500, 20, fill="black")
        canvas.grid(row=self.current_row,column=0,columnspan=2)
        self.current_row += 1
        
        '''Plot Panel'''
        
        #Drop down menu for time units
        Label(window, text="Time Units").grid(row = self.current_row, column = 0, columnspan=2)
        self.current_row += 1
        
        self.timevar = StringVar(window)
        # Use dictionary to map time unit to second -> unit conversion factor
        self.choices = collections.OrderedDict()
        self.choices['Seconds'] = 1
        self.choices['Minutes'] = 1/60
        self.choices['Hours'] = 1/3600
        self.choices['Days'] = 1/86400
        self.choices['Weeks'] = 1/604800

        self.timevar.set('Seconds')
        self.unit = self.timevar.get()
        self.unit_conv = self.choices[self.unit]
        self.option = OptionMenu(window, self.timevar, *self.choices, command = self.updateTimeUnit)
            
        self.option.grid(row = self.current_row, column=0, columnspan=2)
        #self.unit_conv = self.choices[self.timevar.get()]
        self.current_row += 1
 

        #EvT Plot 
        self.EvT_color = Label (window, text= "Data Color:")
        self.EvT_color.grid(row=self.current_row,column=0)  
        self.EvT_size = Label (window, text= "Data Size:")
        self.EvT_size.grid(row=self.current_row,column=1) 
        self.current_row += 1
        
        self.EvT_color_text = StringVar(value='black')
        EvT_Entry = Entry(window,textvariable=self.EvT_color_text).grid(row=self.current_row,column=0)        
        self.EvT_size_text = IntVar(value=10)
        Entry(window, textvariable=self.EvT_size_text).grid(row=self.current_row,column=1)
        self.current_row += 1        
     
        self.EvT_button = Button (window, text="Energy v Time Plot", command=self.EvT_plot)
        self.EvT_button.grid(row=self.current_row,column=0,columnspan=2)
        self.current_row += 1        

        
        
        #Energy Histogram Plot
        self.EHist_color = Label (window, text= "Hist Color:")
        self.EHist_color.grid(row=self.current_row,column=0) 
        self.EHist_bin = Label (window, text= "Hist Binning:")
        self.EHist_bin.grid(row=self.current_row,column=1)    
        self.current_row += 1  
            
        self.EHist_color_text = StringVar(value='black')
        Entry(window, textvariable=self.EHist_color_text).grid(row=self.current_row,column=0)         
        self.EHist_bin_text = IntVar(value=20)
        Entry(window, textvariable=self.EHist_bin_text).grid(row=self.current_row,column=1)
        self.current_row += 1  
        
        self.EHist_button = Button (window, text="Energy Histogram", command=self.EHist_plot)
        self.EHist_button.grid(row=self.current_row,column=0,columnspan=2)
        self.current_row += 1 
        
        
        #Time Histogram Settings (replace with error bar of count v time)
        self.THist_color = Label (window, text= "Hist Color:")
        self.THist_color.grid(row=self.current_row,column=0)  
        self.THist_bin = Label (window, text= "Hist Binning:")
        self.THist_bin.grid(row=self.current_row,column=1)  
        self.current_row += 1 
        
        self.THist_color_text = StringVar(value='black')
        Entry(window, textvariable=self.THist_color_text).grid(row=self.current_row,column=0)
        self.THist_bin_text = IntVar(value=20)
        Entry(window, textvariable=self.THist_bin_text).grid(row=self.current_row,column=1)        
        self.current_row += 1 
        
        self.THist_button = Button (window, text="Time Histogram", command=self.THist_plot)
        self.THist_button.grid(row=self.current_row,column=0,columnspan=2)
        self.current_row += 1  
        
        '''Line to separate segments'''
        canvas = Canvas(master=window, width=500, height=40)
        canvas.create_line(0, 20, 500, 20, fill="black")
        canvas.grid(row=self.current_row,column=0,columnspan=2)
        self.current_row += 1
           
            
        '''Event Selection Panel'''    
        self.select_button = Button (window, text="Select Data from EvT", command=self.EvT_select)
        self.select_button.grid(row=self.current_row,column=0,columnspan=1)
        
        self.cancel_button = Button (window, text="Cancel Selection", command=self.EvT_select_cancel)
        self.cancel_button.grid(row=self.current_row,column=1,columnspan=1)
        self.current_row += 1

        self.delete_button = Button (window, text="Delete Selection", command=self.EvT_delete)
        self.delete_button.grid(row=self.current_row,column=0,columnspan=1)

        self.select_only_button = Button (window, text="Delete All But Selection", command=self.EvT_select_only)
        self.select_only_button.grid(row=self.current_row,column=1,columnspan=1)
        self.current_row += 1
         
            
        '''Line to separate segments'''
        canvas = Canvas(master=window, width=500, height=40)
        canvas.create_line(0, 20, 500, 20, fill="black")
        canvas.grid(row=self.current_row,column=0,columnspan=2)
        self.current_row += 1

            
        '''Saving Panel'''  
        self.save_label = Label (window, text= "Save Location: ")
        self.save_label.grid(row=self.current_row, column=0, columnspan=2)          
        self.current_row += 1 
               
        self.save_text = StringVar()
        Entry(window, textvariable=self.save_text).grid(row=self.current_row,column=0, columnspan=2)
        self.current_row += 1 
        
        self.save_button = Button (window, text="Save Data", command=self.save_to_csv)
        self.save_button.grid(row=self.current_row,column=0,columnspan=2)
        self.current_row += 1
        
            
    '''Plot Setting Functions'''
    #Change timestamp units
    def updateTimeUnit(self,event):
        #Verify that the data frame is not empty
        if self.df.empty:
             messagebox.showinfo('Selection Error', 'Error: Please collect data before selecting time units.')
        else:
            self.unit = self.timevar.get()
            self.unit_conv = self.choices[self.unit]
            
            #Change units on timestamp plots if they exist
            self.df['timestamp']=self.df['orig_timestamp'] * self.unit_conv
            
            #Figure two is the EvT plot:
            if plt.fignum_exists(2):
                self.EvT_plot()
            
            #Figure one is the Time Histogram
            if plt.fignum_exists(1):
                self.THist_plot()
     
    #Set Energy histogram options                
    def EHist_set (self):
        self.EHist_color_set = self.EHist_color_text.get()
        self.EHist_bin_set = self.EHist_bin_text.get() 

    #Set Time histogram options
    def THist_set (self):
        self.THist_color_set = self.THist_color_text.get()
        self.THist_bin_set = self.THist_bin_text.get()

    #Set EvT plotting options
    def EvT_set (self):
        self.EvT_color_set = self.EvT_color_text.get()
        self.EvT_size_set = self.EvT_size_text.get()

    '''Plotting Functions'''   
    #Energy Histogram
    def EHist_plot (self):
        #Verify that the data frame is not empty
        if self.df.empty:
             messagebox.showinfo('Plot Error', 'Error: No data has been collected.  Can not produce empty plot.')
        else:
            self.EHist_color_set = self.EHist_color_text.get()
            self.EHist_bin_set = self.EHist_bin_text.get() 
            
            plt.figure(0)
            plt.clf()
            plt.hist(self.df['amplitude'],self.EHist_bin_set, normed=0, facecolor=self.EHist_color_set)
            plt.title("Energy Histogram", fontsize=16)
            plt.ylabel("Count",fontsize=14)
            plt.xlabel("Energy",fontsize=14)
            plt.show(block=False)
        
    #Time Histogram    
    def THist_plot (self):
        #Verify that the data frame is not empty
        if self.df.empty:
             messagebox.showinfo('Plot Error', 'Error: No data has been collected.  Can not produce empty plot.')
        else:
            self.THist_color_set = self.THist_color_text.get()
            self.THist_bin_set = self.THist_bin_text.get()
        
            #Energy V Time plot
            plt.figure(1)
            plt.clf()
            plt.hist(self.df['timestamp'],self.THist_bin_set, normed=0, facecolor=self.THist_color_set)
            plt.title("Time Histogram", fontsize=16)
            plt.ylabel("Count",fontsize=14)
            plt.xlabel(self.unit,fontsize=14)
            plt.show(block=False)
        
    #Energy V Time plot    
    def EvT_plot (self):
        #Verify that the data frame is not empty
        if self.df.empty:
             messagebox.showinfo('Plot Error', 'Error: No data has been collected.  Can not produce empty plot.')
        else:
            self.EvT_color_set = self.EvT_color_text.get()
            self.EvT_size_set = self.EvT_size_text.get()
            
            #Energy V Time plot
            plt.figure(2)
            plt.clf()
            plt.scatter(self.df['timestamp'],self.df['amplitude'],color=self.EvT_color_set,s=self.EvT_size_set)
            plt.title("Energy v Time", fontsize=16)
            plt.ylabel("Amplitude",fontsize=14)
            plt.xlabel(self.unit,fontsize=14)
            plt.show(block=False)
    
    
    '''Event Selection Functions'''
    def EvT_select (self):
        #Verify that the data frame is not empty
        if self.df.empty:
             messagebox.showinfo('Plot Error', 'Error: No data has been collected.  Can not select absent data.')
        elif not plt.fignum_exists(2):
             messagebox.showinfo('Plot Error', 'Please open the EvT plot before trying to select data.')            
        else:        
            plt.figure(2)
            plt.scatter(self.df['timestamp'],self.df['amplitude'],color=self.EvT_color_set,s=self.EvT_size_set)
            self.EvT_interface = PlotInterface(self.df,self.EvT_color_set,self.EvT_size_set)

    #Cancel event selection
    def EvT_select_cancel (self):
        if self.EvT_interface:
            del self.EvT_interface
            self.EvT_plot()

    def EvT_delete (self):
        #Verify that the data frame is not empty
        if self.df.empty:
             messagebox.showinfo('Plot Error', 'Error: No data has been collected.  Can not delete absent data.')
        else:          
            #If we have selected some data, delete that data
            if self.EvT_interface:
                sub = ['timestamp', 'amplitude']
                mask = self.df[sub].isin(self.EvT_interface.sub_df[sub].to_dict(orient='list')).all(axis=1)
                self.df = self.df[~mask]
                self.EvT_interface = None
                
            #Update the plots      
            #Figure two is the EvT plot:
            if plt.fignum_exists(2):
                self.EvT_plot()
            
            #Figure one is the Time Histogram
            if plt.fignum_exists(1):
                self.THist_plot()
                
            #Figure zero is the Energy Histogram
            if plt.fignum_exists(0):
                self.EHist_plot()     
            
    def EvT_select_only (self):
        #Verify that the data frame is not empty
        if self.df.empty:
             messagebox.showinfo('Plot Error', 'Error: No data has been collected.  Can not delete absent data.')
        else:   
            #If we have selected some data, delete that data
            if self.EvT_interface:
                sub = ['timestamp', 'amplitude']
                mask = self.df[sub].isin(self.EvT_interface.sub_df[sub].to_dict(orient='list')).all(axis=1)
                self.df = self.df[mask]
                self.EvT_interface = None
                
            #Update the plots      
            #Figure two is the EvT plot:
            if plt.fignum_exists(2):
                self.EvT_plot()
            
            #Figure one is the Time Histogram
            if plt.fignum_exists(1):
                self.THist_plot()
                
            #Figure zero is the Energy Histogram
            if plt.fignum_exists(0):
                self.EHist_plot()         
     
    
    '''Saving to CSV'''
    
    #Could add more options for saving if desired... such as save to txt, or changing delimiter
    def save_to_csv(self):
        #Verify that the data frame is not empty
        if self.df.empty:
             messagebox.showinfo('Save Error', 'Error: Can not save empty data frame')

        #If we have data
        else:
            #Get the save location and extract the directory path
            self.save_loc = self.save_text.get() 
            self.save_dir = os.sep.join(self.save_loc.split(os.sep)[:-1])

            #Verify that the directory exists
            if os.path.isdir(self.save_dir):

                #Verify that the filename is more than just '.csv'
                if len(self.save_loc.split(os.sep)[-1]) < 5:
                    messagebox.showinfo('Save Error', 'Error: Path should end in filename.csv')

                #Verify that the filename ends in '.csv'
                elif self.save_loc.split(os.sep)[-1][-4:] != '.csv':
                    messagebox.showinfo('Save Error', 'Error: File should end in filename.csv')

                #Save the data if it passes error checking
                else:
                    self.df.to_csv(self.save_loc)
                    
            #Error if the directory doesn't exist        
            else:
                messagebox.showinfo('Save Error', 'Error: File directory does not exist')
    
    '''Query Functions'''
    def sendQuery(self):
        start_t = self.start_t_text.get()
        end_t = self.end_t_text.get()
        run_tag = self.runtag_text.get()
        limit = self.limit_text.get()
        
        #Check that date is valid
        time_format = re.compile('^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}$')
        start_t_okay=None
        if (start_t) and (time_format.match(start_t)):
            start_year=int(start_t[0:4])
            start_month=int(start_t[5:7])
            start_day=int(start_t[8:10])
            start_hour=int(start_t[11:13])
            start_min=int(start_t[14:16])
            start_sec=int(start_t[17:19])
            try:
                datetime.datetime(start_year,start_month,start_day,start_hour,start_min,start_sec)
                start_t_okay='yes'
            except:
                start_t_okay='no'


        end_t_okay=None
        if (end_t) and (time_format.match(end_t)):
            end_year=int(end_t[0:4])
            end_month=int(end_t[5:7])
            end_day=int(end_t[8:10])
            end_hour=int(end_t[11:13])
            end_min=int(end_t[14:16])
            end_sec=int(end_t[17:19])
            try:
                datetime.datetime(end_year,end_month,end_day,end_hour,end_min,end_sec)
                end_t_okay='yes'
            except:
                end_t_okay='no'        

        '''Error Handling'''
        #Check that the date is in the correct format
        if (start_t) and (not time_format.match(start_t)):
            messagebox.showinfo('Query Error', 'Error: Start time Must have the format "YYYY-MM-DD HH:MM:SS" or be empty')
        elif (end_t) and (not time_format.match(end_t)):
            messagebox.showinfo('Query Error', 'Error: End time Must have the format "YYYY-MM-DD HH:MM:SS" or be empty')
        #Check that the date is a valid date
        elif (start_t) and (start_t_okay == 'no'):        
            messagebox.showinfo('Query Error', 'Error: Start time is not valid')
        elif (end_t) and (end_t_okay == 'no'):
            messagebox.showinfo('Query Error', 'Error: End time is not valid')
        #Check that start_t < end_t
        elif (start_t) and (end_t) and (start_t > end_t):         
            messagebox.showinfo('Query Error', 'Start time must come before end time')    
        else:
            #Default to not including the conditions in the query
            run_tag_string = ''
            time_string = ''
            limit_string = ''
            
            #Toggle query conditions on if the GUI field isn't blank
            
            #Run tag condition
            if run_tag != '':
                run_tag_string = 'AND run_tag = "%s"' %(run_tag)
                
            #Time window condition (in order: both set, only start set, only end set)
            if start_t != '' and end_t != '':
                time_string = 'AND clock_time BETWEEN "%s" AND "%s"' %(start_t,end_t)
            elif start_t != '':
                time_string = 'AND clock_time > "%s"' %(start_t)
            elif end_t != '':
                time_string = 'AND clock_time < "%s"' %(end_t)           
                
            #Limit query
            if limit != '':
                limit_string = 'LIMIT %s' %(limit)
                            
            #Setting the query
            #We don't really need "amplitude IS NOT NULL" condition, it was added so that
            #we can add conditions with simple 'and ...'
            #This way we don't have to worry about whether or not to include a 'where' clause
            self.query = "SELECT clock_time, amplitude FROM summarydata where amplitude IS NOT NULL %s %s %s" % (run_tag_string, time_string, limit_string)
            
            con = mdb.connect(self.ip, self.user, self.pswd, self.db);

            #With will close the connection after the code is done, 
            #regardless of how the code exists. Use as an alternative to 'finally' statement
            with con:

                cur = con.cursor()

                #Send select statement
                #columns are [id, 'run_tag', clock_time, centroid_time, amplitude, rmse]
                cur.execute(self.query)

                #Get all of the MySQL return at once
                #Returns a tuple of tuples, with each inner tupple being one row
                self.df = pd.DataFrame()
                rows = cur.fetchall()
                self.df['timestamp']=[time.mktime(rows[i][0].timetuple()) for i in range(len(rows))]
          
                #Set up default 1 second time stamps
                self.start_timestamp =  min(self.df['timestamp'])
                self.df['timestamp'] = self.df['timestamp'] - self.start_timestamp
                self.df['orig_timestamp']=self.df['timestamp']
                self.df['amplitude'] = [rows[i][1] for i in range(len(rows))]

            self.querystatus.configure(text = "Finished")

            #Check if the query returned data
            if self.df.empty:
                 messagebox.showinfo('Query Warning', 'Warning: Query returned no data')
        

#Use a main wrapper so that the code isn't executed if it is ever imported as a module
#Arguements need to be IP of database, 
def main(argv):
	'''
	Usage: python radonGUI.py databaseIP username password databaseTable
	'''
	ip = argv[1]
	user = argv[2]
	pswd = argv[3]
	db = argv[4]

	window= Tk()
	start= App(window,ip,user,pswd,db)
	window.mainloop()

if __name__ == "__main__":
    main(sys.argv)
