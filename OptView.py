"""

Provides interactive visualization of optimization results created by
pyOptSparse. Figures produced here can be saved as images or pickled
for future customization.

John Jasa 2015-2017

"""

# ======================================================================
# Standard Python modules
# ======================================================================
import os
import argparse
import shelve
import tkFont
import Tkinter as Tk
import re
import warnings

# ======================================================================
# External Python modules
# ====================================================================== 
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,\
    NavigationToolbar2TkAgg
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
warnings.filterwarnings("ignore",category=UserWarning)
import numpy as np
from pyoptsparse import SqliteDict

class Display(object):

    """
    Container for display parameters, properties, and objects.
    This includes a canvas for MPL plots and a bottom area with widgets.
    """

    def __init__(self, histList, outputDir):

        # Initialize the Tkinter object, which will contain all graphical
        # elements.
        self.root = Tk.Tk()
        self.root.wm_title("OptView")

        # Load the OptView icon
        try:
            icon_dir = os.path.dirname(os.path.abspath(__file__))
            icon_name = 'OptViewIcon.gif'
            icon_dir_full = os.path.join(icon_dir, icon_name)
            img = Tk.PhotoImage(file=icon_dir_full)
            self.root.tk.call('wm', 'iconphoto', self.root._w, img)
        except: # bare except because error is not in standard Python
            pass

        # If the screen is bigger than 1080p, use a large window
        if self.root.winfo_screenheight() > 1100:
            figsize = (14, 10)
        else: # Otherwise, use a slightly smaller window
              # so everything fits on the screen
            figsize = (6, 6)

        # Instantiate the MPL figure
        self.f = plt.figure(figsize=figsize, dpi=100, facecolor='white')

        # Link the MPL figure onto the TK canvas and pack it
        self.canvas = FigureCanvasTkAgg(self.f, master=self.root)
        self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        # Add a toolbar to explore the figure like normal MPL behavior
        toolbar = NavigationToolbar2TkAgg(self.canvas, self.root)
        toolbar.update()
        self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        # Increase the font size
        matplotlib.rcParams.update({'font.size': 16})

        # Initialize lists, dicts, and save inputs from user
        self.arr_active = 0
        self.plots = []
        self.annotate = None
        self.histList = histList
        self.outputDir = outputDir
        self.bounds = {}

        # Actually setup and run the GUI
        self.OptimizationHistory()

    def OptimizationHistory(self):
        """
        Reads in database history file and stores contents.
        Function information is stored as a dict in func_data,
        variable information is stored as a dict in var_data,
        and bounds information is stored as a dict in bounds.
        """

        # Initialize dictionaries for design variables and unknowns.
        # The data is saved redundantly in dicts for all iterations and then
        # for major iterations as well.
        self.func_data_all = {}
        self.func_data_major = {}
        self.var_data_all = {}
        self.var_data_major = {}
        db = {}
        self.num_iter = 0

        # Loop over each history file name provided by the user.
        for histIndex, histFileName in enumerate(self.histList):

            # If they only have one history file, we don't change the keys' names
            if len(self.histList) == 1:
                histIndex = ''
            else: # If multiple history files, append letters to the keys,
                  # such that 'key' becomes 'key_A', 'key_B', etc
                histIndex = '_' + chr(histIndex + ord('A'))
            self.histIndex = histIndex

            try: # This is the classic method of storing history files
                db = shelve.open(histFileName, 'r')
                OpenMDAO = False
            except: # Bare except because error is not in standard Python.
                # If the db has the 'iterations' tag, it's an OpenMDAO db.
                db = SqliteDict(histFileName, 'iterations')
                OpenMDAO = True

                # If it has no 'iterations' tag, it's a pyOptSparse db.
                if db.keys() == []:
                    OpenMDAO = False
                    db = SqliteDict(histFileName)

            # Specific instructions for OpenMDAO databases
            if OpenMDAO:

                # Get the number of iterations by looking at the largest number
                # in the split string names for each entry in the db
                string = db.keys()[-1].split('|')
                nkey = int(string[-1])
                self.solver_name = string[0]

                # Initalize a list detailing if the iterations are major or minor
                self.iter_type = np.zeros(nkey)

                # Get the keys of the database where derivatives were evaluated.
                # These correspond to major iterations, while no derivative
                # info is calculated for gradient-free linesearches.
                deriv_keys = SqliteDict(histFileName, 'derivs').keys()
                self.deriv_keys = [int(key.split('|')[-1]) for key in deriv_keys]

                # Save information from the history file for the unknowns.
                self.SaveDBData(db, self.func_data_all, self.func_data_major, OpenMDAO=OpenMDAO, data_str='Unknowns')

                # Save information from the history file for the design variables.
                self.SaveDBData(db, self.var_data_all, self.var_data_major, OpenMDAO=OpenMDAO, data_str='Parameters')

                # Add labels to OpenMDAO variables.
                # Corresponds to constraints, design variables, and objective.
                try:
                    db = SqliteDict(histFileName, 'metadata')
                    self.SaveOpenMDAOData(db)

                except KeyError: # Skip metadata info if not included in OpenMDAO hist file
                    pass

            else:

                # Get the number of iterations
                nkey = int(db['last']) + 1
                self.nkey = nkey

                # Initalize a list detailing if the iterations are major or minor
                self.iter_type = np.zeros(nkey)

                # Check to see if there is bounds information in the db file.
                # If so, add them to self.bounds to plot later.
                try:
                    bounds_dict = dict(db['varBounds'].items() + db['conBounds'].items())
                    for key in bounds_dict.keys():
                        bounds_dict[key + histIndex] = bounds_dict.pop(key)
                    self.bounds.update(bounds_dict)
                except KeyError:
                    pass

                # Check to see if there is proper saved info about iter type
                if 'isMajor' in db['0'].keys():
                    self.storedIters = True
                else:
                    self.storedIters = False

                # Save information from the history file for the funcs.
                self.SaveDBData(db, self.func_data_all, self.func_data_major, OpenMDAO=OpenMDAO, data_str='funcs')

                # Save information from the history file for the design variables.
                self.SaveDBData(db, self.var_data_all, self.var_data_major, OpenMDAO=OpenMDAO, data_str='xuser')

        # Set the initial dictionaries to reference all iterations.
        # Later this can be set to reference only the major iterations.
        self.func_data = self.func_data_all
        self.var_data = self.var_data_all

        # Find the maximum length of any variable in the dictionaries and
        # save this as the number of iterations.
        for data_dict in [self.func_data, self.var_data]:
            for key in data_dict.keys():
                length = len(data_dict[key])
                if length > self.num_iter:
                    self.num_iter = length

    def SaveDBData(self, db, data_all, data_major, OpenMDAO, data_str):
        """ Method to save the information within the database corresponding
            to a certain key to the relevant dictionaries within the Display
            object. This method is called twice, once for the design variables
            and the other for the outputs. """

        # Loop over each optimization iteration
        for i, iter_type in enumerate(self.iter_type):

            # If this is an OpenMDAO file, the keys are of the format
            # 'rank0:SNOPT|1', etc
            if OpenMDAO:
                key = '{}|{}'.format(self.solver_name, i+1) # OpenMDAO uses 1-indexing
            else: # Otherwise the keys are simply a number
                key = '%d' % i

            # If this is the 'funcs' key, perform some special operations,
            # such as determining the iteration type
            if data_str == 'funcs':

                # Only actual optimization iterations have 'funcs' in them.
                # pyOptSparse saves info for two iterations for every
                # actual major iteration. In particular, one has funcs
                # and the next has funcsSens, but they're both part of the
                # same major iteration.
                if any('funcs' == s for s in db[key].keys()):

                    # If the proper history is stored coming out of
                    # pyoptsparse, use that for filtering major iterations.
                    if self.storedIters:
                        self.iter_type[i] = int(db[key]['isMajor'])

                        # If this iteration has 'funcs' within it, but it's not
                        # flagged as major, then it's a minor iteration.
                        if self.iter_type[i] == 0:
                            self.iter_type[i] = 2

                    else: # Otherwise, use a spotty heuristic to see if the
                        # iteration is major or not. NOTE: this is often
                        # inaccurate, especially if the optimization used
                        # gradient-enhanced line searches.
                        try:
                            keyp1 = '%d' % (i + 1)
                            db[keyp1]['funcsSens']
                            self.iter_type[i] = 1 # for 'major' iterations
                        except KeyError:
                            self.iter_type[i] = 2 # for 'minor' iterations

                else:
                    self.iter_type[i] = 0 # this is not a real iteration,
                                          # just the sensitivity evaluation
            elif data_str == 'Unknowns':
                if i in self.deriv_keys:
                    self.iter_type[i] = 1 # for 'major' iterations
                else:
                    self.iter_type[i] = 2 # for 'minor' iterations

            # Do this for both major and minor iterations
            if self.iter_type[i]:

                # Get just the info in the dict for this iteration
                iter_data = db[key][data_str]

                # Loop through each key within this iteration
                for key in sorted(iter_data):

                    # Format a new_key string where we append a modifier
                    # if we have multiple history files
                    new_key = key + '{}'.format(self.histIndex)

                    # If this key is not in the data dictionaries, add it
                    if new_key not in data_all:
                        data_all[new_key] = []
                        data_major[new_key] = []

                    # Process the data from the key. Convert it to a np
                    # array, keep only the real part, squeeze any 1-dim
                    # axes out of it, then flatten it.
                    data = np.squeeze(np.array(iter_data[key]).real).flatten()

                    # Append the data to the entries within the dictionaries.
                    data_all[new_key].append(data)
                    if self.iter_type[i] == 1:
                        data_major[new_key].append(data)

    def SaveOpenMDAOData(self, db):
        """ Examine the OpenMDAO dict and save tags if the variables are
            objectives (o), constraints (c), or design variables (dv). """

        # Loop over each key in the metadata db
        for tag in db:

            # Only look at variables and unknowns
            if tag in ['Unknowns', 'Parameters']:
                for old_item in db[tag]:

                    # We'll rename each item, so we need to get the old item
                    # name and modify it
                    item = old_item + '{}'.format(self.histIndex)

                    # Here we just have an open parenthesis, and then we will
                    # add o, c, or dv. Note that we could add multiple flags
                    # to a single item. That's why we have a sort of convoluted
                    # process of adding the tags.
                    new_key = item + ' ('
                    flag_list = []

                    # Check each flag and see if they have the relevant entries
                    # within the dict; if so, tag them.
                    for flag in db[tag][old_item]:
                        if 'is_objective' in flag:
                            flag_list.append('o')
                        if 'is_desvar' in flag:
                            flag_list.append('dv')
                        if 'is_constraint' in flag:
                            flag_list.append('c')

                    # Create the new_key based on the flags for each variable
                    for flag in flag_list:
                        if flag == flag_list[-1]:
                            new_key += flag + ')'
                        else:
                            new_key += flag + ', '

                    # If there are actually flags to add, pop out the old items
                    # in the dict and re-add them with the new name.
                    if flag_list:
                        try:
                            if 'dv' in flag_list:
                                self.var_data_all[new_key] = self.func_data_all.pop(item)
                                self.var_data_major[new_key] = self.func_data_major.pop(item)

                            else:
                                self.func_data_all[new_key] = self.func_data_all.pop(item)
                                self.func_data_major[new_key] = self.func_data_major.pop(item)

                        except KeyError:
                            pass

    def quit(self):
        """
        Destroy GUI window cleanly if quit button pressed.
        """
        self.root.quit()
        self.root.destroy()

    def error_display(self, string="That option not supported"):
        """
        Display error string on canvas when invalid options selected.
        """
        self.f.clf()
        a = self.f.add_subplot(111)
        a.text(0.05, .9,
            "Error: " + string,
            fontsize=20,
            transform=a.transAxes)
        self.canvas.show()

    def warning_display(self, string="That option not supported"):
        """
        Display warning message on canvas as necessary.
        """
        a = plt.gca()
        a.text(0.05, 1.04,
            "Warning: " + string,
            fontsize=20,
            transform=a.transAxes)
        self.canvas.show()

    def plot_bounds(self, val, a, color):
        """
        Plot the bounds if selected.
        """

        if val not in self.bounds:
            for ii, char in enumerate(reversed(val)):
                if char == '_':
                    split_loc = len(val) - ii
                    break
            val_name = val[:split_loc - 1]
            val_num = int(val[split_loc:])
            lower = [self.bounds[val_name]['lower'][val_num]]
            upper = [self.bounds[val_name]['upper'][val_num]]
        else:
            lower = self.bounds[val]['lower']
            upper = self.bounds[val]['upper']

        lower = list(lower)
        upper = list(upper)

        a.margins(None, .02)
        a.set_color_cycle(color)
        for lower_bound in lower:
            if lower_bound is not None and abs(lower_bound) < 1e18:
                a.plot(
                    [0, self.num_iter - 1], [
                        lower_bound, lower_bound],
                    "--", linewidth=2, clip_on=False
                )

        a.set_color_cycle(color)
        for upper_bound in upper:
            if upper_bound is not None and abs(upper_bound) < 1e18:
                a.plot(
                    [0, self.num_iter - 1], [
                        upper_bound, upper_bound],
                    "--", label=val + ' bounds', linewidth=2, clip_on=False)

    def orig_plot(self, dat, val, values, a, i=0):
        """
        Plots the original data values from the history file.
        """
        cc = (
            matplotlib.rcParams['axes.color_cycle'] * 10
        )
        color = cc[i]

        try:
            array_size = len(dat[val][0])
            if self.var_minmax.get():
                a.set_color_cycle(color)
                minmax_list = []
                for minmax in dat[val]:
                    minmax_list.append(
                        [np.min(minmax), np.max(minmax)])
                plots = a.plot(minmax_list, "o-", label=val,
                    markeredgecolor='none', clip_on=False)

            elif array_size < 20 or self.var_showall.get():
                if i > 0:
                    a.set_color_cycle(color)
                plots = a.plot(dat[val], "o-", label=val,
                    markeredgecolor='none', clip_on=False)

                a.set_ylabel(val)
                self.color_error_flag = 1
            else:
                self.error_display("Too many values to display")

        except TypeError:
            a.set_color_cycle(color)
            if self.var.get() == 0:
                pass
            else:
                a.set_ylabel(val)
            plots = a.plot(dat[val], "o-", label=val,
                markeredgecolor='none', clip_on=False)

        except (KeyError, IndexError):
            self.warning_display("No 'major' iterations")
        try:
            if len(plots) > 1:
                for i, plot in enumerate(plots):
                    self.plots.append([plot, i])
            else:
                self.plots.append([plots[0], -1])
        except UnboundLocalError:
            self.error_display("Too many values to display")
        try:
            if self.var_bounds.get():
                self.plot_bounds(val, a, color)

        except (UnboundLocalError, ValueError):
            if len(values) > 1:
                pass
            else:
                self.error_display("No bounds information")

    def plot_selected(self, values, dat):
        """
        Plot data based on selected keys within listboxes.
        """
        fail = 0
        self.color_error_flag = 0
        self.f.clf()
        self.plots = []
        try:
            if self.var_bounds.get():
                try:
                    self.bounds
                except AttributeError:
                    self.error_display("No bounds information in history file")
                    fail = 1

            # Plot on shared axes
            if self.var.get() == 0 and not fail:
                a = self.f.add_subplot(111)

                # Calculate and plot the delta values if selected
                if self.var_del.get():
                    for idx, val in enumerate(values):
                        newdat = []
                        for i, value in enumerate(dat[val], start=1):
                            newdat.append(abs(value - dat[val][i - 2]))
                        plots = a.plot(
                            range(1, self.num_iter),
                            newdat[1:],
                            "o-",
                            label=val,
                            markeredgecolor='none', clip_on=False)
                        if len(plots) > 1:
                            for i, plot in enumerate(plots):
                                self.plots.append([plot, i])
                        else:
                            self.plots.append([plots[0], -1])

                # Otherwise plot original data
                else:
                    for i, val in enumerate(values):
                        self.orig_plot(dat, val, values, a, i)

                if self.color_error_flag and self.var_bounds.get():
                    self.warning_display(
                        "Line color for bounds may not match data color")

                # Plot using log scale if selected
                if self.var_log.get():
                    a.set_yscale('log')

                if self.var_legend.get():
                    a.legend(loc='best')

                plt.subplots_adjust(right=.95)
                a.set_xlabel('iteration')
                a.set_xlim(0, self.num_iter - 1)
                self.canvas.show()

            # Plot on individual vertical axes
            elif self.var.get() == 1 and not fail:

                # Set window sizing parameters for when additional axes are
                # added
                n = len(values)
                plt.figure(self.f.number)
                par_list = [[] for i in range(n)]  # make into array
                par_list[0] = host_subplot(111, axes_class=AA.Axes)
                size_list = [.95, .95, .93, .83, .73, .63, .53, .43, .33]
                plt.subplots_adjust(right=size_list[n])

                for i in range(1, n):
                    par_list[i] = par_list[0].twinx()

                offset = 60
                for i in range(2, n):
                    new_fixed_axis = par_list[
                        i].get_grid_helper().new_fixed_axis
                    par_list[i].axis["right"] = new_fixed_axis(
                        loc="right", axes=par_list[i], offset=(offset * i ** 1.15, 0))
                    par_list[i].axis["right"].toggle(all=True)

                p_list = [[] for i in range(n)]

                # Compute and plot delta values if selected
                if self.var_del.get():
                    for i, val in enumerate(values):
                        newdat = []
                        for idx, value in enumerate(dat[val], start=1):
                            newdat.append(abs(value - dat[val][idx - 2]))
                        p_list[i], = par_list[i].plot(range(1, self.num_iter),
                                                    newdat[1:], "o-", label=val,
                                                    markeredgecolor='none', clip_on=False)
                        par_list[i].set_ylabel(val)
                # Otherwise plot original data
                else:
                    for i, val in enumerate(values):
                        cc = (matplotlib.rcParams['axes.color_cycle'] * 10)
                        par_list[i].set_color_cycle(cc[i])
                        p_list[i], = par_list[i].plot(
                            dat[val], "o-", label=val, markeredgecolor='none', clip_on=False)
                        par_list[i].set_ylabel(val)

                        try:
                            if self.var_bounds.get():
                                self.plot_bounds(val, par_list[i], cc[i])
                        except (ValueError, UnboundLocalError):
                            if len(values) > 1:
                                pass
                            else:
                                self.error_display("No bounds information")

                # Plot using log scale if selected
                if self.var_log.get():
                    for ax in par_list:
                        ax.set_yscale('log')

                par_list[0].set_xlim(0, self.num_iter - 1)
                par_list[0].set_xlabel('iteration')
                if self.var_legend.get():
                    par_list[0].legend(loc='best')
                for i, plot in enumerate(p_list):
                    self.plots.append([plot, i])

                self.canvas.show()

            # Plot on stacked axes with shared x-axis
            elif self.var.get() == 2 and not fail:
                n = len(values)

                # Compute and plot delta values if selected
                if self.var_del.get():
                    a = []
                    for i, val in enumerate(values):
                        newdat = []
                        for idx, value in enumerate(dat[val], start=1):
                            newdat.append(abs(value - dat[val][idx - 2]))
                        a.append(self.f.add_subplot(n, 1, i + 1))
                        plots = a[i].plot(range(1, self.num_iter), newdat[1:],
                                          "o-", label=val, markeredgecolor='none', clip_on=False)
                        a[i].set_ylabel('delta ' + val)
                        self.plots.append([plots[0], -1])

                # Otherwise plot original data
                else:
                    a = []
                    for i, val in enumerate(values):
                        a.append(self.f.add_subplot(n, 1, i + 1))
                        self.orig_plot(dat, val, values, a[i])

                # Plot using log scale if selected
                if self.var_log.get():
                    for ax in a:
                        ax.set_yscale('log')

                # Turn off horiztonal axes above the bottom plot
                a[-1].set_xlabel('iteration')
                for ax in a:
                    if ax != a[-1]:
                        ax.spines['bottom'].set_visible(False)
                        ax.set_xticklabels([])
                        ax.xaxis.set_major_locator(plt.NullLocator())
                    ax.spines['top'].set_visible(False)
                    for tic in ax.xaxis.get_major_ticks():
                        tic.tick2On = False
                    ax.tick_params(
                        axis='y',
                        which='both',
                        labelleft='off',
                        labelright='on')
                    ax.set_xlim(0, self.num_iter - 1)

                plt.subplots_adjust(right=.95)
                self.canvas.show()
        except ValueError:
            self.error_display()

    def onselect(self, evt, data_name):
        """
        Update current plot with selected data from listboxes.
        Also checks if the data is array-type and provides an
        additional listbox to select data within that array.
        """
        w = evt.widget
        values = [w.get(int(i)) for i in w.curselection()]
        self.update_graph()
        if len(values) == 1:
            try:
                data = data_name[values[0]]
                data[0][0]
                self.v.set(values[0])
                self.lb_arr.delete(0, Tk.END)
                for i, val in enumerate(data[0]):
                    self.lb_arr.insert(Tk.END, values[0] + '_' + str(i))
                self.arr_title.pack(side=Tk.TOP)
                self.scrollbar_arr.pack(side=Tk.RIGHT, fill=Tk.Y)
                self.lb_arr.pack(side=Tk.RIGHT)
                self.arr_active = 1
            except (IndexError, TypeError):
                self.lb_arr.pack_forget()
                self.scrollbar_arr.pack_forget()
                self.arr_title.pack_forget()
                self.arr_active = 0
            except KeyError:
                self.warning_display("No 'major' iterations")

    def onselect_arr(self, evt):
        """
        Obtain the selected plotting values from the array-based variable listbox.
        """
        w = evt.widget
        values = [int(i) for i in w.curselection()]

        # Get the currently selected functions/variables
        func_sel = self.lb_func.curselection()
        var_sel = self.lb_var.curselection()
        if len(func_sel):
            values_orig = [self.lb_func.get(i) for i in func_sel]
            dat = self.func_data[values_orig[0]]
        elif len(var_sel):
            values_orig = [self.lb_var.get(i) for i in var_sel]
            dat = self.var_data[values_orig[0]]

        # Add the array-based information to the listbox for selection
        self.arr_data = {}
        self.val_names = []
        for i, val in enumerate(values):
            self.val_names.append(values_orig[0] + '_{0}'.format(val))
            self.arr_data[self.val_names[i]] = []
            for ind_dat in dat:
                self.arr_data[self.val_names[i]].append(ind_dat[val])
        self.plot_selected(self.val_names, self.arr_data)

    def update_graph(self):
        """
        Produce an updated graph based on user options.
        """

        if self.var_minmax.get() and self.var_showall.get():
            self.error_display("Cannot show all and min/max at same time")

        else:
            func_sel = self.lb_func.curselection()
            var_sel = self.lb_var.curselection()
            arr_sel = self.lb_arr.curselection()
            values = []
            dat = {}

            if len(arr_sel) and self.arr_active:
                self.plot_selected(self.val_names, self.arr_data)

            elif len(func_sel) or len(var_sel):
                values.extend([self.lb_func.get(i) for i in func_sel])
                dat = self.func_data.copy()
                values.extend([self.lb_var.get(i) for i in var_sel])
                dat.update(self.var_data)
                self.plot_selected(values, dat)

    def set_mask(self):

        if self.var_mask.get():
            self.func_data = self.func_data_major
            self.var_data = self.var_data_major
        else:
            self.func_data = self.func_data_all
            self.var_data = self.var_data_all

        self.num_iter = 0

        for key in self.func_data.keys():
            length = len(self.func_data[key])
            if length > self.num_iter:
                self.num_iter = length

        self.update_graph()

    def save_figure(self):
        """
        Save the current figure using the selected variables as the filename.
        """
        func_sel = self.lb_func.curselection()
        var_sel = self.lb_var.curselection()
        arr_sel = self.lb_arr.curselection()
        values = []
        if len(arr_sel) and self.arr_active:
            values = self.val_names
        elif len(func_sel):
            values = [self.lb_func.get(i) for i in func_sel]
        elif len(var_sel):
            values = [self.lb_var.get(i) for i in var_sel]
        groups = ''
        for string in values:
            groups += string + '_'
        fname = groups + '.png'
        fpathname = os.path.join(self.outputDir, fname)
        plt.savefig(fpathname)
        fname = 'saved_figure.pickle'
        fpathname = os.path.join(self.outputDir, fname)
        try:
            import dill
            dill.dump(self.f, file(fpathname, 'wb'))
        except ImportError:
            pass

    def save_all_figues(self):
        """
        Batch save all individual figures from functions and variables.
        """
        for data_name in [self.func_data, self.var_data]:
            for key in data_name:
                fig = plt.figure()
                plt.plot(data_name[key], 'ko-')
                plt.title(key)
                fname = key + '.png'
                fpathname = os.path.join(self.outputDir, fname)
                plt.savefig(fpathname)
                plt.clf()

    def save_tec(self):
        """
        Output selected data to tec file.
        """
        func_sel = self.lb_func.curselection()
        var_sel = self.lb_var.curselection()
        arr_sel = self.lb_arr.curselection()
        dat = {}
        if len(arr_sel) and self.arr_active:
            for name in self.val_names:
                dat[name] = self.arr_data[name]
        elif len(func_sel):
            values = [self.lb_func.get(i) for i in func_sel]
            for name in values:
                dat[name] = self.func_data[name]
        elif len(var_sel):
            values = [self.lb_var.get(i) for i in var_sel]
            for name in values:
                dat[name] = self.var_data[name]

        keys = dat.keys()

        num_vars = len(keys)
        num_iters = len(dat[keys[0]])
        full_data = np.arange(num_iters, dtype=np.float_).reshape(num_iters, 1)
        var_names = ['Iteration']
        for key in keys:
            small_data = np.asarray(dat[key])

            if len(small_data.shape) == 1:
                full_data = np.c_[full_data, small_data]
                var_names.append(key)

            else:
                m = small_data.shape[0]
                n = small_data.shape[1]
                indiv_data = np.empty((m, 1))
                for i in range(n):
                    for j in range(m):
                        indiv_data[j] = small_data[j][i]
                    full_data = np.c_[full_data, indiv_data]
                    var_names.append(key + '_{}'.format(i))

        filename = 'OptView_tec.dat'
        self._file = open(filename, 'w')
        self._file.write('Title = \"OptView data output\"" \n')
        self._file.write('Variables = ')
        for name in var_names:
            self._file.write('\"' + name + '\" ')
        self._file.write('\n')

        self._file.write('Zone T= \"OptView_tec_data\", ' + \
                    'I={}, '.format(num_iters) + 'F=POINT\n')
        np.savetxt(self._file, full_data)
        self._file.close()

    def var_search(self, _):
        """
        Remove listbox entries that do not contain user-inputted string,
        used to search through outputted data.
        """
        self.lb_func.delete(0, Tk.END)
        self.lb_var.delete(0, Tk.END)
        for key in sorted(self.func_data):
            self.lb_func.insert(Tk.END, key)
        for key in sorted(self.var_data):
            self.lb_var.insert(Tk.END, key)

        search_entry = self.entry_search.get()
        func_range = range(len(self.func_data))
        for i in func_range[::-1]:
            if not re.search(search_entry.lower(), self.lb_func.get(i).lower()):
                self.lb_func.delete(i)

        var_range = range(len(self.var_data))
        for i in var_range[::-1]:
            if not re.search(search_entry.lower(), self.lb_var.get(i).lower()):
                self.lb_var.delete(i)

        if not self.lb_var.get(1) and not self.lb_func.get(1):
            if self.lb_var.get(0):
                self.lb_var.select_set(0)
            if self.lb_func.get(0):
                self.lb_func.select_set(0)
            self.update_graph()

    def update_font(self, val):
        """
        Set the font for matplotlib based on slider.
        """
        matplotlib.rcParams.update({'font.size': int(val)})
        self.update_graph()

    def refresh_history(self):
        """
        Refresh opt_his data if the history file has been updated.
        """
        old_funcs = []
        for key in self.func_data:
            old_funcs.append(key)
        old_vars = []
        for key in self.var_data:
            old_vars.append(key)

        self.OptimizationHistory()

        new_funcs = []
        for key in self.func_data:
            new_funcs.append(key)
        new_vars = []
        for key in self.var_data:
            new_vars.append(key)

        if not (old_funcs == new_funcs and old_vars == new_vars):
            self.var_search('dummy')

    def refresh_history_init(self):
        self.refresh_history()
        self.set_mask()

    def auto_ref(self):
        """
        Automatically refreshes the history file, which is
        useful if examining a running optimization.
        """
        if self.var_ref.get():
            self.root.after(1000, self.auto_ref)
            self.refresh_history()
            self.set_mask()

    def clear_selections(self):
        """
        Deselects all currently-selected variables, functions, and array options
        """
        self.lb_func.selection_clear(0, Tk.END)
        self.lb_var.selection_clear(0, Tk.END)
        self.lb_arr.selection_clear(0, Tk.END)

    def on_move(self, event):
        """
        Checks to see if the cursor is over a plot and provides a
        hovering label if necessary.
        """
        try:
            self.annotation.remove()
        except (AttributeError, ValueError):
            pass
        if event.xdata:
            visibility_changed = False
            point_selected = None
            for point in self.plots:
                if point[0].contains(event)[0]:
                    point_selected = point

            # Prevent error message if we move out of bounds while hovering
            # over a point on a line
            if point_selected:
                visibility_changed = True
                ax = point_selected[0].get_axes()
                label = point_selected[0].get_label()
                if point_selected[1] >= 0:
                    label = label + '_' + str(point_selected[1])

                xdat = point_selected[0].get_xdata()
                ydat = point_selected[0].get_ydata()

                iter_count = np.round(event.xdata, 0)
                ind = np.where(xdat == iter_count)[0][0]

                label = label + '\niter: {0:d}\nvalue: {1}'.format(int(iter_count), ydat[ind])
                self.annotation = ax.annotate(label,
                                              xy=(event.xdata,
                                                  event.ydata), xycoords='data',
                                              xytext=(
                                              event.xdata, event.ydata), textcoords='data',
                                              horizontalalignment="left",
                                              bbox=dict(
                                              boxstyle="round", facecolor="w",
                                              edgecolor="0.5", alpha=0.8),
                                              )
        else:
            try:
                self.annotation.remove()
            except (AttributeError, ValueError):
                pass

        self.canvas.show()

    def draw_GUI(self):
        """
        Create the frames and widgets in the bottom section of the canvas.
        """
        font = tkFont.Font(family="Helvetica", size=10)

        sel_frame = Tk.Frame(self.root)
        sel_frame.pack(side=Tk.LEFT)

        # Produce a frame and listbox to contain function information
        func_frame = Tk.Frame(sel_frame)
        func_frame.pack(side=Tk.LEFT, fill=Tk.Y, padx=20)
        func_title = Tk.Label(func_frame, text="Functions", font=font)
        func_title.pack(side=Tk.TOP)
        self.lb_func = Tk.Listbox(
            func_frame,
            name='lb_func',
            selectmode=Tk.EXTENDED,
            font=font,
            width=30,
            exportselection=0)
        self.lb_func.pack(side=Tk.LEFT)
        for key in sorted(self.func_data):
            self.lb_func.insert(Tk.END, key)
        scrollbar_func = Tk.Scrollbar(func_frame)
        scrollbar_func.pack(side=Tk.RIGHT, fill=Tk.Y)
        self.lb_func.config(yscrollcommand=scrollbar_func.set)
        scrollbar_func.config(command=self.lb_func.yview)
        self.lb_func.bind(
            '<<ListboxSelect>>',
            lambda event: self.onselect(
                event,
                self.func_data))

        # Produce a frame and listbox to contain variable information
        var_frame = Tk.Frame(sel_frame)
        var_frame.pack(side=Tk.RIGHT, fill=Tk.Y, padx=20)
        var_title = Tk.Label(var_frame, text="Design Variables", font=font)
        var_title.pack(side=Tk.TOP)
        scrollbar_var = Tk.Scrollbar(var_frame)
        scrollbar_var.pack(side=Tk.RIGHT, fill=Tk.Y)
        self.lb_var = Tk.Listbox(
            var_frame,
            name='lb_var',
            selectmode=Tk.EXTENDED,
            font=font,
            width=30,
            exportselection=0)
        self.lb_var.pack(side=Tk.RIGHT)
        for key in sorted(self.var_data):
            self.lb_var.insert(Tk.END, key)
        self.lb_var.config(yscrollcommand=scrollbar_var.set)
        scrollbar_var.config(command=self.lb_var.yview)
        self.lb_var.bind(
            '<<ListboxSelect>>',
            lambda event: self.onselect(
                event,
                self.var_data))

        # Produce a frame and listbox to contain array-based variable
        # information
        arr_frame = Tk.Frame(self.root)
        arr_frame.pack(side=Tk.RIGHT, fill=Tk.Y, padx=20)
        self.v = Tk.StringVar()
        self.arr_title = Tk.Label(
            arr_frame,
            text="Array variables",
            font=font,
            textvariable=self.v)
        self.scrollbar_arr = Tk.Scrollbar(arr_frame)
        self.lb_arr = Tk.Listbox(
            arr_frame,
            name='lb_arr',
            selectmode=Tk.EXTENDED,
            font=font,
            width=30,
            exportselection=0)
        self.lb_arr.config(yscrollcommand=self.scrollbar_arr.set)
        self.scrollbar_arr.config(command=self.lb_arr.yview)
        self.lb_arr.bind('<<ListboxSelect>>', self.onselect_arr)

        # Create options frame with buttons and checkboxes for user
        options_frame = Tk.Frame(self.root)
        options_frame.pack(side=Tk.LEFT, fill=Tk.Y, pady=10)

        button0 = Tk.Button(
            options_frame,
            text='Clear selections',
            command=self.clear_selections,
            font=font)
        button0.grid(row=0, column=2, padx=5, sticky=Tk.W)

        button1 = Tk.Button(
            options_frame,
            text='Refresh history',
            command=self.refresh_history_init,
            font=font)
        button1.grid(row=1, column=2, padx=5, sticky=Tk.W)

        button2 = Tk.Button(
            options_frame,
            text='Save all figures',
            command=self.save_all_figues,
            font=font)
        button2.grid(row=2, column=2, padx=5, sticky=Tk.W)

        button3 = Tk.Button(
            options_frame,
            text='Save figure',
            command=self.save_figure,
            font=font)
        button3.grid(row=3, column=2, padx=5, sticky=Tk.W)

        button4 = Tk.Button(
            options_frame,
            text='Save tec file',
            command=self.save_tec,
            font=font)
        button4.grid(row=4, column=2, padx=5, sticky=Tk.W)

        button5 = Tk.Button(
            options_frame,
            text='Quit',
            command=self.quit,
            font=font)
        button5.grid(row=5, column=2, padx=5, sticky=Tk.W)

        # Plot options
        self.var = Tk.IntVar()
        c1 = Tk.Radiobutton(
            options_frame, text="Shared axes", variable=self.var,
            command=self.update_graph, font=font, value=0)
        c1.grid(row=0, column=0, sticky=Tk.W)

        c2 = Tk.Radiobutton(
            options_frame, text="Multiple axes", variable=self.var,
            command=self.update_graph, font=font, value=1)
        c2.grid(row=1, column=0, sticky=Tk.W)

        c3 = Tk.Radiobutton(
            options_frame, text="Stacked plots", variable=self.var,
            command=self.update_graph, font=font, value=2)
        c3.grid(row=2, column=0, sticky=Tk.W)

        self.var_del = Tk.IntVar()
        c4 = Tk.Checkbutton(
            options_frame,
            text="Absolute delta values",
            variable=self.var_del,
            command=self.update_graph,
            font=font)
        c4.grid(row=0, column=1, sticky=Tk.W)

        self.var_log = Tk.IntVar()
        c5 = Tk.Checkbutton(
            options_frame,
            text="Log scale",
            variable=self.var_log,
            command=self.update_graph,
            font=font)
        c5.grid(row=1, column=1, sticky=Tk.W)

        # Option to only show the min and max of array-type variables
        self.var_minmax = Tk.IntVar()
        c6 = Tk.Checkbutton(
            options_frame,
            text="Min/max for arrays",
            variable=self.var_minmax,
            command=self.update_graph,
            font=font)
        c6.grid(row=2, column=1, sticky=Tk.W)

        # Option to show all values for array-type variables
        self.var_showall = Tk.IntVar()
        c7 = Tk.Checkbutton(
            options_frame,
            text="Show all for arrays",
            variable=self.var_showall,
            command=self.update_graph,
            font=font)
        c7.grid(row=3, column=1, sticky=Tk.W)

            # Option to show legend
        self.var_legend = Tk.IntVar()
        c8 = Tk.Checkbutton(
            options_frame,
            text="Show legend",
            variable=self.var_legend,
            command=self.update_graph,
            font=font)
        self.var_legend.set(1)
        c8.grid(row=4, column=1, sticky=Tk.W)

        # Option to show bounds
        self.var_bounds = Tk.IntVar()
        c9 = Tk.Checkbutton(
            options_frame,
            text="Show bounds",
            variable=self.var_bounds,
            command=self.update_graph,
            font=font)
        c9.grid(row=5, column=1, sticky=Tk.W)

        # Option to only show major iterations
        self.var_mask = Tk.IntVar()
        c10 = Tk.Checkbutton(
            options_frame,
            text="Show major iterations",
            variable=self.var_mask,
            command=self.set_mask,
            font=font)
        c10.grid(row=6, column=1, sticky=Tk.W, pady=6)

        # Option to automatically refresh history file
        # especially useful for running optimizations
        self.var_ref = Tk.IntVar()
        c11 = Tk.Checkbutton(
            options_frame,
            text="Automatically refresh",
            variable=self.var_ref,
            command=self.auto_ref,
            font=font)
        c11.grid(row=7, column=1, sticky=Tk.W, pady=6)

        lab = Tk.Label(
            options_frame,
            text="Search for a function/variable:",
            font=font)
        lab.grid(row=8, column=0, columnspan=2, pady=10, sticky=Tk.W)

        # Search box to filter displayed functions/variables
        vs = Tk.StringVar()
        vs.trace("w", lambda name, index, mode, vs=vs: self.var_search(vs))
        self.entry_search = Tk.Entry(
            options_frame, text="Search", textvariable=vs,
            font=font)
        self.entry_search.grid(row=8, column=2, pady=10, sticky=Tk.W)

        lab_font = Tk.Label(
            options_frame,
            text="Font size for plots:",
            font=font)
        lab_font.grid(row=9, column=0, sticky=Tk.S)

        w = Tk.Scale(options_frame, from_=6, to=24, orient=Tk.HORIZONTAL,
                     resolution=2, command=self.update_font, font=font)
        w.set(16)
        w.grid(row=9, column=1)

if __name__ == '__main__':
    # Called only if this script is run as main.

    # ======================================================================
    # Input Information
    # ======================================================================
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'histFile', nargs='*', type=str, default='opt_hist.hst',
        help="Specify the history file to be plotted")
    parser.add_argument('--output', nargs='?', type=str, default='./',
                        help="Specify the output directory")
    args = parser.parse_args()
    histList = args.histFile
    outputDir = args.output

    histFileName = histList[0]

    # Check that the output directory is available. Create it if not
    if not os.path.isdir(outputDir):
        os.makedirs(outputDir)
    # Initialize display parameters, obtain history, and draw GUI
    disp = Display(histList, outputDir)
    disp.draw_GUI()
    disp.root.protocol("WM_DELETE_WINDOW", disp.quit)
    on_move_id = disp.f.canvas.mpl_connect('motion_notify_event', disp.on_move)
    Tk.mainloop()
