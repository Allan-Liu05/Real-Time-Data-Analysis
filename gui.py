import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
import threading
import pandas as pd
import os
import re

class GUI:
    def __init__(self, numReferenceWidgets=2, numDataWidgets=2):
        # the plots for each widget should be stores as an element in the array except for ratio
        self.numReferenceWidgets = [None] * numReferenceWidgets
        self.numDataWidgets = [None] * numDataWidgets
        self.ratioPlot = None
        
        # setting up the main window that everything else is contained in
        self.root = tk.Tk()
        self.root.title("RT Transmission Analysis")
        self.root.state("zoomed")
        
        # initilises the main frame and prepares the locations for which the widgets that contain the plots to reside
        self.mainFrame = tk.Frame(self.root)
        self.mainFrame.grid(row=0, column=0, sticky="nsew")

        # configure a grid for the ROOT window, weight = 0 to not allow resiving of grids, weight = 1 to allow for dynamic resizing
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)
        
        self.buttonFrame = tk.Frame(self.root)
        self.buttonFrame.grid(row=1, column=0, columnspan=2, pady=0, padx=0)

        # "Save Plots" button
        self.guiButton = tk.Button(
            self.buttonFrame, 
            text="Save Plots", 
            command=self.savePlots,
            width=20,
            height=1,
            font=("Arial", 20)
        )


 
        self.startButton = tk.Button(
            self.buttonFrame, 
            text="Start", 
            command= self.startDataAcquisition,
            width=20,
            height=1,
            font=("Arial", 20)
        )

        # configuring grid for start and save buttons
        self.guiButton.grid(row=0, column=1, padx=5, pady=0)
        self.startButton.grid(row=0, column=0, padx=5, pady=0)

        self.buttonFrame.grid_columnconfigure(0, weight=1)  # Let the buttons stretch to fill screen
        self.buttonFrame.grid_columnconfigure(1, weight=1) 

        self.numRows = 3 # Defines number of rows and columns in GUI
        self.numColumns = 2 

        self.dataAcquisitionFunction = None 

        # variables for the .csv and .png outputs of the Analysis class.
        # saved to a target folder using the save button
        self.pngOutput = None
        self.csvOutput = None 

        # temporary variable that stores the x,y values of ratio plot for optional save
        self.ratioPlotToCSV = None
        self.ratioPlotToPNG = None
    # this function is used to update the GUI with plots after being run through pulse fitting algorithm
    def updateGUI(self):
        for i in range(len(self.numReferenceWidgets)):
            widget = self.numReferenceWidgets[i].get_tk_widget()
            widget.grid(row=0, column=i, sticky="nsew")  
        
        for i in range(len(self.numDataWidgets)):
            widget = self.numDataWidgets[i].get_tk_widget()
            widget.grid(row=1, column=i, sticky="nsew")  

        if self.ratioPlot:
            self.ratioPlot.get_tk_widget().grid(row=2, column=0, columnspan=2, sticky="nsew")  

    # prepares either reference, data, or ratio frames
    def addWidgets(self, figure, location, index):
        if location.lower() == "reference":
            self.numReferenceWidgets[index] = FigureCanvasTkAgg(figure, self.mainFrame)
        elif location.lower() == "data":
            self.numDataWidgets[index] = FigureCanvasTkAgg(figure, self.mainFrame)
        elif location.lower() == "ratio":
            self.ratioPlot = FigureCanvasTkAgg(figure, self.mainFrame)
            self.ratioPlotToPNG = figure # extra step for optional save of the ratio plot
        else:
            raise Exception("Location must be 'reference', 'data', or 'ratio'.")

    def setGridding(self):
        # Configure grid rows and columns to accommodate the 2-2-1 layout
        self.mainFrame.grid_rowconfigure(0, weight=1)  # Reference 
        self.mainFrame.grid_rowconfigure(1, weight=1)  # Data 
        self.mainFrame.grid_rowconfigure(2, weight=1)  # Ratio

        self.mainFrame.grid_columnconfigure(0, weight=1)  # For Reference and Data only
        self.mainFrame.grid_columnconfigure(1, weight=1)  

    # creates a new thread that runs whatever function is assigned to self.dataAcquisitionFunction
    def startDataAcquisition(self):
        threading.Thread(target=self.dataAcquisitionFunction).start()

    def savePNG(self):
        os.makeDirs('results', exist_ok=True)
    
    # saves ratioPlot as a .csv
    def saveCSV(self):
        # figures out number of file (transmisison ratio1, transmission ratio2, etc.)
        files = [f for f in os.listdir('results/') if f.startswith('transmission ratio') and f.endswith('.csv')]
        fileNumber = None
        for fileName in files:
            if fileName == "transmission ratio":
                fileNumber = 0
            match = re.search(rf'transmission ratio(\d+)\.csv', fileName)
            if match:
                fileNumber = int(match.group(1))
        # defaults save file to a directory called results 
        out = np.array(self.ratioPlotToCSV).T
        np.savetxt('results/transmission ratio'+str(fileNumber), out, delimiter=',')

    def savePlots(self):
        # figures out number of file (transmisison ratio1, transmission ratio2, etc.)
        files = [f for f in os.listdir('results/') if f.startswith('transmission ratio') and f.endswith('.csv')]
        fileNumber = None
        for fileName in files:
            if fileName == "transmission ratio":
                fileNumber = 0
            match = re.search(rf'transmission ratio(\d+)\.csv', fileName)
            if match:
                fileNumber = int(match.group(1))
        
        # defaults save file to a directory called results 
        out = np.array(self.ratioPlotToCSV).T
        np.savetxt('results/transmission ratio'+str(fileNumber), out, delimiter=',') # saving ratio as .csv
        self.ratioPlotToPNG.savefig('transmission ratio'+str(filNumber)+'.png')

    def begin(self):
        # starts the GUI
        self.root.mainloop()


# testing class, save button does not work using the code below as the save button is designed to work with data 
# from RTTransmittivity.Analysis object

if __name__ == "__main__":

    # Example data
    x = np.random.rand(10)
    y = np.random.rand(10)

    # Create figures
    figureRef = Figure(figsize=(5, 4), dpi=100)
    ax1 = figureRef.add_subplot(111)
    ax1.plot(x, y, label="Reference Trace")
    ax1.legend()

    figureData = Figure(figsize=(5, 4), dpi=100)
    ax2 = figureData.add_subplot(111)
    ax2.plot(y, x, label="Data Trace")
    ax2.legend()

    display = GUI()
    display.setGridding()

    # Add reference and data widgets
    for i in range(2):
        display.addWidgets(figureRef, "reference", i)
        display.addWidgets(figureData, "data", i)

    # Add a ratio plot (using reference figure as a placeholder)
    display.addWidgets(figureRef, "ratio", 0)

    # Update the GUI to reflect the added plots
    display.updateGUI()

    # Start the GUI application
    display.begin()