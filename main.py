from RTTransmittivity import Analysis
from gui import GUI
import threading
from lecroypy import LeCroyWaveRunner 

osc = LeCroyWaveRunner('USB0::0x05FF::0x1023::LCRY3704N16293::INSTR') 
print(osc.idn) # confirms connection by pinging the oscilloscope for identification, example: LECROY, WAVESURFER, LCRY3704N16293, 9.1.0


gui = GUI() 

def dataProcessing():
    analysis = Analysis('C1trace0', 'C3trace0',thresh=0.01, dataStart=343) 
    osc.wait_for_single_trigger()
    dataC1 = osc.get_waveform(n_channel=1) 
    dataC3 = osc.get_waveform(n_channel=3) 
    # take data from dataC1 and dataC3 and passes them into the pulse fitting algorithm
    analysis.runAnalysis(data1Array=[dataC1['waveforms'][0]['Time (s)'], dataC1['waveforms'][0]['Amplitude (V)']],
                         data2Array=[dataC3['waveforms'][0]['Time (s)'], dataC3['waveforms'][0]['Amplitude (V)']])

    gui.ratioPlotToCSV = analysis.outputPlots[4]
    for j in range(2):
        gui.addWidgets(analysis.outputPlots[j], "reference", j)
        gui.addWidgets(analysis.outputPlots[j+2], "data", j)
    gui.addWidgets(analysis.outputPlots[4], "ratio", 0)

    gui.root.after(0,gui.updateGUI()) # update the GUI
gui.dataAcquisitionFunction = dataProcessing # assigns the dataProcessing() function to the start button of the GUI

gui.setGridding() 
gui.begin() # starts GUI event loop