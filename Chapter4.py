import math
import sys

import PyQt5.QtWidgets as QtWidgets

import pyqtgraph.opengl as gl
from PyQt5.QtGui import QSurfaceFormat

# Force OpenGL Compatibility Profile
fmt = QSurfaceFormat()
fmt.setRenderableType(QSurfaceFormat.OpenGL)
fmt.setProfile(QSurfaceFormat.CompatibilityProfile)  # Allows legacy OpenGL
fmt.setVersion(2, 1)  # Use OpenGL 2.1 (widely compatible)
QSurfaceFormat.setDefaultFormat(fmt)

import ece163.Display.baseInterface as baseInterface
import ece163.Containers.Inputs as Inputs
import ece163.Display.GridVariablePlotter
import ece163.Display.SliderWithValue
import ece163.Simulation.Chapter4Simulate
import ece163.Display.DataExport
import ece163.Display.WindControl as WindControl

from ece163.Utilities.Joystick import Joystick
from ece163.Constants import JoystickConstants as JSC


stateNamesofInterest = ['pn', 'pe', 'pd', 'yaw', 'pitch', 'roll', 'u', 'v', 'w', 'p', 'q', 'r', 'alpha', 'beta']
systemInputs = [('Throttle', 0, 1, 0),
				('Aileron', -0.3, 0.3, 0),
				('Elevator', -0.3, 0.3, 0),
				('Rudder', -0.3, 0.3, 0)]

positionRange = 200

class Chapter4(baseInterface.baseInterface):
	def __init__(self, parent=None):
		self.simulateInstance = ece163.Simulation.Chapter4Simulate.Chapter4Simulate()
		super().__init__(parent)
		self.setWindowTitle("ECE163 Chapter 4")
		plotElements = [[x] for x in stateNamesofInterest]
		plotElements.append(['Va', 'Vg'])
		titleNames = list(stateNamesofInterest)
		titleNames.append('Va & Vg')
		legends = [False] * len(stateNamesofInterest) + [True]
		self.stateGrid = ece163.Display.GridVariablePlotter.GridVariablePlotter(5, 3, plotElements, titles=titleNames, useLegends=legends)

		self.outPutTabs.addTab(self.stateGrid, "States")
		self.outPutTabs.setCurrentIndex(2)
		self.stateUpdateDefList.append(self.updateStatePlots)

		self.exportWidget = ece163.Display.DataExport.DataExport(self.simulateInstance,'Chapter4')
		self.outPutTabs.addTab(self.exportWidget, "Export Data")

		self.inputControlsWidget = QtWidgets.QWidget()
		gridSquish = QtWidgets.QVBoxLayout()
		self.inputGrid = QtWidgets.QGridLayout()
		gridSquish.addLayout(self.inputGrid)
		gridSquish.addStretch()
		self.inputControlsWidget.setLayout(gridSquish)
		self.inputTabs.addTab(self.inputControlsWidget, "Control Inputs")
		# self.inputLayout.addLayout(self.inputGrid)
		self.windControl = WindControl.WindControl(self.simulateInstance.underlyingModel)
		self.inputTabs.addTab(self.windControl, WindControl.widgetName)

		resetSlidersButton = QtWidgets.QPushButton("Reset Sliders")
		gridSquish.addWidget(resetSlidersButton)
		resetSlidersButton.clicked.connect(self.resetSliders)

		self.inputLayout.addStretch()
		self.inputSliders = list()

		for row in range(2):
			for col in range(2):
				index = col+row*2
				name, minValue, maxValue, startValue = systemInputs[index]
				newSlider = ece163.Display.SliderWithValue.SliderWithValue(name, minValue, maxValue, startValue)
				self.inputSliders.append(newSlider)
				self.inputGrid.addWidget(newSlider, row, col)

		# self.playButton.setDisabled(True)
		self.showMaximized()

		self.joystick = Joystick()

		#For convenience, if a controller is active, go directly to the wind tab
		if self.joystick.active:
			self.inputTabs.setCurrentIndex(1)

		####Simulation Update code###
		# # Updates the simulation when tab is being changed
		self.outPutTabs.currentChanged.connect(self.newTabClicked)
		self.outPutTabs.setCurrentIndex(0)
		self.plotWidgets = [self.stateGrid]
		# Default for all graphs to be turned off
		self.updatePlotsOn()
		self.updatePlotsOff()
		# Overwrite simulationTimedThread function with modified sliderChangeResponse
		self.simulationTimedThread.timeout.connect(self.UpdateSimulationPlots)



		return

	def resetSliders(self):
		for slider in self.inputSliders:
			slider.resetSlider()
		return

	def updateStatePlots(self, newState):
		self.updatePlotsOff()
		stateList = list()
		for key in stateNamesofInterest:
			newVal = getattr(newState, key)
			if key in ['yaw', 'pitch', 'roll', 'p', 'q', 'r', 'alpha', 'beta']:
				newVal = math.degrees(newVal)
			stateList.append([newVal])
		stateList.append([newState.Va, math.hypot(newState.u, newState.v, newState.w)])

		self.stateGrid.addNewAllData(stateList, [self.simulateInstance.time]*(len(stateNamesofInterest) + 1))
		self.updatePlotsOn()
		return

	def getVehicleState(self):
		return self.simulateInstance.underlyingModel.getVehicleState()

	def runUpdate(self):
		inputControls = Inputs.controlInputs()
		
		#If a controller was initialized, use the values for input
		if self.joystick.active:
			joystick_vals = self.joystick.get_joystick_values().control_axes
			inputControls.Aileron = joystick_vals.Aileron * JSC.CHAPTER4_MAX_THROW
			inputControls.Elevator = joystick_vals.Elevator * JSC.CHAPTER4_MAX_THROW
			inputControls.Rudder = -joystick_vals.Rudder * JSC.CHAPTER4_MAX_THROW
			inputControls.Throttle = joystick_vals.Throttle
			for slider in self.inputSliders:
				slider.setSlider(getattr(inputControls,slider.name))
		
		else:
			for control in self.inputSliders:
				setattr(inputControls, control.name, control.curValue)
		
		self.simulateInstance.takeStep(inputControls)

		return

	# def sliderChangeResponse(self, newValue, name):
	# 	if name in ['yaw', 'pitch', 'roll']:
	# 		setattr(self.vehicleState, name, math.radians(newValue))
	# 	else:
	# 		if name == 'z':
	# 			newValue = -1*newValue
	# 		setattr(self.vehicleState, name, newValue)
	# 	self.runSimulation()
	# 	return

	def resetSimulationActions(self):
		self.simulateInstance.reset()
		self.stateGrid.clearDataPointsAll()
		self.vehicleInstance.reset(self.simulateInstance.underlyingModel.getVehicleState())
		self.outPutTabs.setCurrentIndex(0)

	#### Simulation Update Code ##########

	# Updates a simulation widget when new tab clicked
	def UpdateSimulationPlots(self):

		currentWidget = self.outPutTabs.currentWidget()
		# Ensure that that the timer is only enabled for states, sensors, and control response widgets
		if (currentWidget in self.plotWidgets):
			#self.runUpdate()
			self.updatePlotsOn()
			self.updatePlotsOff()
		return
	def newTabClicked(self):
		self.updatePlotsOn()
		self.updatePlotsOff()
		return

	# toggles the state grid widget
	def togglestateGridPlot(self, toggleIn):
		self.stateGrid.setUpdatesEnabled(toggleIn)
		return

	# Turns on all simulation plots
	def updatePlotsOn(self):
		# print("Turning on plot update")
		self.togglestateGridPlot(True)
		return

	# Turns off all simulation plots
	def updatePlotsOff(self):
		# print("Turning off plot update")
		self.togglestateGridPlot(False)
		return

sys._excepthook = sys.excepthook

def my_exception_hook(exctype, value, tracevalue):
	# Print the error and traceback
	import traceback
	with open("LastCrash.txt", 'w') as f:
		traceback.print_exception(exctype, value, tracevalue, file=f)
		# traceback.print_tb(tracevalue, file=f)
	print(exctype, value, tracevalue)
	# Call the normal Exception hook after
	sys._excepthook(exctype, value, tracevalue)
	sys.exit(0)

# Set the exception hook to our wrapping function
sys.excepthook = my_exception_hook



qtApp = QtWidgets.QApplication(sys.argv)
ourWindow = Chapter4()
ourWindow.show()
qtApp.exec()