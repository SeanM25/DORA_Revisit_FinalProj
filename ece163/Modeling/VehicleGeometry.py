"""
Holds the vehicle graphics, only operation on it is to return a set of points and meshes with the appropriate rotation/translation
currently just returns the modified points, does not update the base ones. Module uses its baseUnit variable to self.SCale the model to an
arbitrary size for good rendering in the display window.
"""

from ..Utilities import MatrixMath
from ..Utilities import Rotations
from ..Constants import VehiclePhysicalConstants as VPC
import math

baseUnit = 1.0

class VehicleGeometry():
	def __init__(self):
		"""
		defines the vehicle in NED coordinates around the local body frame origin. Rotations and translations will be
		around the [0,0,0] point of this local frame. Has to be in NED or the rotation matrices will not work. Vehicle is
		self.SCaled to match the wing span of the actual vehicle, this all the points are in meters.

		"vertices" is an [n x 3] matrix of xyz points for each vertex of the vehicle;
		"faces" is an [m x3] index matrix of which vertices connect to which face (only triangles allowed for faces);
		"colors" is an [m x 4] matrix where each row is a CMYK definition of that face color
		"""


		red = [1., 0., 0., 1]
		green = [0., 1., 0., 1]
		blue = [0., 0., 1., 1]
		yellow = [1., 1., 0., 1]
		white = [1.,1.,1.,1]
		black = [0., 0., 0., 0.]

		# Add more colors ?

		# squat squarish vehicle as a more complicated practice

		# self.vertices = [[baseUnit,0,0], # [0] nose
		# 				 [baseUnit/2,-baseUnit/2,-baseUnit/4],   # [1] front top left
		# 				 [-baseUnit/2,-baseUnit/2,-baseUnit/4],  # [2] rear top left
		# 				 [-baseUnit/2,baseUnit/2,-baseUnit/4],   # [3] rear top right
		# 				 [baseUnit/2,baseUnit/2,-baseUnit/4],     # [4] front top right
		# 				 [baseUnit/2,baseUnit/2,baseUnit/4],    # [5] front bottom right
		# 				 [baseUnit/2,-baseUnit/2,baseUnit/4],   # [6] front bottom left
		# 				 [-baseUnit/2,-baseUnit/2,baseUnit/4],   # [7] rear bottom left
		# 				 [-baseUnit/2,baseUnit/2,baseUnit/4],    # [8] rear bottom right
		# 				 [-baseUnit/2,0,-0.75*baseUnit],			 # [9] fin top point
		# 				 [-baseUnit/2,0,0],                      # [10] fin rear point
		# 				 [0,0,0]]                      # [11] fin front point
		#
		# self.faces = [[0,1,4],         # [0] nose top
		# 			  [1,3,4],[1,2,3], # [1],[2] body top
		# 			  [0,1,6],         # [3] nose left
		# 			  [1,7,6],[1,7,2], # [4],[5] body left
		# 			  [2,3,8],[2,8,7], # [6],[7] body rear
		# 			  [3,4,8],[8,4,5], # [8],[9] body right
		# 			  [5,0,4],         # [10] nose right
		# 			  [0,5,6],         # [11] nose bottom
		# 			  [6,8,5],[6,8,7], # [12],[13] body bottom
		# 			  [9,10,11]]       # [14] fin
		#
		# self.colors=[yellow,yellow,yellow, # top
		# 			 blue,blue,blue, # left
		# 			 red,red, # back
		#              green,green,green, # right
		# 			 white,white,white, # bottom
		# 			 blue] # tail

		# actual MAV model from Beard Chapter 2

		# define MAV body parameters
		self.SC = 0.001 # self.SCalingUnit # VPC.b / 6.0	# self.SCaling determined by the wingspan of the aircraft in VehiclePhysicalConstants

		# fuse_h = self.SCalingUnit
		# fuse_w = self.SCalingUnit
		# fuse_l1 = 2 * self.SCalingUnit
		# fuse_l2 = self.SCalingUnit
		# fuse_l3 = 4 * self.SCalingUnit
		# wing_l = self.SCalingUnit
		# wing_w = 6 * self.SCalingUnit	# will match the wingspan (VPC.b) in meters exactly
		# tail_h = self.SCalingUnit
		# tail_l = self.SCalingUnit
		# tail_w = 2 * self.SCalingUnit

		wing_length = 700 # mm

		wing_depth = 100 # mm

		# Editable flap coordinates
		self.aelronFlapRadius = 60*self.SC # CONST: for calculating the angle based coordinates
		self.aeleronFlapRAngle = -30 # adjusts based on controller: for calculating the angle based coordinates
		self.aeleronFlapLAngle = 30 # adjusts based on controller: for calculating the angle based coordinates
		
		# Where R = radius of aeleron; math: x = x_0 - (R*cos(theta) - R), y = y, z = R*sin(theta)
		self.aeleron_R_right = [-160*self.SC - (self.aelronFlapRadius*math.cos(self.aeleronFlapRAngle*math.pi/180) - self.aelronFlapRadius), 600*self.SC, self.aelronFlapRadius*math.sin(self.aeleronFlapRAngle*math.pi/180)] # back right of aeleronn flap R | default: [-160*self.SC, 600*self.SC, 0.01]
		self.aeleron_R_left = [-160*self.SC - (self.aelronFlapRadius*math.cos(self.aeleronFlapRAngle*math.pi/180) - self.aelronFlapRadius), 200*self.SC, self.aelronFlapRadius*math.sin(self.aeleronFlapRAngle*math.pi/180)] # back left of aeleronn flap R | default: [-160*self.SC, 200*self.SC, 0.01]
		self.aeleron_L_right = [-160*self.SC - (self.aelronFlapRadius*math.cos(self.aeleronFlapRAngle*math.pi/180) - self.aelronFlapRadius), -200*self.SC, self.aelronFlapRadius*math.sin(self.aeleronFlapLAngle*math.pi/180)] # back right of aeleron flap L | default: [-160*self.SC, -200*self.SC, 0.01]
		self.aeleron_L_left = [-160*self.SC - (self.aelronFlapRadius*math.cos(self.aeleronFlapRAngle*math.pi/180) - self.aelronFlapRadius), -600*self.SC, self.aelronFlapRadius*math.sin(self.aeleronFlapLAngle*math.pi/180)] # back left of aeleron flap L | default: [-160*self.SC, -600*self.SC, 0.01]

		self.vertices = [
						
						# Main Wing Body
						
						[0, -700*self.SC, 0], # [0] Front Left Corner of Wing 
				   		[0, 700*self.SC, 0], # [1] Front Right Corner of Wing 
						[-100*self.SC, 700*self.SC, 0], # [2] Back Right Corner of Wing / front right of solid Wing R
						[-100*self.SC, -700*self.SC, 0], # [3] Back Left Corner of Wing / front left of solid Wing L

						# Solid Wing R of Wing

						[-160*self.SC, 700*self.SC, 0], # [4] back right corner of Solid Wing R
						[-160*self.SC, 600*self.SC, 0], # [5] back left corner of Solid Wing R
						[-100*self.SC, 600*self.SC, 0], # [6] front left corner of Solid Wing R / front right of aeleronn flap R

						# aeleron Flap R: adjustable

						self.aeleron_R_right, # [7] back right of aeleronn flap R
						self.aeleron_R_left, # [8] back left of aeleronn flap R

						# Solid Wing C of Wing

						[-100*self.SC, 200*self.SC, 0], # [9] front right of Solid Wing C / front left of aeleronn flap R
						[-160*self.SC, 200*self.SC, 0], # [10] back right of Solid Wing C
						[-160*self.SC, -200*self.SC, 0], # [11] back left of Solid Wing C
						[-100*self.SC, -200*self.SC, 0], # [12] front left of Solid Wing C / front right of aeleron flap L

						# aeleron Flap L: adjustable
						
						self.aeleron_L_right, # [13] back right of aeleron flap L
						self.aeleron_L_left, # [14] back left of aeleron flap L

						# Solid Wing L (of wing)

						[-100*self.SC, -600*self.SC, 0], # [15] front right of Solid Wing L / front left of aeleron flap L
						[-160*self.SC, -600*self.SC, 0], # [16] back right of Solid Wing L
						[-160*self.SC, -700*self.SC, 0], # [17] back left of Solid Wing L

						# Body Top

						[160*self.SC, 41*self.SC, 0], # [18]

						[160*self.SC, -41*self.SC, 0], # [19]

						[-578*self.SC, -41*self.SC, 0], # [20]

						[-578*self.SC, -41*self.SC, 0], # [21]

						[-578*self.SC, 41*self.SC, 0], # [22]

						[160*self.SC, 41*self.SC, 0], # [23]

						# Body Bottom

						[95*self.SC, 41*self.SC, 102*self.SC], # [24]

						[95*self.SC, -41*self.SC, 102*self.SC], # [25]

						[-220*self.SC, -41*self.SC, 102*self.SC], # [26]


						[-220*self.SC, -41*self.SC, 102*self.SC], # [27]

						[-220*self.SC, 41*self.SC, 102*self.SC], # [28]

						[95*self.SC, 41*self.SC, 102*self.SC], # [29]

						[-220*self.SC, -41*self.SC, 102*self.SC], # [30]

						[-220*self.SC, 41*self.SC, 102*self.SC], # [31]

						[-578*self.SC, -41*self.SC, 37*self.SC], # [32]


						[-578*self.SC, -41*self.SC, 37*self.SC], # [33]

						[-578*self.SC, 41*self.SC, 37*self.SC], # [34]

						[-220*self.SC, 41*self.SC, 102*self.SC], # [35]


						[-578*self.SC, -41*self.SC, 0], # [36]

						[-578*self.SC, -41*self.SC, 37*self.SC], # [37]

						[-578*self.SC, 41*self.SC, 37*self.SC], # [38]


						[-578*self.SC, 41*self.SC, 37*self.SC], # [39]

						[-578*self.SC, 41*self.SC, 0], # [40]

						[-578*self.SC, -41*self.SC, 0], # [41]

						# Body Left Side

						[160*self.SC, -41*self.SC, 37*self.SC], # [42]

						[95*self.SC, -41*self.SC, 102*self.SC], # [43]

						[95*self.SC, -41*self.SC, 37*self.SC], # [44]


					  [95*SC, -41*SC, 102*SC], # [45]

					  [95*SC, -41*SC, 37*SC], # [46]

					  [-220*SC, -41*SC, 37*SC], # [47]


					  [-220*SC, -41*SC, 37*SC], # [48]

					  [-220*SC, -41*SC, 102*SC], # [49]

					  [95*SC, -41*SC, 102*SC], # [50]


					  [-220*SC, -41*SC, 102*SC], # [51]

					  [-220*SC, -41*SC, 37*SC], # [52]

					  [-578*SC, -41*SC, 37*SC] # [53]


















				   
				   ]
		

		self.faces = [

					# Main Wing

					  [0, 1, 2], # Connect Front Left, Front Right, and Back Right Wing vertices
				
					  [0, 3, 2], # Connect Front Left, Back Left, and Back Right Wing vertices

					  # Solid Wing R - green

					  [6, 2, 4], # Connect Front Left, Front Right, and Back Right Solid Wing R vertices

					  [6, 5, 4], # Connect Front Left, Back Left, and Back Right Solid Wing R vertices

					  # Flap R

					  [9, 6, 7], # Connect Front Left, Front Right, and Back Right Flap R vertices

					  [9, 8, 7], # Connect Front Left, Back Left, and Back Right Flap vertices

					  # Solid Wing C

					  [12, 9, 10], # Connect Front Left, Front Right, and Back Right Solid Wing C vertices

					  [12, 11, 10], # Connect Front Left, Back Left, and Back Right Solid Wing C vertices

					  # Flap L

					  [15, 12, 13], # Connect Front Left, Front Right, and Back Right Flap L vertices

					  [15, 14, 13], # Connect Front Left, Back Left, and Back Right Flap L vertices

					  # Solid Wing L

					  [3, 15, 16], # Connect Front Left, Front Right, and Back Right Solid Wing L vertices

					  [3, 17, 16], # Connect Front Left, Back Left, and Back Right Solid Wing L vertices


					  # Body Top 

					  [18, 19, 20],

					  [21, 22, 23],

					  # Body Bottom

					  [24, 25, 26],

					  [27, 28, 29],

					  [30, 31, 32],

					  [33, 34, 35],

					  [36, 37, 38],

					  [39, 40, 41],

					  # Body Left Side

					  [42, 43, 44],

					  [45, 46, 47],

					  [48, 49, 50],

					  [51, 52, 53]



					]
		

		self.colors = [green, green, green, green, red, red, green, green, red, red, green, green, white, white, blue, blue, blue, blue, blue, blue, blue, blue, blue, blue]

		'''

		self.vertices = [[fuse_l1, 0, 0],  # point 1 [0]
						 [fuse_l2, fuse_w / 2.0, -fuse_h / 2.0],  # point 2 [1]
						 [fuse_l2, -fuse_w / 2.0, -fuse_h / 2.0],  # point 3 [2]
						 [fuse_l2, -fuse_w / 2.0, fuse_h / 2.0],  # point 4 [3]
						 [fuse_l2, fuse_w / 2.0, fuse_h / 2.0],  # point 5 [4]
						 [-fuse_l3, 0, 0],  # point 6 [5]
						 [0, wing_w / 2.0, 0],  # point 7 [6]
						 [-wing_l, wing_w / 2.0, 0],  # point 8 [7]
						 [-wing_l, -wing_w / 2.0, 0],  # point 9 [8]
						 [0, -wing_w / 2.0, 0],  # point 10 [9]
						 [-fuse_l3 + tail_l, tail_w / 2.0, 0],  # point 11 [10]
						 [-fuse_l3, tail_w / 2.0, 0],  # point 12 [11]
						 [-fuse_l3, -tail_w / 2.0, 0],  # point 13 [12]
						 [-fuse_l3 + tail_l, -tail_w / 2.0, 0],  # point 14 [13]
						 [-fuse_l3 + tail_l, 0, 0],  # point 15 [14]
						 [-fuse_l3, 0, -tail_h]]  # point 16 [15]

		self.faces = [[0, 1, 2],
					  [0, 1, 4],
					  [0, 3, 4],
					  [0, 3, 2],
					  [5, 2, 3],
					  [5, 1, 2],
					  [5, 1, 4],
					  [5, 3, 4],
					  [6, 7, 9], # Wing Orig 6 7 9
					  [7, 8, 9], # Wing Orig 7 8 9
					  [10, 11, 12],
					  [10, 12, 13],
					  [5, 14, 15]]
					  

		self.colors = [yellow, yellow, yellow, yellow,
					   blue, blue, blue,
					   red,
					   green, green, green, green,
					   blue]
'''
		return

	def updateAngles(self, Throttle=0.0, Aileron=0.0, Elevator=0.0, Rudder=0.0):
		
		self.aeleronFlapRAngle = -Aileron # adjusts based on controller: for calculating the angle based coordinates
		self.aeleronFlapLAngle = Aileron # adjusts based on controller: for calculating the angle based coordinates
		print(self.aeleron_R_right)
		# Where R = radius of aeleron; math: x = x_0 - (R*cos(theta) - R), y = y, z = R*sin(theta)
		self.aeleron_R_right = [-160*self.SC - (self.aelronFlapRadius*math.cos(self.aeleronFlapRAngle*math.pi/180) - self.aelronFlapRadius), 600*self.SC, self.aelronFlapRadius*math.sin(self.aeleronFlapRAngle*math.pi/180)] # back right of aeleronn flap R | default: [-160*self.SC, 600*self.SC, 0.01]
		self.aeleron_R_left = [-160*self.SC - (self.aelronFlapRadius*math.cos(self.aeleronFlapRAngle*math.pi/180) - self.aelronFlapRadius), 200*self.SC, self.aelronFlapRadius*math.sin(self.aeleronFlapRAngle*math.pi/180)] # back left of aeleronn flap R | default: [-160*self.SC, 200*self.SC, 0.01]
		self.aeleron_L_right = [-160*self.SC - (self.aelronFlapRadius*math.cos(self.aeleronFlapRAngle*math.pi/180) - self.aelronFlapRadius), -200*self.SC, self.aelronFlapRadius*math.sin(self.aeleronFlapLAngle*math.pi/180)] # back right of aeleron flap L | default: [-160*self.SC, -200*self.SC, 0.01]
		self.aeleron_L_left = [-160*self.SC - (self.aelronFlapRadius*math.cos(self.aeleronFlapRAngle*math.pi/180) - self.aelronFlapRadius), -600*self.SC, self.aelronFlapRadius*math.sin(self.aeleronFlapLAngle*math.pi/180)] # back left of aeleron flap L | default: [-160*self.SC, -600*self.SC, 0.01]
		print(self.aeleron_R_right)

	def getNewPoints(self, x, y, z, yaw, pitch, roll, Throttle=0.0, Aileron=0.0, Elevator=0.0, Rudder=0.0):
		"""
		Function to get new ENU points of the vehicle in inertial space from Euler angles, NED displacements, and base
		drawing contained within the __init__ function. That is, points to be remapped are contained within self.vertices

		:param x: North Displacement (Pn) in [m]
		:param y: East Displacement (Pe) in [m]
		:param z: Down Displacement (Pd) in [m]
		:param yaw: rotation about inertial down [rad]
		:param pitch: rotation about intermediate y-axis [rad]
		:param roll: rotation about body x-axis [rad]
		:return: Points in inertial EAST-NORTH-UP frame (for plotting)
		"""
		# EXPERIMENTAL Test!!!
		self.updateAngles(Throttle=Throttle, Aileron=Aileron, Elevator=Elevator, Rudder=Rudder)
		self.vertices[7] = self.aeleron_R_right
		self.vertices[8] = self.aeleron_R_left
		self.vertices[13] = self.aeleron_L_right
		self.vertices[14] = self.aeleron_L_left

		newPoints = self.vertices

		rmatrix = Rotations.euler2DCM(yaw, pitch, roll) # Get the rotation matrix from the euler angles using euler2dcm()

		coords_rotated = MatrixMath.multiply(newPoints,rmatrix) # multiply aircraft points by the rotation matrix to apply the appropriate rotation in NED

		displacements = [x,y,z] # Given x, y, z NED displacements

		for row in range(len(coords_rotated[0])): # Go through the rows of the vehicle points matrix

			for col in range((len(coords_rotated))): # Go through the cols of the vehicle points matrix

				coords_rotated[col][row] += displacements[row] # Add the x displacement to col 1, y to 2, and z to 3 for each row in nx3
		
		newPoints = Rotations.ned2enu(coords_rotated) # Convert dispalced NED matrix to ENU coords

		return newPoints # Return new ENU matrix of vehicle points
