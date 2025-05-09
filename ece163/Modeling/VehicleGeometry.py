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


class VehicleGeometry:
    def __init__(self):
        """
        defines the vehicle in NED coordinates around the local body frame origin. Rotations and translations will be
        around the [0,0,0] point of this local frame. Has to be in NED or the rotation matrices will not work. Vehicle is
        self.SCaled to match the wing span of the actual vehicle, this all the points are in meters.

        "vertices" is an [n x 3] matrix of xyz points for each vertex of the vehicle;
        "faces" is an [m x3] index matrix of which vertices connect to which face (only triangles allowed for faces);
        "colors" is an [m x 4] matrix where each row is a CMYK definition of that face color
        """

        red = [1.0, 0.0, 0.0, 1]
        green = [0.0, 1.0, 0.0, 1]
        blue = [0.0, 0.0, 1.0, 1]
        yellow = [1.0, 1.0, 0.0, 1]
        white = [1.0, 1.0, 1.0, 1]
        black = [0.0, 0.0, 0.0, 0.0]

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
        self.SC = 0.001  # self.SCalingUnit # VPC.b / 6.0	# self.SCaling determined by the wingspan of the aircraft in VehiclePhysicalConstants

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

        wing_length = 700  # mm

        wing_depth = 100  # mm

        # Editable flap coordinates
        self.aelronFlapRadius = 60 * self.SC  # CONST: for calculating the angle based coordinates
        self.aeleronFlapRAngle = -15 # adjusts based on controller: for calculating the angle based coordinates
        self.aeleronFlapLAngle = 15  # adjusts based on controller: for calculating the angle based coordinates

        # Where R = radius of aeleron; math: x = x_0 - (R*cos(theta) - R), y = y, z = R*sin(theta)
        self.aeleron_R_right = [
            -160 * self.SC - (self.aelronFlapRadius * math.cos(self.aeleronFlapRAngle * math.pi / 180) - self.aelronFlapRadius),
            600 * self.SC,
            self.aelronFlapRadius * math.sin(self.aeleronFlapRAngle * math.pi / 180)
        ]  # back right of aeleronn flap R | default: [-160*self.SC, 600*self.SC, 0.01]
        self.aeleron_R_left = [
            -160 * self.SC - (self.aelronFlapRadius * math.cos(self.aeleronFlapRAngle * math.pi / 180) - self.aelronFlapRadius),
            200 * self.SC,
            self.aelronFlapRadius * math.sin(self.aeleronFlapRAngle * math.pi / 180)
        ]  # back left of aeleronn flap R | default: [-160*self.SC, 200*self.SC, 0.01]
        self.aeleron_L_right = [
            -160 * self.SC - (self.aelronFlapRadius * math.cos(self.aeleronFlapRAngle * math.pi / 180) - self.aelronFlapRadius),
            -200 * self.SC,
            self.aelronFlapRadius * math.sin(self.aeleronFlapLAngle * math.pi / 180)
        ]  # back right of aeleron flap L | default: [-160*self.SC, -200*self.SC, 0.01]
        self.aeleron_L_left = [
            -160 * self.SC - (self.aelronFlapRadius * math.cos(self.aeleronFlapRAngle * math.pi / 180) - self.aelronFlapRadius),
            -600 * self.SC,
            self.aelronFlapRadius * math.sin(self.aeleronFlapLAngle * math.pi / 180)
        ]  # back left of aeleron flap L | default: [-160*self.SC, -600*self.SC, 0.01]
        
        # Propellars
        self.PropellarRadius = 113 * self.SC  # CONST: for calculating the angle based coordinates
        self.propellarIncrementer = 20 # multiplier that scales propellar speed
        self.propellarAngle = 30  # adjusts based on controller: for calculating the angle based coordinates
        self.propellar_R_top = [
            161 * self.SC,
            self.PropellarRadius * math.cos(self.propellarAngle * math.pi / 180 + 0.07),
            30 * self.SC + self.PropellarRadius * math.sin(self.propellarAngle * math.pi / 180 + 0.07)
        ]
        self.propellar_R_bottom = [
            161 * self.SC,
            self.PropellarRadius * math.cos(self.propellarAngle * math.pi / 180 - 0.07),
            30 * self.SC + self.PropellarRadius * math.sin(self.propellarAngle * math.pi / 180 - 0.07)
        ]
        self.propellar_L_top = [
            161 * self.SC,
            self.PropellarRadius * math.cos(self.propellarAngle * math.pi / 180 + 3.07),
            30 * self.SC + self.PropellarRadius * math.sin(self.propellarAngle * math.pi / 180 + 3.07)
        ]
        self.propellar_L_bottom = [
            161 * self.SC,
            self.PropellarRadius * math.cos(self.propellarAngle * math.pi / 180 + 3.21),
            30 * self.SC + self.PropellarRadius * math.sin(self.propellarAngle * math.pi / 180 + 3.21)
        ]
        
		# Elevator
        self.ElevatorFlapRadius = 38 * self.SC  # CONST: for calculating the angle based coordinates
        self.elevatorFlapAngle = 30  # adjusts based on controller: for calculating the angle based coordinates
        self.elevator_L_rearflap = [
            -616 * self.SC - (self.ElevatorFlapRadius * math.cos(self.elevatorFlapAngle * math.pi / 180) - self.ElevatorFlapRadius),
            -175*self.SC,
            37*self.SC + self.ElevatorFlapRadius * math.sin(self.elevatorFlapAngle * math.pi / 180)
        ]
        self.elevator_R_rearflap = [
            -616 * self.SC - (self.ElevatorFlapRadius * math.cos(self.elevatorFlapAngle * math.pi / 180) - self.ElevatorFlapRadius),
            175*self.SC,
            37*self.SC + self.ElevatorFlapRadius * math.sin(self.elevatorFlapAngle * math.pi / 180)
        ]
        
		# Rudder
        self.RudderFlapRadius = 38 * self.SC  # CONST: for calculating the angle based coordinates
        self.rudderFlapAngle = 30  # adjusts based on controller: for calculating the angle based coordinates
        self.rudder_top_rearflap = [
            -616 * self.SC - (self.RudderFlapRadius * math.cos(self.rudderFlapAngle * math.pi / 180) - self.RudderFlapRadius),
            self.RudderFlapRadius * math.sin(self.rudderFlapAngle * math.pi / 180),
            -122*self.SC
        ]
        self.rudder_bottom_rearflap = [
            -616 * self.SC - (self.RudderFlapRadius * math.cos(self.rudderFlapAngle * math.pi / 180) - self.RudderFlapRadius),
            self.RudderFlapRadius * math.sin(self.rudderFlapAngle * math.pi / 180),
            0
        ]
		# self.rudder_top_rearflap = [-616*self.SC, 0, -122*self.SC]
        # self.rudder_bottom_rearflap = [-616*self.SC, 0, 0]
        
        
        self.vertices = [
            # Main Wing Body
            [0, -700 * self.SC, -self.SC],  # [0] Front Left Corner of Wing
            [0, 700 * self.SC, -self.SC],  # [1] Front Right Corner of Wing
            [-100 * self.SC, 700 * self.SC, -self.SC],  # [2] Back Right Corner of Wing / front right of solid Wing R
            [-100 * self.SC,-700 * self.SC, -self.SC],  # [3] Back Left Corner of Wing / front left of solid Wing L
            # Solid Wing R of Wing
            [-160 * self.SC, 700 * self.SC, -self.SC],  # [4] back right corner of Solid Wing R
            [-160 * self.SC, 600 * self.SC, -self.SC],  # [5] back left corner of Solid Wing R
            [-100 * self.SC, 600 * self.SC, -self.SC],  # [6] front left corner of Solid Wing R / front right of aeleronn flap R
            # aeleron Flap R: adjustable
            self.aeleron_R_right,  # [7] back right of aeleronn flap R
            self.aeleron_R_left,  # [8] back left of aeleronn flap R
            # Solid Wing C of Wing
            [-100 * self.SC,200 * self.SC, -self.SC],  # [9] front right of Solid Wing C / front left of aeleronn flap R
            [-160 * self.SC, 200 * self.SC, -self.SC],  # [10] back right of Solid Wing C
            [-160 * self.SC, -200 * self.SC, -self.SC],  # [11] back left of Solid Wing C
            [-100 * self.SC,-200 * self.SC, -self.SC],  # [12] front left of Solid Wing C / front right of aeleron flap L
            # aeleron Flap L: adjustable
            self.aeleron_L_right,  # [13] back right of aeleron flap L
            self.aeleron_L_left,  # [14] back left of aeleron flap L
            # Solid Wing L (of wing)
            [-100 * self.SC,-600 * self.SC, -self.SC],  # [15] front right of Solid Wing L / front left of aeleron flap L
            [-160 * self.SC, -600 * self.SC, -self.SC],  # [16] back right of Solid Wing L
            [-160 * self.SC, -700 * self.SC, -self.SC],  # [17] back left of Solid Wing L
            # Body Top
            [160 * self.SC, 41 * self.SC, 0],  # [18]
            [160 * self.SC, -41 * self.SC, 0],  # [19]
            [-578 * self.SC, -41 * self.SC, 0],  # [20]
            [-578 * self.SC, -41 * self.SC, 0],  # [21]
            [-578 * self.SC, 41 * self.SC, 0],  # [22]
            [160 * self.SC, 41 * self.SC, 0],  # [23]
            # Body Bottom
            [95 * self.SC, 41 * self.SC, 102 * self.SC],  # [24]
            [95 * self.SC, -41 * self.SC, 102 * self.SC],  # [25]
            [-220 * self.SC, -41 * self.SC, 102 * self.SC],  # [26]
            [-220 * self.SC, -41 * self.SC, 102 * self.SC],  # [27]
            [-220 * self.SC, 41 * self.SC, 102 * self.SC],  # [28]
            [95 * self.SC, 41 * self.SC, 102 * self.SC],  # [29]
            [-220 * self.SC, -41 * self.SC, 102 * self.SC],  # [30]
            [-220 * self.SC, 41 * self.SC, 102 * self.SC],  # [31]
            [-578 * self.SC, -41 * self.SC, 37 * self.SC],  # [32]
            [-578 * self.SC, -41 * self.SC, 37 * self.SC],  # [33]
            [-578 * self.SC, 41 * self.SC, 37 * self.SC],  # [34]
            [-220 * self.SC, 41 * self.SC, 102 * self.SC],  # [35]
            [-578 * self.SC, -41 * self.SC, 0],  # [36]
            [-578 * self.SC, -41 * self.SC, 37 * self.SC],  # [37]
            [-578 * self.SC, 41 * self.SC, 37 * self.SC],  # [38]
            [-578 * self.SC, 41 * self.SC, 37 * self.SC],  # [39]
            [-578 * self.SC, 41 * self.SC, 0],  # [40]
            [-578 * self.SC, -41 * self.SC, 0],  # [41]
            # Body Left Side
            [160 * self.SC, -41 * self.SC, 37 * self.SC],  # [42]
            [95 * self.SC, -41 * self.SC, 102 * self.SC],  # [43]
            [95 * self.SC, -41 * self.SC, 37 * self.SC],  # [44]
            [95 * self.SC, -41 * self.SC, 102 * self.SC],  # [45]
            [95 * self.SC, -41 * self.SC, 37 * self.SC],  # [46]
            [-220 * self.SC, -41 * self.SC, 37 * self.SC],  # [47]
            [-220 * self.SC, -41 * self.SC, 37 * self.SC],  # [48]
            [-220 * self.SC, -41 * self.SC, 102 * self.SC],  # [49]
            [95 * self.SC, -41 * self.SC, 102 * self.SC],  # [50]
            [-220 * self.SC, -41 * self.SC, 102 * self.SC],  # [51]
            [-220 * self.SC, -41 * self.SC, 37 * self.SC],  # [52]
            [-578 * self.SC, -41 * self.SC, 37 * self.SC],  # [53]
            [-578 * self.SC, -41 * self.SC, 0],  # [54]
            [-578 * self.SC, -41 * self.SC, 37 * self.SC],  # [55]
            [160 * self.SC, -41 * self.SC, 0],  # [56]
            [160 * self.SC, -41 * self.SC, 37 * self.SC],  # [57]
            [160 * self.SC, -41 * self.SC, 0],  # [58]
            [-578 * self.SC, -41 * self.SC, 37 * self.SC],  # [59]
            # Body Right Side same as left but with opposite Y
            [160 * self.SC, 41 * self.SC, 37 * self.SC],  # [60]
            [95 * self.SC, 41 * self.SC, 102 * self.SC],  # [61]
            [95 * self.SC, 41 * self.SC, 37 * self.SC],  # [62]
            [95 * self.SC, 41 * self.SC, 102 * self.SC],  # [63]
            [95 * self.SC, 41 * self.SC, 37 * self.SC],  # [64]
            [-220 * self.SC, 41 * self.SC, 37 * self.SC],  # [65]
            [-220 * self.SC, 41 * self.SC, 37 * self.SC],  # [66]
            [-220 * self.SC, 41 * self.SC, 102 * self.SC],  # [67]
            [95 * self.SC, 41 * self.SC, 102 * self.SC],  # [68]
            [-220 * self.SC, 41 * self.SC, 102 * self.SC],  # [69]
            [-220 * self.SC, 41 * self.SC, 37 * self.SC],  # [70]
            [-578 * self.SC, 41 * self.SC, 37 * self.SC],  # [71]
            [-578 * self.SC, 41 * self.SC, 0],  # [72]
            [-578 * self.SC, 41 * self.SC, 37 * self.SC],  # [73]
            [160 * self.SC, 41 * self.SC, 0],  # [74]
            [160 * self.SC, 41 * self.SC, 37 * self.SC],  # [75]
            [160 * self.SC, 41 * self.SC, 0],  # [76]
            [-578 * self.SC, 41 * self.SC, 37 * self.SC],  # [77]
            # Body Front
            [160 * self.SC, 41 * self.SC, 0],  # [78]
            [160 * self.SC, 41 * self.SC, 37 * self.SC],  # [79]
            [160 * self.SC, -41 * self.SC, 37 * self.SC],  # [80]
            [160 * self.SC, -41 * self.SC, 37 * self.SC],  # [81]
            [160 * self.SC, -41 * self.SC, 0],  # [82]
            [160 * self.SC, 41 * self.SC, 0],  # [83]
            [160 * self.SC, 41 * self.SC, 37 * self.SC],  # [84]
            [160 * self.SC, -41 * self.SC, 37 * self.SC],  # [85]
            [95 * self.SC, -41 * self.SC, 102 * self.SC],  # [86]
            [95 * self.SC, -41 * self.SC, 102 * self.SC],  # [87]
            [95 * self.SC, 41 * self.SC, 102 * self.SC],  # [88]
            [161 * self.SC, 41 * self.SC, 37 * self.SC],  # [89]
            # Propellor
            self.propellar_R_top,  # [161*self.SC, 113*self.SC, 22*self.SC], # [90] propellar right top
            self.propellar_R_bottom,  # [161*self.SC, 113*self.SC, 38*self.SC], # [91] propellar right bottom
            [161 * self.SC, 0, 30 * self.SC],  # [92]
            self.propellar_L_top,  # [161*self.SC, -113*self.SC, 22*self.SC], # [93] propellar left top
            self.propellar_L_bottom,  # [161*self.SC, -113*self.SC, 38*self.SC], # [94] propelar left bottom
            [161 * self.SC, 0, 30 * self.SC],  # [95]
            
            # Elevator
            self.elevator_L_rearflap, # [96] elevator Left Rearflap
			self.elevator_R_rearflap, #[97] elevator Right Rearflap
            [-578 * self.SC, 175 * self.SC, 37 * self.SC], # [98] elevator Left Rear
            [-578 * self.SC, -175 * self.SC, 37 * self.SC], # [99] elevator Right Rear
            [-498 * self.SC, 175 * self.SC, 37 * self.SC], # [100] elevator Left Front
            [-498 * self.SC, -175 * self.SC, 37 * self.SC], # [101] elevator Right Front
            

			# Rudder
			self.rudder_top_rearflap, # [102] rudder top rearflap
			self.rudder_bottom_rearflap, # [103] rudder bottom rearflap
            [-578* self.SC, 0, 37* self.SC], # [104] rudder bottom rear flap
            [-578* self.SC, 0, -122* self.SC], # [105] rudder top rear hinge
            [-578* self.SC, 0, 0], # [106] rudder bottom rear hinge
            [-522* self.SC, 0, 0], # [107] rudder bottom middle
            [-522* self.SC, 0, -122* self.SC], # [108] rudder top middle
            [-443* self.SC, 0, 0], # [109] rudder bottom front     
        ]
        

        self.faces = [
            # Main Wing
            [0, 1, 2],  # Connect Front Left, Front Right, and Back Right Wing vertices
            [0, 3, 2],  # Connect Front Left, Back Left, and Back Right Wing vertices
            # Solid Wing R - green
            [6,2,4],  # Connect Front Left, Front Right, and Back Right Solid Wing R vertices
            [6,5,4],  # Connect Front Left, Back Left, and Back Right Solid Wing R vertices
            # Flap R
            [9,6,7],  # Connect Front Left, Front Right, and Back Right Flap R vertices
            [9, 8, 7],  # Connect Front Left, Back Left, and Back Right Flap vertices
            # Solid Wing C
            [12,9,10],  # Connect Front Left, Front Right, and Back Right Solid Wing C vertices
            [12,11,10],  # Connect Front Left, Back Left, and Back Right Solid Wing C vertices
            # Flap L
            [15,12,13],  # Connect Front Left, Front Right, and Back Right Flap L vertices
            [15,14,13],  # Connect Front Left, Back Left, and Back Right Flap L vertices
            # Solid Wing L
            [3,15,16],  # Connect Front Left, Front Right, and Back Right Solid Wing L vertices
            [3,17,16],  # Connect Front Left, Back Left, and Back Right Solid Wing L vertices
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
            [51, 52, 53],
            [54, 55, 56],
            [57, 58, 59],
            # Body Right Side
            [60, 61, 62],
            [63, 64, 65],
            [66, 67, 68],
            [69, 70, 71],
            [72, 73, 74],
            [75, 76, 77],
            # Body Front
            [78, 79, 80],
            [81, 82, 83],
            [84, 85, 86],
            [87, 88, 89],
            # Propeller
            [90, 91, 92],
            [93, 94, 95],
            # Elevator
			[96, 97, 98],
			[96, 99, 98],
			[98, 99, 100],
			[99, 100, 101],
            # Rudder
            [102, 103, 104],
            [102, 105, 104],
            [106, 107, 108],
            [106, 105, 108],
            [108, 109, 107],
        ]

        self.colors = [
            green,
            green,
            green,
            green,
            red,
            red,
            green,
            green,
            red,
            red,
            green,
            green,
            blue, # boddy top
            blue, # boddy top
            blue,
            blue,
            blue,
            blue,
            blue,
            blue,
            blue,
            blue,
            blue,
            blue,
            blue,
            blue,
            blue,
            blue,
            blue,
            blue,
            blue,
            blue,
            yellow,
            yellow,
            yellow,
            yellow,
            red,
            red,
            red,
            red,
            green,
            green,
            yellow,
            yellow,
            green,
            green,
            green,
        ]

        """

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
"""
        return

    def updateAngles(self, Throttle=0.0, Aileron=0.0, Elevator=0.0, Rudder=0.0):
		### Ailerons
        self.aeleronFlapRAngle = (-Aileron * 2)  # adjusts based on controller: for calculating the angle based coordinates
        self.aeleronFlapLAngle = (Aileron * 2)  # adjusts based on controller: for calculating the angle based coordinates
        # print(self.aeleron_R_right)
        ## Where R = radius of aeleron; math: x = x_0 - (R*cos(theta) - R), y = y, z = R*sin(theta)
        self.aeleron_R_right = [
            -160 * self.SC - (self.aelronFlapRadius * math.cos(self.aeleronFlapRAngle * math.pi / 180) - self.aelronFlapRadius),
            600 * self.SC,
            self.aelronFlapRadius * math.sin(self.aeleronFlapRAngle * math.pi / 180),
        ]  # back right of aeleronn flap R | default: [-160*self.SC, 600*self.SC, 0.01]
        self.aeleron_R_left = [
            -160 * self.SC - (self.aelronFlapRadius * math.cos(self.aeleronFlapRAngle * math.pi / 180) - self.aelronFlapRadius),
            200 * self.SC,
            self.aelronFlapRadius * math.sin(self.aeleronFlapRAngle * math.pi / 180),
        ]  # back left of aeleronn flap R | default: [-160*self.SC, 200*self.SC, 0.01]
        self.aeleron_L_right = [
            -160 * self.SC - (self.aelronFlapRadius * math.cos(self.aeleronFlapRAngle * math.pi / 180) - self.aelronFlapRadius),
            -200 * self.SC,
            self.aelronFlapRadius * math.sin(self.aeleronFlapLAngle * math.pi / 180)
        ]  # back right of aeleron flap L | default: [-160*self.SC, -200*self.SC, 0.01]
        self.aeleron_L_left = [
            -160 * self.SC - (self.aelronFlapRadius * math.cos(self.aeleronFlapRAngle * math.pi / 180) - self.aelronFlapRadius),
            -600 * self.SC,
            self.aelronFlapRadius * math.sin(self.aeleronFlapLAngle * math.pi / 180)
        ]  # back left of aeleron flap L | default: [-160*self.SC, -600*self.SC, 0.01]
        # print(self.aeleron_R_right)


        ### Propellars
        # print(self.propellarAngle)
        self.PropellarRadius = (
            113 * self.SC
        )  # CONST: for calculating the angle based coordinates
        self.propellarAngle = (
            self.propellarAngle + Throttle * self.propellarIncrementer
        )  # adjusts based on controller: for calculating the angle based coordinates
        if self.propellarAngle > 360:
            self.propellarAngle = self.propellarAngle - 360  # wrap propellar angle
        self.propellar_R_top = [
            161 * self.SC,
            self.PropellarRadius * math.cos(self.propellarAngle * math.pi / 180 + 0.07),
            30 * self.SC + self.PropellarRadius* math.sin(self.propellarAngle * math.pi / 180 + 0.07)
        ]
        self.propellar_R_bottom = [
            161 * self.SC,
            self.PropellarRadius * math.cos(self.propellarAngle * math.pi / 180 - 0.07), 
            30 * self.SC + self.PropellarRadius* math.sin(self.propellarAngle * math.pi / 180 - 0.07)
        ]
        self.propellar_L_top = [
            161 * self.SC,
            self.PropellarRadius * math.cos(self.propellarAngle * math.pi / 180 + 3.07),
            30 * self.SC + self.PropellarRadius* math.sin(self.propellarAngle * math.pi / 180 + 3.07)
        ]
        self.propellar_L_bottom = [
            161 * self.SC,
            self.PropellarRadius * math.cos(self.propellarAngle * math.pi / 180 + 3.21),
            30 * self.SC + self.PropellarRadius * math.sin(self.propellarAngle * math.pi / 180 + 3.21)
        ]
        # Elevator
        self.elevatorFlapAngle = Elevator * 2  # adjusts based on controller: for calculating the angle based coordinates
        self.elevator_L_rearflap = [
            -616 * self.SC - (self.ElevatorFlapRadius * math.cos(self.elevatorFlapAngle * math.pi / 180) - self.ElevatorFlapRadius),
            -175*self.SC,
            37*self.SC + self.ElevatorFlapRadius * math.sin(self.elevatorFlapAngle * math.pi / 180)
        ]
        self.elevator_R_rearflap = [
            -616 * self.SC - (self.ElevatorFlapRadius * math.cos(self.elevatorFlapAngle * math.pi / 180) - self.ElevatorFlapRadius),
            175*self.SC,
            37*self.SC + self.ElevatorFlapRadius * math.sin(self.elevatorFlapAngle * math.pi / 180)
        ]
        # Rudder
        self.rudderFlapAngle = -Rudder * 2  # adjusts based on controller: for calculating the angle based coordinates
        self.rudder_top_rearflap = [
            -616 * self.SC - (self.RudderFlapRadius * math.cos(self.rudderFlapAngle * math.pi / 180) - self.RudderFlapRadius),
            self.RudderFlapRadius * math.sin(self.rudderFlapAngle * math.pi / 180),
            -122*self.SC
        ]
        self.rudder_bottom_rearflap = [
            -616 * self.SC - (self.RudderFlapRadius * math.cos(self.rudderFlapAngle * math.pi / 180) - self.RudderFlapRadius),
            self.RudderFlapRadius * math.sin(self.rudderFlapAngle * math.pi / 180),
            0
        ]
         

    def getNewPoints(self,x,y,z,yaw,pitch,roll,Throttle=0.0,Aileron=0.0,Elevator=0.0,Rudder=0.0):
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
        # aileron
        self.vertices[7] = self.aeleron_R_right
        self.vertices[8] = self.aeleron_R_left
        self.vertices[13] = self.aeleron_L_right
        self.vertices[14] = self.aeleron_L_left
        # propellar:
        self.vertices[90] = self.propellar_R_top
        self.vertices[91] = self.propellar_R_bottom
        self.vertices[93] = self.propellar_L_top
        self.vertices[94] = self.propellar_L_bottom
        # Elevator:
        self.vertices[96] = self.elevator_L_rearflap
        self.vertices[97] = self.elevator_R_rearflap
        # Rudder:
        self.vertices[102] = self.rudder_top_rearflap
        self.vertices[103] = self.rudder_bottom_rearflap

        newPoints = self.vertices

        rmatrix = Rotations.euler2DCM(yaw, pitch, roll)  # Get the rotation matrix from the euler angles using euler2dcm()

        coords_rotated = MatrixMath.multiply(newPoints, rmatrix)  # multiply aircraft points by the rotation matrix to apply the appropriate rotation in NED

        displacements = [x, y, z]  # Given x, y, z NED displacements

        for row in range(len(coords_rotated[0])):  # Go through the rows of the vehicle points matrix
            for col in range((len(coords_rotated))):  # Go through the cols of the vehicle points matrix
                coords_rotated[col][row] += displacements[row]  # Add the x displacement to col 1, y to 2, and z to 3 for each row in nx3
        newPoints = Rotations.ned2enu(coords_rotated)  # Convert dispalced NED matrix to ENU coords

        return newPoints  # Return new ENU matrix of vehicle points
