from dataclasses import dataclass
import pyx, bac
import numpy as np

## 坐标系
@dataclass
class quadrature_coordinate(object):
    xlength: float = 10
    ylength: float = 10

    def draw(self, pu, xmirror=False, ymirror=False, settings=[]):
        pu.pyx_arrow(                  [ self.xlength, 0], settings=settings)
        if xmirror: pu.pyx_line([0,0], [-self.xlength, 0], settings=settings)
        pu.pyx_arrow(                  [0,  self.ylength], settings=settings)
        if ymirror: pu.pyx_line([0,0], [0, -self.ylength], settings=settings)


## 绘制任意函数波形
@dataclass
class function_as_graph(object):
    def draw(self, pu, list_of_points, settings=[]):
        last_point = list_of_points[0]
        for point in list_of_points[1:]:
            pu.pyx_line(last_point, point, settings=settings)
            last_point = point

global_size = 1.5
global_label_bias_scale = 1.5

## 电路：电压源
@dataclass
class voltage_source(object):
    radius: float = global_size
    center: tuple = (0,0)
    bool_vertical: bool = False

    def __post_init__(self):
        center = self.center
        radius = self.radius
        self.anchors = [
                        (center[0]+radius, center[1]),
                        (center[0]-radius, center[1]),
                        (center[0], center[1]+radius),
                        (center[0], center[1]-radius)
                        ]

    def draw(self, pu, label=None):

        # 圆
        pu.pyx_circle(radius=self.radius, center=self.center)

        # 加号
        sign_bias     = self.radius/2
        stroke_length = self.radius/4
        if self.bool_vertical:
            center_bias1 = [self.center[0], self.center[1]+sign_bias-stroke_length]
            center_bias2 = [self.center[0], self.center[1]+sign_bias+stroke_length]
            center_bias3 = [self.center[0]-stroke_length, self.center[1]+sign_bias]
            center_bias4 = [self.center[0]+stroke_length, self.center[1]+sign_bias]
        else:
            center_bias1 = [self.center[0], self.center[1]+sign_bias-stroke_length]
            center_bias2 = [self.center[0], self.center[1]+sign_bias+stroke_length]
            center_bias3 = [self.center[0]-stroke_length, self.center[1]+sign_bias]
            center_bias4 = [self.center[0]+stroke_length, self.center[1]+sign_bias]
        pu.pyx_line(center_bias1, center_bias2)
        pu.pyx_line(center_bias3, center_bias4)

        # 减号
        sign_bias = - sign_bias
        if self.bool_vertical:
            center_bias3 = [self.center[0]-stroke_length, self.center[1]+sign_bias]
            center_bias4 = [self.center[0]+stroke_length, self.center[1]+sign_bias]
        else:
            center_bias3 = [self.center[0]-stroke_length, self.center[1]+sign_bias]
            center_bias4 = [self.center[0]+stroke_length, self.center[1]+sign_bias]
        pu.pyx_line(center_bias3, center_bias4)

        if label is not None:
            self.add_label(pu, label)

    def add_label(self, pu, label):
        if self.bool_vertical:
            pu.pyx_text([0.5*(self.anchors[0][0]+self.anchors[1][0])-global_label_bias_scale*self.radius, 
                         0.5*(self.anchors[0][1]+self.anchors[1][1])+global_label_bias_scale*self.radius], 
                         label, settings=['night-blue'])
        else:
            pu.pyx_text([0.5*(self.anchors[0][0]+self.anchors[1][0])-global_label_bias_scale*self.radius, 
                         0.5*(self.anchors[0][1]+self.anchors[1][1])+global_label_bias_scale*self.radius], 
                         label, settings=['night-blue'])


## 电路：被动元件父类
@dataclass
class passive_components(object):
    size: float = global_size*2
    location: tuple = (0,0)
    bool_vertical: bool = False
    angle: float = 0.0

    def __post_init__(self):
        location = self.location
        size = self.size
        if self.bool_vertical:
            self.anchors = [
                        (location[0], location[1]),
                        (location[0], location[1]+size)
                            ]
        else:
            self.anchors = [
                        (location[0], location[1]),
                        (location[0], location[1]+size)
                            ]

    @abc.abstractmethod
    def draw(self, pu):
        pass

    def add_label(self, pu, label):
        if not self.bool_vertical:
            pu.pyx_text([0.5*(self.anchors[0][0]+self.anchors[1][0]), 
                         0.5*(self.anchors[0][1]+self.anchors[1][1])+global_label_bias_scale*self.size], 
                         label, settings=['night-blue'])
        else:
            pu.pyx_text([0.5*(self.anchors[0][0]+self.anchors[1][0]-global_label_bias_scale*self.size), 
                         0.5*(self.anchors[0][1]+self.anchors[1][1])], 
                         label, settings=['night-blue'])

## 电路：电感
@dataclass
class inductance(passive_components):
    def draw(self, pu, label=None):
        graph = function_as_graph()
        def f1(t):
            return  np.abs(np.sin(2*np.pi* ( (0.5)*t) ))
        if self.bool_vertical:
            graph.draw( pu, [ (-f1(t)+self.location[0], 
                                    t+self.location[1]) \
                              for t in np.arange(0, self.size+.05, 0.05)] )
        else:
            graph.draw( pu, [ (t +self.location[0],
                             f1(t)+self.location[1]) \
                             for t in np.arange(0, self.size+.05, 0.05)] )

        if label is not None:
            self.add_label(pu, label)

## 电路：电感
@dataclass
class resistance(passive_components):
    def draw(self, pu, label=None):
        if self.bool_vertical:
            pu.pyx_line([self.location[0],             self.location[1]],
                        [self.location[0]-self.size/4, self.location[1]+self.size*1/8],)
            pu.pyx_line([self.location[0]-self.size/4, self.location[1]+self.size*1/8],
                        [self.location[0],             self.location[1]+self.size*2/8],)

            pu.pyx_line([self.location[0],             self.location[1]+self.size*2/8],
                        [self.location[0]-self.size/4, self.location[1]+self.size*3/8])
            pu.pyx_line([self.location[0]-self.size/4, self.location[1]+self.size*3/8],
                        [self.location[0],             self.location[1]+self.size*4/8],)

            pu.pyx_line([self.location[0],             self.location[1]+self.size*4/8],
                        [self.location[0]-self.size/4, self.location[1]+self.size*5/8])
            pu.pyx_line([self.location[0]-self.size/4, self.location[1]+self.size*5/8],
                        [self.location[0],             self.location[1]+self.size*6/8],)

            pu.pyx_line([self.location[0],             self.location[1]+self.size*6/8],
                        [self.location[0]-self.size/4, self.location[1]+self.size*7/8])
            pu.pyx_line([self.location[0]-self.size/4, self.location[1]+self.size*7/8],
                        [self.location[0],             self.location[1]+self.size*8/8],)

        if label is not None:
            self.add_label(pu, label)

## 电路：简单连接
def distance(p1, p2):
    return np.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)
def easy_connect(pu, obj1, obj2):
    dict_of_distances = {(p1,p2): distance(p1,p2) for p1 in obj1.anchors for p2 in obj2.anchors}

    # debug
    # for k,v in dict_of_distances.items():
    #     print(k, ':', v)

    # 找到距离最小的一对点
    p1, p2 = min(dict_of_distances, key=dict_of_distances.get)
    pu.pyx_line(p1, p2)

## 电路：连接点
@dataclass
class connection_point(object):
    location: tuple = (0,0)
    size: float = 0.15

    def __post_init__(self):
        location = self.location
        self.anchors = [
                    (location[0], location[1]),
                            ]
    def draw(self, pu):
        pu.pyx_marker(self.location, size=self.size)
