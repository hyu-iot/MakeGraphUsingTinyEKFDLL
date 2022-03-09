# MakeGraphUsingTinyEKFDLL

This uses DLL with python ctypes module to draw a 3D graph.
DLL is compiled by MS Visual Studio and the cpp file before the compilation includes the TinyEKF codes
which the License.md is for.
You can check the original code at https://github.com/simondlevy/TinyEKF 

I used the TinyEKF to use Extended Kalman Filter while deriving the position of a imu sensor in 3-axis from 3 types of data(accelerometer-a,angular velocity-w, geomagnetic field-h) the sensor is measuring.
Python file(graph_maker.py) includes how to use ctypes module in order to use the function in DLL file.
graph_maker.py first makes 3 arrays(a, w, h), after reading txt files, which will be the parameters for the function I call later.
After the function call it gets an array s which is the position of the sensor, and draws a 3D graph expressing the position.
