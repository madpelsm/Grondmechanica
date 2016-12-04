<h1> Grondmechanica: Zettingsberekening</h1>
<h2>Compile</h2>
<p>Using cmake, adjust cmake to the proper include directories for nanogui. </br>
Clone nanogui to wherever you want (recursive), build it (shared). place the lib file in /dep.</p>
</br>On windows: place the nanogui dll in the same folder as the executable</br>
<h2>Usage</h2>
<h3>Load Type</h3>
<p>First create a load type</br>
<h5>Plate load</h5>
<p>This has now completely replaced the old uniform strip load</p>
Start and end point represent respectively the width and height of the load. It will spawn a rectangle with diagonal (0,0)->(beginPoint,EndPoint).</br>
Define the position using X and Y for the calculation point (X,Y). This represents the relative position to (0,0) of the rectangle. </br>
The point is within the load surface as long as (0<=X<=width) and (0<=Y<=height).
</br>
<h5>Uniform strip load</h5>
<b>Deprecated</b>
You just have to define the width of the load.(via starting point and endpoint -> end-start = width). 
any calculations will be done using the X position relative to the begin and end point of the stripload, since any change in Y does not affect the results if X is constant. Since the assumption of infinitely long uniform strip load was made. 
</br> 
<h6>Then define a point relative to the position of the load, in which you want to calculate the consolidation.</h6>
<h3>Solid types</h3>
Positions of the solid layers are relative, but <b>the order in which entries are made are important.</b> Layers do not overlap automaticly. E.g.
if the entry says upper limit Z1, lower limit Z2, then the thickness will be calculated and that is what will be used in the calculations. Also note that the position of the water table has to be correct relative to the ground layer.
There is no autosort. 
</br>
<b>It is expected that one starts with the top layer and works his way down to the bottom layer</b>
</br>
<h3>Calculation</h3>
After all entries are made, choose the gridsize (in meters) and press calculate. You can also save your configuration for later usage.</br>
The savefile contains the load type, position of the calculation point and the configuration of the solid layers. 
Since the savefile contains the load type and the load type is proper to the consolidation calculation, it is possible to open (load) several save files, press calculate once and get the results for all the layers you loaded.
</br>
<h3>Export</h3>
The export button will generate a report of the calculations, summarize the total consolidation, layer configuration, load type and position
all in one text file. It contains the results of all the calculated consolidation configurations.
</br>
<h3>Review</h3>
After the calculations are done, you can change the number of 'Zettingspunt' under 'Sonderingsinformatie'. If that number represents the number of a consolidation calculation, 
it will be used to retain certain information of the calculation. Such as to get the stresses on a certain depth, the time dependcy of the consolidation, the graph rendering etc.
</br>
<h3>Load files</h3>
The save files are stored in json format, so you can actually edit certain properties and parameters very easily.
