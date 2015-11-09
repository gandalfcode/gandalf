#==============================================================================
# example16.py
# Retrieve data and save it inside a variable
#==============================================================================
from gandalf.analysis.facade import *

# Load in simulation from disk (from example10.py)
sim = loadsim("KHI1")

# Get the data in variable x
x=get_data('x')
print x

# Do a rendered plot and save the image
image=get_render_data('x','y','rho')

# Use the matplotlib to plot the image
import matplotlib.pyplot as plt
plt.imshow(image,interpolation='nearest')
plt.show()

# You can also the do the following when doing a plot:
data=plot('x','y')

# Let's do a matplotlib plot in a new figure to show how you can access the data
plt.figure()
plt.plot(data.x_data,data.y_data,'.')

# It works with render too
plt.figure()
data=render('x','y','rho')
plt.figure()
plt.imshow(data.render_data,interpolation='nearest')
plt.show()