HOMEWORK 2: CLOTH & FLUID SIMULATION

NAME:  < zhiyangxu >


TOTAL TIME SPENT:  < 60 >
Please estimate the number of hours you spent on this assignment.


COLLABORATORS AND OTHER RESOURCES: List the names of everyone you
talked to about this assignment and all of the resources (books,
online reference material, etc.) you consulted in completing this
assignment.

< insert collaborators / resources >

Remember: Your implementation for this assignment must be done on your
own, as described in "Academic Integrity for Homework" handout.


DISCUSS STABILITY OF YOUR CLOTH SIMULATION:
My cloth simulation is very stable and with Provot correction, it looks pretty good.

DESCRIBE YOUR NEW CLOTH TEST SCENE & 
THE COMMAND LINE TO RUN YOUR EXAMPLE:
I change the k_structural,k_shear and k_bend so the cloth looks softer. 

./simulation -cloth ../src/my_test_case.txt -timestep 0.001


DISCUSS THE ACCURACY & STABILITY OF YOUR FLUID SIMULATION:
My fluid simulation with full cells is very stable and the pressures are all 0 across all cells.

KNOWN BUGS IN YOUR CODE:
Please be concise!
The only part that doesn’t work is the free surface part. Other parts work well.

NEW FEATURES OR EXTENSIONS FOR EXTRA CREDIT:
Include instructions for use and test cases and sample output as appropriate.
