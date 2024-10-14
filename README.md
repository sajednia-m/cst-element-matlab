# cst-element-matlab
CST elements analyzing using finite element method in MATLAB

= Analyzing CST elements using MATLAB =
= Developed by Mohammad Sajednia =

Inputs:

1. nodes.txt
Column 1: node number
Column 2: X coordinate of the node
Column 3: Y coordinate of the node

2. elements.txt
Column 1: the number of the element
Column 2: the number of the first node
Column 3: the number of the second node
Column 4: the number of the third node

3. boundaries.txt
Column 1: the number of the constrained node
Column 2: the direction of the constraint (1 for X, 2 for Y)
Column 3: the value of the constraint (0 if the constraint is a support and any value for support settlements)

4. nloads.txt
Column 1: the number of the node
Column 2: the magnitude of the nodal force
Column 3: the angle between the force vector and the positive direction of the X axis

5. tractions.txt
Column 1: the number of the first node
Column 2: the number of the second node
Column 3: the magnitude of the traction between two nodes (on the edge)
Column 4: the angle between the traction vector and the positive direction of the X axis

6. vloads.txt
Column 1: the number of the element that includes the body force
Column 2: the magnitude of the body force
Column 3: the angle between the body force vector and the positive direction of the X axis
