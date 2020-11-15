import numpy as np
import ipdb
import sys




num_elem = int(sys.argv[1]) 
left = float(sys.argv[2])
right = float(sys.argv[3])
dx = (right-left)/num_elem
temp = np.random.rand(num_elem+1)

#a = 0.5
#temp = np.random.power(a, num_elem + 1)

temp[0] = 0
x = np.cumsum(temp)
x = left + (right-left) * x / x[-1]
h = x[1:] - x[:-1]
print("Minimum volume fraction: " , np.min(h/dx) )

merge_type = 0
TOL = dx

print("\n###### MERGING META DATA ######\n")
print("min volume fraction %.16e" % np.min(h/dx) )
if merge_type == 0:
    print("merging type = left and right merging (not periodic)")
elif merge_type == 1:
    print("merging type = left and right merging (periodic)")
elif merge_type == 2:
    print("merging type = left merging (periodic)")
elif merge_type == 3:
    print("merging type = right merging (periodic)")
else:
    print("merging type not recognized...\n")
    quit()
print("merging tolerance %lf\n" % (TOL/dx) )
print("###### ----------------- #######\n\n")











m = np.zeros(num_elem)
M = np.zeros(num_elem)
overlaps = np.ones(num_elem)
for elem in range(num_elem):
    left_length = 0
    right_length = 0

    m[elem] = elem
    M[elem] = elem

    if h[elem] > TOL:
        continue

#    ipdb.set_trace(context=21)

    if merge_type == 0: # merge to left and right (nonperiodic)
        curr = elem-1
        while left_length < TOL and curr > 0:
            left_length = left_length + h[curr]
            m[elem] = curr
            overlaps[curr] = overlaps[curr] + 1
            curr = curr - 1
        
        curr = elem+1
        while right_length < TOL and curr < num_elem:
            right_length = right_length + h[curr]
            M[elem] = curr
            overlaps[curr] = overlaps[curr] + 1
            curr = curr + 1
    elif merge_type == 1: # merge to left and right (periodic)
        curr = (elem-1)%num_elem
        while left_length < TOL:
            left_length = left_length + h[curr]
            m[elem] = curr
            overlaps[curr] = overlaps[curr] + 1
            curr = (curr - 1)%num_elem
        
        curr = (elem+1)%num_elem
        while right_length < TOL:
            right_length = right_length + h[curr]
            M[elem] = curr
            overlaps[curr] = overlaps[curr] + 1
            curr = (curr + 1)%num_elem
    elif merge_type == 2: # merge to left (periodic)
        curr = (elem-1)%num_elem
        while left_length < TOL:
            left_length = left_length + h[curr]
            m[elem] = curr
            overlaps[curr] = overlaps[curr] + 1
            curr = (curr - 1)%num_elem

    elif merge_type == 3: # merge to right (periodic)
        curr = (elem+1)%num_elem
        while right_length < TOL:
            right_length = right_length + h[curr]
            M[elem] = curr
            overlaps[curr] = overlaps[curr] + 1
            curr = (curr + 1)%num_elem
       





f=open("grid_"+str(num_elem)+".mdat", 'w')
f.write(str(dx)+"\n")
f.write(str(merge_type)+"\n")
f.write(str(TOL)+"\n")
f.close()
np.savetxt("grid_"+str(num_elem)+".dat", x)

pdata = np.hstack( (m.reshape( (-1,1) ), M.reshape( (-1,1) ), overlaps.reshape( (-1,1) ) )  ).astype(int)
np.savetxt("grid_"+str(num_elem)+".pdat", pdata, fmt="%i")


