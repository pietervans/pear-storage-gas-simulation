import matplotlib.pyplot as plt
import numpy as np

errors = np.array([4.08087e-01, 5.84866e-02, 9.54651e-03, 8.89867e-04, 1.04003e-05, 1.52169e-09, 1.87006e-14])

y = np.emath.logn(1.6, np.abs(np.log(errors)))
plt.plot(y, 'o-')
plt.plot(range(len(y)), '--')
plt.legend(['Observed', 'exp($-1.6^k$)'])
plt.xlabel('Iteration')
plt.ylabel('log$_a$(|log(backward error)|)')
plt.show()

# for i in range(1,len(errors)):
#     print(np.log(errors[i])/np.log(errors[i-1]))