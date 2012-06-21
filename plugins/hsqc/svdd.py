import numpy as np
from openopt import QP
#from sklearn import svm
#import pylab
from scipy.spatial.distance import pdist, squareform
import scipy
def svdd(X, K, C):
	n = K.shape[0] 	# num. of samples
	
	# for QP sover
	H = 2 * K					# n x n
	f = -np.diag(K)				# n x 1
	l=np.size(X[0,:])
	Aeq = np.ones(n)			# p x n
	beq = 1.0					# p x 1
	
	lb = np.zeros(n)
	ub = np.zeros(n) + C
	
	p = QP(H, f, Aeq = Aeq, beq = beq, lb = lb, ub = ub)
	r = p.solve('nlp:ralg',iprint = -2)#, r = p.solve('nlp:algencan')
	f_opt, x_opt = r.ff, r.xf
	alpha = np.array(x_opt).reshape(-1)
	
	epsilon = C * 1e-4
	svi = np.arange(l)[alpha > epsilon]
	
	nsv = len(svi)

	sv = X[:, svi]
	alpha_sv = alpha[svi].reshape(-1, 1)
	
	#pylab.scatter(sv[0, :], sv[1, :], s = 50, c = 'b')
	#pylab.scatter(X[0, :], X[1, :], c = 'r')
	#pylab.show()
	
	#print 'Support Vectors : %d (%3.1f%%)' % (nsv, 100.0*float(nsv)/float(l))
	
	svii = np.arange(l)[(alpha > epsilon) & (alpha < (C-epsilon))]

	l_svii = len(svii)
	
	#print '# of points in boundary', l_svii

	if l_svii > 0:
		Kt = K[:, svi]
		
		Ksv = Kt[svi, :]
		Kbound = Kt[svii, :]
		Kbound_diag = np.diag(K)[svii]
		
		b = np.dot(np.dot(alpha_sv.T, Ksv), alpha_sv).reshape(-1)
		R2_all = Kbound_diag.reshape(-1) - 2 * np.dot(alpha_sv.T, Kbound).reshape(-1) + b
		R2 = R2_all.mean()
		
		bias = R2 - b
	else:
		print 'No support vectors on the decision boundary'
		
	svdd_model = {'alpha': alpha_sv, 'sv': sv, 'R2': R2, 'bias': bias}
	return svdd_model
	
def svdd_test(model, z):
	Kz = np.vdot(z, z)
	Kzx = np.vdot(model['alpha'], np.dot(z, model['sv']))
	
	r = Kz - 2 * Kzx - model['bias']
	r1 = model['bias']/(Kz - 2 * Kzx)
	return r, r <= 0.0, r1,model['bias'],(Kz - 2 * Kzx)
	

	#print s1, s2
	#print is_in1, is_in2
	
