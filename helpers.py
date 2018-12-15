# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from scipy.optimize import fsolve
from decimal import Decimal
import pandas as pd
from uncertainties import ufloat, unumpy
from uncertainties.umath import *
from subprocess import call, Popen, PIPE
import os
import sys
import re
from decimal import Decimal
import math

'''
Special datatype for floats with uncertainties and significant figures
'''
class physical:
	'''
	@param {physical} self      - reference to the instance
	@param {float|int|ufloat} n - nominal value
	@param {float|int} s    	- standard deviation of the nominal value 
	@param {int} sf  			- significant figures of the nominal value
	@param {int} sfs  			- significant figures of the standard deviation
	@param {string} tag 		- string to tag the number
	@param {int} sfpm 			- significant figures after the comma
	@return {physical}			- physical instance of the given quantities

	Can be called by:

	x = physical(ufloat(n, s), sf, sfs, tag, sfpm)
	x = physical(n, s, sf, sfs, tag, sfpm)
	'''
	#def __init__(self, n, s=0.0, sf=None, tag=None):
	def __init__(self, n, s=None, sf=None, sfs=None, tag=None, sfpm=None):
		if hasattr(n, 'n'):
			self.d = n
			self.sf = np.inf if s == None else s
			self.sfs = np.inf if sf == None else sf
			self.tag = sfs
			self.sfpm = tag if tag != None else self.significant_figures_after_comma(n.n)
		else:
			self.s = 0 if s == None else s
			self.d = ufloat(n, self.s, tag)
			self.sf = np.inf if sf == None else sf
			self.sfs = np.inf if sfs == None else sfs
			self.tag = tag
			self.sfpm = sfpm if sfpm != None else self.significant_figures_after_comma(n)

		self.n = self.d.n
		self.nominal_value = self.d.n
		self.s = self.d.s
		self.std_dev = self.d.s

	def __add__(self, operand):
		if hasattr(operand, "__len__"):
			return (np.vectorize( lambda se, op: se + op)(self, operand) )
		elif hasattr(operand, 'd'):
			d = self.d + operand.d
			sfpm = self.sfpm if self.sfpm <= operand.sfpm else operand.sfpm
			sf = sfpm + len(str(int(d.n)))
			sfs = 0 #TODO
		else:
			d = self.d + operand
			sfpm = self.sfpm
			sf = self.sf
			sfs = 0 #TODO

		return physical(d, sf, sfs, None, sfpm)

	def __sub__(self, operand):
		return self.sub(operand)

	def __mul__(self, operand):
		if hasattr(operand, "__len__"):
			return (np.vectorize( lambda se, op: se * op)(self, operand) )
		elif hasattr(operand, 'd'):
			d = self.d * operand.d
			sf = self.sf if self.sf <= operand.sf else operand.sf
			sfpm = 0 if sf < len(str(int(d.n))) else sf - len(str(int(d.n)))
			sfs = 0 #TODO
		else:
			d = self.d * operand
			sf = self.sf
			sfpm = self.sfpm
			sfs = 0 #TODO

		return physical(d, sf, sfs, None, sfpm)

	def __div__(self, operand):
		return self.div(operand)

	# basic operations from the right
	def __radd__(self, operand):
		return self.__add__(operand)
	def __rsub__(self, operand):
		return self.sub(operand, False)
	def __rmul__(self, operand):
		return self.__mul__(operand)
	def __rdiv__(self, operand):
		return self.div(operand, False)

	def __abs__(self):
		return abs(self.d)

	def __mod__(self, operand):
		return self.d % operand

	def __neg__(self):
		self.d = -self.d
		self.n = -self.n
		self.nominal_value = -self.nominal_value
		return self

	def __pos__(self):
		return self

	def __eq__(self, operand):
		return self.d == operand
	def __lt__(self, operand):
		return self.d < operand
	def __le__(self, operand):
		return self.d <= operand
	def __gt__(self, operand):
		return self.d > operand
	def __ge__(self, operand):
		return self.d >= operand
	def __ne__(self, operand):
		return self.d != operand

	# for Python3
	def __truediv__(self, operand):
		return self.div(operand)

	def __rtruediv__(self, operand):
		return self.div(operand, False)

	def __str__(self):
		return str(self.d)

	def __repr__(self):
		return repr(self.d)

	def __format__(self, a):
		return self.d.__format__(a)

	@classmethod
	def significant_figures_after_comma(self, nr):
		if '.' in str(nr):
			return len(str(Decimal(str(nr)) % 1)) - (2 if nr >= 0 else 3)
		else:
			return 0

	# subtraction and division are no cummutative
	def sub(self, operand, left = True):
		if hasattr(operand, "__len__"):
			return ( np.vectorize( lambda a, b: a - b if left else b - a )(self, operand, left) )
		elif hasattr(operand, 'd'):
			d = self.d - operand.d if left else operand.d - self.d
			sfpm = self.sfpm if self.sfpm <= operand.sfpm else operand.sfpm
			sf = sfpm + len(str(int(d.n)))
			sfs = 0 #TODO
		else:
			d = self.d - operand if left else operand - self.d
			sfpm = self.sfpm
			sf = self.sf
			sfs = 0 #TODO

		return physical(d, sf, sfs, None, sfpm)

	def div(self, operand, left = True):
		if hasattr(operand, "__len__"):
			return ( np.vectorize( lambda a, b, r: a / b if left else b / a )(self, operand, left) )
		elif hasattr(operand, 'd'):
			d = self.d / operand.d if left else operand.d / self.d
			sf = self.sf if self.sf <= operand.sf else operand.sf
			sfpm = 0 if sf < len(str(int(d.n))) else sf - len(str(int(d.n)))
			sfs = 0 #TODO
		else:
			d = self.d / operand if left else operand / self.d
			sf = self.sf
			sfpm = self.sfpm
			sfs = 0 #TODO

		return physical(d, sf, sfs, None, sfpm)

	#tan = lambda self: physical(unumpy.tan(self.d), self.sf, self.sfs, self.tag, self.sfpm)
	#sin = lambda self: physical(unumpy.sin(self.d), self.sf, self.sfs, self.tag, self.sfpm)

class pnumpy:
	sin = lambda x: run_function(x, unumpy.sin)
	cos = lambda x: run_function(x, unumpy.cos)
	tan = lambda x: run_function(x, unumpy.tan)
	arcsin = lambda x: run_function(x, unumpy.arcsin)
	arccos = lambda x: run_function(x, unumpy.arccos)
	arctan = lambda x: run_function(x, unumpy.arctan)

def run_function(x, func):
	if hasattr(x, "__len__"):
		u = func(np.vectorize( lambda a: a.d if hasattr(a, "d") else a )(x))
		return ( np.vectorize(
			lambda a, b: physical(a.n, a.s, b.sf, b.sfs, b.tag, b.sfpm) if hasattr(a, "n") else a,
			otypes=[object])
			(u, x) )
	else:
		u = func(np.array([ x.d if hasattr(x, "d") else x ]))[0]
		return physical(u.n, u.s, x.sf, x.sfs, x.tag, x.sfpm) if hasattr(u, "n") else u


def pharray(nom, errs, sfs, sfss=None, tags=None, sfpms=None):
	return (np.vectorize(
		lambda a, b, c, d, e, f: physical(a, b, c, d, e, f), otypes=[object])
		(nom, errs, sfs, sfss, tags, sfpms))

def sf(nr):
	if re.match(r'0\.[0]+', str(nr)):
		return int(len(str(nr)) - 1)
	elif str(nr) == "0":
		return 1
	else:
		a = str(nr)
		a = re.sub('[eE]+.*$', '', a)
		a = re.sub('[-\.,]', '', a)
		a = re.sub('^[0]+', '', a)
		return int(len(a))

def sf_ac(nr):
	if '.' in str(nr):
		return len(str(Decimal(str(nr)) % 1)) - (2 if float(nr) >= 0 else 3)
	else:
		return 0

'''
Fetch data from a csv/excel file
@param {string} file        - csv/excel file containing data
@param {string} col         - column header name
@param {float|array} err    - error of the data set 
@param {int|string} sheet   - sheet in case of an excel file
@return {numpy.ndarray}		- numpy array containing the data (and their errors if given)
'''
def fetch (file, col, err = 0, sheet = "Sheet1"):
	ext = os.path.splitext(file)[1]
	if ext == ".csv":
		df = pd.read_csv(file)
	elif ext == ".xlsx":
		df = pd.read_excel(open(file,'rb'), sheetname=sheet, converters={col: parse})

	#for a in df[col]:
	#	print(str(a), type(a))
	#for i, a in enumerate(df[col]*10**3):
	#	print(i, a, type(a))
	#print(df[col])
	#exit()

	if err == 0:
		return np.array(df[col])
	else:
		res = unumpy.uarray(np.zeros(df[col].size), np.zeros(df[col].size))
		errarr = np.zeros(df[col].size) + err
		for i, item in enumerate(df[col]):
			sign = len(re.sub('[\.,-]', '', item)) # TODO: calculation not yet correct
			#print(float(item))
			res[i] = ufloat(float(item), errarr[i], tag=str(sign))
			#x = ufloat(float(item), errarr[i], tag=str(sign))
			#print(x.tag)
		return res
		#return unumpy.uarray(df[col], np.zeros(df[col].size) + err)

def fetch2 (file, col, err = 0, sheet = 0):
	ext = os.path.splitext(file)[1]
	if ext == ".csv":
		df = pd.read_csv(file)
	elif ext == ".xlsx":
		df = pd.read_excel(open(file,'rb'), sheetname=sheet, converters={col: parse})

	if not hasattr(err, "__len__") and err == 0:
		return np.array(df[col]).astype(float)
	else:
		
		if isinstance(err, str) and "%" in err:
			rel_err = float(re.sub('[\%]', '', err))/100
			earr = np.zeros(df[col].size) + rel_err*np.abs(np.array(df[col]).astype(float))
		else:
			earr = np.zeros(df[col].size) + err

		sfarr = np.zeros(df[col].size).astype(int) + (np.vectorize(
			lambda a: sf(a), otypes=[int])
			(df[col]))
		sfpmarr = np.zeros(df[col].size).astype(int) + (np.vectorize(
			lambda a: sf_ac(a), otypes=[int])
			(df[col]))

		return pharray(np.array(df[col]).astype(float), earr, sfarr, None, None, sfpmarr)

def parse (item):
	return str(item)

def addError (arr, err):
	for i, el in enumerate(arr):
		s = arr[i].s + err
		arr[i] = physical(arr[i].n, s, arr[i].sf, arr[i].sfs, arr[i].tag, arr[i].sfpm)
	return arr

def nominal (a):
	if hasattr(a, "__len__"):
		ret = []
		for i, e in enumerate(a):
			ret.append(nominal(e))
		ret = np.array(ret)
	else:
		ret = a.n if hasattr(a, 'n') else a
	return ret

def stddev (a):
	if hasattr(a, "__len__"):
		ret = []
		for i, e in enumerate(a):
			ret.append(stddev(e))
		ret = np.array(ret)
	else:
		ret = a.s if hasattr(a, 's') else a
	return ret

# TODO: propagation not correct, there are rules for */+-
def sig_figures (a):
	r = np.inf
	if hasattr(a, 'd'):
		return '.' + str(a.sf) + 'g'
	elif hasattr(a, 'n'):
		for (var, error) in a.error_components().items():
			tag = len(str(a.n)) if (var.tag == None) else var.tag
			r = int(tag) if int(tag) < r else r

	return '.' + str(r) + 'g'

def sig_figures2 (a):
	r = np.inf
	if hasattr(a, 'd'):
		r = a.sf
	elif hasattr(a, 'n'):
		for (var, error) in a.error_components().items():
			tag = len(str(a.n)) if (var.tag == None) else var.tag
			r = int(tag) if int(tag) < r else r
	else:
		return 3

	return r

def sig_figures_err (a):
	r = np.inf
	if hasattr(a, 'd'):
		r = a.sfs
	else:
		r = sig_figures2(a)

	return r

'''
Compare two arrays of numbers with an absolute tolerance
@param {array} A	- array A
@param {array} B	- array B
@param {float} tol  - absolte tolerance between values
@return {bool}		- whether A ~= B or not
'''
def equal (A, B, tol):
	for d in np.abs(A - B) - tol:
		if d > 0:
			return False
	return True


'''
Format a number (with its error if given) in latex
@param {int|float|ufloat} nr	- number
@param {string} sign			- number of significant figures
@return {string}				- formatted number in latex
'''
def fmt_number (nr, sign = None):

	sign_default = 3
	s = "NaN"

	n = nr.n if hasattr(nr, 'n') else nr
	if (np.abs(n) > 9999 or np.abs(n) < 0.01) and n != 0:
		if hasattr(nr, 'n'):
			sf = sign if isinstance(sign, int) else sig_figures2(nr)
			float_str = to_precision(nr.n, sf, True)
			base, exponent = float_str.split("e")
			sign_err = sf_ac(base)
			base, exponent = float(base), int(exponent)
			err = nr.s / 10**exponent
			ac = str(sf_ac(to_precision(base, sf, False)))
			fmtNull = (r"{0:." + ac + r"f}").format(0.0)
			base = to_precision(base, sf, False)

			if (nr.n == 10**exponent and sf == 1):
				s = (r"10^{{{0}}} \pm {1}").format(exponent, fmt_number(nr.s, sf))
			elif (r"{0:." + ac + r"f}").format(err) == fmtNull:
				s = (r"{0} \times 10^{{{1}}} \pm {2}").format(base, exponent, fmt_number(nr.s, sf))
			else:
				# correction of the sig. figures of the error if the error and the nominal
				# value of very close to each other
				sign_err = sign_err + 1 if np.abs(err - float(base)) < 9 else sign_err
				s = (r"\left( {0} \pm {1:." + ac + r"f} \right) \times 10^{{{2}}}").format(base, float(to_precision(err, sign_err, False)), exponent)

		else:
			sf = sign if isinstance(sign, int) else sign_default
			float_str = to_precision(nr, sf, True)
			base, exponent = float_str.split("e")
			exponent = int(exponent)
			if (nr == 10**exponent and sign == 1):
				s = (r"10^{{{0}}}").format(exponent)
			else:
				s = (r"{0} \times 10^{{{1}}}").format(base, exponent)
	else:
		if hasattr(nr, 'n'):
			sf = sign if isinstance(sign, int) else sig_figures2(nr)

			ac = str(sf_ac(to_precision(nr.n, sf, False)))
			fmtNull = (r"{0:." + ac + r"f}").format(0.0)

			if (nr.s == 0.0):
				s = to_precision(nr.n, sf, False)
			elif (r"{0:." + ac + r"f}").format(nr.s) == fmtNull:
				s = (r"{0} \pm {1}").format(to_precision(nr.n, sf, False), fmt_number(nr.s, sf))
			else:
				s = (r"{0} \pm {1:." + ac + r"f}").format(to_precision(nr.n, sf, False), nr.s)
		elif isinstance(nr, int):
			s = str(nr)
		else:
			sf = sign if isinstance(sign, int) else sign_default
			s = to_precision(nr, sf, False)

	return s


'''
Format a 2d-array to a latex table
@param {array} m	- 2d-array with rows and columns
@return {string}	- formatted table in latex
'''
def fmt_table (m):
	#m = format(np.atleast_2d(m)) if fmt else np.atleast_2d(m)
	m = np.atleast_2d(m)
	t = []
	for i, v in enumerate(m):
		t.append(' & '.join(str(x) for x in v))
	
	tbl = ''
	tbl += t[0] + '\\\\\n\\hline\n\hline\n'
	t.pop(0)
	tbl += '\\\\\n'.join(str(x) for x in t)
	return tbl

'''
Format a string, array or a number in latex
@param {string|array|int|float|ufloat} s 	- element to format
@param {bool} fmt							- whether to put '$'-signs around numbers/etc...
@return {string}							- formatted element in latex
'''
def format(s, fmt=False):
	if isinstance(s, str): # it is a string
		return s

	elif hasattr(s, "__len__"): # it is an array
		s = np.atleast_2d(s)
		for i, a in enumerate(s):
			for j, b in enumerate(a):
				s[i][j] = format(b, True)
		return fmt_table(s)		
	elif (hasattr(s, 'n')): # it has uncertainties
		return "$" + fmt_number(s) + "$" if fmt else fmt_number(s)
		#return "$" + '{:L}'.format(s) + "$" if fmt else '{:L}'.format(s)
	elif isinstance(s, int) or isinstance(s, float):
		return "$" + fmt_number(s) + "$" if fmt else fmt_number(s)
	else: # it is a regular number (float, int)
		return str(s)

'''
Create latex equation for given polynomial coefficients
@param {array} C 	- array containig the coefficients (highest power first)
@param {int} deg	- degree of the resulting latex function
@return {string}	- the formatted latex equation
'''
def fmt_fit(C, deg=None, var='x'):
	if (deg == None):
		deg = len(C)-1

	pf = ""
	for i, c in (enumerate(C)):
		e = len(C) - (i+1) # exponent
		if e > deg:
			continue

		if c != 0:
			coeff = fmt_number(c) if c > 0 else (fmt_number(-c) if c < 0 else "")
			pf += " + " if (c > 0 and deg != e) else ""
			pf += " - " if c < 0 else ""
			pf += r"\left( " if ( ( hasattr(c, 'n') or r'\times' in coeff ) and e != 0) else ""
			#pf += fmt_number(c) if c > 0 else ""
			#pf += fmt_number(-c) if c < 0 else ""
			pf += coeff
			pf += r"\right) " if ( (hasattr(c, 'n') or r'\times' in coeff ) and e != 0) else ""
			pf += r" " + str(var) if e == 1 else ""
			pf += r" " + str(var) + r"^{" + str(e) +"}" if e > 1 else ""

	return pf


'''
Replace an identifier in the tex document with the replace value
@param {string} s 			- indentifier to search the document for
@param {mixed} r			- replace value
@param {bool} fmt 			- whether to put '$'-signs around numbers/etc... defaults to False
@param {string} texfile		- the tex-file to search in
@return {void}
'''
first = True
def replace(s, r, fmt=False, texfile = "main.tex"):
	global first

	file = texfile if first else 'log/intermediate.tex'

	with open(file, 'rU') as file :
		filedata = file.read()

	# Replace the target string
	filedata = filedata.replace(str("[[" + s + "]]"), format(r, fmt))

	# create the log directory if it does not already exist
	if not os.path.exists('log'):
		os.makedirs('log')

	# Write the file out again
	with open('log/intermediate.tex', 'w') as file:
		file.write(filedata)

	first = False

def compile (pdffile = "main", converter = "pdflatex"):
	#call([converter, '-jobname=' + pdffile, 'log/intermediate.tex'], stdin=None, stdout=None, stderr=None)

	# replacing all TODOs with yellow boxes
	pattern = re.compile("\[\[(TODO[^\]]*)\]\]")
	for i, line in enumerate(open('main.tex')):
		for match in re.finditer(pattern, line):
			replace(str(match.groups()[0]), todo(str(match.groups()[0])))

	if not os.path.exists('log'):
		os.makedirs('log')
	p = Popen([converter, '-halt-on-error', '-jobname=' + pdffile, 'log/intermediate.tex'], stdin=None, stdout=PIPE, stderr=PIPE)
	output, err = p.communicate()
	rc = p.returncode

	if rc != 0:
		print("Compilation of main.tex failed:")
		print(output)
	else:
		os.system("evince main.pdf &")
		#os.remove(pdffile + ".aux")
		#os.rename(pdffile + ".aux", "log/" + pdffile + ".aux")
		#os.rename(pdffile + ".log", "log/" + pdffile + ".log")

def todo(st):
	return r'\colorbox{yellow!30}{' + str(st) + r'}'

'''
Wrapper for np.polyfit to cope with the ufloat datatype
@param {array} x 		- array with x values (can be an array of ufloats)
@param {array} y		- array with y values (can be an array of ufloats)
@param {int} deg 		- degree of the polynomial fit
@param {mixed} kwargs	- arguments for np.polyfit
@return {uarray}		- the coefficients as an array of ufloats (highest powert first)
'''
def phpolyfit (x, y, deg, **kwargs):
	# numpy polyfit needs float type arrays to process
	x = np.array(nominal(x), dtype=float)
	y = np.array(nominal(y), dtype=float)
	c, cov = np.polyfit(x, y, deg, cov=True, **kwargs)
	
	if True in np.isnan(np.sqrt(np.diagonal(cov))):
		print("Something went wrong in calculating the stddev of the coeffs:")
		print(x, y)
		cov = cov_matrix(x, y, deg, **kwargs)

	stderr = np.sqrt(np.diagonal(cov))
	sf = np.zeros(c.size).astype(int) + 3

	return pharray(c, stderr, sf)


def cov_matrix (x, y, deg, **kwargs):
	x_new = np.concatenate((x, x[-1:]))
	y_new = np.concatenate((y, y[-1:]))
	weights = np.concatenate((np.ones(x.size, dtype=float), np.array([sys.float_info.epsilon])))
	_, cov = np.polyfit(x_new, y_new, deg, w=weights, cov=True, **kwargs)
	return cov

# From http://randlet.com/blog/python-significant-figures-format/
def to_precision(x, p, force=None):
	"""
	returns a string representation of x formatted with a precision of p

	Based on the webkit javascript implementation taken from here:
	https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
	"""

	x = float(x)
	out = []

	if p <= 0:
		return 0

	if x == 0.:
		out.append("0")
		if (p - 1 > 0):
			out.append(".")
			out.append("0"*(p - 1))
		return "".join(out)

	if x < 0:
		out.append("-")
		x = -x

	e = int(math.log10(x))
	tens = math.pow(10, e - p + 1)
	n = math.floor(x/tens)

	if n < math.pow(10, p - 1):
		e = e - 1
		tens = math.pow(10, e - p + 1)
		n = math.floor(x / tens)

	if abs((n + 1.) * tens - x) <= abs(n * tens - x):
		n = n + 1

	if n >= math.pow(10,p):
		n = n / 10.
		e = e + 1

	m = "%.*g" % (p, n)

	if (e < -2 or e >= p or force == True) and force != False:
		out.append(m[0])
		if p > 1:
			out.append(".")
			out.extend(m[1:p])
		out.append('e')
		if e > 0:
			out.append("+")
		out.append(str(e))
	elif e == (p - 1):
		out.append(m)
	elif e >= 0:
		out.append(m[:e + 1])
		if e+1 < len(m):
			out.append(".")
			out.extend(m[e + 1:])
	else:
		out.append("0.")
		out.extend(["0"]*-(e + 1))
		out.append(m)

	return "".join(out)