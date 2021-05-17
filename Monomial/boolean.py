import PySimpleGUI as sg
import numpy as np
from math import log, ceil
from beautifultable import BeautifulTable
from itertools import permutations, combinations

def derivative(func, u):
	res = ''
	for i in range(len(func)):
		res += str((int(func[i]) + int(func[u^i]))%2) 
	return to_vector(polynom([int(k) for k in res]))

def non_degenerate(func):
	if degree(polynom(func)) == 1:
		return 'Нет, аффинная функция.'
	if func.count(1)%2 == 1:
		return "Да"
	max_deg = degree(polynom(func))-1
	for i in range(1,len(func)):
		if degree(polynom(derivative(func, i))) < max_deg:
			return "Нет, deg(Duf) = "+str(degree(polynom(derivative(func, i))))+' при u = '+str(format(i,'b').zfill(int(log(len(func),2))))
	return "Да"	

def is_bent(func):
	for i in range(len(func)):
		if not(abs(W_H_transform(func, i)) == 2**(int(log(len(func),2))/2)):
			return 'Нет'
	return 'Да'

def weight(func):
	res = 0
	for i in func:
		if int(i) == 1:
			res += 1
	return res

def M_n_d(n,d,func):
	var = []
	for i in range(n):
		var.append('x'+str(i+1))
	row = [1]
	for i in var:
		row.append(i)
	for i in range(2,d+1):
		combs = list(combinations(var,i))
		for j in range(len(combs)):
			buf = ''
			for k in range(len(combs[j])):
				buf+=combs[j][k]
			combs[j] = buf
		for j in combs:
			row.append(j)
	supp = []
	for k in range(len(func)):
		if func[k] == 1:
			supp.append(k)
	mat = [[0 for k in range(len(row))] for j in range(len(supp))]
	for k in range(len(supp)):
		for l in range(len(row)):
			z = []
			if l == 0:
				mat[k][l] = 1
				continue
			for a in row[l]:
				if '0'<=a<='9':
					z.append(a)
			flag = 0
			for a in z:
				if format(supp[k],'b').zfill(n)[int(a)-1] == '0':
					flag = 1
			if flag:
				mat[k][l] = 0
			else:
				mat[k][l] = 1
	return mat

def a_d(func,n):
	if func.count(1) == 0:
		return 0
	rev_func = [0 for k in range(len(func))]
	for i in range(len(rev_func)):
		rev_func[i] = (func[i]+1)%2
	for i in range(1, ceil(n/2)):
		z = np.linalg.matrix_rank(np.array(M_n_d(n,i,func)))
		zz = np.linalg.matrix_rank(np.array(M_n_d(n,i,rev_func)))
		maxv = len(M_n_d(n,i,rev_func)[0])
		if z<maxv or zz<maxv:
			return i
	return ceil(n/2)

def correlation_immunity(func):
	for m in range(1,len(func)):
		ms = []
		for i in range(len(func)):
			u = format(i,'b').zfill(int(log(len(func),2)))
			if 1 <= weight(u) <= m:
				ms.append(i)
		for i in ms:
			if not(W_H_transform(func, i) == 0):
				return m-1
	return m

def W_H_transform(func, u):
	res = 0
	uu = format(u, 'b').zfill(int(log(len(func),2)))
	prod = 0
	x = 0
	for i in range(len(func)):
		x = format(i,'b').zfill(int(log(len(func),2)))
		for j in range(len(x)):
			prod += (int(x[j])*int(uu[j]))
		res += (-1)**((func[i]+(prod%2))%2)
		prod = 0
	return res

def spectre(func):
	res =[]
	for i in range(len(func)):
		res.append(W_H_transform(func,i))
	return res

def degree(pol):
	counter = 0
	deg = 0
	if pol == "1" or pol == "0":
		return 0
	for i in pol:
		if i == "+":
			if counter > deg:
				deg = counter
			counter = 0
		else:
			if i > '9':
				counter += 1
	if counter > deg:
		deg = counter
	return deg

def jevons_eq(func, n):
	s = set()
	perms = list(permutations([i for i in range(n)]))
	for i in list(perms):
		for j in range(2**n):
			a = format(j,'b').zfill(n)
			ff=[]
			for k in range(2**n):
				kk = format(k,'b').zfill(n)
				kkk=''
				for l in range(n):
					kkk += str((int(kk[i[l]]) + int(a[l]))%2)
				ff.append(func[int(kkk,2)])
			s.add(tuple(ff))
	return s, len(s)

def tabl(polynom, func, n):
	s = set()
	for i in range(len(polynom)):
		if 'a' <= polynom[i] <= 'z' or 'A' <= polynom[i] <= 'Z' :
			ss = polynom[i]
			j = 1
			while (i+j) <= len(polynom)-1 and "0" <= polynom[i+j] <= "9":
				ss += polynom[i+j]
				j += 1
			s.add(ss)
	if len(s) < n:
		s = ['x'+str(i) for i in range(1,n+1)]
	table = BeautifulTable()
	table.columns.header = [i for i in sorted(s)] + ['f(x)']
	for i in range(2**n):
		table.rows.append([int(k) for k in format(i,'b').zfill(n)]+[func[i]])
	return table

def to_vector(polynom):
	s = set()
	ss=''
	for i in range(len(polynom)):
		if '0'<=polynom[i]<='9' and polynom[i-1] == '+':
			return 'Функция введена неправильно!'
		if 'a' <= polynom[i] <= 'z' or 'A' <= polynom[i] <= 'Z' :
			ss = polynom[i]
			j = 1
			while (i+j) <= len(polynom)-1 and "0" <= polynom[i+j] <= "9":
				ss += polynom[i+j]
				j += 1
			s.add(ss)
			ss =''
	res_func = []
	if len(s) == 0:
		return [int(polynom) for k in range(8)]
	for i in range(2**len(s)):
		res = polynom
		u = format(i,'b').zfill(len(s))
		for j in range(len(s)-1, -1, -1):
			res = res.replace(sorted(s)[j], u[j])
		sub_res = 1
		func_res = 0
		for k in res:
			if not k == "+":
				sub_res = (sub_res * int(k))
			else:
				func_res = (func_res + sub_res)%2
				sub_res = 1
		func_res = (func_res+sub_res)%2
		res_func.append(func_res)
	return res_func

def is_hex(s):
	try:
		int(s, 16)
		return True
	except ValueError:
		return False

def is_bin(func):
	res = True
	for i in func:
		if not i == '1' and not i == '0':
			res = False
	return res

def polynom(func):
	buf = func
	if not(1 in func):
		return "0"
	buf_2 = []
	res = [func[0]]
	while not(len(buf_2) == 1):
		buf_2 = []
		for i in range(len(buf)-1):
			buf_2.append((buf[i]+buf[i+1])%2)
		res.append(buf_2[0])
		buf = buf_2
	res_2 = ''
	vars = ["x"+str(i+1) for i in range(int(log(len(func),2)))]
	for i in range(len(func)):
		if i == 0 :
			res_2 += res[0]*"1"
			if len(res_2) > 0:
				res_2 += "+"
		else:
			if res[i] > 0:
				var = ''
				for j in range(int(log(len(func),2))):
					var += vars[j]*int(format(i,'b').zfill(int(log(len(func),2)))[j])
				res_2 += var+'+'  


	return res_2[:len(res_2)-1]

def func_read(func):
	if func[0:2] == '0x' and is_hex(func):
		n = int(ceil(log(len(format(int(func[2:], 16), 'b')),2)))
		if  n == 0:
			return 'n = 0!'
		res = [int(k) for k in format(int(func[2:], 16), 'b').zfill(2**n)]
		return res, polynom(res), n
	elif is_bin(func):
		if log(len(func),2).is_integer() and len(func)>1:
			res = [int(k) for k in func]
			return res, polynom(res), int(log(len(func),2))
		else:
			return 'Вектор должен состоять из 2^n значений!' 
	elif func[0] > '9' or func[0] < '0':
		res = to_vector(func)
		if res == 'Функция введена неправильно!':
			return res
		return res, func, int(log(len(res), 2))
	else:
		return 'Функция введена неправильно!'
layout = [
	[sg.Text('f(x) = '), sg.InputText()],
	[sg.Checkbox('Таблица истинности'), sg.Checkbox('Эквивалентные функции')],
	[sg.Output(size=(100, 30), key = 'out')],
	[sg.Button('Запуск'), sg.Button('Очистить'), sg.Button('Выход')]
]
window = sg.Window('Параметры булевой функции', layout)

while True:                           
	event, values = window.read()
	if event in (None, sg.WIN_CLOSED, 'Выход'):
		break
	if event == 'Очистить':
		window.FindElement('out').Update('')
	if event == 'Запуск':
		if not values[0]:
			print('Введите функцию!')
			continue
		func = func_read(values[0])
		if len(func) > 3:
			print(func)
			continue
		print('Введенная функция:', values[0])
		print('Вектор значений:', ''.join(str(k) for k in func[0]))
		print('АНФ:', func[1])
		print('Параметры функции:')
		print('deg(f) =', degree(func[1]))
		print('wt(f) =', func[0].count(1))
		print('nl(f) =', int(2**(func[2]-1)-1/2*(max([abs(k) for k in spectre(func[0])]))))
		print('cor(f) =', correlation_immunity(func[0]))
		print('AI(f) =', a_d(func[0], func[2]))
		print('Бент:', is_bent(func[0]))
		print('Невырождена:', non_degenerate(func[0]))
		if values[1]:
			print('Таблица истинности:')
			print(tabl(func[1], func[0], func[2]))
		if values[2]:
			print('Эквивалентные функции:')
			s = jevons_eq(func[0], func[2])
			c = 1
			for i in s[0]:
				print(c,'. f(x) = ', polynom(list(i)), sep ='')
				c+=1
		print(168*'-')
window.close()