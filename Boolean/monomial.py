import PySimpleGUI as sg
import numpy as np
from math import log, ceil
from beautifultable import BeautifulTable
from itertools import permutations, combinations
from finitefield import *

def trace(g,n,l,k):
	GF2_n = FiniteField(2, g)
	lambd = FiniteFieldElt(GF2_n, l)
	res = []
	for i in range(2**n):
		b = FiniteFieldElt(GF2_n, [int(z) for z in format(i,'b').zfill(n)])
		summa = FiniteFieldElt(GF2_n, [0]*n)
		for j in range(n):
			summa+=(lambd*b**k)**(2**j)
		res.append(int(str(summa)[1]))
	return res

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

def is_bin(func):
	res = True
	for i in func:
		if not i == '1' and not i == '0':
			res = False
	return res

def is_int(s):
	try:
		int(s)
		return True
	except ValueError:
		return False

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

def check(g, n, l, k):
	if not is_int(n):
		return "n должно быть натуральным числом!"
	if not is_int(k):
		return "k должно быть целым числом!"
	elif not is_bin(g):
		return 'Введите двоичный вектор в поле g(x)!'
	elif not is_bin(l):
		return 'Введите двоичный вектор в поле l!'
	elif not len(l) == int(n) or not len(g) == int(n) :
		return "l и g(x) должны иметь длину n!"
	else:
		return [[int(z) for z in g], [int(z) for z in l], int(n), int(k)]


layout = [
	[sg.Text('g = '), sg.InputText(), sg.Text('n = '), sg.InputText()],
	[sg.Text('l = '), sg.InputText(), sg.Text('k = '), sg.InputText()],
	[sg.Output(size=(100, 30), key = 'out')],
	[sg.Button('Запуск'), sg.Button('Очистить'), sg.Button('Выход')]
]

window = sg.Window('Мономиальные булевы функции', layout)

while True:                           
	event, values = window.read()
	if event in (None, sg.WIN_CLOSED, 'Выход'):
		break
	if event == 'Очистить':
		window.FindElement('out').Update('')
	if event == 'Запуск':
		if not values[0] or not values[1] or not values[2] or not values[3]:
			print('Заполните все поля!')
			continue
		inp = check(values[0], values[1], values[2], values[3])
		if len(inp) > 4:
			print(inp)
			continue
		func = trace(inp[0],inp[2],inp[1], inp[3])
		print('Вектор значений:', ''.join(str(k) for k in func))
		print('Бент:', is_bent(func))
		print('Невырождена: ', non_degenerate(func))
		print(168*'-')
window.close()



























