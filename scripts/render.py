#!/usr/bin/env python3 

from mako.template import Template

SKIP_FUNCS = [
   'rando',
   'airy_ai',
   'airy_bi',
]


items = [

   #prefix  name           func_name
   (1,      'log1p',       None),
   (1,      'expm1',       None),
   (2,      'hypot',       None),
   (3,      'hypot3',      None),
   (1,      'acosh',       None),
   (1,      'asinh',       None),
   (1,      'atanh',       None),
   (2,      'ldexp',       None),
   (2,      'pow_int',     None),
   (1,      'pow_2',       None),
   (1,      'pow_3',       None),
   (1,      'pow_4',       None),
   (1,      'pow_5',       None),
   (1,      'pow_6',       None),
   (1,      'pow_7',       None),
   (1,      'pow_8',       None),
   (1,      'pow_9',       None),
   (2,      'rando',       None),
   (1,      'airy_ai',     'sf_airy_Ai'),
   (1,      'airy_bi',     'sf_airy_Bi'),
   (1,      'bessel',      'sf_bessel_J0'),
   (1,      'clausen',     'sf_clausen'),
   (1,      'dawson',      'sf_dawson'),
   (1,      'debye_1',     'sf_debye_1'),
   (1,      'debye_2',     'sf_debye_2'),
   (1,      'debye_3',     'sf_debye_3'),
   (1,      'debye_4',     'sf_debye_4'),
]



# def hash(str):
#    h = 0
#    for c in str:
#        h += ord(c)
#    return h

def hash(str):
   h = 5381
   for c in str:
      h = ((h << 5) + h) + ord(c)
   return h





class Func:
   def __init__(self, nargs, name, func_name):
      self.nargs = nargs
      self.name = name
      self.func_name = func_name

   @property
   def hashed(self):
      h = 5381
      for c in self.name:
         h = ((h << 5) + h) + ord(c)
      return h      

   @property
   def ftype(self):
      return {
         1: 'ufunc',
         2: 'bfunc',
         3: 'tfunc',
      }[self.nargs]

   @property
   def slots(self):
      _slots = ['A_DEFFLOAT'] * self.nargs
      return ", ".join(_slots)

   @property
   def fullname(self):
      if self.func_name:
         return self.func_name
      else:
         return self.name


def render(**kwds):
    tmpl = Template(filename='scripts/psl.c.mako')
    rendered = tmpl.render(**kwds)
    return rendered

funcs = [Func(i[0],i[1], i[2]) for i in items]

with open('psl.c', 'w') as f:
   f.write(render(funcs=funcs, skip=SKIP_FUNCS))
