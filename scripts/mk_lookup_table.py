from mako.template import Template

names = [
   'rando',
   'airy_ai',
   'bessel',
   'clausen',
   'dawson',
   'debye_1',
   'debye_2',
   'debye_3',
   'debye_4',
]

def hash(str):
   h = 0
   for c in str:
       h += ord(c)
   return h

def render(template, **kwds):
    templ = Template(text=template)
    rendered = templ.render(**kwds)
    return rendered


tmpl = """

enum FUNC {
   % for name in names:
   ${name.upper()} = ${hash(name)},
   % endfor
};


// set default function from symbol
void select_default_function(t_symbol *s) {

   switch (hash(s->s_name)) {
      % for name in names:
      case ${name.upper()}:
         x->ufunc = &psl_${name};
         break;
      % endfor
      default:
         break;
   }
}
"""

print(render(tmpl, names=names, hash=hash))