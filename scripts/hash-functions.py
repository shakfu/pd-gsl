def hash1(str):
   h = 0
   for c in str:
       h += ord(c)
   return h

def hash2(str):
   h = 0
   for c in str:
      h += (h << 1) + ord(c)
   return h

def hash3(str):
   h = 5381
   for c in str:
      h = ((h << 5) + h) + ord(c)
   return h



