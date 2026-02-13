with open('Cmat.txt','r') as f:
    cmat = f.read()
f.close()

cmat = cmat.replace(' ','')
cmat = cmat.replace('[','')
cmat = cmat.replace(']','')
cmat = cmat.replace(',\n\n','\n')
cmat = cmat.replace('\n\n','\n')

with open('Cmat.txt','w') as g:
    g.write(cmat)