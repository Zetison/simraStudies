import numpy as np
import os
import click
import matplotlib.pyplot as plt


def parseLogFile(filename):
	with open(filename) as f:
		content = f.readlines()
	
	#df = pd.DataFrame()
	df = {}
	for i in range(0,len(content)):
		sep_eq = content[i].find('=') != -1
		sep_col = content[i].find(':') != -1
		if sep_eq or sep_col:
			if sep_eq:
				temp = content[i].split('=')
			else:
				temp = content[i].split(':')

			LHS = temp[0].split(',')
			RHS = " ".join(temp[1].split()).split(' ')
			for j in range(0,len(LHS)):
				LHS_j = LHS[j]
				LHS_j = LHS_j.replace('residuals', '').replace('residual', '').replace('(', '').replace(')', '').replace(' ', '').replace('.','_').replace('dP','').replace('_pT','').strip()
				try:
					RHS_j = float(RHS[j])
				except:
					RHS_j = RHS[j]

				try:
					df[LHS_j] = np.append(df[LHS_j],RHS_j)
				except:
					df[LHS_j] = np.array(RHS_j)

	return df


@click.command()
@click.option('--infile', type=str)
def main(infile):
	df = parseLogFile(infile)
	U = df['L2-normu1']
	dU = (U-U[-1])/U[-1]
	V = df['L2-normu2']
	dV = (V-V[-1])/V[-1]
	W = df['L2-normu3']
	dW = (W-W[-1])/W[-1]
	P = df['L2-normpd']
	dP = (P-P[-1])/P[-1]
	k = df['L2-normtk']
	dk = (k-k[-1])/k[-1]
	d = df['L2-normtd']
	dd = (d-d[-1])/d[-1]
	v = df['L2-normvt']
	dv = (v-v[-1])/v[-1]
	plt.semilogy(dU, label = "L2-norm of U")
	plt.semilogy(dV, label = "L2-norm of V")
	plt.semilogy(dW, label = "L2-norm of W")
	plt.semilogy(dP, label = "L2-norm of pd")
	plt.semilogy(dk, label = "L2-norm of tk")
	plt.semilogy(dd, label = "L2-norm of td")
	plt.semilogy(dv, label = "L2-norm of vt")
	plt.legend()
	plt.xlabel('iteration')
	plt.ylabel('Relative error compared to final iteration')
	plt.savefig(os.path.dirname(infile)+'/L2errors.png', dpi=300)
#	plt.ylim([1e-15,1])
	plt.show() 

if __name__ == '__main__':
    main()

