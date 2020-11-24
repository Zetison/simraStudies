import numpy as np
import click
import matplotlib.pyplot as plt


def parseLogFile(filename):
	with open(filename) as f:
		content = f.readlines()
	
	#df = pd.DataFrame()
	df = {}
	for i in range(0,len(content)):
		if content[i].find('=') != -1:
			temp = content[i].split('=')
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
	plt.semilogy(df['dU'], label = "dU")
	plt.semilogy(df['dV'], label = "dV")
	plt.semilogy(df['dW'], label = "dW")
	plt.semilogy(df['K'], label = "K")
	plt.semilogy(df['epsilon'], label = "epsilon")
	plt.legend()
	plt.xlabel('iteration')
	plt.ylabel('Residual')
	plt.show() 

if __name__ == '__main__':
    main()

